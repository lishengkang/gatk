package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import htsjdk.samtools.SAMFlag;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVFastqUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndexCache;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembler;
import org.broadinstitute.hellbender.utils.fermi.FermiLiteAssembly;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.*;

/** This LocalAssemblyHandler aligns assembly contigs with BWA, along with some optional writing of intermediate results. */
public final class FermiLiteAssemblyHandler implements FindBreakpointEvidenceSpark.LocalAssemblyHandler {
    private static final long serialVersionUID = 1L;
    private final String alignerIndexFile;
    private final int maxFastqSize;
    private final String fastqDir;
    private final boolean writeGFAs;

    public FermiLiteAssemblyHandler( final String alignerIndexFile, final int maxFastqSize,
                                     final String fastqDir, final boolean writeGFAs ) {
        this.alignerIndexFile = alignerIndexFile;
        this.maxFastqSize = maxFastqSize;
        this.fastqDir = fastqDir;
        this.writeGFAs = writeGFAs;
    }

    @Override
    public AlignedAssemblyOrExcuse apply( final Tuple2<Integer, List<SVFastqUtils.FastqRead>> intervalAndReads ) {
        final int intervalID = intervalAndReads._1();
        final String assemblyName = AlignedAssemblyOrExcuse.formatAssemblyID(intervalID);
        final List<SVFastqUtils.FastqRead> readsList = intervalAndReads._2();

        final int fastqSize = readsList.stream().mapToInt(FastqRead -> FastqRead.getBases().length).sum();
        if ( fastqSize > maxFastqSize ) {
            return new AlignedAssemblyOrExcuse(intervalID, "no assembly -- too big (" + fastqSize + " bytes).");
        }

        if ( fastqDir != null ) {
            final String fastqName = String.format("%s/%s.fastq", fastqDir, assemblyName);
            final ArrayList<SVFastqUtils.FastqRead> sortedReads = new ArrayList<>(readsList);
            sortedReads.sort(Comparator.comparing(SVFastqUtils.FastqRead::getHeader));
            SVFastqUtils.writeFastqFile(fastqName, sortedReads.iterator());
        }

        final long timeStart = System.currentTimeMillis();
        final FermiLiteAssembly assembly = new FermiLiteAssembler().createAssembly(readsList);
        final int secondsInAssembly = (int)((System.currentTimeMillis() - timeStart + 500)/1000);

        if ( fastqDir != null && writeGFAs ) {
            final String gfaName =  String.format("%s/%s.gfa", fastqDir, assemblyName);
            try ( final OutputStream os = BucketUtils.createFile(gfaName) ) {
                assembly.writeGFA(os);
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Can't write "+gfaName, ioe);
            }
        }

        final List<byte[]> tigSeqs =
                assembly.getContigs().stream()
                        .map(FermiLiteAssembly.Contig::getSequence)
                        .collect(SVUtils.arrayListCollector(assembly.getNContigs()));
        final AlignedAssemblyOrExcuse result;
        try ( final BwaMemAligner aligner = new BwaMemAligner(BwaMemIndexCache.getInstance(alignerIndexFile)) ) {
            aligner.setIntraCtgOptions();
            final List<List<BwaMemAlignment>> alignments = aligner.alignSeqs(tigSeqs);
            result = new AlignedAssemblyOrExcuse(intervalID, assembly, secondsInAssembly, alignments);
        }
        final int nTigs = tigSeqs.size();
        if ( nTigs == 0 ) {
            return result;
        }

        final String tmpDir = "/tmp/fbes";
        try {
            IOUtils.createDirectory(tmpDir);
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't create " + tmpDir, ioe);
        }
        final String fastaFile = String.format("%s/%s.fasta", tmpDir, assemblyName);
        try ( final BufferedOutputStream os = new BufferedOutputStream(BucketUtils.createFile(fastaFile)) ) {
            for (int tigId = 0; tigId != nTigs; ++tigId) {
                final String seqName = assemblyName + "." + tigId;
                final byte[] line1 = (">" + seqName + "\n").getBytes();
                final byte[] seq = tigSeqs.get(tigId);
                os.write(line1);
                os.write(seq);
                os.write('\n');
            }
        }
        catch ( final IOException ioe ) {
            throw new GATKException("Unable to write fasta file of contigs from assembly " + assemblyName, ioe);
        }

        final String imageFile = String.format("%s/%s.img", tmpDir, assemblyName);
        BwaMemIndex.createIndexImageFromFastaFile(fastaFile, imageFile);
        try { BucketUtils.deleteFile(fastaFile); }
        catch ( final IOException ioe ) { throw new GATKException("unable to delete " + fastaFile, ioe); }
        try ( final BwaMemIndex assemblyIndex = new BwaMemIndex(imageFile);
              final BwaMemAligner aligner = new BwaMemAligner(assemblyIndex) ) {
            aligner.alignPairs();
            final List<List<BwaMemAlignment>> alignments =
                    aligner.alignSeqs(readsList, SVFastqUtils.FastqRead::getBases);
            final Map<Link, Integer> linkCounts = new HashMap<>();
            final int nReads = readsList.size();
            for ( int idx = 0; idx < nReads; idx += 2 ) {
                final List<BwaMemAlignment> alignList1 = alignments.get(idx);
                final List<BwaMemAlignment> alignList2 = alignments.get(idx + 1);
                for ( final BwaMemAlignment alignment1 : alignList1 ) {
                    for ( final BwaMemAlignment alignment2 : alignList2 ) {
                        if ( alignment1.getRefId() != alignment2.getRefId() ) {
                            final Link link = new Link(alignment1.getRefId(),
                                                        SAMFlag.READ_REVERSE_STRAND.isSet(alignment1.getSamFlag()),
                                                        alignment2.getRefId(),
                                                        !SAMFlag.READ_REVERSE_STRAND.isSet(alignment2.getSamFlag()));
                            linkCounts.merge(link, 1, Integer::sum);
                        }
                    }
                }
            }
        }
        try { BucketUtils.deleteFile(imageFile); }
        catch ( final IOException ioe ) { throw new GATKException("unable to delete " + imageFile, ioe); }
        return result;
    }

    private static final class Link {
        private final int contigID1;
        private final boolean isRC1;
        private final int contigID2;
        private final boolean isRC2;

        public Link( final int contigID1, boolean isRC1, final int contigID2, boolean isRC2 ) {
            if ( contigID1 < contigID2 ) {
                this.contigID1 = contigID1;
                this.isRC1 = isRC1;
                this.contigID2 = contigID2;
                this.isRC2 = isRC2;
            } else {
                this.contigID1 = contigID2;
                this.isRC1 = !isRC2;
                this.contigID2 = contigID1;
                this.isRC2 = !isRC1;
            }
        }

        @Override public boolean equals( final Object obj ) {
            return obj instanceof Link && equals((Link)obj);
        }

        public boolean equals( final Link link ) {
            return this == link ||
                    (contigID1 == link.contigID1 && isRC1 == link.isRC1 &&
                     contigID2 == link.contigID2 && isRC2 == link.isRC2);
        }

        @Override public int hashCode() {
            int hashVal = 113;
            hashVal = 47 * (hashVal + contigID1);
            hashVal = 47 * (hashVal + (isRC1 ? 0 : 1));
            hashVal = 47 * (hashVal + contigID2);
            return 47 * (hashVal + (isRC2 ? 0 : 1));
        }
    }
}
