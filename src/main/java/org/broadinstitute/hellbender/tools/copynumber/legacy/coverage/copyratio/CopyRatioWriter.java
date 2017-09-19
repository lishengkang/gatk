package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.NamedSampleFile;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;

import java.io.File;
import java.io.IOException;

final class CopyRatioWriter extends TableWriter<CopyRatio> {

    private final String sampleName;

    CopyRatioWriter(final File file,
                           final String sampleName) throws IOException {
        super(file, CopyRatioTableColumn.COLUMNS);
        this.sampleName = Utils.nonNull(sampleName);
    }

    void writeSampleName() {
        try {
            writeComment(NamedSampleFile.SAMPLE_NAME_COMMENT_PREFIX + sampleName);
        } catch (final IOException e) {
            throw new UserException("Could not write sample name.");
        }
    }

    @Override
    protected void composeLine(final CopyRatio record, final DataLine dataLine) {
        Utils.nonNull(record);
        Utils.nonNull(dataLine);
        dataLine.append(record.getInterval().getContig())
                .append(record.getInterval().getStart())
                .append(record.getInterval().getEnd())
                .append(record.getLog2CopyRatioValue());
    }
}