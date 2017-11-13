package org.broadinstitute.hellbender.utils.nio;

import com.google.common.jimfs.Configuration;
import com.google.common.jimfs.Jimfs;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.FileSystem;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

public class PathLineIteratorUnitTest extends GATKBaseTest {

    final String[] opus = {"Hello world", "What's new?"};

    @Test
    public void testLineIteratorInForEach() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path path = jimfs.getPath("test.txt");
            try (BufferedWriter bufferedWriter = Files.newBufferedWriter(path)) {
                bufferedWriter.write(String.join("\n", opus));
            }
            ArrayList<String> got = new ArrayList<>();
            try (PathLineIterator lines = new PathLineIterator(path)) {
                for (String s : lines) {
                    got.add(s);
                }
            }
            Assert.assertEquals(got.toArray(), opus);
        }
    }

    @Test
    public void testLineIteratorAsResource() throws IOException {
        try (FileSystem jimfs = Jimfs.newFileSystem(Configuration.unix())) {
            final Path path = jimfs.getPath("test.txt");
            try (BufferedWriter bufferedWriter = Files.newBufferedWriter(path)) {
                bufferedWriter.write(String.join("\n", opus));
            }
            ArrayList<String> got = new ArrayList<>();
            try (PathLineIterator lines = new PathLineIterator(path)) {
                lines.forEach(s -> got.add(s));
            }
            Assert.assertEquals(got.toArray(), opus);
        }
    }

}