package sketch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.hasItems;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.husonlab.sketch.FracMinHashSketch;
import org.husonlab.sketch.SequenceGrouper;
import org.junit.Test;

import jloda.thirdparty.MurmurHash;
import jloda.util.FileLineIterator;
import jloda.util.progress.ProgressSilent;

public class FracMinHashSketchTests {
    /**
     * @throws IOException
     */
    @Test
    public void shouldCalculateFracMinSketch() throws IOException {
        try (FileLineIterator it = new FileLineIterator("src/test/resources/virus1.fasta")) {
			byte[] sequence =  it.stream().filter(line -> !line.startsWith(">")).map(line -> line.replaceAll("\\s+", "")).collect(Collectors.joining()).getBytes();
            FracMinHashSketch sketch = FracMinHashSketch.compute("test", Collections.singleton(sequence), true, 10, 21, 42, true, true, new ProgressSilent());
            
            // the test file contains the kmer "TTGGATGAAACGCACCCGCTAT". For
            // this, the reverse complement is "ATAGCGGGTGCGTTTCATCCA", which
            // is smaller.
            long expectedHash = MurmurHash.hash64("ATAGCGGGTGCGTTTCATCCA".getBytes(), 0, 21, 42);
            
            // The expected hash is -8556973683090215967, which is lesser than 
            // Long.MIN + (Long.MAX *  1/s * 2) with s = 10 as used above.
            // Thus, we expect it to be included in the result.

            List<Long> values = new ArrayList<>();
            for (long value : sketch.getValues()) {
                values.add(value);
            }
            assertThat(values, hasItems(
                expectedHash
            ));
            System.out.println();
        }
    }

    @Test
    public void testSegments() throws IOException {
        String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/225/685/GCA_026225685.1_Pinf1306_UCR_1/GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz";
        long initialHeap = Runtime.getRuntime().totalMemory();
        long start = System.currentTimeMillis();
        SequenceGrouper grp = new SequenceGrouper(url);
        long groupingHeap = Runtime.getRuntime().totalMemory();
        FracMinHashSketch sketch = FracMinHashSketch.compute("test", grp.getGroups(), true, 1000, 21, 42, false, false, new ProgressSilent());
        long finalHeap = Runtime.getRuntime().totalMemory();
        long end = System.currentTimeMillis();
        System.out.println(String.format("sketch size: %d", sketch.getValues().length));
        System.out.println(String.format("runtime: %d", (end-start)/1000));
        System.out.println(String.format("initial heap: %d\ngrouping heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, groupingHeap / 1024 / 1024, finalHeap / 1024 /1024));
    }
}
