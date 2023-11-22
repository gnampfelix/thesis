package sketch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.hasItems;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.husonlab.sketch.FracMinHashSketch;
import org.husonlab.sketch.SequenceGrouper;
import org.husonlab.util.KMerIterator;
import org.junit.Ignore;
import org.junit.Test;

import jloda.thirdparty.MurmurHash;
import jloda.util.FileLineBytesIterator;
import jloda.util.progress.ProgressSilent;

public class FracMinHashSketchTests {
    /**
     * @throws IOException
     */
    @Test
    public void shouldCalculateFracMinSketch() throws IOException {
        try (FileLineBytesIterator it = new FileLineBytesIterator("src/test/resources/virus1.fasta")) {
			KMerIterator kmers = new KMerIterator(it, 21);
            FracMinHashSketch sketch = FracMinHashSketch.compute("test", kmers, 420, true, 21, 42, true, true, new ProgressSilent());
            
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
    public void compareSketches() throws IOException {
        String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/225/685/GCA_026225685.1_Pinf1306_UCR_1/GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz";
        int k = 21;
        int s = 1000;
        
        SequenceGrouper group = new SequenceGrouper(url);
        FracMinHashSketch s1 = FracMinHashSketch.compute(
            "test", 
            group.getGroups(), 
            true,
            s,
            k,
            42,
            false,
            false,
            new ProgressSilent()
        );

        FracMinHashSketch s2;
        try (FileLineBytesIterator it = new FileLineBytesIterator(url)) {
            KMerIterator kmers = new KMerIterator(it, k);
            s2 = FracMinHashSketch.compute(
                "test",
                kmers,
                246895792,
                true,
                s,
                42,
                false,
                false,
                new ProgressSilent()
            );
        }
        assertEquals(s1.getValues().length, s2.getValues().length);
        // order should also be the same!
        assertArrayEquals(s1.getValues(), s2.getValues());
    }

    @Test
    @Ignore // Not a real test, but for understanding the performance
    public void benchmarkBloomPerformance() throws IOException {
        String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/225/685/GCA_026225685.1_Pinf1306_UCR_1/GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz";
        int k = 21;
        int s = 1000;
        long initialHeap = Runtime.getRuntime().totalMemory();
        long start = System.currentTimeMillis();
        
        try (FileLineBytesIterator it = new FileLineBytesIterator(url)) {
            KMerIterator kmers = new KMerIterator(it, k);
            FracMinHashSketch.compute(
                "test",
                kmers,
                246895792,
                true,
                s,
                42,
                true,
                false,
                new ProgressSilent()
            );
            long finalHeap = Runtime.getRuntime().totalMemory();
            long end = System.currentTimeMillis();
            System.out.println(String.format("runtime: %d", (end-start)/1000));
            System.out.println(String.format("initial heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, finalHeap / 1024/1024));
        }
    }

    @Test
    @Ignore // Not a real test, but for understanding the performance
    public void testSegments() throws IOException {
        String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/225/685/GCA_026225685.1_Pinf1306_UCR_1/GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz";
        long initialHeap = Runtime.getRuntime().totalMemory();
        long start = System.currentTimeMillis();
        try (FileLineBytesIterator it = new FileLineBytesIterator(url)) {
			KMerIterator kmers = new KMerIterator(it, 21);
            FracMinHashSketch sketch = FracMinHashSketch.compute("test", kmers, 246895792, true, 1000, 42, false, false, new ProgressSilent());
            long finalHeap = Runtime.getRuntime().totalMemory();
            long end = System.currentTimeMillis();
            System.out.println(String.format("sketch size: %d", sketch.getValues().length));
            System.out.println(String.format("runtime: %d", (end-start)/1000));
            System.out.println(String.format("initial heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, finalHeap / 1024/1024));
        }        
    }
}
