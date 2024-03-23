package org.husonlab.fmhdist.sketch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.hasItems;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.husonlab.fmhdist.util.FastKMerIterator;
import org.junit.Ignore;
import org.junit.Test;

import jloda.thirdparty.MurmurHash;

public class FracMinHashSketchTests {
    /**
     * @throws IOException
     */
    @Test
    public void shouldCalculateFracMinSketch() throws IOException {
        try (FastKMerIterator kmers = new FastKMerIterator(21, "src/test/resources/virus1.fasta", true)) {
			FracMinHashSketch sketch = FracMinHashSketch.compute("test", kmers, true, 21, 42);
            
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
    @Ignore // Not a real test, but for understanding the performance
    public void benchmarkBloomPerformance() throws IOException {
        String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/026/225/685/GCA_026225685.1_Pinf1306_UCR_1/GCA_026225685.1_Pinf1306_UCR_1_genomic.fna.gz";
        int k = 21;
        int s = 1000;
        long initialHeap = Runtime.getRuntime().totalMemory();
        long start = System.currentTimeMillis();
        
        try (FastKMerIterator kmers = new FastKMerIterator(k, url, true)) {
            FracMinHashSketch.compute(
                "test",
                kmers,
                true,
                s,
                42
            );
            long finalHeap = Runtime.getRuntime().totalMemory();
            long end = System.currentTimeMillis();
            System.out.println(String.format("runtime: %d", (end-start)/1000));
            System.out.println(String.format("initial heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, finalHeap / 1024/1024));
        }
    }

    // @Ignore // Not a real test, but for understanding the performance
    @Test
    public void testSegments() throws IOException {
        String url = "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/test_genomes/PhyInf1.fna.gz";
        long initialHeap = Runtime.getRuntime().totalMemory();
        long start = System.currentTimeMillis();
        try (FastKMerIterator kmers = new FastKMerIterator(21, url, true)) {
			FracMinHashSketch sketch = FracMinHashSketch.compute("test", kmers, true, 1000, 42);
            long finalHeap = Runtime.getRuntime().totalMemory();
            long end = System.currentTimeMillis();
            System.out.println(String.format("sketch size: %d", sketch.getValues().length));
            System.out.println(String.format("runtime: %d", (end-start)/1000));
            System.out.println(String.format("initial heap: %d\nfinal heap: %d", initialHeap / 1024 / 1024, finalHeap / 1024/1024));
        }        
    }

}
