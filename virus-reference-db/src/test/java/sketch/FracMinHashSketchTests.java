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
        String url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/731/815/GCF_009731815.1_ASM973181v1/GCF_009731815.1_ASM973181v1_genomic.fna.gz";
        try (FileLineIterator it = new FileLineIterator(url)) {
			byte[] sequence =  it.stream().filter(line -> !line.startsWith(">")).map(line -> line.replaceAll("\\s+", "")).collect(Collectors.joining()).getBytes();
            FracMinHashSketch sketch = FracMinHashSketch.compute("test", Collections.singleton(sequence), true, 10, 21, 42, false, true, new ProgressSilent());
            System.out.println(sketch.getValues().length);
        }
        SequenceGrouper grp = new SequenceGrouper(url);
        FracMinHashSketch sketch = FracMinHashSketch.compute("test", grp.getGroups(), true, 10, 21, 42, false, true, new ProgressSilent());
        System.out.println(sketch.getValues().length);
    }
}
