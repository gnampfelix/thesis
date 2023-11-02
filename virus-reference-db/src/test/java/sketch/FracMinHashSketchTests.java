package sketch;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.hasItems;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import org.husonlab.sketch.FracMinHashSketch;
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
            FracMinHashSketch sketch = FracMinHashSketch.compute("test", Collections.singleton(sequence), true, 1.0f / 10.0f, 21, 42, true, true, new ProgressSilent());
            
            // the test file contains the kmer "TTGGATGAAACGCACCCGCTAT". For
            // this, the reverse complement is "ATAGCGGGTGCGTTTCATCCA", which
            // is smaller.
            long expectedHash = MurmurHash.hash64("ATAGCGGGTGCGTTTCATCCA".getBytes(), 0, 21, 42);
            
            // The expected hash is -8556973683090215967, which is lesser than 
            // Long.MIN + (Long.MAX * s * 2) with s = 1/10 as used above.
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
}
