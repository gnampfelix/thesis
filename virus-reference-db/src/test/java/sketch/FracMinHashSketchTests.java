package sketch;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.Collections;
import java.util.stream.Collectors;

import org.husonlab.sketch.FracMinHashSketch;
import org.junit.Test;

import jloda.util.FileLineIterator;
import jloda.util.progress.ProgressSilent;

public class FracMinHashSketchTests {
    @Test
    public void shouldCalculateFracMinSketch() throws IOException {
        try (FileLineIterator it = new FileLineIterator("src/test/resources/virus1.fasta")) {
			byte[] sequence =  it.stream().filter(line -> !line.startsWith(">")).map(line -> line.replaceAll("\\s+", "")).collect(Collectors.joining()).getBytes();
            FracMinHashSketch sketch = FracMinHashSketch.compute("test", Collections.singleton(sequence), true, 1000, 21, 42, true, false, new ProgressSilent());
		}
    }
}
