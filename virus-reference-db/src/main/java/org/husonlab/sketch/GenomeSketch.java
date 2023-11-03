package org.husonlab.sketch;

import java.io.IOException;
import java.util.Collections;
import java.util.logging.Logger;
import java.util.stream.Collectors;

import org.husonlab.ncbi.Genome;

import jloda.util.FileLineIterator;
import jloda.util.progress.ProgressSilent;

public class GenomeSketch {
    private Genome genome;
    private FracMinHashSketch sketch;
    private static Logger logger = Logger.getLogger(GenomeSketch.class.getName());

    public static GenomeSketch sketch(Genome genome, int kSize, int sParam, int seed) throws IOException {
        logger.info("Calculating sketch for " + genome.getAccession());
        final GenomeSketch result = new GenomeSketch(genome);
        try (FileLineIterator it = new FileLineIterator(genome.getDownloadLink())) {
            // TODO: change grouping behaviour
			byte[] sequence =  it.stream().filter(line -> !line.startsWith(">")).map(line -> line.replaceAll("\\s+", "")).collect(Collectors.joining()).getBytes();
            result.sketch = FracMinHashSketch.compute(genome.getAccession(), Collections.singleton(sequence), true, sParam, kSize, seed, false, true, new ProgressSilent());
        }
        if (result.sketch.getValues().length == 0) {
            logger.warning("empty sketch for " + genome.getAccession());
        }
        return result;
    }

    private GenomeSketch(Genome genome) {
        this.genome = genome;
    }

    public Genome getGenome() {
        return this.genome;
    }

    public FracMinHashSketch getSketch() {
        return this.sketch;
    }
}
