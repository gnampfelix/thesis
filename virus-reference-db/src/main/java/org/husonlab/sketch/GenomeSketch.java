package org.husonlab.sketch;

import java.io.IOException;
import java.util.List;
import java.util.logging.Logger;

import org.husonlab.ncbi.Genome;

public class GenomeSketch {
    private Genome genome;
    private FracMinHashSketch sketch;
    private static Logger logger = Logger.getLogger(GenomeSketch.class.getName());

    public static GenomeSketch sketch(Genome genome, int kSize, int sParam, int seed) throws IOException {
        logger.fine("Calculating sketch for " + genome.getAccession());
        final GenomeSketch result = new GenomeSketch(genome);
        List<byte[]> sequences = new SequenceGrouper(genome.getFastaUrl()).getGroups();
        // result.sketch = FracMinHashSketch.compute(genome.getAccession(), sequences, true, sParam, kSize, seed, true, false, new ProgressSilent());
        // if (result.sketch.getValues().length == 0) {
        //     logger.warning("empty sketch for " + genome.getAccession());
        // }
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
