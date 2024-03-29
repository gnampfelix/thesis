package org.husonlab.fmhdist.sketch;

import java.io.IOException;
import java.util.logging.Logger;

import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.util.FastKMerIterator;
import org.husonlab.fmhdist.util.KMerIterator;
import org.husonlab.fmhdist.util.LineKMerIterator;

public class GenomeSketch {
    private Genome genome;
    private FracMinHashSketch sketch;
    private static Logger logger = Logger.getLogger(GenomeSketch.class.getName());

    public static GenomeSketch sketch(Genome genome, int kSize, int sParam, int seed) throws IOException {
        logger.fine("Calculating sketch for " + genome.getAccession());
        final GenomeSketch result = new GenomeSketch(genome);
        try (KMerIterator kmers = new LineKMerIterator(kSize, genome.getFastaUrl(), true)) {
            result.sketch = FracMinHashSketch.compute(genome.getAccession(), kmers, true, sParam, seed);
        }
        return result;
    }

    private GenomeSketch(Genome genome) {
        this.genome = genome;
    }

    public GenomeSketch(Genome genome, FracMinHashSketch sketch) {
        this.genome = genome;
        this.sketch = sketch;
    }

    public Genome getGenome() {
        return this.genome;
    }

    public FracMinHashSketch getSketch() {
        return this.sketch;
    }
}
