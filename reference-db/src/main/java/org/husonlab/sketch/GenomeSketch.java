package org.husonlab.sketch;

import java.io.IOException;
import java.util.logging.Logger;

import org.husonlab.ncbi.Genome;
import org.husonlab.util.FastKMerIterator;

import jloda.util.progress.ProgressListener;

public class GenomeSketch {
    private Genome genome;
    private FracMinHashSketch sketch;
    private static Logger logger = Logger.getLogger(GenomeSketch.class.getName());

    public static GenomeSketch sketch(Genome genome, int kSize, int sParam, int seed, boolean filterUniqueKMers, boolean saveKMers, ProgressListener progress) throws IOException {
        logger.fine("Calculating sketch for " + genome.getAccession());
        final GenomeSketch result = new GenomeSketch(genome);
        try (FastKMerIterator kmers = new FastKMerIterator(kSize, genome.getFastaUrl())) {
            result.sketch = FracMinHashSketch.compute(genome.getAccession(), kmers, genome.getGenomeSize(), true, sParam, seed, filterUniqueKMers, saveKMers, progress);
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
