package org.husonlab.fmhdist.sketch;

import java.io.IOException;
import java.util.logging.Logger;

import org.husonlab.fmhdist.ncbi.Genome;
import org.husonlab.fmhdist.util.KMerIterator;
import org.husonlab.fmhdist.util.LineKMerIteratorWithCoordinates;

import net.openhft.hashing.LongHashFunction;

public class GenomeSketch {
    private Genome genome;
    private FracMinHashSketch sketch;
    private static Logger logger = Logger.getLogger(GenomeSketch.class.getName());

    public static GenomeSketch sketch(Genome genome, int kSize, int sParam, LongHashFunction hashFunction, int seed, boolean prepareCoordinates) throws IOException {
        logger.fine("Calculating sketch for " + genome.getAccession());
        final GenomeSketch result = new GenomeSketch(genome);
        
        KMerIterator kmers;
        kmers = new LineKMerIteratorWithCoordinates(kSize, genome.getFastaUrl(), true);
        result.sketch = FracMinHashSketch.compute(genome.getAccession(), kmers, sParam, hashFunction, seed, prepareCoordinates);
        kmers.close();
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
