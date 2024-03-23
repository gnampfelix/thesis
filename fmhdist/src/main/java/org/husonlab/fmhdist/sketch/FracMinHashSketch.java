package org.husonlab.fmhdist.sketch;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

import org.husonlab.fmhdist.util.FastKMerIterator;
import org.husonlab.fmhdist.util.KMerCoordinates;

import jloda.seq.SequenceUtils;
import jloda.thirdparty.MurmurHash;
import jloda.util.ByteInputBuffer;
import jloda.util.ByteOutputBuffer;


/**
 * Class to store the FracMinHash sketch as defined by Irber et al.. The size of
 * the sketch depends on the scaling parameter s and the underlying sequences.
 * To get an insight into the composition of the sketch, users may wish to
 * obtain the set of KMerCoordinates with getCoordinates() which describes the
 * sketch in more detail (what KMers are included, where are they coming from?)
 */
public class FracMinHashSketch {
    public static final int MAGIC_INT = 1213415759; // for starters, I've just increased the number

    private static final byte[] complementTable = new byte[128];
    static {
        complementTable['A'] = 'T';
        complementTable['T'] = 'A';
        complementTable['G'] = 'C';
        complementTable['C'] = 'G';
        complementTable['a'] = 't';
        complementTable['t'] = 'a';
        complementTable['g'] = 'c';
        complementTable['c'] = 'g';
    }

    private final int sParam;
    private final int kSize;

    private long[] hashValues;

    private String name;
    private int seed;

    private List<KMerCoordinates> coordinates;

    private FracMinHashSketch(int sParam, int kSize, String name, boolean isNucleotides, int seed) {
        this.sParam = sParam;
        this.kSize = kSize;
        this.name = name;
        this.seed = seed;
        this.coordinates = new ArrayList<>();
    }

    /**
     * Creates a new FracMinHash sketch of the given k-mers according to Irber
     * et al.
     *
     * For this, the method calculates the canonical k-mer by selecting the
     * lexicographically smaller of the k-mer and its reverse complement. This
     * canonical k-mer is then hashed using the MurMur64 hash function and the
     * given random seed. If the sequence is not a nucleotide sequence, the
     * original k-mer is automatically the canonical k-mer.
     *
     * If the hash is below the threshold that is given by H/s where H is the
     * maximum output of the HashFunction, then it is part of the sketch.
     * @param name A name to describe the sketch, this has no impact on the
     * algorithm
     * @param kmers An instance of the FastKMerIterator that provides all k-mers
     * that should be considered for the sketch
     * @param isNucleotides flag to indicate if the sequence consists of
     * nucleotides. If so, a canonical k-mer is chosen by selecting the
     * lexicographically smaller of the k-mer and its reverse complement.
     * @param sParam The scaling param s of the algorithm
     * @param seed A random seed that should be used for hashing
     * @return A new FracMinHashSketch
     */
    public static FracMinHashSketch compute(
        String name, 
        FastKMerIterator kmers,
        boolean isNucleotides, 
        int sParam, 
        int seed
    ) {
        final FracMinHashSketch sketch = new FracMinHashSketch(sParam, kmers.getK(), name, isNucleotides, seed);
        final TreeSet<Long> sortedSet = new TreeSet<>();

        // Irber et al define the hash function as h: o -> [0, H]. However, in
        // the case of our Java Long hashes, the range is h: o -> [-H, H-1].
        // Thus, we need to shift the threshold accordingly.
        final double fraction = Long.MAX_VALUE * (1/(double)sParam) * 2; //the complete range is 2H, thus a fraction is 2Hs
        final double threshold = Long.MIN_VALUE + fraction;
        
        final byte[] kMer = new byte[kmers.getK()];
        final byte[] kMerReverseComplement = new byte[kmers.getK()];
        while (kmers.hasNext()) {
            byte[] next = kmers.next();
            System.arraycopy(next, 0, kMer, 0, sketch.kSize);
            final byte[] kMerUse;
            if (isNucleotides) {
                System.arraycopy(kmers.getReverseComplement(), 0, kMerReverseComplement, 0, kmers.getK());
                if (SequenceUtils.compare(kMer, kMerReverseComplement) <= 0) {
                    kMerUse = kMer;
                } else {
                    kMerUse = kMerReverseComplement;
                }
            } else {
                kMerUse = kMer;
            }

            final long hash = MurmurHash.hash64(kMerUse, 0, sketch.kSize, seed);

            if (hash < threshold) {
                KMerCoordinates coords = kmers.getCoordinates();
                coords.setHash(hash);
                sortedSet.add(hash);
                sketch.coordinates.add(coords);
            }
        }
        sketch.hashValues = new long[sortedSet.size()];
        int pos = 0;
        for (Long value : sortedSet) {
            sketch.hashValues[pos++] = value;
        }
        
        return sketch;        
    }

    /**
     * Returns all hash values that are part of the sketch.
     * @return
     */
    public long[] getValues() {
        return this.hashValues;
    }

    /**
     * Returns the k-mers size k that the FastKmerIterator used for the k-mer
     * decomposition.
     * @return
     */
    public int getKSize() {
        return this.kSize;
    }

    /**
     * Returns the scaling param s of the algorithm.
     * @return
     */
    public int getSParam() {
        return this.sParam;
    }

    /**
     * Returns the random seed that was used for hashing.
     * @return
     */
    public int getSeed() {
        return this.seed;
    }

    /**
     * Converts the sketch into a byte sequence. All values are encoded in
     * little endian in the following order:
     * 
     * 1. Magic Number (4B)
     * 2. S Param (4B)
     * 3. K Size (4B)
     * 4. Random seed (4B)
     * 5. Sketch size (4B)
     * 6. Hash values (sketch size x 8B)
     * @return
     */
    public byte[] getBytes() {
        ByteOutputBuffer bytes = new ByteOutputBuffer();
        bytes.writeIntLittleEndian(MAGIC_INT);
        bytes.writeIntLittleEndian(this.sParam);
        bytes.writeIntLittleEndian(this.kSize);
        bytes.writeIntLittleEndian(this.seed);
        bytes.writeIntLittleEndian(this.hashValues.length);
        for (int i = 0; i < this.hashValues.length; i++) {
            bytes.writeLongLittleEndian(this.hashValues[i]);
        }
        return bytes.copyBytes();
    }

    /**
     * Returns the descriptive name of the sketch that was given during
     * computation.
     * @return
     */
    public String getName() {
        return this.name;
    }

    /**
     * Sets a descriptive name for the sketch.
     * @param name
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * Parses a byte sequence as a FracMinHashSketch. The byte sequence is
     * described in the getBytes() documentation.
     * @param bytes
     * @return
     * @throws IOException
     */
    public static FracMinHashSketch parse(byte[] bytes) throws IOException {
        final ByteInputBuffer buffer = new ByteInputBuffer(bytes);

        if (buffer.readIntLittleEndian() != MAGIC_INT)
            throw new IOException("Incorrect magic number");
        int sParam = buffer.readIntLittleEndian();
        int kMerSize = buffer.readIntLittleEndian();
        int seed = buffer.readIntLittleEndian();
        int sketchSize = buffer.readIntLittleEndian();

        final FracMinHashSketch sketch = new FracMinHashSketch(sParam, kMerSize, "", true, seed);
        sketch.hashValues = new long[sketchSize];
        for (int i = 0; i < sketchSize; i++) {
            sketch.hashValues[i] = buffer.readLongLittleEndian();
        }
        return sketch;
    }

    /**
     * Returns the coordinates of a freshly computed sketch. This is not
     * applicable if the sketch was reconstructed using the parse() method, i.e.
     * the returned list will be empty.
     * @return
     */
    public List<KMerCoordinates> getCoordinates() {
        return new ArrayList<>(this.coordinates);
    }
}
