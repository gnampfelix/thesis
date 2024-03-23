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
    private final int lastKmerIndex;

    private long[] hashValues;

    private String name;
    private int seed;

    private List<KMerCoordinates> coordinates;

    private FracMinHashSketch(int sParam, int kSize, String name, boolean isNucleotides, int seed) {
        this.sParam = sParam;
        this.kSize = kSize;
        this.lastKmerIndex = kSize-1;
        this.name = name;
        this.seed = seed;
        this.coordinates = new ArrayList<>();
    }

    /**
     * 
     * @param name
     * @param sequences
     * @param isNucleotides
     * @param sParam scaling factor, 0 <= 1/s <= 1; 1/s*H defines the threshold
     * @param kSize
     * @param seed
     * @param filterUniqueKMers
     * @param saveKMers
     * @param progress
     * @return
     */
    public static FracMinHashSketch compute(
        String name, 
        FastKMerIterator kmers,
        long genomeSize, 
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

    public long[] getValues() {
        return this.hashValues;
    }

    public int getKSize() {
        return this.kSize;
    }

    public int getSParam() {
        return this.sParam;
    }

    public int getSeed() {
        return this.seed;
    }

    public void fastReverseComplement(byte[]kmer, byte[]result) {
        for (int i = 0; i < this.kSize; i++) {
            result[i] = complementTable[kmer[this.lastKmerIndex - i]];
        }
    }

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

    public String getName() {
        return this.name;
    }

    public void setName(String name) {
        this.name = name;
    }

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

    public List<KMerCoordinates> getCoordinates() {
        return new ArrayList<>(this.coordinates);
    }
}
