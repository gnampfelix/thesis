package org.husonlab.sketch;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import jloda.kmers.bloomfilter.BloomFilter;
import jloda.seq.SequenceUtils;
import jloda.thirdparty.MurmurHash;
import jloda.util.CanceledException;
import jloda.util.StringUtils;
import jloda.util.progress.ProgressListener;

public class FracMinHashSketch {
    public static final int MAGIC_INT = 1213415758; // for starters, I've just increased the number

    private final double sParam;
    private final int kSize;
    private final String name;
    private final boolean isNucleotides;

    private long[] hashValues;
    private byte[][] kmers;

    public FracMinHashSketch(double sParam, int kSize, String name, boolean isNucleotides) {
        this.sParam = sParam;
        this.kSize = kSize;
        this.name = name;
        this.isNucleotides = isNucleotides;
    }

    /**
     * 
     * @param name
     * @param sequences
     * @param isNucleotides
     * @param sParam scaling factor, 0 <= s <= 1; s*H defines the threshold
     * @param kSize
     * @param seed
     * @param filterUniqueKMers
     * @param saveKMers
     * @param progress
     * @return
     */
    public static FracMinHashSketch compute(
        String name, 
        Collection<byte[]> sequences, 
        boolean isNucleotides, 
        double sParam, 
        int kSize, 
        int seed, 
        boolean filterUniqueKMers, 
        boolean saveKMers, 
        ProgressListener progress
    ) {
        final FracMinHashSketch sketch = new FracMinHashSketch(sParam, kSize, name, isNucleotides);
        final TreeSet<Long> sortedSet = new TreeSet<>();

        final Map<Long, byte[]> hash2kmer = saveKMers ? new HashMap<>() : null;

        // s = H / s', 0 <= s <= H
        final double threshold = Long.MAX_VALUE * sParam;

        final BloomFilter bloomFilter;
        if (filterUniqueKMers) {
            bloomFilter = new BloomFilter(sequences.stream().mapToInt(s -> s.length).sum(), 500000000);
        } else {
            bloomFilter = null;
        }

        try {
            final byte[] kMer = new byte[kSize];
            final byte[] kMerReverseComplement = new byte[kSize];

            for (byte[] sequence : sequences) {
                final int top = sequence.length - kSize;
                for (int offset = 0; offset < top; offset++) {
                    // check if kMer contains ambiguous nucleotide "N"
                    if (isNucleotides) {
                        final int ambiguousPos = StringUtils.lastIndexOf(sequence, offset, kSize, 'N');
                        if (ambiguousPos != -1) {
                            offset = ambiguousPos;
                            continue;
                        }
                    }

                    SequenceUtils.getSegment(sequence, offset, kSize, kMer);
                    final byte[] kMerUse;
                    if (isNucleotides) {
                        SequenceUtils.getReverseComplement(sequence, offset, kSize, kMerReverseComplement);

                        if (SequenceUtils.compare(kMer, kMerReverseComplement) <= 0) {
                            kMerUse = kMer;
                        } else {
                            kMerUse = kMerReverseComplement;
                        }
                    } else {
                        kMerUse = kMer;
                    }

                    if (bloomFilter != null && bloomFilter.add(kMerUse)) {
                        continue; // we have just seen this k-mer for the first time
                    }

                    final long hash = MurmurHash.hash64(kMerUse, 0, kSize, seed);

                    if (hash < threshold) {
                        sortedSet.add(hash);
                        if (hash2kmer != null) {
                            hash2kmer.put(hash, kMerUse.clone());
                        }
                    }

                    progress.checkForCancel();
                }
            }
            sketch.hashValues = new long[sortedSet.size()];
            int pos = 0;
            for (Long value : sortedSet) {
                sketch.hashValues[pos++] = value;
            }
            progress.incrementProgress();
        } catch (CanceledException ignored) {}

        if (saveKMers && hash2kmer != null) {
            sketch.kmers = new byte[hash2kmer.size()][];
            int i = 0;
            for (byte[] kmer : hash2kmer.values()) {
                sketch.kmers[i++] = kmer;
            }
        }
        progress.reportTaskCompleted();
        return sketch;        
    }
}
