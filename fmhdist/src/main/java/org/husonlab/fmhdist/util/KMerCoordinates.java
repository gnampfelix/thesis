package org.husonlab.fmhdist.util;

public class KMerCoordinates {
    private int recordIndexInFile;
    private int sequenceIndexInRecord;
    private int sequenceIndexInRecordIncludingAmbiguous;
    private int sequenceIndexInFile;
    private int sequenceIndexInFileIncludingAmbiguous;
    private byte[] kmer;

    /**
     * Creates a new instance of a kmer.
     * @param recordIndexInFile - the index of the record in the FASTA file from
     * which the k-mer is
     * @param sequenceIndexInFile - the index of the k-mer in the file (assuming
     * the first k-mer is 0), excluding all ambiguous k-mers when counting
     * @param sequenceIndexInRecord - the index of the k-mer in the record
     * (assuming the first k-mer is 0), excluding all ambiguous k-mers when
     * counting
     * @param sequenceIndexInFileIncludingAmbiguous - the index of the k-mer in
     * the file (assuming the first k-mer is 0), including  all ambiguous k-mers
     * when counting
     * @param sequenceIndexInRecordIncludingAmbiguous - the index of the k-mer in
     * the record (assuming the first k-mer is 0), including  all ambiguous k-mers
     * when counting
     * @param kmer A reference to the k-mer content. Will be clonedm on initialization.
     */
    public KMerCoordinates(
        int recordIndexInFile, 
        int sequenceIndexInFile, 
        int sequenceIndexInRecord, 
        int sequenceIndexInFileIncludingAmbiguous, 
        int sequenceIndexInRecordIncludingAmbiguous, 
        byte[] kmer
    ) {
        this.recordIndexInFile = recordIndexInFile;
        this.sequenceIndexInFile = sequenceIndexInFile;
        this.sequenceIndexInRecord = sequenceIndexInRecord;
        this.sequenceIndexInFileIncludingAmbiguous = sequenceIndexInFileIncludingAmbiguous;
        this.sequenceIndexInRecordIncludingAmbiguous = sequenceIndexInRecordIncludingAmbiguous;
        this.kmer = kmer.clone();
    } 

    public int getRecordIndexInFile() {
        return this.recordIndexInFile;
    }

    public int getSequenceIndexInRecord() {
        return this.sequenceIndexInRecord;
    }

    public int getSequenceIndexInRecordIncludingAmbiguous() {
        return this.sequenceIndexInRecordIncludingAmbiguous;
    }

    public int getSequenceIndexInFile() {
        return this.sequenceIndexInFile;
    }

    public int getSequenceIndexInFileIncludingAmbiguous() {
        return this.sequenceIndexInFileIncludingAmbiguous;
    }

    public byte[] getKmer() {
        return this.kmer;
    }

    public String toString() {
        return String.format(
            "%d,%d,%d,%d,%d,%s",
            this.recordIndexInFile, 
            this.sequenceIndexInFile,
            this.sequenceIndexInRecord,
            this.sequenceIndexInFileIncludingAmbiguous,
            this.sequenceIndexInRecordIncludingAmbiguous,
            new String(this.kmer));
    }

    public static KMerCoordinates fromString(String input) {
        String[] parts = input.split(",");
        int recordIndexInFile = Integer.parseInt(parts[0]);
        int sequenceIndexInFile = Integer.parseInt(parts[1]);
        int sequenceIndexInRecord = Integer.parseInt(parts[2]);
        int sequenceIndexInFileIncludingAmbiguous = Integer.parseInt(parts[3]);
        int sequenceIndexInRecordIncludingAmbiguous = Integer.parseInt(parts[4]);
        String kmer = parts[5].strip();
        return new KMerCoordinates(
            recordIndexInFile, 
            sequenceIndexInFile, 
            sequenceIndexInRecord, 
            sequenceIndexInFileIncludingAmbiguous, 
            sequenceIndexInRecordIncludingAmbiguous, 
            kmer.getBytes()
        );
    }
}
