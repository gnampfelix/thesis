package org.husonlab.fmhdist.util;

public class KMerCoordinates {
    private int recordIndexInFile;
    private int sequenceIndexInRecord;
    private int sequenceIndexInRecordIncludingAmbiguous;
    private int sequenceIndexInFile;
    private int sequenceIndexInFileIncludingAmbiguous;
    private byte[] kmer;

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
            "(%d, %d, %d)/(%d, %d, %d): %s", 
            this.recordIndexInFile, 
            this.sequenceIndexInFile,
            this.sequenceIndexInRecord,
            this.recordIndexInFile,
            this.sequenceIndexInFileIncludingAmbiguous,
            this.sequenceIndexInRecordIncludingAmbiguous,
            new String(this.kmer));
    }
}
