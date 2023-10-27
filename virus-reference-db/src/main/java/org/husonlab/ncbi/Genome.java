package org.husonlab.ncbi;

import java.nio.file.Paths;

public class Genome {
    private String organismName;
    private String accession;
    private int taxonId;
    private String assemblyName;
    private String downloadLink;
    private int genomeSize;

    public Genome(String organismName, String accession, int taxonId, String assemblyName, String downloadLink, int genomeSize) {
        this.organismName = organismName;
        this.accession = accession;
        this.taxonId = taxonId;
        this.assemblyName = assemblyName;
        this.downloadLink = downloadLink;
        this.genomeSize = genomeSize;
    }

    public String getOrganismName() {
        return this.organismName;
    }

    public void setOrganismName(String organismName) {
        this.organismName = organismName;
    }

    public String getAccession() {
        return this.accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public int getTaxonId() {
        return this.taxonId;
    }

    public void setTaxonId(int taxonId) {
        this.taxonId = taxonId;
    }

    public String getAssemblyName() {
        return this.assemblyName;
    }

    public void setAssemblyName(String assemblyName) {
        this.assemblyName = assemblyName;
    }

    public String getDownloadLink() {
        return this.downloadLink;
    }

    public void setDownloadLink(String downloadLink) {
        this.downloadLink = downloadLink;
    }

    public String getFastaUrl() {
        return this.downloadLink + "/" + this.accession + "_" + this.assemblyName + "_genomic.fna.gz";
    }

    public void setGenomeSize(int genomeSize) {
        this.genomeSize = genomeSize;
    }

    public int getGenomeSize() {
        return this.genomeSize;
    }
    
}
