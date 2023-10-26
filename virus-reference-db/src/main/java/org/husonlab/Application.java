package org.husonlab;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.husonlab.ncbi.Genome;
import org.husonlab.ncbi.NcbiApi;

public class Application {

    public static void main(String[] args) {
        NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha");        
        try {
            List<String> accessionCodes = new ArrayList<>(Arrays.asList("GCF_000819615.1", "GCF_000836805.1"));
            List<Genome> genomes = api.getGenomes(accessionCodes);
            for (Genome genome : genomes) {
                System.out.println(genome.getOrganismName());
                System.out.println(genome.getDownloadLink());
            }
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
