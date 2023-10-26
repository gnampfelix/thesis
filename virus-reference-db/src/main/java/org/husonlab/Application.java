package org.husonlab;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.husonlab.ncbi.Genome;
import org.husonlab.ncbi.NcbiApi;
import org.husonlab.ncbi.TaxonomyTree;

public class Application {

    public static void main(String[] args) {
        NcbiApi api = new NcbiApi("https://api.ncbi.nlm.nih.gov/datasets/v2alpha"); 

        try {
            List<String> accessionCodes = new ArrayList<>(Arrays.asList("GCF_000819615.1", "GCF_000836805.1"));
            List<Genome> genomes = api.getGenomes(accessionCodes);
            TaxonomyTree tree = api.getTaxonomyTreeForGenomes(genomes);
            System.out.println(tree.getTree());
        } catch (Exception e) {
            System.out.println("well, f****");
            e.printStackTrace();
        }
    }
}
