package org.husonlab.ncbi;

import java.util.HashMap;
import java.util.Map;

import com.google.common.graph.Graph;
import com.google.common.graph.Graphs;

public class TaxonomyTree {
    private Graph<Taxon> tree;
    private Map<Integer, Taxon> taxa;
    
    public TaxonomyTree(Graph<Taxon> tree, Map<Integer,Taxon> taxa) {
        this.tree = tree;
        this.taxa = taxa;
    }
    
    public Graph<Taxon> getTree() {
        return Graphs.copyOf(this.tree);
    }

    public Map<Integer,Taxon> getTaxa() {
        return new HashMap<>(this.taxa);
    }

    public Taxon getRoot() {
        // This is the root node in the NCBI taxonomy
        return this.taxa.get(1);
    }
}
