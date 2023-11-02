package org.husonlab.ncbi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.logging.Logger;

import org.openapitools.client.ApiClient;
import org.openapitools.client.ApiException;
import org.openapitools.client.api.GenomeApi;
import org.openapitools.client.api.TaxonomyApi;
import org.openapitools.client.model.V2AssemblyLinksReply;
import org.openapitools.client.model.V2AssemblyLinksReplyAssemblyLink;
import org.openapitools.client.model.V2AssemblyLinksReplyAssemblyLinkType;
import org.openapitools.client.model.V2TaxonomyMatch;
import org.openapitools.client.model.V2TaxonomyMetadataRequestContentType;
import org.openapitools.client.model.V2TaxonomyMetadataResponse;
import org.openapitools.client.model.V2reportsAssemblyDataReport;
import org.openapitools.client.model.V2reportsAssemblyDataReportPage;

import com.google.common.graph.GraphBuilder;
import com.google.common.graph.MutableGraph;

public class NcbiApi {
    private ApiClient client;
    private final int TAXON_PAGE_SIZE = 100;
    private final Logger logger;
    private final GenomeApi genomes;
    private final TaxonomyApi taxonomy;

    public NcbiApi(String basePath) {
        this.client = new ApiClient();
        this.client.setBasePath(basePath);
        this.logger = Logger.getLogger(NcbiApi.class.getName());
        this.genomes = new GenomeApi(this.client);
        this.taxonomy = new TaxonomyApi(this.client);
    }

    public List<Genome> getGenomes(List<String> accessionCodes) throws ApiException, AmbiguousDataException {
        List<Genome> result = new ArrayList<>();
        Map<String, String> downloadLinks = new HashMap<>();      
        String nextPageToken = null;
        
        // somehow, the links API does not require/provide pagination
        V2AssemblyLinksReply linksReply = this.genomes.genomeLinksByAccession(accessionCodes);
        List<V2AssemblyLinksReplyAssemblyLink> links = linksReply.getAssemblyLinks();
        if (links != null) {
            for (V2AssemblyLinksReplyAssemblyLink link : links) {
                if(link.getAssemblyLinkType().equals(V2AssemblyLinksReplyAssemblyLinkType.FTP_LINK)) {
                    if (downloadLinks.containsKey(link.getAccession())) {
                        throw new AmbiguousDataException("multiple valid download links for accession code found");
                    }
                    downloadLinks.put(link.getAccession(), link.getResourceLink());
                }
            }
        }
        
        do {
            V2reportsAssemblyDataReportPage reportPage = 
                    this.genomes.genomeDatasetReport(
                        accessionCodes, 
                        null, 
                        null, 
                        null,
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        null, 
                        nextPageToken, 
                        null, 
                        null, 
                        null
                    );
                nextPageToken = reportPage.getNextPageToken();
                List<V2reportsAssemblyDataReport> reports = reportPage.getReports();
                if (reports != null) {
                    for(V2reportsAssemblyDataReport report : reports) {
                        
                        result.add(
                            new Genome(
                                report.getOrganism().getOrganismName(), 
                                report.getAccession(), 
                                report.getOrganism().getTaxId(), 
                                report.getAssemblyInfo().getAssemblyName(),
                                downloadLinks.get(report.getAccession()),
                                Integer.parseInt(report.getAssemblyStats().getTotalSequenceLength())
                            )
                        );
                    }
                }
            } while (nextPageToken != null);

        return result;
    }

    private void fetchLeafBatch(
        MutableGraph<Taxon> tree, 
        Map<Integer, Taxon> taxa, 
        Queue<List<String>> lineageQueryQueue, 
        List<String> taxonBatch, 
        Map<Integer, Integer> genomeParents
    ) throws ApiException {
        V2TaxonomyMetadataResponse response = this.taxonomy.taxonomyMetadata(taxonBatch, V2TaxonomyMetadataRequestContentType.COMPLETE);
        List<V2TaxonomyMatch> nodes = response.getTaxonomyNodes();
        if (nodes != null) {
            for(V2TaxonomyMatch node : nodes) {
                List<String> lineage = new ArrayList<>();
                int directParent = 0;
                for (int parent : node.getTaxonomy().getLineage()) {
                    lineage.add(String.valueOf(parent));
                    directParent = parent;
                }
                lineageQueryQueue.add(lineage);
                Taxon taxon = new Taxon(node.getTaxonomy().getOrganismName(), node.getTaxonomy().getTaxId());
                taxa.put(taxon.getTaxonId(), taxon);
                tree.addNode(taxon);
                genomeParents.put(taxon.getTaxonId(), directParent);
            }                    
        }
    }
    
    public TaxonomyTree getTaxonomyTreeForGenomes(List<Genome> genomes) throws ApiException {
        logger.info("Fetching taxonomy tree...");
        

        Map<Integer, Taxon> taxa = new HashMap<>();
        Map<Integer, Integer> genomeParents = new HashMap<>();
        MutableGraph<Taxon> tree = GraphBuilder.directed().build();
        // The taxonomy API does not support pagination as of now but the
        // requests could be come very large if we query the complete lineage of
        // all taxa. So, maybe query the lineage for one taxon.
        Queue<List<String>> lineageQueryQueue = new LinkedList<>();
        
        logger.info("Processing leaves...");
        // First, get the lineage of all given genomes and create taxon for genome
        List<String> taxonBatch = new ArrayList<>();
        for (Genome g : genomes) {
            taxonBatch.add(String.valueOf(g.getTaxonId()));
            if (taxonBatch.size() >= TAXON_PAGE_SIZE) {
                logger.info("Leaf batch for " + taxonBatch.toString());
                fetchLeafBatch(tree, taxa, lineageQueryQueue, taxonBatch, genomeParents);
                taxonBatch.clear();
            }
        }
        // Process final batch
        logger.info("Leaf batch for " + taxonBatch.toString());
        fetchLeafBatch(tree, taxa, lineageQueryQueue, taxonBatch, genomeParents);
        
        
        // Now, resolve all lineages, all elements are findable by construction
        logger.info("Processing lineages...");
        while (!lineageQueryQueue.isEmpty()) {
            List<String> lineageQuery = lineageQueryQueue.remove();
            logger.info("Fetching lineage " +lineageQuery.toString());

            V2TaxonomyMetadataResponse response = taxonomy.taxonomyMetadata(lineageQuery, V2TaxonomyMetadataRequestContentType.COMPLETE);
            
            // First, create all nodes (their order is not the same as in the query!)
            List<V2TaxonomyMatch> nodes = response.getTaxonomyNodes();
            if (nodes != null) {
                for(V2TaxonomyMatch node : nodes) {
                    logger.info("processing taxon " + node.getTaxonomy().getTaxId() + "...");
                    Taxon current;
                    // We might have seen this taxon before (shared lineage)
                    if (!taxa.containsKey(node.getTaxonomy().getTaxId())) {
                        logger.info("taxon not known, creating new...");
                        current = new Taxon(node.getTaxonomy().getOrganismName(), node.getTaxonomy().getTaxId());
                        taxa.put(current.getTaxonId(), current);
                        tree.addNode(current);
                    }
                }

                // Then, add relationships as indicated by the lineage order
                logger.info("checking edges...");
                Taxon prev = null;
                for(String lineageItem : lineageQuery) {
                    Taxon current = taxa.get(Integer.parseInt(lineageItem));
                    if (prev != null && !tree.hasEdgeConnecting(prev, current)) {
                        logger.info("inserting edge (" + prev.getTaxonId() + "," + current.getTaxonId() + ")...");
                        tree.putEdge(prev, current);
                    }
                    prev = current;
                }
            }
        }

        // Finally, insert the leaves
        logger.info("Adding genome linkes...");
        for (Entry<Integer, Integer> genomeLink : genomeParents.entrySet()) {
            if (!tree.hasEdgeConnecting(taxa.get(genomeLink.getValue()), taxa.get(genomeLink.getKey()))) {
                logger.info("inserting edge (" + genomeLink.getValue() + "," + genomeLink.getKey() + ")...");
                tree.putEdge(taxa.get(genomeLink.getValue()), taxa.get(genomeLink.getKey()));
            }
        }
        logger.info("Fetched taxonomy tree!");
        return new TaxonomyTree(tree, taxa);
    }
}
