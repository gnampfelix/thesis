package org.husonlab.ncbi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.openapitools.client.ApiClient;
import org.openapitools.client.ApiException;
import org.openapitools.client.api.GenomeApi;
import org.openapitools.client.model.V2AssemblyLinksReply;
import org.openapitools.client.model.V2AssemblyLinksReplyAssemblyLink;
import org.openapitools.client.model.V2AssemblyLinksReplyAssemblyLinkType;
import org.openapitools.client.model.V2reportsAssemblyDataReport;
import org.openapitools.client.model.V2reportsAssemblyDataReportPage;

public class NcbiApi {
    private ApiClient client;

    public NcbiApi(String basePath) {
        this.client = new ApiClient();
        this.client.setBasePath(basePath);
    }

    public List<Genome> getGenomes(List<String> accessionCodes) throws ApiException, AmbiguousDataException {
        List<Genome> result = new ArrayList<>();
        Map<String, String> downloadLinks = new HashMap<>();      
        GenomeApi genomes = new GenomeApi(this.client);
        String nextPageToken = null;
        
        // somehow, the links API does not require/provide pagination
        V2AssemblyLinksReply linksReply = genomes.genomeLinksByAccession(accessionCodes);
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
                    genomes.genomeDatasetReport(
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
                                downloadLinks.get(report.getAccession())
                            )
                        );
                    }
                }
            } while (nextPageToken != null);

        return result;
    }        
}
