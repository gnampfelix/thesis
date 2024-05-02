# Summary

In this directory, I am trying to analyze the claim that FracMinHash works
better for genomes with different sizes than Mash. For this, I've setup a
reference database consisting of all currently (as of April 2024) published
_Phytophthora_ reference genomes.

As query sequences, I am using different _Phytophthora cinnamomi_ genomes as
well as genomes of species that are typically found in soil samples of Avocado
trees (which is a target of _P. cinnamomi_) according to [Solís-García _et
al._](https://pubmed.ncbi.nlm.nih.gov/33510714/).

I've sketched all sequences using $s=2000$, $k=21$ and FarmHash as the hash
function and estimated the phylogenetic distances. For comparison, I have
calculated distances using MashTree and the default parameters.

While the overall placement is similar, MashTree seems to produce more jaccard
similarities of 0, especially when dealing with the bacerial and fungal queries.
This leads to a placement in the final outline that hides an interesting detail:
_Phytophthora_ and the query fungi have a non empty intersection.

I am trying to find a ground truth to compare this against as ANI stops working
beyond genus level comparisons. AAI seems promising, but the results are not yet
in.