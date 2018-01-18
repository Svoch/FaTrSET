# FaTrSET

Abstract

The current eco-toxicological approach to identify sites of biological concern uses measures like median lethal dose and lethal concentration analysis which do not link the effects of toxins to any other level of biological organization. These traditional methods, would not inform of the presence of hazardous chemicals not before a significant environmental damage has already been dealt. Being motivated by the idea that essential pathways and genes are conserved across the taxa, we present a framework that incorporates transcriptomic data of inhabitant species to monitor the environment.
This framework mines chemical-gene interaction data from publicly available databases like the Comparative Toxicology Database, retrieving an initial set of genes in the input organisms that are known to be affected by target chemicals. Next, it retrieves translated protein sequences related to the obtained genes and performs homology-based searches in the genome of target organisms to predict a more thorough list of affected sequences. Afterwards, it clusters these protein sequences to obtain groups of orthologous sequences. Moreover, for genes in the ortholog groups, this pipeline retrieves supplementary data, like related pathways and Gene Ontology terms. Lastly, for a user-defined number of pathways with the most genes present, the pipeline learns a Bayesian Network to demonstrate meaningful correlations between shared pathways within target organisms and input chemicals.
To demonstrate this pipeline, we ran the pipeline on a test data consisting of a group of phylogenetically distant taxa, namely Mus musculus, Danio rerio, Drosophila melanogaster and Caenorhabditis elegans, with input chemical groups of heavy-metals and dioxins. A total number of 8700 and 34102 affected genes were retrieved for heavy-metals and dioxins respectively. Based on retrieved protein sequences for these genes, they were clustered into 1866 and 7740 homologous gene groups in 101 and 125 pathways. The shared highly affected pathways suggested by our pipeline correspond with findings of previous studies on toxic effects of these chemicals.

