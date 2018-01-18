#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# I just copied this little section... It's all around the internet, like a "HelloWorld.R" script.
if (length(args)<2) {
  stop("Initial \"Affected Gene list\" and organism Entrez Taxa ID must be provided.", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "/Users/siavash/Desktop/CTD_geneIDs.csv"
}

library(readr)
library("mygene")
library("biomaRt")



#write.table(df_out, file=args[2], row.names=FALSE)
#unique_ids <- read_csv("~/Desktop/unique_ids.txt",col_names = FALSE)
#View(unique_ids)

ids <- read_csv(args[1], col_names = TRUE)
# make sure unique_ids doesn't have duplicates
unique_ids <- unique(ids[ , c("Gene Symbol","Gene ID") ])
# write.csv(unique_ids,args[3] , row.names=FALSE )

## config before each run
taxid <- args[2]
#ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "celegans_gene_ensembl",  host = "mar2016.archive.ensembl.org" )
## complete this list. I know it's not the best practice to have this list entered manually...
organism_mart_list <- list(
    "6239" = "celegans_gene_ensembl" ,
    "10090" = "mmusculus_gene_ensembl" ,
    "7227" = "dmelanogaster_gene_ensembl" ,
    "7955" = "drerio_gene_ensembl"
    )

ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = organism_mart_list[[args[2]]],  host = "mar2016.archive.ensembl.org" )

## use dotfield=TRUE to get dot fields... and probably nicely formed results
result = queryMany(unique_ids$"Gene Symbol", scopes="symbol", fields=c("ensembl.gene, name"), species=taxid)
temp <- result[complete.cases(result[c("ensembl.gene","name")]),]


genes = getGenes((result$ensembl.gene) , fields = c("ensembl.protein","symbol","ensembl.gene") , species=taxid )
genes_temp = getGenes((temp$ensembl.gene) , fields = c("ensembl.protein","symbol","ensembl.gene") , species=taxid )


p = data.frame(Protein_id=character(), Protein_seq=character(),stringsAsFactors = FALSE)

CTD_geneIDs <- genes_temp[c("ensembl.protein","ensembl.gene","symbol","query")]
write.csv(CTD_geneIDs, args[3] )

# to count total protein IDs before pruning...
protein_id_list = NULL
for( i in 1:nrow(genes) ) {
    for(protein_id in genes$ensembl.protein[[i]]) {
        protein_id_list <- c(protein_id_list,(protein_id))
    }
} 

cat( length(protein_id_list) , " peptide IDs found.\n")


for( j in 1:nrow(genes)) {
    if(!is.null(genes$ensembl.protein[[j]])) {
        maxprotein = getSequence( id = genes$ensembl.protein[[j]][1] , type = "ensembl_peptide_id" , seqType = "peptide" , mart = ensembl)
        for( i in 1:length(genes$ensembl.protein[[j]])) {
            protein = getSequence( id = genes$ensembl.protein[[j]][i]  , type = "ensembl_peptide_id" , seqType = "peptide" , mart = ensembl)
            if( nrow(protein)!=0 ) {
                if( nrow(maxprotein)!=0) {
                    if( nchar(protein$peptide) >= nchar(maxprotein$peptide) ) {
                        maxprotein <- protein
                    }
                }
            }
        # print(genes$ensembl.protein[[1]][i])
        }
    }
    if(is.null(genes$ensembl.protein[[j]])) {
        if(!is.null(genes$ensembl[[j]])) {
            maxprotein = getSequence( id = genes$ensembl[[j]]$protein[1] , type = "ensembl_peptide_id" , seqType = "peptide" , mart = ensembl)
            ## for fucked up Mus Musculus gene/protein non-retrivals...
            ## buggy piece of shit here.
            if(nrow(maxprotein)==0) {
                maxprotein = getSequence( id = genes$ensembl[[j]]$protein[2] , type = "ensembl_peptide_id" , seqType = "peptide" , mart = ensembl)
            }
            for( i in 1:length(genes$ensembl[[j]]$protein)) {
                protein = getSequence( id = genes$ensembl[[j]]$protein[i]  , type = "ensembl_peptide_id" , seqType = "peptide" , mart = ensembl)
                if( nrow(protein)!=0) {
                    if( nchar(protein$peptide) >= nchar(maxprotein$peptide) ) {
                        maxprotein <- protein
                    }
                }
            # print(genes$ensembl.protein[[1]][i])
            }
        }
    }
    if(nrow(maxprotein)!=0) {
        if(nrow(maxprotein)==1) {
                seq = maxprotein$peptide    
                p[nrow(p)+1,] = c(maxprotein$ensembl_peptide_id,seq)
            }
        if(nrow(maxprotein)>=2) {
            maxseq = maxprotein$peptide[1]
            maxid = maxprotein$ensembl_peptide_id[1]
            for( k in 1:nrow(maxprotein)) {
                seq = maxprotein$peptide[k]
                if(nchar(seq) >= nchar(maxseq) ) {
                    maxseq = seq
                    maxid = maxprotein$ensembl_peptide_id[k]
                }
            }
            p[nrow(p)+1,] = c(maxid,maxseq)
            }
        }
        else{
            count <- count+1
        }
        if((j %% 100 )==0) {
            cat((j%/%100)*100 , "protein sequences retrieved...\n")
        }
}

cat( nrow(p) , " protein sequences retrieved in total.\n")

write.csv(p,args[3], row.names=FALSE )

