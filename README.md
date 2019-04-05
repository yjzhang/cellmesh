The data consists of two matrices: gene-mesh cell type and gene-pmid.

the gene-mesh matrix is corpus.mm.
cell_info.json is a dict of index to mesh cell type.
gene.dict.text contains info for each gene, where the indices are the same as corpus.mm. gene\_index(e.g. gene\_index 0 means the 0-th row of gene/cell matrix)   ncbi's taxid,geneid,gene\_symbol,  dfs(doc freq, number of cells the gene co-occurs with)
