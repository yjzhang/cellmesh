# this is heavily based on https://github.com/tanghaibao/goatools/blob/master/notebooks/goea_nbt3102.ipynb
import os

# uses goatools to query given a gene set
from goatools.base import download_go_basic_obo

from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.test_data.genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

PATH = os.path.dirname(__file__)
go_basic_path = os.path.join(PATH, 'data', 'go-basic.obo')
gene2go_path = os.path.join(PATH, 'data', 'gene2go')

try:
    # Get http://geneontology.org/ontology/go-basic.obo
    go_basic_path = download_go_basic_obo(go_basic_path)
    # Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    gene2go_path = download_ncbi_associations(gene2go_path)
except:
    # if directory is not writeable for whatever reason, just save to /tmp
    go_basic_path = os.path.join('/tmp', 'go-basic.obo')
    gene2go_path = os.path.join('/tmp', 'gene2go')
    go_basic_path = download_go_basic_obo(go_basic_path)
    gene2go_path = download_ncbi_associations(gene2go_path)

obodag = GODag(go_basic_path)
# Read NCBI's gene2go. Store annotations in a list of namedtuples
objanno = Gene2GoReader(gene2go_path, taxids=[10090])
ns2assoc = objanno.get_ns2assc()
symbols_to_ids = {val.Symbol : key for key, val in GeneID2nt_mus.items()}
ids_to_symbols = {val : key for key, val in symbols_to_ids.items()}

def gene_set_query(genes, fdr_threshold=0.10, return_header=False):
    """
    Runs a GO enrichment analysis query using goatools.
    The GO dataset here is for mouse, but it might apply to human as well.
    """
    goeaobj = GOEnrichmentStudyNS(
            GeneID2nt_mus.keys(), # List of mouse protein-coding genes
            ns2assoc, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            alpha=fdr_threshold, # default significance cut-off
            methods=['fdr_bh']) # defult multipletest correction method
    genes = [x.capitalize() for x in genes]
    gene_ids = [symbols_to_ids[x] for x in genes if x in symbols_to_ids]
    print('gene_ids:', gene_ids)

    results = goeaobj.run_study(gene_ids)
    results_sig = [r for r in results if r.p_fdr_bh < fdr_threshold]
    results_table = []
    for r in results_sig:
        results_table.append([r.goterm.id, r.goterm.name, r.p_fdr_bh, [ids_to_symbols[gene_id] for gene_id in r.study_items]])
    print(results_table)
    results_table.sort(key=lambda x: x[2])
    if return_header:
        results_table = [['GO ID', 'Name', 'FDR', 'Overlapping Genes']] + results_table
    print('GO results_table:', results_table)
    return results_table

