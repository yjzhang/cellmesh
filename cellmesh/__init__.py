import os
import sqlite3

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'cellmesh.db')
DB_TFIDF_DIR = os.path.join(PATH, 'data', 'cellmesh_tfidf.db')
# number of cell types that pass the threshold
N_CELLS_THRESHOLD = 372
TFIDF_THRESHOLD = 0.034154
N_CELLS_THRESHOLD_TFIDF = 534

def get_all_genes(db_dir=DB_DIR):
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT DISTINCT gene FROM gene_info')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

@lru_cache(maxsize=None)
def get_all_cell_id_names(db_dir=DB_DIR, include_cell_components=True, include_chromosomes=False):
    """
    Returns a list of all unique cell ids + names
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT DISTINCT cellID, cellName FROM cell_name')
    results = C.fetchall()
    if not include_cell_components:
        with open(os.path.join(PATH, 'data', 'cell_component_ids.txt')) as f:
            cell_components = set(x.strip() for x in f.readlines())
            results = [x for x in results if x[0] not in cell_components]
    if not include_chromosomes:
        with open(os.path.join(PATH, 'data', 'chromosome_ids.txt')) as f:
            chromosomes = set(x.strip() for x in f.readlines())
            results = [x for x in results if x[0] not in chromosomes]
    conn.close()
    return results

@lru_cache(maxsize=None)
def get_cell_genes_pmids(cell, threshold=3, db_dir=DB_DIR):
    """
    Given a cell ID, this returns a list of all genes associated with that cell.
    The threshold is the minimum count for the gene to be included.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT gene, pmids, count FROM cell_gene WHERE cellID=?', (cell,))
    results = C.fetchall()
    results = [x[:2] for x in results if x[2] > threshold]
    conn.close()
    return results

@lru_cache(maxsize=None)
def get_cells_threshold(threshold=3, db_dir=DB_DIR, include_cell_components=True, include_chromosomes=False):
    """
    Returns a list of all cell types with their max citation count, as a tuple of (cellID, cellName, count).
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT DISTINCT cell_name.cellID, cellName, MAX(count) FROM cell_name INNER JOIN cell_gene ON cell_name.cellID = cell_gene.cellID GROUP BY cell_name.cellID;')
    results = C.fetchall()
    if not include_cell_components:
        with open(os.path.join(PATH, 'data', 'cell_component_ids.txt')) as f:
            cell_components = set(x.strip() for x in f.readlines())
            results = [x for x in results if x[0] not in cell_components]
    if not include_chromosomes:
        with open(os.path.join(PATH, 'data', 'chromosome_ids.txt')) as f:
            chromosomes = set(x.strip() for x in f.readlines())
            results = [x for x in results if x[0] not in chromosomes]
    return [x for x in results if x[2] > threshold]

def hypergeometric_test(genes, return_header=False, include_cell_components=False, include_chromosomes=False):
    """
    Uses a hypergeometric test to identify the most relevant cell types.

    Returns:
        list of 5-tuples: MeSH ID, cell name, p-value, overlapping genes, pmids, in order
        of ascending p-value.
    """
    from scipy import stats
    genes = [x.upper() for x in genes]
    all_cells = get_all_cell_id_names(include_cell_components=include_cell_components, include_chromosomes=include_chromosomes)
    all_genes = get_all_genes()
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell_id, cell_name in all_cells:
        genes_pmids = set(get_cell_genes_pmids(cell_id))
        cell_genes = [x[0] for x in genes_pmids]
        overlapping_genes = genes.intersection(cell_genes)
        if len(overlapping_genes) == 0:
            continue
        pmids = {}
        for gene, pmid in genes_pmids:
            if gene in overlapping_genes:
                pmids[gene] = pmid.split(',')
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(cell_genes), len(genes))
        overlapping_genes = list(overlapping_genes)
        cell_p_vals[cell_id] = (cell_name, 1 - pv, overlapping_genes, pmids)
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1][1])
    # merge items
    cell_p_vals = [(x[0],) + x[1] for x in cell_p_vals]
    if return_header:
        header = ['MeSH ID', 'Cell Name', 'P-value', 'Overlapping Genes', 'PMIDs']
        cell_p_vals = [header] + cell_p_vals
    return cell_p_vals

def normed_hypergeometric_test(genes, return_header=False, include_cell_components=False, include_chromosomes=False):
    """
    This hypergeometric test is on the tf-idf matrix, with a calibrated threshold.

    Returns:
        list of 5-tuples: MeSH ID, cell name, p-value, overlapping genes, pmids, in order
        of ascending p-value.
    """
    from scipy import stats
    all_cells = get_all_cell_id_names(db_dir=DB_TFIDF_DIR, include_cell_components=include_cell_components, include_chromosomes=include_chromosomes)
    all_genes = get_all_genes(db_dir=DB_TFIDF_DIR)
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell_id, cell_name in all_cells:
        genes_pmids = set(get_cell_genes_pmids(cell_id, db_dir=DB_TFIDF_DIR, threshold=0.034154))
        cell_genes = [x[0] for x in genes_pmids]
        overlapping_genes = genes.intersection(cell_genes)
        if len(overlapping_genes) == 0:
            continue
        pmids = {}
        for gene, pmid in genes_pmids:
            if gene in overlapping_genes:
                pmids[gene] = pmid.split(',')
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(cell_genes), len(genes))
        overlapping_genes = list(overlapping_genes)
        cell_p_vals[cell_id] = (cell_name, 1 - pv, overlapping_genes, pmids)
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1][1])
    # merge items
    cell_p_vals = [(x[0],) + x[1] for x in cell_p_vals]
    if return_header:
        header = ['MeSH ID', 'Cell Name', 'P-value', 'Overlapping Genes', 'PMIDs']
        cell_p_vals = [header] + cell_p_vals
    return cell_p_vals


# TODO: what is a more sophisticated test that accounts for the same genes being present in many different cell types?
