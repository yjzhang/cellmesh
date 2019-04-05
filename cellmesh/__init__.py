import os
import sqlite3

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'cellmesh.db')

def get_all_genes():
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT DISTINCT gene FROM gene_info')
    results = C.fetchall()
    conn.close()
    return [x[0] for x in results]

def get_all_cell_id_names():
    """
    Returns a list of all unique cell ids + names
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT DISTINCT cellID, cellName  FROM cell_name')
    results = C.fetchall()
    conn.close()
    return results

@lru_cache(maxsize=None)
def get_cell_genes(cell):
    """
    Given a cell ID, this returns a list of all genes associated with that cell.
    """
    conn = sqlite3.connect(DB_DIR)
    C = conn.cursor()
    C.execute('SELECT gene, count FROM cell_gene WHERE cellID=?', (cell,))
    results = C.fetchall()
    results = [x[0] for x in results]
    conn.close()
    return results

def hypergeometric_test(genes, return_header=False):
    """
    Uses a hypergeometric test to identify the most relevant cell types.
    """
    from scipy import stats
    all_cells = get_all_cell_id_names()
    all_genes = get_all_genes()
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell_id, cell_name in all_cells:
        cell_genes = set(get_cell_genes(cell_id))
        overlapping_genes = list(genes.intersection(cell_genes))
        if len(overlapping_genes) == 0:
            continue
        #pmids = {}
        #for gene in overlapping_genes:
        #    pmids[gene] = get_papers_cell_gene(cell, gene)
        k = len(overlapping_genes)
        pv = stats.hypergeom.cdf(k - 1, len(all_genes), len(cell_genes), len(genes))
        cell_p_vals[cell_id] = (cell_name, 1 - pv, overlapping_genes)
    cell_p_vals = list(cell_p_vals.items())
    cell_p_vals.sort(key=lambda x: x[1][1])
    # merge items
    cell_p_vals = [(x[0],) + x[1] for x in cell_p_vals]
    if return_header:
        header = ['MeSH ID', 'Cell Name', 'P-value', 'Overlapping Genes']
        cell_p_vals = [header] + cell_p_vals
    return cell_p_vals

# TODO: what is a more sophisticated test that accounts for the same genes being present in many different cell types?
