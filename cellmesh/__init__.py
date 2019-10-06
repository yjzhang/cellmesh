import os
import sqlite3

try:
    from functools import lru_cache
except ImportError:
    from backports.functools_lru_cache import lru_cache

PATH = os.path.dirname(__file__)
DB_DIR = os.path.join(PATH, 'data', 'cellmesh.db')
DB_TFIDF_DIR = os.path.join(PATH, 'data', 'cellmesh_tfidf.db')
ANATOMY_DB_DIR = os.path.join(PATH, 'data', 'anatomy_mesh.db')
CELLMESH_CELLMARKER_MAP = os.path.join(PATH, 'data', 'cell_ontology_mesh_mapping.tsv')
ROOT_MESH_ID_NAMES_DIR = os.path.join(PATH, 'data', 'root_mesh_id_names.txt')
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
def get_descendants(cell_types, db_dir=ANATOMY_DB_DIR):
    """
    cell_types has to be a tuple.
    Returns a list of all descendants of the given cell types, where cell_types is a list of MeSH IDs (D0...)
    Returns a list of cell IDs
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    descendants = []
    for cell in cell_types:
        C.execute('SELECT children FROM cell_children WHERE cellID=?', (cell,))
        results = C.fetchall()
        r = results[0][0].split(',')
        descendants += r
    return descendants

@lru_cache(maxsize=None)
def get_all_cell_id_names(db_dir=DB_DIR, include_cell_components=True, include_chromosomes=False, include_cell_lines=False,
        cell_type_subset=None):
    """
    Returns a list of all unique cell ids + names

    cell_type_subset is a tuple of mesh IDs.
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
    if not include_cell_lines:
        with open(os.path.join(PATH, 'data', 'cell_line_ids.txt')) as f:
            cell_lines = set(x.strip() for x in f.readlines())
            results = [x for x in results if x[0] not in cell_lines]
    if cell_type_subset is not None:
        try:
            descendants = get_descendants(tuple(cell_type_subset))
            results = [x for x in results if x[0] not in descendants]
        except:
            results = [x for x in results if x[0] not in cell_type_subset]
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
def get_cell_genes_pmids_count(cell, threshold=3, db_dir=DB_DIR):
    """
    Given a cell ID, this returns a list of all genes associated with that cell.
    The threshold is the minimum count for the gene to be included.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT gene, pmids, count FROM cell_gene WHERE cellID=?', (cell,))
    results = C.fetchall()
    results = [x for x in results if x[2] > threshold]
    conn.close()
    return results

@lru_cache(maxsize=None)
def get_cells_threshold(threshold=3, db_dir=DB_DIR, include_cell_components=True, include_chromosomes=False, include_cell_lines=False):
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
    if not include_cell_lines:
        with open(os.path.join(PATH, 'data', 'cell_line_ids.txt')) as f:
            cell_lines = set(x.strip() for x in f.readlines())
            results = [x for x in results if x[0] not in cell_lines]
    return [x for x in results if x[2] > threshold]

def hypergeometric_test(genes, return_header=False, include_cell_components=False, include_chromosomes=False,
        include_cell_lines=False, cell_type_subset=None, db_dir=DB_DIR):
    # TODO: cell type subset - only include a subset of cell types
    """
    Uses a hypergeometric test to identify the most relevant cell types.

    Args:
        genes (list): upper-case gene names
        return_header (bool): if True, returns a tuple as the first element.
        include_cell_components, include_chromosomes, include_cell_lines: booleans that should probably always be False.
        cell_type_subset (list or None): list of cell types to include, along with their descendants..

    Returns:
        list of 5-tuples: MeSH ID, cell name, p-value, overlapping genes, pmids, in order
        of ascending p-value.
    """
    from scipy import stats
    genes = [x.upper() for x in genes]
    if isinstance(cell_type_subset, list):
        cell_type_subset = tuple(cell_type_subset)
    all_cells = get_all_cell_id_names(include_cell_components=include_cell_components,
            include_chromosomes=include_chromosomes, include_cell_lines=include_cell_lines,
            db_dir=db_dir, cell_type_subset=cell_type_subset)
    all_genes = get_all_genes(db_dir=db_dir)
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell_id, cell_name in all_cells:
        genes_pmids = set(get_cell_genes_pmids(cell_id, db_dir=db_dir))
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

def normed_hypergeometric_test(genes, return_header=False, include_cell_components=False, include_chromosomes=False, include_cell_lines=False):
    """
    This hypergeometric test is on the tf-idf matrix, with a calibrated threshold.

    Returns:
        list of 5-tuples: MeSH ID, cell name, p-value, overlapping genes, pmids, in order
        of ascending p-value.
    """
    from scipy import stats
    all_cells = get_all_cell_id_names(db_dir=DB_TFIDF_DIR, include_cell_components=include_cell_components, include_chromosomes=include_chromosomes, include_cell_lines=include_cell_lines)
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


def cellmesh_to_cellonto(name_or_id, is_name=True):
    """
    Maps a cellmesh name to cellmarker id + name pairs

    Args:
        name_or_id (str or list): string or list of names
    """
    import pandas as pd
    data = pd.read_table(CELLMESH_CELLMARKER_MAP)
    outputs = []
    if isinstance(name_or_id, str):
        return cellmesh_to_cellonto_single(name_or_id, data, is_name)
    else:
        for name in name_or_id:
            outputs.append(cellmesh_to_cellonto_single(name, data, is_name))
    return outputs

def cellmesh_to_cellonto_single(name_or_id, data=None, is_name=True):
    """
    Returns: a list of (cell ontology id, cell ontology name) pairs
    """
    import pandas as pd
    if data is None:
        data = pd.read_table(CELLMESH_CELLMARKER_MAP)
    if is_name:
        data_subset = data[data['MeSH Name(s)'] == name_or_id]
    else:
        data_subset = data[data['MeSH UID'] == name_or_id]
    if len(data_subset) == 0:
        return []
    else:
        return [(value[1], value[2]) for value in data_subset[['Cell Ontology ID', 'Cell Ontology Name']].itertuples()]


def cellonto_to_cellmesh(name_or_id, is_name=True):
    """
    Maps a cellmarker name to cellmesh id + name pairs

    Args:
        name_or_id (str or list): string or list of names
    """
    import pandas as pd
    data = pd.read_table(CELLMESH_CELLMARKER_MAP)
    outputs = []
    if isinstance(name_or_id, str):
        return cellonto_to_cellmesh_single(name_or_id, data, is_name)
    else:
        for name in name_or_id:
            outputs.append(cellonto_to_cellmesh_single(name, data, is_name))
    return outputs


def cellonto_to_cellmesh_single(name_or_id, data=None, is_name=True):
    """
    Returns: a list of (cellmesh_id, cellmesh_name) pairs
    """
    import pandas as pd
    if data is None:
        data = pd.read_table(CELLMESH_CELLMARKER_MAP)
    if is_name:
        data_subset = data[data['Cell Ontology Name'] == name_or_id]
    else:
        data_subset = data[data['Cell Ontology ID'] == name_or_id]
    if len(data_subset) == 0:
        return []
    else:
        return [(value[1], value[2]) for value in data_subset[['MeSH UID', 'MeSH Name(s)']].itertuples()]
