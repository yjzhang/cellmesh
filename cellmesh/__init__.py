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
COUNT_THRESHOLD = 3
N_CELLS_THRESHOLD_TFIDF = 534
SPECIES_MAP = {'mus_musculus': 10090, 'homo_sapiens': 9606, 'human': 9606, 'mouse': 10090, 'worm': 6239, 'c_elegans': 6239}
SPECIES_TFIDF_THRESHOLDS = {'mus_musculus': 0.03856, 'mouse': 0.03856, 'homo_sapiens': 0.03342, 'human': 0.03342, 'worm': 0, 'c_elegans': 0}
SPECIES_TFIDF_THRESHOLDS_ANATOMY = {'mus_musculus': 0.04891, 'mouse': 0.04891, 'homo_sapiens': 0.04911, 'human': 0.04911, 'worm': 0, 'c_elegans': 0}
BOTH_TAXIDS = [9606, 10090]

@lru_cache(maxsize=None)
def get_all_genes(db_dir=DB_DIR, species='human', uppercase_names=True):
    """
    Returns a list of all unique gene symbols.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    if species == 'both':
        results = []
        for taxid in BOTH_TAXIDS:
            C.execute('SELECT DISTINCT gene FROM gene_info WHERE taxid=?', (taxid,))
            results += C.fetchall()
        conn.close()
        if uppercase_names:
            results = list(set(x[0].upper() for x in results))
        else:
            results = list(set(x[0] for x in results))
        return results
    else:
        taxid = SPECIES_MAP[species]
        C.execute('SELECT DISTINCT gene FROM gene_info WHERE taxid=?', (taxid,))
        results = C.fetchall()
        conn.close()
        if uppercase_names:
            return [x[0].upper() for x in results]
        else:
            return [x[0] for x in results]

@lru_cache(maxsize=None)
def get_immediate_descendants(cell_type, db_dir=ANATOMY_DB_DIR):
    """
    Returns a list of cell IDs representing direct descendents of the given cell type.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    C.execute('SELECT directChildren FROM cell_children WHERE cellID=?', (cell_type,))
    results = C.fetchall()
    r = results[0][0].split(',')
    if len(r) == 1 and r[0] == '':
        return []
    return r

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
    Returns a list of all unique cell ids + names - tuple (cell_id, cell_name)

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
def get_cell_genes_pmids(cell, threshold=3, db_dir=DB_DIR, species='homo_sapiens', return_count=False, use_tfidf=False,
        uppercase_gene_names=True):
    """
    Given a cell ID, this returns a list of all genes associated with that cell.
    The threshold is the minimum count for the gene to be included.
    """
    conn = sqlite3.connect(db_dir)
    C = conn.cursor()
    # TODO: deal with 'both'
    statement = 'SELECT gene, pmids, count FROM cell_gene WHERE cellID=? AND taxid=?'
    if use_tfidf:
        statement = 'SELECT gene, pmids, tfidf FROM cell_gene WHERE cellID=? AND taxid=?'
        if threshold == 3:
            if db_dir == DB_DIR:
                threshold == SPECIES_TFIDF_THRESHOLDS[species]
            elif db_dir == ANATOMY_DB_DIR:
                threshold == SPECIES_TFIDF_THRESHOLDS_ANATOMY[species]
    if species == 'both':
        # merge results by gene
        gene_results = {}
        for taxid in BOTH_TAXIDS:
            C.execute(statement, (cell, taxid))
            results = C.fetchall()
            for gene, pmids, count in results:
                gene = gene.upper()
                if gene in gene_results:
                    _, old_pmids, old_counts = gene_results[gene]
                    gene_results[gene] = (gene, old_pmids + ',' + pmids, count + old_counts)
                else:
                    gene_results[gene] = (gene, pmids, count)
        results = [x for x in gene_results.values() if x[2] > threshold]
    else:
        taxid = SPECIES_MAP[species]
        C.execute(statement, (cell, taxid))
        results = C.fetchall()
        results = [x for x in results if x[2] > threshold]
    conn.close()
    if uppercase_gene_names:
        results = [(x[0].upper(),) + x[1:] for x in results]
    if not return_count:
        results = [x[:2] for x in results]
    return results


def hypergeometric_test(genes, return_header=False, include_cell_components=False, include_chromosomes=False,
        include_cell_lines=False, cell_type_subset=None, db_dir=DB_DIR, additional_table=None, species='homo_sapiens',
        use_tfidf=False):
    """
    Uses a hypergeometric test to identify the most relevant cell types.

    Args:
        genes (list): upper-case gene names
        return_header (bool): if True, returns a tuple as the first element.
        include_cell_components, include_chromosomes, include_cell_lines: booleans that should probably always be False.
        cell_type_subset (list or None): list of cell types to include, along with their descendants..
        db_dir (str): either DB_DIR, ANATOMY_DIR
        additional_table (DataFrame or None): table with cell types as columns and genes as names.

    Returns:
        list of 5-tuples: MeSH ID, cell name, p-value, overlapping genes, pmids, in order
        of ascending p-value.
    """
    from scipy import stats
    # both query genes and database genes are converted to upper case.
    genes = [x.upper() for x in genes]
    if isinstance(cell_type_subset, list):
        cell_type_subset = tuple(cell_type_subset)
    all_cells = get_all_cell_id_names(include_cell_components=include_cell_components,
            include_chromosomes=include_chromosomes, include_cell_lines=include_cell_lines,
            db_dir=db_dir, cell_type_subset=cell_type_subset)
    all_genes = get_all_genes(db_dir=db_dir, species=species, uppercase_names=True)
    cell_p_vals = {}
    genes = set(genes)
    # each cell should have 4 items: cell type, p-value, overlapping genes, PMIDs
    for cell_id, cell_name in all_cells:
        genes_pmids = set(get_cell_genes_pmids(cell_id, db_dir=db_dir, species=species,
            use_tfidf=use_tfidf, return_count=False, uppercase_gene_names=True))
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

def normed_hypergeometric_test(genes, **kwargs):
    """
    This hypergeometric test is on the tf-idf matrix, with a calibrated threshold.

    Returns:
        list of 5-tuples: MeSH ID, cell name, p-value, overlapping genes, pmids, in order
        of ascending p-value.
    """
    return hypergeometric_test(genes, use_tfidf=True, **kwargs)


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

def get_cellmesh_anatomy_root_terms():
    """
    Returns a list of pairs (cellmesh_term, anatomy)
    """
    terms = []
    with open(ROOT_MESH_ID_NAMES_DIR) as f:
        for line in f.readlines():
            data = line.split()
            mesh_id = data[0]
            mesh_name = ' '.join(data[1:])
            terms.append((mesh_id, mesh_name))
    return terms

def get_cellmesh_anatomy_tree():
    """
    Returns a tree-structure: a dict of dicts...
    """
    root_terms = get_cellmesh_anatomy_root_terms()
    all_term_ids = get_all_cell_id_names(db_dir=ANATOMY_DB_DIR, include_cell_components=True, include_chromosomes=True, include_cell_lines=True)
    id_to_name = {x[0]: x[1] for x in all_term_ids}
    tree = {}
    to_visit = [(mesh_id, name, tree) for mesh_id, name in root_terms]
    while to_visit:
        mesh_id, name, parent = to_visit.pop()
        descendants = get_immediate_descendants(mesh_id)
        if descendants:
            sub_tree = {}
            parent[mesh_id] = sub_tree
            for descendant in descendants:
                to_visit.append((descendant, id_to_name[descendant], sub_tree))
    return tree, id_to_name
