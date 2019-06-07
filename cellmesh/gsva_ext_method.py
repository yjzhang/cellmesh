import pdb
import numpy as np
import multiprocessing

from cellmesh import \
  DB_TFIDF_DIR,\
  get_all_cell_id_names, \
  get_all_genes, \
  get_cell_genes_pmids_count

def calc_gsva_ext_one_query_one_cell(args):
  '''
  calc gsva_ext ES (enrichment score) for one query and one cell

  Input:
    genes: a ranked list of gene symbols, as query Q
    cell_id: MeSH cell id of a candidate cell C (wrt a column in DB)
    cell_gene_count: list of (gene symbol, cnt wrt C)
    overlapping_genes: set of gene symbols occurring in Q and C
      currently not used
    params: contains config params of gsva_ext_test. see gsva_ext_test() for description
    N_all_genes: total number of genes in DB
  Output:
    a tuple of cell id of MeSH cell C and ES wrt Q (Q for query)
  '''
  genes, cell_id, cell_gene_count, overlapping_genes, params, N_all_genes = args

  print("process cell %s, Kc=%d, k=%d"%(cell_id, len(cell_gene_count), len(overlapping_genes)))

  return (cell_id, 0) 

def gsva_ext_test_default_params():
  params = {}
  params["n_proc"] = 10
  params["db_cnt_thre"] = 0
  return params

def gsva_ext_test(genes, return_header=False, include_cell_components=False, include_chromosomes=False, params=None):
  '''
  This is the GSVA ext (aka modified) query test on the tf-idf matrix
    the function is modified so that it has similar I/O structure as normed_hypergeometric_test

  Input:
    genes: a ranked list of gene symbols
    return_header:
    include_cell_components:
    include_chromosomes:
    params: None or contains config params of gsva_ext_test
      {
        "n_proc": number of processes for parallel processing,
        "db_cnt_thre": a gene g is considered to co-occur with a cell c if db(g, c) > db_cnt_thre
      }
  Output:
    cell_ES_vals: list of 5-tuples: MeSH ID, cell name, ES (enrichment score) val, overlapping genes, pmids, 
      in descending order
  '''

  all_cells = get_all_cell_id_names(
    db_dir=DB_TFIDF_DIR, include_cell_components=include_cell_components, include_chromosomes=include_chromosomes)
  all_genes = get_all_genes(db_dir=DB_TFIDF_DIR)
  N_all_genes = len(all_genes)

  genes_set = set(genes)

  #----- algo configuration
  if params is None:
    params = gsva_ext_test_default_params()

  #----- prepare info for per (cell_id, cell_name, ES=n/a, overlapping_genes, pmids)
  cell_ES_vals = {}
  args_list = []
  for cell_id, cell_name in all_cells:
    #a set of (gene_symbol, its pmids connected by ",", cnt val) wrt candidate cell_id
    genes_pmids_count = set(get_cell_genes_pmids_count(cell_id, db_dir=DB_TFIDF_DIR, threshold=params["db_cnt_thre"]))
    cell_genes = [x[0] for x in genes_pmids_count]
    overlapping_genes = genes_set.intersection(cell_genes)
    if len(overlapping_genes) == 0:
      continue

    pmids = {}
    for gene, pmid, _ in genes_pmids_count:
      if gene in overlapping_genes:
        pmids[gene] = pmid.split(',')

    overlapping_genes = list(overlapping_genes)
    cell_ES_vals[cell_id] = [cell_name, -np.inf, overlapping_genes, pmids]

    cell_gene_count = [(x[0], x[2]) for x in genes_pmids_count]
    args = (
      genes,
      cell_id,
      cell_gene_count,
      overlapping_genes,
      params,
      N_all_genes
      )
    args_list.append(args)

  #----- calc ES scores, in parallel if possible
  n_proc = params["n_proc"]
  res = []
  if n_proc==1:
    for args in args_list:
      r = calc_gsva_ext_one_query_one_cell(args)
      res.append(r)
  else:
    p = multiprocessing.Pool(n_proc)
    res = p.map(calc_gsva_ext_one_query_one_cell, args_list)
  for cell_id, ES in res:
    cell_ES_vals[cell_id][1] = ES
    cell_ES_vals[cell_id] = tuple(cell_ES_vals[cell_id])

  #----- wrap up outputs -----

  cell_ES_vals = list(cell_ES_vals.items())
  cell_ES_vals.sort(key=lambda x: -x[1][1]) #descending order
  # merge items
  cell_ES_vals = [(x[0],) + x[1] for x in cell_ES_vals]
  if return_header is True:
    header = ['MeSH ID', 'Cell Name', 'Prob-value', 'Overlapping Genes', 'PMIDs']
    cell_ES_vals = [header] + cell_ES_vals

  return cell_ES_vals

def test_gsva_ext_test():
  #tabula-muris-dropseq, top ~20 genes (taxid,geneid,symbol) for B cell
  #based on scLit, the top retrieval should be:
  #  D001402, B-Lymphocytes, -89.49026841217798
  genes = [
    '9606,973,CD79A',
    '9606,931,MS4A1',
    '9606,974,CD79B',
    '9606,115650,TNFRSF13C',
    '9606,55024,BANK1',
    '9606,1380,CR2',
    '9606,930,CD19',
    '9606,951,CD37',
    '9606,933,CD22',
    '9606,115350,FCRL1',
    '9606,84824,FCRLA',
    '9606,972,CD74',
    '9606,4050,LTB',
    '9606,640,BLK',
    '9606,5450,POU2AF1'
  ]
  genes = [g.split(",")[2] for g in genes]
  params = gsva_ext_test_default_params()

  cell_ES_vals = gsva_ext_test(genes=genes, params=params)
  
  for i in range(min(len(cell_ES_vals), 10)):
    t = cell_ES_vals[i]
    print("i=%d, id=%s, name=%s, ES=%f"%(i, t[0], t[1], t[2]))
  return

if __name__ == '__main__':
  test_gsva_ext_test()