import pdb
import numpy as np
import multiprocessing
from math import pow

from cellmesh import \
  DB_TFIDF_DIR,\
  get_all_cell_id_names, \
  get_all_genes, \
  get_cell_genes_pmids_count

def calc_cell_gene_ranks(cell_gene_count, all_genes):
  '''
  Input:
    cell_gene_count: list of (gene symbol, cnt)
    all_genes: set of all genes in DB
  Output:
    cell_gene_rank_sorted: list of (gene symbol, rank)
      including all genes

  TODO:
    The same cell has calc_cell_gene_ranks called repeatedly for different queries,
      this is not efficient. Needs to incorporate this in DB.
  '''
  cell_genes = set([x[0] for x in cell_gene_count])
  non_cell_genes = all_genes.difference(cell_genes)
  cell_gene_sorted = cell_gene_count + [(g, 0) for g in non_cell_genes]
  cell_gene_sorted = sorted(cell_gene_sorted, key=lambda x: -x[1])
  
  N = len(all_genes)
  p = int(N/2)
  cell_gene_rank_sorted = [(cell_gene_sorted[i][0], p-i) for i in range(N)]

  return cell_gene_rank_sorted

def calc_gsva_ext_one_query_one_cell(args):
  '''
  calc gsva_ext ES (enrichment score) for one query and one cell

  Input:
    genes_set: a set of query gene symbols
    cell_id: MeSH cell id of a candidate cell C (wrt a column in DB)
    cell_gene_count: list of (gene symbol, cnt wrt C)
    params: contains config params of gsva_ext_test. see gsva_ext_test() for description
    all_genes: set of all genes in DB
  Output:
    a tuple of cell id of MeSH cell C and ES wrt Q (Q for query)
  '''
  genes_set, cell_id, cell_gene_count, params, all_genes = args

  print("process cell %s, Kc=%d, k=%d"%(cell_id, len(cell_gene_count), len(overlapping_genes)))

  cell_gene_rank_sorted = calc_cell_gene_ranks(cell_gene_count, all_genes)

  walk_step_list = []
  for g, r in cell_gene_rank_sorted:
    if g in genes_set:
      walk_step = pow(np.abs(r), params.get("tau", 1))
    else:
      walk_step = -1
    walk_step_list.append(walk_step)

  denom1 = sum([w for w in walk_step_list if w>=0])
  denom2 = N - len(genes_set)

  if denom1<=0 or denom2<=0:
    return (cell_id, 0)

  ES = -np.inf 
  v = 0.0

  for i in range(N):
    inc_val = float(walk_step_list[i])/denom1 if walk_step_list[i]>=0 else (-1.0/denom2)
    v = v + inc_val
    ES = max(v, ES)  

  if ES<0: ES=0

  return (cell_id, ES) 

def gsva_ext_test_default_params():
  params = {}
  params["n_proc"] = 10
  params["db_cnt_thre"] = 0
  params["tau"] = 1
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
        "tau": used to calc walk step
      }
  Output:
    cell_ES_vals: list of 5-tuples: MeSH ID, cell name, ES (enrichment score) val, overlapping genes, pmids, 
      in descending order
  '''

  all_cells = get_all_cell_id_names(
    db_dir=DB_TFIDF_DIR, include_cell_components=include_cell_components, include_chromosomes=include_chromosomes)
  all_genes = get_all_genes(db_dir=DB_TFIDF_DIR)
  all_genes = set(all_genes)

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
      genes_set,
      cell_id,
      cell_gene_count,
      params,
      all_genes
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