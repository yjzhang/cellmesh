from collections import defaultdict
import json
import sqlite3
import pandas as pd
from scipy import sparse
from scipy.io import mmread

# load corpus
corpus = mmread('data/corpus.mm')
corpus = sparse.csr_matrix(corpus)
cells, genes = corpus.shape

# load cell types
with open('data/cell_info.json') as f:
    cell_info = json.load(f)

# load gene names
gene_names = pd.read_table('data/gene.dict.text', header=None)
gene_names['taxid'], gene_names['gene_id'], gene_names['gene_symbol'] = gene_names[1].str.split(',').str

# create sqlite3 db
conn = sqlite3.connect('cellmesh/data/cellmesh.db')
c = conn.cursor()

# create a table representing a cell type - gene mapping.
try:
    c.execute('CREATE TABLE cell_gene(cellID text, gene text, count integer)')
    # TODO: iterate over nonzero entries of corpus
    for i1, i2 in zip(*corpus.nonzero()):
        count = corpus[i1, i2]
        cell_record = cell_info[str(i1+1)]
        gene_record = gene_names.iloc[i2]
        cell = cell_record['id']
        gene = gene_record.gene_symbol
        c.execute('INSERT INTO cell_gene VALUES (?, ?, ?)', (cell, gene, count))
except Exception as e:
    print(str(e))
try:
    c.execute('CREATE INDEX cell_gene_id_index ON cell_gene(cellID)')
except:
    pass

# create a table representing a cellID - cell name mapping
try:
    c.execute('CREATE TABLE cell_name(cellID text, cellName text, cellIndex integer)')
    for i, cell_record in cell_info.items():
        cell_id = cell_record['id']
        cell_name = cell_record['name']
        index = int(i)
        c.execute('INSERT INTO cell_name VALUES (?, ?, ?)', (cell_id, cell_name, index))
except Exception as e:
    print(str(e))
try:
    c.execute('CREATE INDEX cell_name_id_index ON cell_name(cellID)')
except:
    pass

# create a table representing a gene symbol - gene info mapping
try:
    c.execute('CREATE TABLE gene_info(gene text, geneID integer, totalCounts integer)')
    for i in range(genes):
        gene_record = gene_names.iloc[i]
        gene = gene_record.gene_symbol
        gene_id = gene_record.gene_id
        total_counts = gene_record[2]
        c.execute('INSERT INTO gene_info VALUES (?, ?, ?)', (gene, gene_id, total_counts))
except Exception as e:
    print(str(e))
try:
    c.execute('CREATE INDEX gene_info_id_index ON gene_info(gene)')
except:
    pass

conn.commit()
conn.close()
