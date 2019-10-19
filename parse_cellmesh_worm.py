from collections import Counter
import json
import sqlite3
import os

import traceback

path = 'data_c_elegans'

# 1. open db
conn_db = sqlite3.connect(os.path.join(path, 'corpus.db'))
cursor_db = conn_db.cursor()

# 1.5. load gene and mesh info from corpus.db
cursor_db.execute('SELECT * FROM gene_dict')
gene_list = cursor_db.fetchall()
gene_dict = {x[2]: x for x in gene_list}
cursor_db.execute('SELECT * FROM cell_dict')
cell_list = cursor_db.fetchall()
cell_dict = {x[1]: x for x in cell_list}
# 2. open json file
with open(os.path.join(path, 'CoCount_gene2pubmed_(C_elegans_protein_coding)_MeSH_cell.json')) as f:
    gene_mesh_pubmed = json.load(f)

# dict of gene id : total count
total_counts_dict = Counter()

# 3. create new db
output_db = sqlite3.connect('cellmesh/data/cellmesh.db')
output_cursor = output_db.cursor()

# create a table representing a MeSH ID - gene mapping.
try:
    # pmids is a comma-separated string of ints
    #output_cursor.execute('CREATE TABLE cell_gene(cellID text, gene text, count integer, pmids text, taxid text)')
    # iterate over nonzero entries of corpus
    for key, item in gene_mesh_pubmed.items():
        s1 = key.split(',')
        taxid = s1[0]
        s2 = s1[1].split(':')
        gene_id, mesh_id = s2

        count = item['cnt']

        cell_record = cell_dict[mesh_id]
        gene_record = gene_dict[gene_id]
        cell_id = cell_record[1]

        gene = gene_record[3]
        total_counts_dict[gene] += count
        taxid = gene_record[1]
        pmids = ','.join(item['place'])
        output_cursor.execute('INSERT INTO cell_gene VALUES (?, ?, ?, ?, ?)', (cell_id, gene, count, pmids, taxid))
except Exception as e:
    text = traceback.format_exc()
    print(text)

# update the gene symbol - gene info mapping
try:
    existing_genes = output_cursor.execute('SELECT geneID FROM gene_info')
    existing_genes = set([x[0] for x in existing_genes])
    for i, gene_record in gene_dict.items():
        gene = gene_record[3]
        taxid = gene_record[1]
        gene_id = gene_record[2]
        if gene_id not in existing_genes:
            total_counts = total_counts_dict[i]
            output_cursor.execute('INSERT INTO gene_info VALUES (?, ?, ?, ?)', (gene, gene_id, total_counts, taxid))
except Exception as e:
    text = traceback.format_exc()
    print(text)

output_cursor.execute('VACUUM')

output_db.commit()
output_db.close()

##################################################################################3

# add data to cellmesh-anatomy db

"""
# 1. open db
conn_db = sqlite3.connect(os.path.join(path, 'corpus.db'))
cursor_db = conn_db.cursor()

# 1.5. load gene and mesh info from corpus.db
cursor_db.execute('SELECT * FROM gene_dict')
gene_list = cursor_db.fetchall()
gene_dict = {x[2]: x for x in gene_list}
cursor_db.execute('SELECT * FROM cell_dict')
cell_list = cursor_db.fetchall()
cell_dict = {x[1]: x for x in cell_list}
# 2. open json file

# dict of gene id : total count
total_counts_dict = Counter()


with open(os.path.join(path, 'CoCount_gene2pubmed_(C_elegans_protein_coding)_MeSH_cell.json')) as f:
    gene_mesh_pubmed = json.load(f)

output_db = sqlite3.connect('cellmesh/data/anatomy_mesh.db')
output_cursor = output_db.cursor()

# create a table representing a MeSH ID - gene mapping.
try:
    # pmids is a comma-separated string of ints
    #output_cursor.execute('CREATE TABLE cell_gene(cellID text, gene text, count integer, pmids text, taxid text)')
    # iterate over nonzero entries of corpus
    for key, item in gene_mesh_pubmed.items():
        s1 = key.split(',')
        taxid = s1[0]
        s2 = s1[1].split(':')
        gene_id, mesh_id = s2

        count = item['cnt']

        cell_record = cell_dict[mesh_id]
        gene_record = gene_dict[gene_id]
        cell_id = cell_record[1]

        gene = gene_record[3]
        total_counts_dict[gene] += count
        taxid = gene_record[1]
        pmids = ','.join(item['place'])
        output_cursor.execute('INSERT INTO cell_gene VALUES (?, ?, ?, ?, ?)', (cell_id, gene, count, pmids, taxid))
except Exception as e:
    text = traceback.format_exc()
    print(text)

# update the gene symbol - gene info mapping
try:
    existing_genes = output_cursor.execute('SELECT geneID FROM gene_info')
    existing_genes = set([x[0] for x in existing_genes])
    for i, gene_record in gene_dict.items():
        gene = gene_record[3]
        taxid = gene_record[1]
        gene_id = gene_record[2]
        if gene_id not in existing_genes:
            total_counts = total_counts_dict[i]
            output_cursor.execute('INSERT INTO gene_info VALUES (?, ?, ?, ?)', (gene, gene_id, total_counts, taxid))
except Exception as e:
    text = traceback.format_exc()
    print(text)

output_db.commit()
output_db.close()
"""


