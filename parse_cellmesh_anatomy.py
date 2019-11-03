from collections import Counter
import json
import sqlite3
import os

import traceback

path = 'cnt_AllAnatomy_T0_RowNormNone_ColNormNone'

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

# TODO: construct hierarchies
tree_ids_dict = {}
# children_dict includes both direct and indirect children.
children_dict = {}
direct_children_dict = {}
parents_dict = {}
direct_parents_dict = {}
root_items = []
for x in cell_list:
    cell_id = x[1]
    tree_id = x[3]
    tree_ids_dict[tree_id] = cell_id
    if cell_id not in children_dict:
        children_dict[cell_id] = set()
        direct_children_dict[cell_id] = set()
        parents_dict[cell_id] = set()
        direct_parents_dict[cell_id] = set()
for tree_id, cell_id in tree_ids_dict.items():
    tree_sub_components = tree_id.split('.')
    if len(tree_sub_components) == 1:
        root_items.append(cell_id)
    for i in range(1, len(tree_sub_components) + 1):
        sub_tree_id = '.'.join(tree_sub_components[:i])
        if sub_tree_id in tree_ids_dict:
            cell_id_parent = tree_ids_dict[sub_tree_id]
            children_dict[cell_id_parent].add(cell_id)
            parents_dict[cell_id].add(cell_id_parent)
            if i == len(tree_sub_components) - 1:
                direct_children_dict[cell_id_parent].add(cell_id)
                direct_parents_dict[cell_id].add(cell_id_parent)
import numpy as np
root_item_names = []
for cell_id in root_items:
    cell_record = cell_dict[cell_id]
    name = cell_record[2]
    root_item_names.append(cell_id + '\t' + name)
np.savetxt('cellmesh/data/root_mesh_id_names.txt', root_item_names, fmt='%s')

# 2. open json file
with open(os.path.join(path, 'CoCount_gene2pubmed_(homo_sapiens_protein_coding)_MeSH_anatomy.json')) as f:
    gene_mesh_pubmed = json.load(f)

# calculate tf-idf
# create count matrix of terms x genes
import numpy as np
data_dict = np.zeros((len(cell_list), len(gene_list)))
total_counts_dict = Counter()
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
    data_dict[cell_dict[cell_id][0], gene_dict[gene_id][0]] += count
from sklearn.feature_extraction.text import TfidfTransformer
tfidf = TfidfTransformer()
data_tfidf = tfidf.fit_transform(data_dict).toarray()

# 3. create new db
output_db = sqlite3.connect('cellmesh/data/anatomy_mesh.db')
output_cursor = output_db.cursor()

# create a table representing a MeSH ID - gene mapping.
try:
    # pmids is a comma-separated string of ints
    output_cursor.execute('CREATE TABLE cell_gene(cellID text, gene text, count integer, tfidf real, pmids text, taxid text)')
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
        taxid = gene_record[1]
        pmids = ','.join(item['place'])
        tfidf_val = float(data_tfidf[cell_dict[cell_id][0], gene_dict[gene_id][0]])
        output_cursor.execute('INSERT INTO cell_gene VALUES (?, ?, ?, ?, ?, ?)', (cell_id, gene, count, tfidf_val, pmids, taxid))
except Exception as e:
    text = traceback.format_exc()
    print(text)

try:
    output_cursor.execute('CREATE INDEX cell_gene_id_index ON cell_gene(cellID, gene, taxid)')
except:
    pass


# create a table representing a cellID - cellName mapping
try:
    output_cursor.execute('CREATE TABLE cell_name(cellID text, cellName text)')
    for i, cell_record in cell_dict.items():
        cell_id = i
        cell_name = cell_record[2]
        output_cursor.execute('INSERT INTO cell_name VALUES (?, ?)', (cell_id, cell_name))
except Exception as e:
    text = traceback.format_exc()
    print(text)
try:
    output_cursor.execute('CREATE INDEX cell_name_id_index ON cell_name(cellID)')
except:
    pass

# create a table representing a gene symbol - gene info mapping
try:
    output_cursor.execute('CREATE TABLE gene_info(gene text, geneID integer, totalCounts integer, taxid text)')
    for i, gene_record in gene_dict.items():
        gene = gene_record[3]
        gene_id = gene_record[2]
        taxid = gene_record[1]
        total_counts = total_counts_dict[i]
        output_cursor.execute('INSERT INTO gene_info VALUES (?, ?, ?, ?)', (gene, gene_id, total_counts, taxid))
except Exception as e:
    text = traceback.format_exc()
    print(text)
try:
    output_cursor.execute('CREATE INDEX gene_info_id_index ON gene_info(gene, taxid)')
except:
    pass

# create a table representing a cell id - children mapping
try:
    output_cursor.execute('CREATE TABLE cell_children(cellID text, children text, directChildren text)')
    for cell_id, children in children_dict.items():
        children_text = ','.join(children)
        direct_children_text = ','.join(direct_children_dict[cell_id])
        output_cursor.execute('INSERT INTO cell_children VALUES (?, ?, ?)', (cell_id, children_text, direct_children_text))
except Exception as e:
    text = traceback.format_exc()
    print(text)
try:
    output_cursor.execute('CREATE INDEX cell_children_id_index ON cell_children(cellID)')
except:
    pass


output_cursor.execute('VACUUM')

output_db.commit()
output_db.close()

