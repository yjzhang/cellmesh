# remove data that belong to the "Cellular Structures" subtree

import json

with open('data/mesh_desc2019_cell.json') as f:
    data = json.load(f)

cell_component_tree = 'A11.284'
chromosome_tree = 'A11.284.187'
cell_components = []
chromosomes = []
for cell_id, desc in data.items():
    tree_nums = desc['TreeNumberList']
    is_cell_component = False
    is_chromosome = False
    for t in tree_nums:
        if t.startswith(cell_component_tree):
            is_cell_component = True
        if t.startswith(chromosome_tree):
            is_chromosome = True
    if is_cell_component:
        cell_components.append(cell_id)
        print('cell component: ', cell_id, desc['DescriptorName'])
    if is_chromosome:
        chromosomes.append(cell_id)
        print('chromosome: ', cell_id, desc['DescriptorName'])

import numpy as np
np.savetxt('cell_component_ids.txt', cell_components, fmt='%s')
np.savetxt('chromosome_ids.txt', chromosomes, fmt='%s')

# TODO: append to db??
# or not...
"""
import sqlite3
conn = sqlite3.connect('cellmesh/data/cellmesh.db')
c = conn.cursor()
c.execute('CREATE TABLE cell_components(cellID text, is_chromosome_or_cell_component text)')
for x in cell_components:
    c.execute('INSERT INTO cell_components VALUES (?, ?)', (x, 'cell_component'))
for x in cell_components:
    c.execute('INSERT INTO cell_components VALUES (?, ?)', (x, 'chromosome'))
c.commit()
c.execute('CREATE INDEX gene_info_id_index ON gene_info(gene)')
"""
