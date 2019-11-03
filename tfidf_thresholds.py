# TODO: calculate tfidf thresholds
import sqlite3
import numpy as np

conn = sqlite3.connect('cellmesh/data/cellmesh.db')
cursor = conn.cursor()

# calculate mouse threshold
cursor.execute('SELECT tfidf FROM cell_gene WHERE taxid=10090')
results = cursor.fetchall()
results = [x[0] for x in results]
results = np.array(results)
results_sorted = np.sort(results)
percentile_90 = results_sorted[int(results_sorted.shape[0]*0.90)]
print('90th percentile tfidf for mouse:', percentile_90)
percentile_93 = results_sorted[int(results_sorted.shape[0]*0.93)]
print('93th percentile tfidf for mouse:', percentile_93)

# calculate human threshold
cursor.execute('SELECT tfidf FROM cell_gene WHERE taxid=9606')
results = cursor.fetchall()
results = [x[0] for x in results]
results = np.array(results)
results_sorted = np.sort(results)
percentile_90 = results_sorted[int(results_sorted.shape[0]*0.90)]
print('90th percentile tfidf for human:', percentile_90)
percentile_93 = results_sorted[int(results_sorted.shape[0]*0.93)]
print('93th percentile tfidf for human:', percentile_93)

# calculate worm threshold
cursor.execute('SELECT tfidf FROM cell_gene WHERE taxid=6239')
results = cursor.fetchall()
results = [x[0] for x in results]
results = np.array(results)
results_sorted = np.sort(results)
percentile_90 = results_sorted[int(results_sorted.shape[0]*0.90)]
print('90th percentile tfidf for worm:', percentile_90)
percentile_93 = results_sorted[int(results_sorted.shape[0]*0.93)]
print('93th percentile tfidf for worm:', percentile_93)

conn.close()

conn = sqlite3.connect('cellmesh/data/anatomy_mesh.db')
cursor = conn.cursor()

print('ANATOMY RESULTS')
# calculate mouse threshold
cursor.execute('SELECT tfidf FROM cell_gene WHERE taxid=10090')
results = cursor.fetchall()
results = [x[0] for x in results]
results = np.array(results)
results_sorted = np.sort(results)
percentile_90 = results_sorted[int(results_sorted.shape[0]*0.90)]
print('90th percentile tfidf for mouse:', percentile_90)
percentile_93 = results_sorted[int(results_sorted.shape[0]*0.93)]
print('93th percentile tfidf for mouse:', percentile_93)

# calculate human threshold
cursor.execute('SELECT tfidf FROM cell_gene WHERE taxid=9606')
results = cursor.fetchall()
results = [x[0] for x in results]
results = np.array(results)
results_sorted = np.sort(results)
percentile_90 = results_sorted[int(results_sorted.shape[0]*0.90)]
print('90th percentile tfidf for human:', percentile_90)
percentile_93 = results_sorted[int(results_sorted.shape[0]*0.93)]
print('93th percentile tfidf for human:', percentile_93)

# calculate worm threshold
cursor.execute('SELECT tfidf FROM cell_gene WHERE taxid=6239')
results = cursor.fetchall()
results = [x[0] for x in results]
results = np.array(results)
results_sorted = np.sort(results)
percentile_90 = results_sorted[int(results_sorted.shape[0]*0.90)]
print('90th percentile tfidf for worm:', percentile_90)
percentile_93 = results_sorted[int(results_sorted.shape[0]*0.93)]
print('93th percentile tfidf for worm:', percentile_93)

conn.close()
