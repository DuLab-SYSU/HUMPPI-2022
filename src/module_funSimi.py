import pandas as pd
import numpy as np
import re
import sys
import argparse
import itertools
import multiprocessing as mp

gobps = {}
with gzip.open(r'../data/MSigDB7.4/c5.bp.v7.4.symbols_GO.ID.txt', 'rt) as f:
    for row in f:
        col = row.strip().split("\t")
        gobps.setdefault(col[0], []).append(col[2])

def f_sim(g1, g2):
    size = []
    for i, g in gobps.items():
        if (g1 in g) & (g2 in g):
            size.append(len(g))
            
    sim = 0 if len(size) == 0 else 2/np.min(size)
    return sim

def ave_fs(gene_list):
    gps = list(itertools.combinations(gene_list, 2))
    with mp.Pool(mp.cpu_count()) as pool:
        res = pool.starmap(f_sim, gps)
        pool.close()
        pool.join()
        return np.mean(res)

def read_(I):
    _out = "../intermediate/batch_MCL_out/out.HUMPPI2022.I" + str(I)
    Clusters = {}
    f = open(_out)
    for i, v in enumerate(f):
        cols = v.strip().split('\t')
        if len(cols) > 2:
            Clusters[str(i+1)] = cols
    f.close()
    return Clusters

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outputfile", action="store", dest="out_filename",
                        required=True, help="non-redundancy complex file")
    args = parser.parse_args()

    out = args.out_filename

    out_f = open(out, "a")
    out_f.write("Inflation" + "\t" + "Cluster_ID" + "\t" + "Mean_similarity" + "\n")
    for i in range(150, 1501):
        clusters = read_(i)
        for c in clusters:
            average_fs = ave_fs(clusters[c])
            out_f.write(str(i/100) + "\t" + c + "\t" + str(average_fs) + "\n")

if __name__ == "__main__":
    main(sys.argv[1:])
