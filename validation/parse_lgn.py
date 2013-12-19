import os
import csv
import pdb

import helper

def read_lgn(filepath):
    edges_in_lgn = []
    with open(filepath) as lgn_file:
        lgn_reader = csv.reader(lgn_file)
        lgn_reader.next() # skip the header row

        edges_in_lgn = [row for row in lgn_reader]

    genes_in_lgn = helper.genes_from_edges(edges_in_lgn)

    return genes_in_lgn, edges_in_lgn

