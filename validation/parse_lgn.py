import os
import csv
import pdb
import logging

import helper

def read_lgn(filepath):
    edges_in_lgn = []
    with open(filepath) as lgn_file:
        lgn_reader = csv.reader(lgn_file, delimiter=';')
        lgn_reader.next() # skip the header row

        edges_in_lgn = []
        for row in lgn_reader:
            try:
                gene_1, gene_2 = helper.swap(row[0], row[1])
                edges_in_lgn.append([gene_1, gene_2])
            except IndexError, e:
                logging.warning("%s is not in the correct format" % row)
                pass

    genes_in_lgn = helper.genes_from_edges(edges_in_lgn)

    return genes_in_lgn, edges_in_lgn

