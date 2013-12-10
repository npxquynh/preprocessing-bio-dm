<<<<<<< HEAD
=======
# internal_validation_2.py
# Refactor code
>>>>>>> b0f20e9bcd1324182d6ecc5e1be4200fc6c710b6

import csv
import os
import fnmatch
import pdb
import operator
import math

from collections import defaultdict

numOfGenes = 0

lgn = defaultdict(list)
#SLGN(i)= subnetwork obtained by removing the i-th gene from the LGN (i=0,.....,n-1)
slgn = []
t = 23
blocks = []

def read_lgn():
    global lgn
    global numOfGenes
    with open('At_LGN_genes.csv', 'rb') as f:
        reader = csv.DictReader(f)
        for row in reader:
            for (k,v) in row.items():
                lgn[k].append(v)

    numOfGenes = len (lgn['probe'])


def create_slgn():
    global slgn

    for x in range(0, len (lgn['probe'])):
         slgn.append([])
         for y in range(0, len (lgn['probe'])):
            slgn[x].append( lgn['probe'][y])

         del slgn[x][x]





def read_blocks():
    print 'read_blocks'
    global blocks
    #i = 0
    for dirpath, dirs, files in os.walk('pc1'):
        for filename in fnmatch.filter(files, '*.csv'):
            with open(os.path.join(dirpath, filename)) as f:
                #blocks.append([])
                datacsv = csv.reader(f, delimiter=',')
                blocks.append (list(datacsv))



                '''
                blocks[i] = defaultdict(list)

                reader = csv.DictReader(f)
                for row in reader:
                    for (k,v) in row.items():
                        print 'key: ', k , '----', v
                        blocks[i][k].append(v)
                        '''

           # i = i + 1

def find_genes_connected_with_LGN(block, LGN):
    # connected_nodes[i] = set() = nodes connected with LGN[i]
    connected_nodes = []
    for gene in LGN:
        nodes = set()
        for edge in block:
            # print edge
            try:
                if edge[0] == gene:
                    nodes.add(edge[1])
                if edge[1] == gene:
                    nodes.add(edge[0])
            except IndexError:
                print "Index error"
        connected_nodes.append(nodes)
    return connected_nodes

# Input:
#   [(1, 10), (1, 2, 3), (2)]
#   0
# Output: (1, 2, 3)
def find_connected_nodes_in_SubLGN_1(connected_nodes, index):
    # connected_nodes = list[set()]
    result = set()
    for (idx, nodes) in enumerate(connected_nodes):
        if idx != index:
            result = result.union(nodes)

    return result

def refine_list_intra(list_intra):
    # TODO: constant number of blocks
    global blocks

    number_of_blocks = float(len(blocks))

    for key in list_intra:
        list_intra[key] /= number_of_blocks

    return list_intra

# Input:
#   [(1, 10), (1, 2, 3), (2)]
#   0
# Output: [1, 2, 3, 2]
def find_connected_nodes_in_SubLGN_2(connected_nodes, index):
    # connected_nodes = list[set()]
    result = list()
    for (idx, nodes) in enumerate(connected_nodes):
        if idx != index:
            for gene in nodes:
                result.append(gene)

    return result

def list_of_sets_to_list(l):
    result = []
    for s in l:
        for element in s:
            result.append(element)
    return result

"""
calculate c[0]
lgn is exactly lgn['probe'], not a dict anymore
"""
def connection_in_extra_genes_step_0(block, lgn):
    control = dict()

    number_of_LGN_genes = len(lgn) + 2

    for edge in block:
        try:
            # number of times it appears in the block
            if edge[0] in control:
                control[edge[0]][0] += 1
            else:
                control[edge[0]] = [0] * number_of_LGN_genes
                control[edge[0]][0] = 1

            if edge[1] in control:
                control[edge[1]][0] += 1
            else:
                control[edge[1]] = [0] * number_of_LGN_genes
                control[edge[1]][0] = 1

            # number of times it's connected with the LGN
            # number of times it's connected with the subLGN
            # will be calculated later

        except IndexError:
            print "Index error 2"
        except KeyError:
            pass

    return control

"""
calculate c[1]: number of times it's directly connected to the LGN

connected_nodes = list[set()] ==> remember to 'flatten' them first
"""

# TODO: connected_nodes = list[set()] so it might remove duplicated node that
# are connected to the same gene in LGN.
# But another thinking suggests to me that this case will never happen,
# since if it happens, we will have duplicated connection
# Like At1 <-> At2, and At1 <-> At2
# So I wonder whether my thinking is correct?
def connection_in_extra_genes_step_1(control, connected_nodes):
    default = [0] * len(lgn['probe'])

    nodes = list_of_sets_to_list(connected_nodes)
    for gene in nodes:
        try:
            control.get(gene, default)[1] += 1
        except IndexError:
            print "Index error 3"

    return control

"""
calculate c[2] ... c[n+1]: number of times it's connected to SLGN(0), ...,
number of times it's connected to SLGN(i), ... ,
number of times it's connected to SLGN(n-1))

connected_nodes_in_subLGN = list[]
"""
def connection_in_extra_genes_step_2(control, connected_nodes_in_subLGN, index):
    default = [0] * len(lgn)

    for gene in connected_nodes_in_subLGN:
        try:
            control.get(gene, default)[index + 2] += 1
        except IndexError:
            print "Index error 4"
    return control

def merge_list_extra(list_extra, control):
    for key in control:
        # TODO: refactor with list_extra.get( , default_value)
        if key in list_extra:
            list_extra[key] = map(operator.add, list_extra[key], control[key])
        else:
            list_extra[key] = control[key]

    return list_extra

def refine_list_extra(list_extra):
    global lgn
    max_length = len(lgn['probe']) + 2

    keys_to_be_deleted = []

    # TODO: refactor constant number of LGN genes
    for key in list_extra:
        if list_extra[key][1] == 0:
            keys_to_be_deleted.append(key)
        else:
            for j in range(2, max_length):
                list_extra[key][j] /= float(list_extra[key][1])

    for key in keys_to_be_deleted:
        list_extra.__delitem__(key)

    return list_extra

def create_list_intra_extra():
    print 'create_list_intra_extra'

    global lgn
    global blocks

    list_intra = dict()
    list_extra = dict()

    # construct list_intra
    for gene in lgn['probe']:
        list_intra[gene] = 0


    lgn_set = set(lgn['probe'])
    lgn_list = list(lgn['probe'])

    for block in blocks:
        connected_nodes = find_genes_connected_with_LGN(block, lgn_set)

        # calculate list_intra
        for (index, gene) in enumerate(lgn['probe']):
            connected_nodes_in_subLGN = find_connected_nodes_in_SubLGN_1(connected_nodes, index)
            if connected_nodes_in_subLGN.__contains__(gene):
                list_intra[gene] += 1

        list_intra = refine_list_intra(list_intra)

        # calculate list_extra
        flattened_connected_nodes = list_of_sets_to_list(connected_nodes)
        """ calculate list_extra, we just simply do it for every genes in the block, and
        then we can remove those genes in LGN
        """
        # number_of_LGN_genes = len(lgn['probe'])
        control = connection_in_extra_genes_step_0(block, lgn['probe'])
        control = connection_in_extra_genes_step_1(control, connected_nodes)
        for (index, gene) in enumerate(lgn['probe']):
            connected_nodes_in_subLGN = find_connected_nodes_in_SubLGN_2(connected_nodes, index)
            control = connection_in_extra_genes_step_2(control, connected_nodes_in_subLGN, index)

        list_extra = merge_list_extra(list_extra, control)
        list_extra = refine_list_extra(list_extra)

    return list_intra, list_extra

def calculate_TPFN(list_intra, list_extra):
    global lgn

    f = [x / 100.0 for x in range(100)]
    ncols = 100
    nrows = len(lgn['probe'])
    TP = [ [0] * ncols for i in range(nrows) ]
    FN = [ [1] * ncols for i in range(nrows) ]
    FP = [ [0] * ncols for i in range(nrows) ]

    # True Positive & False Negative
    for (i, key) in enumerate(list_intra):
        # finding an index in which a remaining part
        # of a row is filled with 1
        # Ex:
        # 0 0 0 1
        # 0 1 1 1
        # 0 0 0 1
        equal_point = int(math.floor(list_intra[key] * 100))
        for j in range(equal_point, ncols):
            TP[i][j] = 1
            FN[i][j] = 0

    # False Positive
    # TODO: number or LGN genes
    for i in range(len(lgn['probe'])):
        c = 0
        try:
            for key in list_extra:
                c += list_extra[key][i + 2]
        except IndexError:
            print "Index error 5"

    return TP, FN, FP
        # FP[i,]

def add_two_dim_array(a1, a2):
    nrows = len(a1)
    ncols = len(a1[0])

    for row in range(nrows):
        for col in range(ncols):
            a1[row][col] += a2[row][col]
    return a1

def divide_two_dim_array(a1, a2):
    nrows = len(a1)
    ncols = len(a1[0])

    for row in range(nrows):
        for col in range(ncols):
            a1[row][col] /= float(a2[row][col])
    return a1

def statistical_result(TP, FN, FP):
    # pdb.set_trace()
    SE = divide_two_dim_array(TP, add_two_dim_array(TP, FN))
    return SE

if __name__ == '__main__':
    read_lgn()
    create_slgn()
    read_blocks()
    list_intra, list_extra = create_list_intra_extra()
    TP, FN, FP = calculate_TPFN(list_intra, list_extra)
    SE = statistical_result(TP, FN, FP)
    pdb.set_trace()

    a = 2
