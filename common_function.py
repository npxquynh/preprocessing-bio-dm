import os
from collections import defaultdict, OrderedDict
import re
from random import randrange

def swap(a,b):
    if a > b:
        return b,a
    return a,b

def read_edgelist_unweighted(filename,delimiter=None,nodetype=str):
    """reads two-column edgelist, returns dictionary
    mapping node -> set of neighbors and a list of edges
    """
    adj = defaultdict(set) # node to set of neighbors
    edges = set()
    nodes = set()
    for line in open(filename, 'U'):
        L = line.strip().split(delimiter)
        ni,nj = nodetype(L[0]),nodetype(L[1]) # other columns ignored
        if ni != nj and ni != "" and nj != "": # skip any self-loops...
            edges.add( swap(ni,nj) )
            nodes.add(ni)
            nodes.add(nj)
            adj[ni].add(nj)
            adj[nj].add(ni) # since undirected
    return dict(adj), edges, nodes

def read_edgelist_get_node_unweighted(filename,delimiter=None,nodetype=str):
    """reads two-column edgelist, returns dictionary
    mapping node -> set of neighbors and a list of edges
    """
    nodes = set()
    for line in open(filename, 'U'):
        L = line.strip().split(delimiter)
        ni,nj = nodetype(L[0]),nodetype(L[1]) # other columns ignored
        nodes.add(ni)
        nodes.add(nj)
    return nodes

def read_node_unweighted(filename, nodetype=str):
    """reads two-column edgelist, returns dictionary
    mapping node -> set of neighbors and a list of edges
    """
    nodes = set()
    for line in open(filename, 'U'):
        L = line.strip()
        nodes.add(nodetype(L))
    return nodes

def write_edges_to_file(edges, filename):
    output_file = open(filename, "w")
    for item in edges:
        output_file.write( "%s,%s\n" % (item[0], item[1]) )

def degree_of_node_set(adj, nodes):
    degree = defaultdict(set)

    for node in nodes:
        number_of_connection = len(adj[node])
        degree[number_of_connection].add(node)

    return degree

def print_node_degree(degree, message = None):
    if message is None:
        message = "- Frequency of node degree:"

    print message

    degree = OrderedDict(sorted(degree.items()))
    for key in degree:
        print("%d\t\t%d" % (key, len(degree[key])))

def find_big_hubs(adj, MAX_DEGREE):
    # filter nodes with large connection
    max_nodes = []

    for node in adj:
        if len(adj[node]) >= MAX_DEGREE:
            max_nodes.append(node)
    return max_nodes

def filter_nodes_with_smaller_degree(adj, nodes, MAX_DEGREE):
    filtered_nodes = set()

    for node in nodes:
        if len(adj[node]) < MAX_DEGREE:

            filtered_nodes.add(node)

    return filtered_nodes

tokenize = re.compile(r'(\d+)|(\D+)').findall

def natural_sortkey(string):
    return tuple(int(num) if num else alpha for num, alpha in tokenize(string))

def find_csv_filenames(path_to_dir, suffix=".csv"):
    filenames = os.listdir(path_to_dir)
    # 0.csv, 1.csv, 2.csv, ..... , 10.csv
    filenames = sorted(filenames, key=natural_sortkey)
    return [ filename for filename in filenames if filename.endswith( suffix ) ]

def random_combination(set_n, k):
    combination = set()
    random_position = set()
    for j in range( len(set_n) - k, len(set_n)):
        t = randrange(0, j)
        if t not in random_position:
            random_position.add(t)
        else:
            random_position.add(j)

    return [list(set_n)[i] for i in random_position]