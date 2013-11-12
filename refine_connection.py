import sys, os
from collections import defaultdict, OrderedDict
from optparse import OptionParser

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
        if ni != nj: # skip any self-loops...
            edges.add( swap(ni,nj) )
            nodes.add(ni)
            adj[ni].add(nj)
            adj[nj].add(ni) # since undirected
    return dict(adj), edges, nodes

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

def node_degree(adj, nodes):
    """
    1st: degree for all nodes
    """
    degree = degree_of_node_set(adj, nodes)
    # l = len(adj)
    # for key in adj:
    #     adj_length = len(adj[key])
    #     degree[adj_length].add(key)

    print_node_degree(degree)

    """
    2nd: degree of those nodes that are neighbor of HUGE node
    """
    print "- Node with more than 3000 connections:"
    MAX_DEGREE = 3000
    max_nodes = []

    # filter nodes with large connection
    for key in degree:
        if key >= MAX_DEGREE:
            max_nodes += list(degree[key])

    # calculate connection for neighbors of those max_nodes
    for node in max_nodes:
        degree = degree_of_node_set(adj, adj[node])
        print_node_degree(degree)

def analysis(adj, edges, nodes):
    print("- Number of nodes: %d" % len(adj))
    print("- Number of edges: %d" % len(edges))
    node_degree(adj, nodes)

if __name__ == '__main__':
  # build option parser:
    class MyParser(OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    usage = "usage: python %prog [options] filename"
    description = """Analysis of Gene Regulatory Network
    """

    epilog = """
    This program gave general information for graph network
    """

    # parse options:
    parser = MyParser(usage, description=description,epilog=epilog)
    parser.add_option("-d", "--delimiter", dest="delimiter", default="\t",
                      help="delimiter of input & output files [default: tab]")

    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    delimiter = options.delimiter
    if delimiter == '\\t':
        delimiter = '\t'


    print "# loading network from edgelist..."
    basename = os.path.splitext(args[0])[0]
    adj,edges,nodes = read_edgelist_unweighted(args[0], delimiter=delimiter)

    print "# start anaylysis..."
    analysis(adj, edges, nodes)

