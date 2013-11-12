import sys, os
from collections import defaultdict, OrderedDict
from optparse import OptionParser

import common_function


def node_degree(adj, nodes):
    """
    1st: degree for all nodes
    """
    degree = degree_of_node_set(adj, nodes)
    print_node_degree(degree)

    """
    2nd: degree of those nodes that are neighbor of HUGE node
    """
    print "- Node with more than 3000 connections:"
    MAX_DEGREE = 3000
    max_nodes = find_big_hubs(adj, MAX_DEGREE)

    # calculate connection for neighbors of those max_nodes
    print "Calculating the degree of neighbor nodes:"
    for node in max_nodes:
        print("%s with %d connections" % (node, len(adj[node])))
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
