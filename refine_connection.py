import sys, os
from collections import defaultdict, OrderedDict
from optparse import OptionParser

import common_function

def edges_from_adj(adj):
    edges = set()
    for origin_node in adj:
        for destination_node in adj[origin_node]:
            edges.add( swap(origin_node, destination_node) )

            # Remove edge from the destination_node to the origin node
            adj[destination_node].remove(origin_node)

    return edges

"""
There are 2 big hubs in the networks, each have more around 4000 connections.
And there are lots of node that has a single link with one of two big hubs.
"""
def refine_connections(adj, nodes):
    MAX_DEGREE = 3000
    max_nodes = find_big_hubs(adj, MAX_DEGREE)

    nodes_with_degree_1 = set()
    for max_node in max_nodes:
        nodes_with_degree_1 = nodes_with_degree_1.union(filter_nodes_with_smaller_degree(adj, adj[max_node], 2))

    # Remove nodes
    new_nodes = nodes.difference(nodes_with_degree_1)

    # Remove connections. We have to do it in both direction
    for node in nodes_with_degree_1:
        # Remove edge from the big hub to nodes.
        for adj_node in adj[node]:
            adj[adj_node].remove(node)
        # Remove edge from the node to the hubs
            adj.pop(node)

    # Generating edges
    new_edges = edges_from_adj(adj)

    # Write the refined connection to file.
    write_edges_to_file(new_edges, "./analysis/refined_grn2.csv")

    print "Number of nodes removed: %d" % len(nodes_with_degree_1)

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

    print "# refine connection by removing node..."
    refine_connections(adj, nodes)

