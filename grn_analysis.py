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
            nodes.add(nj)
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

def find_big_hubs(adj, MAX_DEGREE):
    # filter nodes with large connection
    max_nodes = []

    for node in adj:
        if len(adj[node]) >= MAX_DEGREE:
            max_nodes.append(node)

    return max_nodes

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
    max_nodes = find_big_hubs(degree, MAX_DEGREE)

    # calculate connection for neighbors of those max_nodes
    print "Calculating the degree of neighbor nodes:"
    for node in max_nodes:
        print("%s with %d connections" % (node, len(adj[node])))
        degree = degree_of_node_set(adj, adj[node])
        print_node_degree(degree)

def filter_nodes_with_smaller_degree(adj, nodes, MAX_DEGREE):
    filtered_nodes = set()

    for node in nodes:
        if len(adj[node]) < MAX_DEGREE:

            filtered_nodes.add(node)

    return filtered_nodes

def edges_from_adj(adj):
    edges = set()
    for origin_node in adj:
        for destination_node in adj[origin_node]:
            edges.add( swap(origin_node, destination_node) )

            # Remove edge from the destination_node to the origin node
            adj[destination_node].remove(origin_node)

    return edges

def write_edges_to_file(edges, filename):
    output_file = open(filename, "w")
    for item in edges:
        output_file.write( "%s,%s\n" % (item[0], item[1]) )

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
    write_edges_to_file(new_edges, "./analysis/refined_grn.csv")

    print "Number of nodes removed: %d" % len(nodes_with_degree_1)

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

    print "# refine connection by removing node..."
    refine_connections(adj, nodes)
