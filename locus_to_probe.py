import sys, os
import csv
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

def read_locus_probe_mapping(filename, delimiter=None, nodetype=str):
    mapping = dict() # dict with key = locus, value = probe
    for line in open(filename, 'U'):
        print line
        L = line.strip('"').split(delimiter)
        print L

        probe, locus = nodetype(L[0]), nodetype(L[1])
        if mapping.has_key(locus):
            print "1 locus is not mapping to 1 probe"
        mapping[locus] = probe

    # print mapping
    return mapping

def locus_to_probe_analysis(nodes, mapping):
    """reads two-column that map locus with probe, returns dictionary
    mapping locus -> probe
    """
    locus_without_associative_probe = list()
    for node in nodes:
        if mapping.has_key(node) == False:
            locus_without_associative_probe.append(node)

    """ Output to file
    """
    output_file = open("./analysis/locus_without_probe.csv", "w")
    for item in locus_without_associative_probe:
        output_file.write("%s\n" % item)
    output_file.close()

    print("Locus without associate probe: %d" % len(locus_without_associative_probe))

def write_edges_to_file(edges, filename):
    output_file = open(filename, "w")
    for item in edges:
        output_file.write( "%s,%s\n" % (item[0], item[1]) )

def grn_with_probe(edges, mapping):
    new_edges = set()
    edges_without_associative_probe = set()

    for item in edges:
        x, y = item
        # print("%s - %s - %s -%s" % (x, y, str(mapping.has_key(x)), str(mapping.has_key(y))))
        if(mapping.has_key(x) and mapping.has_key(y)):
            new_edges.add( swap(x,y) )
        else:
            edges_without_associative_probe.add( swap(x,y) )

    # Analysis
    print("Total number of edges: %d", len(edges))
    print("Number of edges without matching probe: %d" % len(edges_without_associative_probe) )

    # Generate new network file
    write_edges_to_file(new_edges, "./analysis/edges_in_probe.csv")

    # Write out those connections that do not have the associative probe.
    write_edges_to_file(edges_without_associative_probe, "./analysis/edges_without_associative_probe.csv")

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
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    delimiter = options.delimiter
    if delimiter == '\\t':
        delimiter = '\t'


    print "# loading network from edgelist..."
    basename = os.path.splitext(args[0])[0]
    adj,edges,nodes = read_edgelist_unweighted(args[0], delimiter=delimiter)

    locus_probe_mapping = read_locus_probe_mapping(args[1], delimiter=delimiter)

    print "# start anaylysis..."
    locus_to_probe_analysis(nodes, locus_probe_mapping)

    print "# create new network file, replacing locus with probe"
    grn_with_probe(edges, locus_probe_mapping)

