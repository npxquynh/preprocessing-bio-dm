import os
import glob
from optparse import OptionParser
from itertools import combinations
import common_function as cf

D = 10
SIZE_OF_SUBLGN = 12

def read_result_from_one_iteration(folder_path, delimiter):
    expansion_lists = list()
    csv_files = cf.find_csv_filenames(folder_path)
    # csv_files = glob.glob(os.path.join(folder_path, "*.csv"))
    # print os.path.join(folder_path, "*.csv")
    # print csv_files
    # import pdb; pdb.set_trace()
    for csv_file in csv_files:
        file_path = os.path.join(folder_path, csv_file)
        nodes = cf.read_edgelist_get_node_unweighted(file_path, delimiter=delimiter)
        expansion_lists.append(nodes)

    return expansion_lists

def merge_nodes_from_expansion_lists(expansion_lists):
    nodes = set()
    for i in range(0, len(expansion_lists)):
        nodes = nodes.union(expansion_lists[i])
    return nodes

def frequency_computation(expansion_lists):
    G = dict()
    for expansion in expansion_lists:
        for probe in expansion:
            if not G.has_key(probe):
                G[probe] = 1
            else:
                G[probe] += 1
    return G

# frequency = { probe: frequency }
def write_frequency(frequency):
    output_file = open("frequency.csv", "w")
    for probe in frequency:
        output_file.write( "%s,%s\n" % (probe, frequency[probe]) )
    output_file.close()

# D: number of the subnetworks of LGN (called subLGN)
# d: size of subLGN
def create_list_subLGN(D, d, LGN):
    subLGNs = list()
    for i in range(0, D):
        subLGNs.append( set(cf.random_combination(LGN, d)) )

    return subLGNs

def internal_assessment(probes_from_expansion, LGN, subLGNs, G, D=D):
    z = 0
    for z in range(0, D):
        intra_probes = LGN.difference( subLGNs[z] )
        TP = len(probes_from_expansion.intersection(intra_probes))
        FP = len(intra_probes.intersection())

        import pdb; pdb.set_trace()

        a = 2

"""
args[0]: Known gene network (LGN)
args[1]: path to folder of output file
"""
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

    LGN = cf.read_node_unweighted(args[0])
    expansion_lists = read_result_from_one_iteration(args[1], delimiter=delimiter)
    probes_from_expansion = merge_nodes_from_expansion_lists(expansion_lists)
    G = frequency_computation(expansion_lists)
    write_frequency(G)
    subLGNs = create_list_subLGN(D, SIZE_OF_SUBLGN, LGN)

    internal_assessment(probes_from_expansion, LGN, subLGNs, G, D)

