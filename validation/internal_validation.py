import math
import operator

# TODO: remove pdb
import pdb

import helper

class InternalValidation:
    def __init__(self, genes_in_tile, blocks, genes_in_lgn, edges_in_lgn):
        self.genes_in_tile = genes_in_tile
        self.blocks = blocks
        self.genes_in_lgn = genes_in_lgn
        self.edges_in_lgn = edges_in_lgn
        self.list_intra = dict()
        self.list_extra = dict()
        self.__generate_zero_matrices()

    def __find_genes_connected_with_LGN(self, block):
        # connected_nodes[i] = set() = nodes connected with LGN[i]
        connected_nodes = []
        for gene in self.genes_in_lgn:
            nodes = set()
            for edge in block:
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
    # LGN = [1, 2, 3]
    #   [(1, 10), (1, 2, 3), (2)]
    #   0
    # Output: (1, 2, 3)
    def __find_connected_nodes_in_SubLGN_1(self, connected_nodes, index):
        # connected_nodes = list[set()]
        result = set()
        for (idx, nodes) in enumerate(connected_nodes):
            if idx != index:
                result = result.union(nodes)

        return result

    def __refine_list_intra(self):
        number_of_blocks = float(len(self.blocks))

        for key in self.list_intra:
            self.list_intra[key] /= number_of_blocks

    # Input:
    # LGN = [1, 2, 3]
    #   [(1, 10), (1, 2, 3), (2)]
    #   0
    # Output: [1, 2, 3, 2]
    def __find_connected_nodes_in_SubLGN_2(self, connected_nodes, index):
        # connected_nodes = list[set()]
        result = list()
        for (idx, nodes) in enumerate(connected_nodes):
            if idx != index:
                for gene in nodes:
                    result.append(gene)

        return result

    def __connection_in_extra_genes_step_0(self, block):
        """
        calculate c[0]: whether or not this 'extra gene' appear in a specific block
        """
        control = dict()

        number_of_LGN_genes = len(self.genes_in_lgn) + 2

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

    # TODO: connected_nodes = list[set()] so it might remove duplicated node that
    # are connected to the same gene in LGN.
    # But another thinking suggests to me that this case will never happen,
    # since if it happens, we will have duplicated connection
    # Like At1 <-> At2, and At1 <-> At2
    # So I wonder whether my thinking is correct?
    def __connection_in_extra_genes_step_1(self, control, connected_nodes):
        """
        calculate c[1]: number of times it's directly connected to the LGN

        connected_nodes = list[set()] ==> remember to 'flatten' them first
        """
        default = [0] * len(self.genes_in_lgn)

        nodes = helper.list_of_sets_to_list(connected_nodes)
        for gene in nodes:
            try:
                control.get(gene, default)[1] += 1
            except IndexError:
                print "Index error 3"

        return control

    def __connection_in_extra_genes_step_2(self, control, connected_nodes_in_subLGN, index):
        """
        calculate c[2] ... c[n+1]: number of times it's connected to SLGN(0), ...,
        number of times it's connected to SLGN(i), ... ,
        number of times it's connected to SLGN(n-1))

        connected_nodes_in_subLGN = list[]
        """
        default = [0] * len(self.genes_in_lgn)

        for gene in connected_nodes_in_subLGN:
            try:
                control.get(gene, default)[index + 2] += 1
            except IndexError:
                print "Index error 4"
        return control

    def __merge_list_extra(self, control):
        for key in control:
            # TODO: refactor with list_extra.get( , default_value)
            if key in self.list_extra:
                self.list_extra[key] = map(operator.add, self.list_extra[key], control[key])
            else:
                self.list_extra[key] = control[key]

        # return list_extra

    def __refine_list_extra(self):
        max_length = len(self.genes_in_lgn) + 2

        keys_to_be_deleted = []

        # TODO: refactor constant number of LGN genes
        for key in self.list_extra:
            if self.list_extra[key][1] == 0:
                keys_to_be_deleted.append(key)
            else:
                for j in range(2, max_length):
                    self.list_extra[key][j] /= float(self.list_extra[key][1])

        for key in keys_to_be_deleted:
            self.list_extra.__delitem__(key)

        # return list_extra

    def create_list_intra_extra(self):
        print 'create_list_intra_extra'

        # construct list_intra
        for gene in self.genes_in_lgn:
            self.list_intra[gene] = 0

        for block in self.blocks:
            connected_nodes = self.__find_genes_connected_with_LGN(block)

            # calculate list_intra
            for (index, gene) in enumerate(self.genes_in_lgn):
                connected_nodes_in_subLGN = self.__find_connected_nodes_in_SubLGN_1(connected_nodes, index)
                if connected_nodes_in_subLGN.__contains__(gene):
                    self.list_intra[gene] += 1

            self.__refine_list_intra()

            # calculate list_extra
            flattened_connected_nodes = helper.list_of_sets_to_list(connected_nodes)
            """ calculate list_extra, we just simply do it for every genes in the block, and
            then we can remove those genes in LGN
            """
            # number_of_LGN_genes = len(lgn['probe'])
            control = self.__connection_in_extra_genes_step_0(block)
            control = self.__connection_in_extra_genes_step_1(control, connected_nodes)
            for (index, gene) in enumerate(self.genes_in_lgn):
                connected_nodes_in_subLGN = self.__find_connected_nodes_in_SubLGN_2(connected_nodes, index)
                control = self.__connection_in_extra_genes_step_2(control, connected_nodes_in_subLGN, index)

            self.__merge_list_extra(control)
            self.__refine_list_extra()

    def __generate_zero_matrices():
        f = [x / 100.0 for x in range(100)]
        ncols = 100
        nrows = len(self.genes_in_lgn)
        self.TP = [ [0] * ncols for i in range(nrows) ]
        self.FN = [ [1] * ncols for i in range(nrows) ]
        self.FP = [ [0] * ncols for i in range(nrows) ]

    def calculate_TPFN(self):
        # True Positive & False Negative
        for (i, key) in enumerate(self.list_intra):
            # finding an index in which a remaining part
            # of a row is filled with 1
            # Ex:
            # 0 0 0 1
            # 0 1 1 1
            # 0 0 0 1
            equal_point = int(math.floor(self.list_intra[key] * 100))
            for j in range(equal_point, ncols):
                self.TP[i][j] = 1
                self.FN[i][j] = 0

        # False Positive
        # TODO: number or LGN genes
        for i in range(self.genes_in_lgn):
            c = 0
            try:
                for key in self.list_extra:
                    c += self.list_extra[key][i + 2]
            except IndexError:
                print "Index error 5"

    def statistical_result(TP, FN, FP):
        SE = divide_two_dim_array(TP, add_two_dim_array(TP, FN))
        return SE