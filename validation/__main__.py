import os
import helper
import parse_xml as px
import parse_lgn as pl
import parse_expanded_network as pen
import write_expansion as we
from internal_validation import *

import pdb

"""
1. Read xml --> find xml_mapping
2. Create blocks
"""

XML_MAPPING_FOLDER = "../xml_mapping"
PC_RESULT_FOLDER = "../expanded_network"
LGN_FOLDER = "../lgns"
PCIM_FOLDER = '../pcim_result'

xml_list = px.xml_list(XML_MAPPING_FOLDER)

for xml_filepath in xml_list:
    xml_filename = os.path.basename(xml_filepath)
    obs_filename, lgn_filename = px.read_single_xml(xml_filepath)

    lgn_filepath = helper.generate_filepath(LGN_FOLDER, lgn_filename)
    genes_in_lgn, edges_in_lgn = pl.read_lgn(lgn_filepath)

    list_of_expanded_network_file = pen.merge_different_pc_run(xml_filename,
        PC_RESULT_FOLDER)
    genes_in_tiles, blocks = pen.read_expanded_network(list_of_expanded_network_file)

    IV = InternalValidation(genes_in_tiles, blocks, genes_in_lgn, edges_in_lgn)

    expansion_list = IV.expansion_list()

    expansion_filepath = helper.generate_filepath(PCIM_FOLDER, xml_filename)
    we.write_expansion_list(expansion_filepath, expansion_list)


