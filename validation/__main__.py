import os
import helper
import parse_xml as px
import parse_lgn as pl
import parse_expanded_network as pen
import write_expansion as we
from internal_validation import *
from internal_validation_2 import *

import pdb

"""
1. Read xml --> find xml_mapping
2. Create blocks
"""

XML_MAPPING_FOLDER = '../xml_mapping'
PC_RESULT_FOLDER = '../expanded_network'
LGN_FOLDER = '../lgns'
PCIM_FOLDER_1 = '../pcim_result_1'
PCIM_FOLDER_2 = '../pcim_result_2'

xml_list = px.xml_list(XML_MAPPING_FOLDER)

for xml_filepath in xml_list:
    # common things used in both type of internal validation
    xml_filename = os.path.basename(xml_filepath)

    obs_filename, lgn_filename = px.read_single_xml(xml_filepath)

    lgn_filepath = helper.generate_filepath(LGN_FOLDER, lgn_filename)
    genes_in_lgn, edges_in_lgn = pl.read_lgn(lgn_filepath)

    list_of_expanded_network_file = pen.merge_different_pc_run(xml_filename,
        PC_RESULT_FOLDER)

    # no expanded network found for this xml file
    # maybe the result will be run at a later time
    if len(list_of_expanded_network_file) == 0:
        continue

    genes_in_tiles, blocks = pen.read_expanded_network(list_of_expanded_network_file)

    print "Internal Validation for %s" % xml_filename
    # 1st internal validation
    IV = InternalValidation(genes_in_tiles, blocks, genes_in_lgn, edges_in_lgn)

    expansion_list = IV.expansion_list()

    expansion_filepath = helper.generate_filepath(PCIM_FOLDER_1, xml_filename)
    we.write_expansion_list(expansion_filepath, expansion_list)

    # 2nd internal validation
    # IV2 = InternalValidationRls(genes_in_tiles, blocks, genes_in_lgn, edges_in_lgn)

    # expansion_list_2 = IV2.expansion_list()
    # expansion_filepath = helper.generate_filepath(PCIM_FOLDER_2, xml_filename)
    # we.write_expansion_list(expansion_filepath, expansion_list_2)

