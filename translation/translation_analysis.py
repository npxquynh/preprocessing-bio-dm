import os, sys, re
import pdb

from itertools import chain

# import file from different folder
sys.path.insert(0, '../')
import common_function as cf

def translation_analysis(not_translated_locus, file_name):
    locus_description = dict()

    locus_probe_file = open(file_name, "U")
    index = 1
    for line in locus_probe_file:
        temp_locus_description = look_for_locus(not_translated_locus, str(line))
        locus_description = merge_locus_description(locus_description, temp_locus_description)
        index += 1

    write_locus_description(locus_description)

def look_for_locus(not_translated_locus, text_line):
    locus_description = dict()
    text_line = text_line.lower()
    for locus in not_translated_locus:
        if text_line.find(locus.lower()) != -1:
            if locus_description.has_key(locus):
                locus_description[locus].append(text_line)
            else:
                locus_description[locus] = [text_line]

    return locus_description

def merge_locus_description(locus_description, temp_locus_description):
    for key in temp_locus_description:
        if locus_description.has_key(key):
            locus_description[key] = chain(locus_description[key], temp_locus_description[key])
        else:
            locus_description[key] = temp_locus_description[key]

    return locus_description

def write_locus_description(locus_description):
    output_file = open("translation_analysis.txt", "w")
    for key in locus_description:
        for description in locus_description[key]:
            line = chain([key], description.strip().split("\t"))
            # pdb.set_trace()
            output_file.write('\t'.join(map(str, line)))
            output_file.write('\n')

    output_file.close()

"""
args[0]: Known gene network (LGN)
args[1]: path to folder of output file
"""
if __name__ == '__main__':
    print "# loading not translated locus..."
    file_name = "not_translated_locus.csv"
    not_translated_locus = cf.read_node_unweighted(file_name)

    print "# locus - probe translationg analysis"
    translation_analysis(list(not_translated_locus), "affy_ATH1_array_elements-2010-12-20.txt")
