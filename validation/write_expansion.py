import os

def write_expansion_list(filepath, file_content):
    """
    file_content = [
        [gene_name, frequency_1],
        [gene_name_2, frequency_2]
    ]
    """
    with open(filepath, 'w') as output_file:
        for line in file_content:
            output_file.write('%s,%s\n' % (line[0], line[1]))
