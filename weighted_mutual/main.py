# main file
import sys, getopt
import numpy as np
import pandas as pd
from weighted_mutual.weight import *
from weighted_mutual.entropy import *
from itertools import product, combinations
from multiprocessing import Pool

# args

def usage():
    print("""
    This is a script to obtain domain weighted mutual information for a bunch of protein

    # Usage
    python main.py --input <filename> --output <filename> [--threads 1]

    input: file has to be comma seperated, with column names as domain names and row names as protein id. see 'test_file/small_test' as an example
    output: select the desired location. will be stored as csv: protein1, protein2, weighted_mutual_information
    threads: how many threads do you want
    """)
def wrap_mutual_info(args):
    protein_one, protein_two = args

    global p
    global entropy_dict

    global chunk
    global other_chunk

    global mode
    global outfile

    # calculate mutual info
    if mode == 'one_chunk':
        H_xy = weighted_joint_entropy(chunk.loc[protein_one].values, chunk.loc[protein_two], p)
    if mode == 'two_chunks':
        H_xy = weighted_joint_entropy(chunk.loc[protein_one].values, other_chunk.loc[protein_two], p)

    m = mutual_info(entropy_dict[protein_one], entropy_dict[protein_two], H_xy)

    # write to file
    with open(outfile, 'a') as f:
        f.write(','.join([protein_one, protein_two, str(m)])+ '\n')

# main
argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv, "hi:o:", ['input=', 'output=', 'threads='])
except getopt.GetoptError as err:
    # print help message and exit:
    print(err)
    usage()
    sys.exit(2)

outfile = None
infile = None
thread = 1
verbose = False
for o, a in opts:
    if o == '-v':
        verbose = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-o", "--output"):
        outfile = a
    elif o in ("-i", "--input"):
        infile = a
    elif o in ("-t", "--threads"):
        thread = int(a)
    else:
        assert False, "unhandled option"
# calculate weight
print("""
Calculating weights and probability for each domain
===================================================
""")
p = probability(weight(infile))

# read file in chunks
chunks = pd.read_csv(infile, chunksize = 50, header = 0, index_col = 0)

# write header
with open(outfile, 'w') as f:
    f.write('gene_one,gene_two,weighted_mutual\n')

entropy_dict = {}
# calculate by chunk
chunk_no = 0
for chunk in chunks:

    if chunk_no == 0: # only in first round there will be no entropy
        print("""
        Calculating weighted entropy for each protein in chunk {0}
        """.format(chunk_no))

        for protein_id in chunk.index:
            e = weighted_entropy(chunk.loc[protein_id].values, p)
            entropy_dict[protein_id] = e

    # calculate joint entropy and mutual_info
    pairs = list(combinations(chunk.index, 2))
    mode = 'one_chunk'
    with Pool(thread) as po:
        po.map(wrap_mutual_info, pairs)



    # combinations regarding to other chunks
    other_chunk_no = 0
    other_chunks = pd.read_csv(infile, chunksize = 50, header = 0, index_col = 0)

    for other_chunk in other_chunks:
        if other_chunk_no <= chunk_no:
            # we have already done that
            other_chunk_no += 1
        else:
            print("""
            Processing combinations in {0} {1}
            """.format(chunk_no, other_chunk_no))

            # calculate entropy for proteins in the other chunk, if it's the first round
            if chunk_no == 0:
                for protein_id in other_chunk.index:
                    e = weighted_entropy(other_chunk.loc[protein_id].values, p)
                    entropy_dict[protein_id] = e


            pairs = list(product(chunk.index, other_chunk.index))

            mode = 'two_chunks'
            with Pool(thread) as po:
                po.map(wrap_mutual_info, pairs)

            other_chunk_no += 1

    chunk_no += 1
