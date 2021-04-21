from scipy.sparse import csr_matrix, lil_matrix, save_npz
import pandas as pd
import os
import numpy as np
from itertools import chain
from optparse import OptionParser

def domain_to_pivot(gold_anno, col = 'domain', index = 'gene_id', term_mapper = None):
    '''
    Convert doamin in gold_anno into np sparse matrix
    '''
    with_data = gold_anno.loc[gold_anno[col].notnull()]
    all_terms = set(list(chain.from_iterable(with_data[col].tolist())))
    
    # mapping term to integer
    if term_mapper:
        pass
    else:
        term_mapper = {}
        for i, term in enumerate(all_terms):
            term_mapper[term]=i
    
    # initialize sparse matrix
    ma = lil_matrix((gold_anno[index].max()+1, len(term_mapper)))
    
       
    for i, terms in zip(with_data[index], with_data[col]):
        gene_id = [i]
        target_id = [term_mapper[t] for t in terms if t in term_mapper.keys()]
        ma[gene_id, target_id] = [1]*len(target_id)
        
    
    return ma, term_mapper

def option_parser():
    
    usage = """
        THIS IS WEIGHTED MUTUAL INFORMATIon FOR DOMAIN NET 1.0.0
        python generate_pivot_table.py --goldanno subset.gold_anno_df --global whole.gold_anno_df --use_global -o weighted_mutual/
        python generate_pivot_table.py -g ../test_file/gold_anno_df --global ../test_file/global_anno_df --use_global -o ../results/ 
        """
    description = """This is a helper script to obtain domain weighted mutual information for a bunch of protein
        python generate_pivot_table.py --goldanno subset.gold_anno_df --global whole.gold_anno_df --use_global -o weighted_mutual/

        Weighted mutual information analysis substantially improves domain-based functional network models https://academic.oup.com/bioinformatics/article/32/18/2824/1744011
        """

    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("-g", "--goldanno",dest="gold_anno",
                  help="path to gold_anno pickle dataframe")
    parser.add_option("--global",dest="global_gold_anno",default = None,
                  help="gold_anno of the whole pan-genome, or proteome space")
    parser.add_option("-o", "--outdir",dest="outdir", default = '/tmp', type = "string",
                  help="output directory")
    parser.add_option("--use_global", dest="use_global", action='store_true',
                  help="generate global pivot table or not")
    
    (options, args) = parser.parse_args()   

    return options

if __name__=='__main__':
    options = option_parser()
    # generate local pivot table
    sub_gold_anno = pd.read_pickle(options.gold_anno)
    ma, term_mapper = domain_to_pivot(sub_gold_anno)
    save_npz(os.path.join(options.outdir, 'domain_binary.npz'), ma.tocsr())   
    print('Saving to {}'.format(os.path.join(options.outdir, 'domain_binary.npz')))

    # Global domain matrix
    if options.use_global:
        gold_anno = pd.read_pickle(options.global_gold_anno)
        gold_anno['gene_id'] = np.arange(gold_anno.shape[0]) # make fake gene id
        ma_glob, term_mapper = domain_to_pivot(gold_anno, term_mapper = term_mapper) # need to use the same term mapper
        save_npz(os.path.join(options.outdir, 'domain_binary_global.npz'), ma.tocsr())   
        print('Saving to {}'.format(os.path.join(options.outdir, 'domain_binary_global.npz')))
        # make sure two pivot table has the same number of domains
        assert ma.shape[1] == ma_glob.shape[1]
