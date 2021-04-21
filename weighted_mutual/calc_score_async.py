from scipy.sparse import load_npz
from itertools import combinations, zip_longest
from sklearn.metrics import mutual_info_score, normalized_mutual_info_score
from weight import *
from entropy import *
import numpy as np
from multiprocessing import Pool, Manager
import multiprocessing
import time
import pandas as pd
import os


def get_weighted_mutual_info(index_1, index_2, **kwargs):
    ''' Calculate weighted mutual information for two protein index '''
    
    weight = np.asarray(kwargs['weight'])
    weighted_entropy = kwargs['weighted_entropy']

    row1 = binned[index_1, :].toarray().ravel()
    row2 = binned[index_2, :].toarray().ravel()

    H_xy = weighted_joint_entropy(row1, row2, weight)

    m = mutual_info(weighted_entropy[index_1], weighted_entropy[index_2], H_xy)
        
    return [index_1, index_2, m]

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def option_parser():
    from optparse import OptionParser
    usage = """
        THIS IS CO-INHERITANCE 1.0.0
        python cal_score_async.py --table ~/data0118/binned.pivot.npz -o ~/data0118
        """
    description = """This is a script to obtain domain weighted mutual information for a bunch of protein
        python main.py --input <filename> --output <filename> [--threads 1]

        Weighted mutual information analysis substantially improves domain-based functional network models https://academic.oup.com/bioinformatics/article/32/18/2824/1744011
        """

    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-p", "--pool", dest="pool", type = "int", default = 8,
                  help="number of multiprocessing threads")
    parser.add_option("-i", "--input",dest="pivot",
                  help="path to pivot table (binned.pivot.npz)")
    parser.add_option("-g", "--global",dest="globalpivot",default = None,
                  help="path to global pivot table (binned.pivot.npz)")
    parser.add_option("-o", "--outdir",dest="outdir", default = '/tmp', type = "string",
                  help="output directory")
    parser.add_option("-c", "--chunk",dest="n_chunk", default = 1000000, type = "int",
                  help="select chunk size to")
    
    
    (options, args) = parser.parse_args()   

    return options

if __name__=='__main__':
    
    options = option_parser()
    print('loading matrix')
    binned = load_npz(options.pivot)
    non_zero_genes = np.where(np.sum(binned, axis = 1)>0)[0]
    n_gene = len(non_zero_genes)

    print('total {} non-zero genes'.format(n_gene))
    n_inference = (n_gene-1)*n_gene/2
    print('Plan to make {} inferences'.format(n_inference))

    # listing out all the tasks
    all_index_combine = combinations(list(non_zero_genes), 2)
    
    temp_dir = os.path.join(options.outdir, 'temp')
    try:
        os.mkdir(temp_dir)
    except:
        print('exists {}'.format(temp_dir))
    
    print('calculate domiain weight and protein weighted entropy')
    
    # domain_weight should come from the global proteome
    if options.globalpivot:
        print('loading global pivot table')
        global_matrix = load_npz(options.globalpivot)
        domain_weight = weight(global_matrix)
    else:
        print('inferring domain distribution locally. If your gene set is not of diverse function, please use global pivot option.')
        domain_weight = weight(binned)
        

    weighted_entropy_per_protein = weighted_entropy(binned, domain_weight) # np.array, len(protein)

    # solve by chunk
    chunk_size = options.n_chunk

    results_computed = 0
    
    for task_id, task_chunk in enumerate(grouper(all_index_combine, chunk_size)):
        outfile=os.path.join(temp_dir, 'weighted_mutual'+'_{}.csv'.format(task_id))

        if os.path.isfile(outfile):
            # already calculated
            print('passing task {}, calucalted'.format(task_id))
        else:

            print('processing task {}'.format(task_id))

            # start a brand new list
            with Pool(options.pool) as pool: 
                
                # start a brand new list
                t1 = time.time()

                
                keywords = {'weight':domain_weight, 'weighted_entropy': weighted_entropy_per_protein}
                sol = [pool.apply_async(get_weighted_mutual_info, t, keywords) for t in task_chunk]
    
                t2 = time.time()
                print('Done computing scores, extracting results to dataframe. Total time: {} secs'.format(t2-t1))
                with open(outfile, 'w') as f:
        
                    for s in sol:
                        #r = s.get()
                        #print(r)
                        try:
                            r = s.get()
                        #    print(r)
                            f.write(','.join([str(s) for s in r])+'\n')
                            results_computed += 1
                        except Exception as e:
                        #    print(e)
                            pass
                
                t3 = time.time()
                print('Extracting takes time: {} secs'.format(t3-t2))
    
    print('Appending all networks to OUTDIR/network.csv')
    print('{} out of {} results saved'.format(results_computed, n_inference))
    os.system('cat {}/*.csv > {}'.format(temp_dir, os.path.join(options.outdir, 'network.csv')))
            

	
                
                

