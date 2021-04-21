from optparse import OptionParser
import os
def option_parser():
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
    parser.add_option("-g", "--goldanno",dest="gold_anno",
                  help="path to gold_anno pickle dataframe")
    parser.add_option("--global",dest="global_gold_anno",default = None,
                  help="gold_anno of the whole pan-genome, or proteome space")
    parser.add_option("-o", "--outdir",dest="outdir", default = '/tmp', type = "string",
                  help="output directory")
    parser.add_option("-c", "--chunk",dest="n_chunk", default = 1000000, type = "int",
                  help="select chunk size to")
    parser.add_option("--use_global", dest="use_global", action='store_true',
                  help="generate global pivot table or not")
    
    
    (options, args) = parser.parse_args()   

    return options

if __name__=='__main__':
    options = option_parser()

    dir_path = os.path.dirname(os.path.realpath(__file__))
    script_path = os.path.join(dir_path, 'generate_pivot_table.py')
    if options.use_global:
        os.system('python {} --goldanno {} --global {}--use_global -o {}'.format(script_path,
                                                                            options.gold_anno,
                                                                            options.global_gold_anno,
                                                                            options.outdir))
        global_pivot = os.path.join(options.outdir, 'domain_binary_global.npz')
    else:
        os.system('python {} --goldanno {} -o {}'.format(script_path,
                                                        options.gold_anno,
                                                        options.outdir))
    
    pivot_path = os.path.join(options.outdir, 'domain_binary.npz')

    script_path = os.path.join(dir_path, 'calc_score_async.py')
    if options.use_global: 
        os.system('python {} -p {} --input {} --global {} -o {}'.format(script_path,
                                                                        options.pool,
                                                                        pivot_path,
                                                                        global_pivot,
                                                                        options.outdir))
    else: 
        os.system('python {} -p {} --input {} -o {}'.format(script_path,
                                                            options.pool,
                                                            pivot_path,
                                                            options.outdir))



    
    