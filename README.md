weighted_mutual
---------------

To use, simply do this in command line:
    >>> python weighted_mutual/main.py -h

Reference
    Weighted mutual information analysis substantially improves domain-based functional network models
    https://academic.oup.com/bioinformatics/article/32/18/2824/1744011


Options:
  -h, --help            show this help message and exit
  -p POOL, --pool=POOL  number of multiprocessing threads
  -g GOLD_ANNO, --goldanno=GOLD_ANNO
                        path to gold_anno pickle dataframe
  --global=GLOBAL_GOLD_ANNO
                        gold_anno of the whole pan-genome, or proteome space
  -o OUTDIR, --outdir=OUTDIR
                        output directory
  -c N_CHUNK, --chunk=N_CHUNK
                        select chunk size to
  --use_global          generate global pivot table or not