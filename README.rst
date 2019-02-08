weighted_mutual
---------------

To use, simply do this in command line:
    >>> python main.py --input <filename> --output <filename> --threads <default 1>

Reference
    Weighted mutual information analysis substantially improves domain-based functional network models
    https://academic.oup.com/bioinformatics/article/32/18/2824/1744011


Input file format
    .csv with columns = protein domain, rows = proteins
    see /test_file/small_test as example

Output file format
    .csv with columns: gene_one, gene_two, weighted_mutual_info
