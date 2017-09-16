from gtf import *

def create_gtf_of_longest_txpt_per_gene(
        in_filename='', out_filename='/opt/lib/gtf_one_txtp_per_gene.txt'):
    g = gtf(in_filename, '/opt/lib/chr_sizes.txt')
    print "loaded %s" % in_filename
    g.subset_to_longest_txpt()
    g.to_csv(out_filename)


if __name__ == '__main__':
    create_gtf_of_longest_txpt_per_gene(
        in_filename='/opt/lib/gtf_with_names_column.txt',
        out_filename='/opt/lib/gtf_one_txtp_per_gene.txt')
