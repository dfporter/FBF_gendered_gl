import collections, glob, time, re, os, sys, HTSeq, argparse
import configparser



# Following:
# http://www.bioconductor.org/help/workflows/rnaseqGene/
#samfile = sys.argv[1]

# HTSeq segment from:
# http://www-huber.embl.de/users/anders/HTSeq/doc/counting.html#counting


def get_gtf(gtf_path=None):
    
    print("Loading gtf {0}...".format(gtf_path))
    
    gtf_file = HTSeq.GFF_Reader(gtf_path)
    features = HTSeq.GenomicArrayOfSets('auto', stranded=True)
    
    for feature in gtf_file:
        if feature.type == 'exon':
            features[feature.iv] += feature.attr['gene_name']
            
    return features


def read_bed(bed_file):
    switch_strand = {'-': '+', '+': '-'}
    
    with open(bed_file, 'r') as f:
        for li in f:
            s = li.rstrip('\n').split('\t')
            yield (HTSeq.GenomicInterval(s[0], int(s[1]), int(s[2]), s[5]), 1)


def map_bed(bed_filename, features):
    
    start_time = time.clock()
    
#    counts = collections.collections.Counter()
    counts = collections.defaultdict(int)
    
    total_reads = 0
    for line_num, (read_iv, number_of_reads), in enumerate(read_bed(bed_filename)):
        
        if line_num and (not line_num % 100000):
            print(line_num)
            elapsed = time.clock() - start_time
            per_ten_million = 1e6 * elapsed/max([1, float(line_num)])
            print("Time elapsed: %f. Seconds per million reads %s." % (
                elapsed, per_ten_million))
#        if not almnt.aligned:
#            counts['_unmapped'] += 1
#            continue

        gene_ids = set()
        for iv, val in features[read_iv].steps():
            gene_ids |= val
            
        if len(gene_ids) == 1:
            gene_id = list(gene_ids)[0]
            counts[gene_id] += number_of_reads
        elif len(gene_ids) == 0:
            counts['_no_feature'] += number_of_reads            
        else:       
            counts['_ambiguous'] += number_of_reads
    
        total_reads += number_of_reads

    print('sum %s' % bed_filename)
    print(sum([counts[x] for x in counts]))
    
    try:
        print('no feature %i ambiguous %i' % (counts['_no_feature'], counts['_ambiguous']))
    except:
        print('either no ambiguous or no "no_feature" reads.')
        
    output_filename = 'counts/' + os.path.basename(bed_filename).partition('.bed')[0]
    output_filename += '_counts.txt'
    print(output_filename)
    
    if not os.path.exists('counts/'):
        os.system('mkdir counts')
        
#    if os.path.exists(output_filename):
#        continue

    with open(output_filename, 'w') as f:
        for gene_id in counts:
            f.write("{gid}\t{cts}\n".format(gid=gene_id, cts=counts[gene_id]))
            
    os.system('cp {a} counts/'.format(a=output_filename))

    # Write the total read count in the bed file.
    if not os.path.exists('total_read_numbers.txt'):
        with open('total_read_numbers.txt', 'w') as f:
            f.write('{}\t{}\n'.format(output_filename, total_reads))
    else:
        with open('total_read_numbers.txt', 'a') as f:
            f.write('{}\t{}\n'.format(output_filename, total_reads))

def fill_in_gaps(folder_name, gtf_filename):
    
    wb_ids = gene_names_in_gtf(gtf_filename)
    
    counts = {}
    known_genes = set()
    
    for fname in glob.glob(folder_name + '/*txt'):
        print(fname)
    
        if re.match('combined_.*', os.path.basename(fname)) is not None:
            continue
            
        counts[fname] = collections.defaultdict(int)
        
        with open(fname, 'r') as f:
            
            for li in f:
                s = li.rstrip('\n').split('\t')
                known_genes.add(s[0])
                counts[fname][s[0]] += int(s[1])
                
    # for fname in counts:
    #     for gene in list(wb_ids - set(counts[fname].keys())):
    #         counts[fname][gene] = 0
    
    gene_output_order = list(known_genes)#counts[counts.keys()[0]].keys()
    
    for fname in counts:
        
        for gene in known_genes:
            if gene not in counts[fname]:
                counts[fname][gene] = 0
                
        with open(fname, 'w') as f:
            for gene in gene_output_order:
                f.write('{i}\t{n}\n'.format(i=gene, n=counts[fname][gene]))


def create_combined_file_of_counts(
        counts_folder, output_filename='combined_counts.txt'):

    _all = collections.defaultdict(dict)
    output_order = []
    
    for fname in glob.glob(counts_folder + '/*.txt'):
        
        f = open(fname, 'r')
        
        for line in f:
            s = line.rstrip('\n').split('\t')
            _all[s[0]][os.path.basename(fname)] = s[1]
            
        output_order.append(os.path.basename(fname))
        
    li = 'gene\t' + '\t'.join(output_order) + '\n'
    
    for gene in _all:
        
        li += '%s\t' % gene
        
        for k in output_order:
            if k in _all[gene]:
                li += '%s\t' % _all[gene][k]
            else:
                li += '0\t'
                
        li += '\n'
        
    with open(output_filename, 'w') as f:
        f.write(li)
        # s[fname] = dict([line.rstrip('\n').split('\t') for line in f])
        # genes.add(set(s.keys()))

def gene_names_in_gtf(gtf_fname):
    wb_ids = set()
    
    with open(gtf_fname, 'r') as f:
        for li in f:
            m = re.match('.*gene_id "([^"]+)";', li)
            if m is not None:
                wb_ids.add(m.groups()[0])
                
    return wb_ids


def wb_to_public_name(gtf_fname):
    public_names = {}
    
    with open(gtf_fname, 'r') as f:
        
        for li in f:
            
            m = re.match('.*gene_id "([^"]+)";', li)
            
            if m is not None:
                name_m = re.match('.*gene_name "([^"]+)";', li)
                if name_m is not None:
                    public_names[m.groups()[0]] = name_m.groups()[0]
                    
    return public_names


def set_of_all_gene_names_in_folder(folder_name):
    gene_names = set()
    
    for fname in glob.glob(folder_name + '/*'):
        with open(fname, 'r') as f:
            for li in f:
                gene_names.add(li.rstrip('\n').split('\t')[0])
                
    return gene_names


def get_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b', '--bed', help="""Folder of bed files.""")
    parser.add_argument('-c', '--config', help="config.ini file, only needed for the path to a gtf file.",
                        default='auto.ini')
    parser.add_argument('-g', '--gtf',
        help='Optional - path to a gtf file (alternative to -c option).',
        default='/opt/lib/Caenorhabditis_elegans.WBcel235.78.noheader.gtf')
    
    args = parser.parse_args()
    
    return args


def run(bed_folder=None, gtf_filename='', lib=None):
    
        
    if bed_folder is None:
        bed_folder = lib['read_beds']
            
    features = get_gtf(gtf_path=gtf_filename)#'./genomes/Saccharomyces_cerevisiae.EF4.70.gtf')
    
    if os.path.exists('total_read_numbers.txt'):
        os.system('rm total_read_numbers.txt')
        
    for bed_filename in glob.glob(bed_folder + '/*.bed'):
        print("Assigning reads from %s..." % bed_filename)
        
        #if not re.search('control_n2.bed', bed_filename):
        #    continue
            
        map_bed(bed_filename, features)
        
    fill_in_gaps('counts', gtf_filename)
    create_combined_file_of_counts(
        'counts', output_filename='combined_counts.txt')

    
def config(
    filepath="/groups/Kimble/Aman Prasad/redo_fbf/analysis/src/config.ini"):
    Config = configparser.ConfigParser()
    Config.read(filepath)
    return ConfigSectionMap('library', Config)


def ConfigSectionMap(section, Config):
    settings = {}
    options = Config.options(section)
    
    for option in options:
        
        try:
            settings[option] = Config.get(section, option)
                
        except:
            print(("Could not parse option %s." % option))
            settings[option] = None
            
    return settings


if __name__ == '__main__':
    """Create a combined_counts.txt file of the reads/gene for a 
    folder of .bed files. Requires onaly one value from lib, lib['gtf_raw'].
    """
    
    args = get_args()
    
    try:
        lib = config(args.config)
        gtf_filename = lib['gtf_raw']
    except:
        gtf_filename = args.gtf
        
    run(bed_folder=args.bed, gtf_filename=gtf_filename, lib=lib)
