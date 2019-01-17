import pandas
import numpy
import re
from collections import defaultdict

def open_clustering_file(path):
    '''
    Reads a clustering file represented in CSV format.
    Assumes the file has a header, so skips first line.
    '''
    clustering = pandas.read_csv(path,header=0,index_col=0,sep='\t') 
    return clustering

def open_distance_file(filename):
    '''
    Reads distance matrix file represented in CSV format.
    Returns distance matrix as a pandas DataFrame matrix.
    '''
    distance = pandas.read_csv(filename,header=0,index_col=0,sep='\t')
    assert( distance.values.shape[0] == distance.values.shape[1] ),\
        "Distance matrix isn't square."
    return distance

def read_mlst_calls(calls_paths):
    '''
    Read MentaLiST MLST calls.
    '''
    # If calls_paths is a string, we assume that it's a path to a file containing the calls paths,
    # and we read the paths from the file.
    if isinstance(calls_paths,str):
        calls_paths_path = calls_paths
        calls_paths = []
        with open(calls_paths_path,'r') as calls_paths_file:
            for line in calls_paths_file:
                calls_paths.append(line.rstrip())
    calls = {}
    #print(calls_paths)
    # Go through the calls file for each sample, and record the MLST calls
    for calls_path in calls_paths:
        calls_path = calls_path.rstrip().split('=')[0]
        with open(calls_path,'r') as calls_file:
            # Skip the header
            calls_file.readline()
            for line in calls_file:
                columns = line.rstrip().split('\t')
                sample = columns[0]
                calls[sample] = numpy.array(columns[1:],dtype=numpy.string_)

    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of MLST calls."

    return calls
    
def read_snp_calls(calls_paths):
    '''
    Read SNP calls from a text file with list of Snippy tsv or snippy-core output.
    '''
    calls = {}
    # If calls_paths is a list of paths or a text file containing the paths, do this. 
    if isinstance(calls_paths,list) or calls_paths.endswith("txt"):
        # collect the paths to the calls into a list
        if isinstance(calls_paths,str):
            calls_paths_path = calls_paths
            calls_paths = []
            with open(calls_paths_path,'r') as calls_paths_file:
                for line in calls_paths_file:
                    calls_path = line.rstrip().split('=')[0]
                    calls_paths.append(calls_path)
        # now get the SNP calls for each of the samples
        ref = {}
        sample = {}
        pos_count = {}
        #initialize dictionaries
        pos_count = defaultdict(int)
        sample = defaultdict(dict)
        for calls_path in calls_paths:
            with open(calls_path,'r') as call_file:
                sample_name=call_file.readline().rstrip()
                sample[sample_name][-1] = -1
                for line in call_file:
                    columns = line.rstrip().split('\t')
                    #ensure there is one entry in a sample with no snps calls 
                    sample[sample_name][-1] = -1
                    if columns[2] == "snp":
                        chrom_pos = columns[0]+"_"+columns[1]
                        pos_count[chrom_pos] += 1
                        sample[sample_name][chrom_pos] = columns[4]
                        ref[chrom_pos] = columns[3]
        pos_count = {k:v for (k,v) in pos_count.items() if (v > 1 and v != len(sample.keys()) ) }
        for sample_name in sample.keys():
            snps = []
            sample_name_keys = sample[sample_name].keys()
            for pos in pos_count.keys():
                if pos in sample_name_keys:
                    snps.append(sample[sample_name][pos])
                else:
                    snps.append(ref[pos])
            calls[sample_name] = numpy.array(snps ,dtype="S1")    

    # Otherwise, it must be snippy-core output 
    elif calls_paths.endswith("tab"):
        snps_union = pandas.read_csv(calls_path, sep='\t')
        snps_union = snps_union.drop(['CHR','POS','Reference','LOCUS_TAG' ,'GENE', 'PRODUCT', 'EFFECT'], axis=1)
        for column in snps_union:
            if column == "CHROM" or column == "POS" or column == "TYPE" or column == "REF":
                continue
            calls[column] = numpy.array(snps_union[column], dtype="S20")
    return calls

def read_snp_calls_with_bed(calls_path, bed_path):
    '''
    Read SNP calls from a text file with list of Snippy tsv or snippy-core output.
    '''
    calls = {}
    bed_filter = {}
#bed_filter = defaultdict(int)
    #bed_path="/home/klw17/filter_test.bed"
    with open(bed_path,'r') as bed_file:
        columns=bed_file.readline().rstrip().split('\t')
        if columns[0].lower() != "chrom":
            for i in range(int(columns[1]),int(columns[2])):
                filter_pos = columns[0]+"_"+str(i)
                bed_filter[filter_pos] = 1        
        for line in bed_file:
            columns = line.rstrip().split('\t')
            for i in range(int(columns[1]),int(columns[2])):
                filter_pos = columns[0]+"_"+str(i)
                bed_filter[filter_pos] = 1        
            
    if calls_path.endswith("txt"):
        ref = {}
        sample = {}
        pos_count = {}
        #initialize dictionaries
        pos_count = defaultdict(int)
        sample = defaultdict(dict)
        with open(calls_path,'r') as calls_file:
            for line in calls_file:
                # We assume that each line of the calls_file is in the form
                # call_path=sample_name
                call_path = line.rstrip().split('=')[0]
                '''
                call_path = line.rstrip()
                sample_name = re.sub('/[A-Za-z.]*.tab', '', call_path)
                sample_name = re.sub('.*/', '', sample_name)        
                '''
                with open(call_path,'r') as call_file:
                    sample_name=call_file.readline().rstrip()
                    sample[sample_name][-1] = -1
                    for line in call_file:
                        columns = line.rstrip().split('\t')
                        #ensure there is one entry in a sample with no snps calls 
                        sample[sample_name][-1] = -1
                        if columns[2] == "snp":
                            chrom_pos = columns[0]+"_"+columns[1]
                            pos_count[chrom_pos] += 1
                            sample[sample_name][chrom_pos] = columns[4]
                            ref[chrom_pos] = columns[3]
        pos_count = {k:v for (k,v) in pos_count.items() if (v > 1 and v != len(sample.keys()) and k not in bed_filter ) }
        for sample_name in sample.keys():
            snps = []
            sample_name_keys = sample[sample_name].keys()
            for pos in pos_count.keys():
                snps.append(sample[sample_name][pos]) if pos in sample_name_keys else snps.append(ref[pos])
            calls[sample_name] = numpy.array(snps ,dtype="S1")    
    '''
        with open(calls_path,'r') as calls_file:
            snps_union = pandas.DataFrame()
            #print(snps_union.shape)
            list_of_df = []
            for line in calls_file:
                call_path = line.rstrip()
                sample = re.sub('/[A-Za-z.]*tab', '', call_path)
                sample = re.sub('.*/', '', sample)
                df=pandas.read_table(call_path)
                snps = df[['CHROM', 'POS', 'TYPE', 'REF', 'ALT']]
                snps = snps.rename(columns={'ALT': sample})
                #snps = snps[snps.TYPE == "snp"]
                snps['POS']=snps['POS'].astype(int)
                list_of_df.append(snps)
            while len(list_of_df) != 1 :
                new_list_of_df = []
                for i in range(0,len(list_of_df),2):
                    if i == len(list_of_df) - 1 and len(list_of_df) % 2 == 1:
                        new_list_of_df[len(new_list_of_df) - 1] = pandas.merge(list_of_df[i], new_list_of_df[len(new_list_of_df) - 1], how='outer', on=['CHROM', 'POS', 'TYPE','REF'])
                        continue
                    new_list_of_df.append(pandas.merge(list_of_df[i], list_of_df[i+1], how='outer', on=['CHROM', 'POS', 'TYPE','REF']))
                list_of_df = new_list_of_df.copy()
            snps_union = list_of_df[0]
            for i in range(4,snps_union.shape[1]):
                snps_union.iloc[:,i].fillna(snps_union.REF, inplace=True)
    '''
    # Using snippy-core output as input
    if calls_path.endswith("tab"):
        snps_union = pandas.read_csv(calls_path, sep='\t')
        snps_union = snps_union.drop(['CHR','POS','Reference','LOCUS_TAG' ,'GENE', 'PRODUCT', 'EFFECT'], axis=1)
        for column in snps_union:
            if column == "CHROM" or column == "POS" or column == "TYPE" or column == "REF":
                continue
            calls[column] = numpy.array(snps_union[column], dtype="S20")
    #modify this later
    #assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
    #    "Samples do not have the same number of MLST calls."
    return calls


def read_cnv_calls(calls_path):
    '''
    Read PRINCE CNV calls.
    @param calls_path: a string to the PRINCE calls file
    @rvalue calls: a dictionary whose keys are sample names, and values are CNV calls represented by 
                   a numpy array
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip().split('=')[0]
            with open(call_path,'r') as call_file:
                # Skip the header
                call_file.readline()
                for line in call_file:
                    columns = line.rstrip().split('\t')
                    sample = columns[0]
                    calls[sample] = numpy.array(columns[1:],dtype=float)
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of CNV calls."
    return calls

def read_spotype_calls(calls_paths):
    '''
    Read SpoTyping calls.
    @param calls_paths: a list containing paths to SpoTyping calls files, OR a path to a file
                        containing the paths to the calls 
    @rvalue calls: a dictionary whose keys are sample names, and values are spoligotypes represented
                   by a numpy array 
    '''
    # if calls_paths is a string, we assume its a file containing the paths to the SpoTyping calls
    if isinstance(calls_paths,str):
        calls_paths_path = calls_path
        calls_paths = []
        with open(calls_paths_path,'r') as calls_paths_file:
            for line in calls_paths_file:
                calls_path = line.rstrip().split('=')[0]
                calls_paths.append(calls_path)
    calls = {}
    for call_path in calls_paths:
        with open(call_path,'r') as call_file:
            # Skip the header
            call_file.readline()
            for line in call_file:
                columns = line.rstrip().split('\t')
                sample = columns[0]
                calls[sample] = numpy.array(columns[1:],dtype=float)
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of SpoType calls."
    return calls

def output_clustering(clustering,output_path):
    '''
    Writes a clustering to file in TSV format.
    '''
    clustering.to_csv(output_path,index=True,sep='\t')

def write_distance_matrix(distance_matrix,output_path):
    '''
    Writes a distance matrix to file in TSV format.
    '''
    distance_matrix.to_csv(output_path,sep='\t')
