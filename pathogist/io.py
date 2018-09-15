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

def read_mlst_calls(calls_path):
    '''
    Read MentaLiST MLST calls.
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip()
            with open(call_path,'r') as call_file:
                # Skip the header
                call_file.readline()
                for line in call_file:
                    columns = line.rstrip().split('\t')
                    sample = columns[0]
                    calls[sample] = numpy.array(columns[1:],dtype=numpy.string_)
    assert( len(set([len(calls[sample]) for sample in calls.keys()])) == 1 ), \
        "Samples do not have the same number of MLST calls."
    return calls
    
def read_snp_calls(calls_path):
    '''
    Read SNP calls from a text file with list of Snippy tsv or snippy-core output.
    '''
    calls = {}
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
                # sample_name=call_path
                sample_name,call_path = line.rstrip().split('=')
                '''
                call_path = line.rstrip()
                sample_name = re.sub('/[A-Za-z.]*.tab', '', call_path)
                sample_name = re.sub('.*/', '', sample_name)        
                '''
                #x=x+1
                with open(call_path,'r') as call_file:
                    #print(x)
                        # Skip the header
                    #call_file.readline()
                    #x=1
                    for line in call_file:
                        columns = line.rstrip().split('\t')
                        #samples=line[2]
                        #ensure there is one entry in a sample with no snps calls 
                        sample[sample_name][-1] = -1
                        if columns[2] == "snp":
                        #if (columns[2] == "snp" or columns[2] == "ins" or columns[2] == "complex" or columns[2] == "mnp" or columns[2] == "del"):
                            chrom_pos = columns[0]+"_"+columns[1]
                            pos_count[chrom_pos] += 1
                            sample[sample_name][chrom_pos] = columns[4]
                            ref[chrom_pos] = columns[3]
        pos_count = {k:v for (k,v) in pos_count.items() if (v > 1 and v != len(sample.keys()) ) }
        for sample_name in sample.keys():
            snps = []
            sample_name_keys = sample[sample_name].keys()
            for pos in pos_count.keys():
                #continue
                snps.append(sample[sample_name][pos]) if pos in sample_name_keys else snps.append(ref[pos])
            #print(1)
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
    if calls_path.endswith("tab"):
        snps_union = pandas.read_csv(calls_path, sep='\t')
        snps_union = snps_union.drop(['CHR','POS','Reference','LOCUS_TAG' ,'GENE', 'PRODUCT', 'EFFECT'], axis=1)
        #print(snps_union.head(30))
        #print(snps_union.shape)
        for column in snps_union:
        #print(column)
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
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip()
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

def read_spotype_calls(calls_path):
    '''
    Read SpoTyping calls.
    '''
    calls = {}
    with open(calls_path,'r') as calls_file:
        for line in calls_file:
            call_path = line.rstrip()
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
