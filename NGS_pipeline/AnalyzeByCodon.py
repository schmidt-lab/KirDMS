# Preprocess Read Files
import os
import argparse
import time
import multiprocessing
import pandas as pd
import tqdm
import math
from itertools import islice, chain
import re
import time

parser = argparse.ArgumentParser(description="Analyze reads for specific mutations")
parser.add_argument('-r1', required=True, help='Read 1 File')
parser.add_argument('-r2', required=True, help='Read 2 File')
parser.add_argument('-adapt3', default='CCAGAAGAGTATAGTCGATTCGAAAT', help='Adapter sequence on 5 prime end')
parser.add_argument('-adapt5', default='GAAGACATTCCCCGGGAA', help='Adapter sequence on 3 prime end')
parser.add_argument('-mutFile', help='File containing mutations found in programmed oligos. List of 19 mutagenic codons for each codon')
parser.add_argument('-numMutations', default=1, type=int, help='Number of mutations in library')
parser.add_argument('-outFile', default='empty')
parser.add_argument('-wtFile', required=True, help='fasta sequence of wild type gene')
args = parser.parse_args()
print("\n-------------------------------------------------")
tmpname = args.r1.split('/')[-1].split('_')[0:3]
tmpname = '-'.join(tmpname)
NGSfile = 'results/' + tmpname + '_mapped.sam'
print("Processing " + tmpname)
print("-------------------------------------------------")
# Combine Reads
if not os.path.exists(NGSfile):
    print(r'filtering reads')
    ret = os.system(r'bbduk.sh in1=' + args.r1 + ' in2=' + args.r2 + ' out1=results/R1_trim_1.fastq out2=results/R2_trim_1.fastq ref=/bin/bbmap/resources/adapters.fa ktrim=r mink=8') #trim illumina adapters
    ret = os.system(r'bbduk.sh in1=results/R1_trim_1.fastq in2=results/R2_trim_1.fastq out1=results/R1_trim_2.fastq out2=results/R2_trim_2.fastq literal=CAGATCCGGCCACC k=8 ktrim=l') # trim kozak/promoter sequence
    ret = os.system(r'bbmerge.sh in1=results/R1_trim_2.fastq in2=results/R2_trim_2.fastq out1=results/R1_trim.fastq out2=results/R2_trim.fastq ecco mix')
    print(r'Aligning Reads - bbmap.sh ref=' + args.wtFile + ' in1=results/R1_trim.fastq in2=results/R2_trim.fastq outm=' + NGSfile + ' local')
    ret = os.system(r'bbmap.sh ref=' + args.wtFile + ' in1=results/R1_trim.fastq in2=results/R2_trim.fastq outm=' + NGSfile + ' local')
    assert ret == 0, 'aligning reads failed'
    print("-------------------------")
else:
    print("Aligned File Exists. Continuing with that file")
    print("-------------------------")

# Analyze Mutational Sequences by Codon
with open(args.wtFile,'r') as reader:
    wt = reader.read().upper()

quality = 20

samFlags = [[]]
def loop_batch(loop_read):
    # Analyze Codons
    totalcounts = []
    while True:  # loop through each read from ngs
        try:
            tmp_r1 = next(loop_read).split('\t')
        except:
            break
        if '@' in tmp_r1[0]:
            continue
        tmp_r2 = next(loop_read).split('\t')
        if tmp_r1[0] == tmp_r2[0]: # paired reads have the same name
            # identify forward strand
            if bin(int(tmp_r1[1]))[::-1][5]:
                r1 = tmp_r2
                r2 = tmp_r1
            else:
                r1 = tmp_r1
                r2 = tmp_r2
            tmpcounts = [r1[0]] # initialize name
            for r3 in [r1,r2]:
                if (bin(int(r3[1]))[::-1][1] or bin(int(r3[1]))[::-1][3]) and int(r3[4]) > quality: # if forward read is aligned and higher than quality threshold
                    # trim the read
                    if 'S' in r3[5]:
                        cut = r3[5].split('S')[0]
                        if all(map(str.isdigit,cut)): #trim beginning
                            r3[9] = r3[9][int(cut):]
                            r3[10] = r3[10][int(cut):]
                        else: #trim end
                            r3[9] = r3[9][:-int(re.findall(r'\d+',cut)[-1])]
                            r3[10] = r3[10][:-int(re.findall(r'\d+',cut)[-1])]
                    # find first codon
                    codon_remainder = int(3-(int(r3[3])-1)%3)
                    ref = int(r3[3]) + codon_remainder -1
                    read = codon_remainder
                    while ref+3<len(wt)-1 and read+3<len(r3[9]):  # cycle through each codon (amino acid in pdz sequence):
                        if r3[9][read:read+3].upper() != wt[ref:ref+3].upper() and 'N' not in r3[9]:  # look for codon
                            if all([ord(r3[10][read+b])-33 > 16 for b in range(0,3)]):
                                tmpcounts.append(str(int((ref)/3+1))+'_'+r3[9][read:read+3])  # record mutation name if found
                        ref += 3
                        read += 3
            if len(tmpcounts) > 1: # only record if mutations found
                # append data
                totalcounts.append(tmpcounts)
        else:
            raise ValueError('Pairs were not found. ' + tmp_r1[0])
    return totalcounts

def batch_iterator(iterable, batch_size):
    """Returns iterator of length batch_size."""
    iterator = iter(iterable)
    for first in iterator:
        yield chain([first],islice(iterator,batch_size-1))


batchsize = 10000
num_lines = sum(1 for line in open(NGSfile))
print("Processed Files have %s reads" % str((num_lines-3)/2))

# Parallel Processing
#if __name__ == '__main__':
#    with multiprocessing.Pool(8) as pool:
#        record_iter = open(NGSfile, 'r')  # should be fast
#        batches = batch_iterator(record_iter, batchsize)
#        pool_counts = []
#        for batch in batches:
#            pool_counts.append(loop_batch(batch))
#            pool_counts = list(tqdm.tqdm(pool.imap(loop_batch,batches),total=math.ceil(num_lines / batchsize)))

start = time.time()
record_iter = open(NGSfile, 'r')
pool_counts = loop_batch(record_iter)
end = time.time()
print('total_time ' + str(int(end-start)) + ' seconds')
#print(pool_counts)

codon_code = {
            'Cys': ['TGT', 'TGC'],
            'Asp': ['GAT', 'GAC'],
            'Ser': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'],
            'Gln': ['CAA', 'CAG'],
            'Met': ['ATG'],
            'Asn': ['AAC', 'AAT'],
            'Pro': ['CCT', 'CCG', 'CCA', 'CCC'],
            'Lys': ['AAG', 'AAA'],
            'Stop': ['TAG', 'TGA', 'TAA'],
            'Thr': ['ACC', 'ACA', 'ACG', 'ACT'],
            'Phe': ['TTT', 'TTC'],
            'Ala': ['GCA', 'GCC', 'GCG', 'GCT'],
            'Gly': ['GGT', 'GGG', 'GGA', 'GGC'],
            'Ile': ['ATC', 'ATA', 'ATT'],
            'Leu': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'],
            'His': ['CAT', 'CAC'],
            'Arg': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'],
            'Trp': ['TGG'],
            'Val': ['GTA', 'GTC', 'GTG', 'GTT'],
            'Glu': ['GAG', 'GAA'],
            'Tyr': ['TAT', 'TAC']}

all_codons = ['TGT', 'TGC', 'GAT', 'GAC', 'TCT', 'TCG', 'TCA',
'TCC', 'AGC', 'AGT', 'CAA', 'CAG', 'ATG', 'AAC', 'AAT', 'CCT',
'CCG', 'CCA', 'CCC', 'AAG', 'AAA', 'TAG', 'TGA', 'TAA', 'ACC',
'ACA', 'ACG', 'ACT', 'TTT', 'TTC', 'GCA', 'GCC', 'GCG', 'GCT',
'GGT', 'GGG', 'GGA', 'GGC', 'ATC', 'ATA', 'ATT', 'TTA', 'TTG',
'CTC', 'CTT', 'CTG', 'CTA', 'CAT', 'CAC', 'CGA', 'CGC', 'CGG',
'CGT', 'AGG', 'AGA', 'TGG', 'GTA', 'GTC', 'GTG', 'GTT', 'GAG',
'GAA', 'TAT', 'TAC']

all_file = open(NGSfile.split('_')[0]+'_all_mutations.csv','w')

#pool_counts = []
#with open('all_mutations.csv','r') as reader:
#    for read in reader:
#        tmp_3 = read.rstrip().split(',')
#        pool_counts.append(tmp_3)

counts = pd.DataFrame(0, index=all_codons, columns = list(range(1,int(len(wt)/3))))
#for d in pool_counts:  # Combine output of parallel pool
for tmp in pool_counts:
    all_file.write(",".join(tmp)+'\n')
    # only record reads with 1 mutation
    # remove duplicates
    tmp_2 = []
    for x in tmp:
        if x not in tmp_2:
            tmp_2.append(x)
    if len(tmp_2)==2:
        col,row = tmp_2[1].split('_')
        counts.loc[row,int(col)] += 1

print('Reads with single Mutations:' + str(sum(counts.sum())))
print(str(sum(counts.sum())* 100 / num_lines/2)[:4] + '% Positive Reads from Filtered Reads')
print("Analysis Done. Writing Files")

# Write to csv
counts.to_csv(NGSfile.split('_')[0]+'_single_mutations.csv')
