import pandas as pd
import os

def reverse_look(codon_code,codon):
    for key, values in codon_code.items():
        for value in values:
            if value == codon:
                AA = key
    return AA

os.mkdir('results/hgvs')

# Analyze Mutational Sequences by Codon
with open('kir_dna.fa','r') as reader:
    wt = reader.read().upper()

usage = pd.read_csv('mutations_trimmed_chart.csv',header=0, index_col=0)

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


for file in os.listdir('results/'):
    if file.endswith('single_mutations.csv'):
        counts_tmp = pd.read_csv('results/'+file,header=0, index_col=0)
        counts = counts_tmp*usage
        with open('results/hgvs/'+file,'w') as hgvs_file:
            hgvs_file.write('variant\tcount\n')
            wt_c = 0
            for mut_codon, rows in counts.iterrows():
                for index2,col in enumerate(rows):
                    if col>0:
                        wt_codon = wt[index2*3:index2*3+3]
                        wt_aa = reverse_look(codon_code, wt_codon)
                        mut_aa = reverse_look(codon_code, mut_codon)
                        if wt_aa == mut_aa:
                            #n_str += ' (p.=), '
                            wt_c += col
                        #else:
                        n_str = ''
                        for i in range(0,3):
                            if wt_codon[i] != mut_codon[i]:
                                n_str += 'c.' + str(index2*3+1+i) + wt_codon[i] + '>' + mut_codon[i]
                                n_str += ' (p.' + wt_aa + str(index2+1) + mut_aa +'), '
                        hgvs_file.write(n_str[:-2]+'\t'+str(col)+'\n')
            hgvs_file.write('_wt\t'+str(wt_c))
