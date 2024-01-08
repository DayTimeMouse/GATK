#!/usr/bin/env python
# coding: utf-8

# In[337]:


###########################################
#2023-12-03
#Jiazheng Lin
#通过hg38的CDS序列得到基因nonsense mutation后的氨基酸序列
#hg38采用gatk的/home/ljz/year1/amlwes/mutation_aa/Homo_sapiens_assembly38.fasta
#使用TBtools对Homo_sapiens_assembly38.fasta提取hg38_cds.fa
#或者gffread hg38.refGene.gtf -g Homo_sapiens_assembly38.fasta -x hg38_cds.fa
##########################################
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
import re
import logging
import sys


# In[338]:


def extract_gene_cds(trans_id, cds_file):
    cds_seq = None
    try:
        for record in SeqIO.parse(cds_file, "fasta"):
            if trans_id in record.id:
                cds_seq = str(record.seq)
                break
    except Exception as e:
        print(f"Error in extract_gene_cds(): {e}")
    return cds_seq

def get_nonsense_mut(maf, gene_id):
    #chunksize读取大文件，再用concat合并
    maf_data = pd.read_csv(maf,sep='\t',chunksize=1000)
    maf_data = pd.concat(maf_data)
    need_colunms = ['Hugo_Symbol','Start_Position','End_Position','Variant_Classification','Tumor_Sample_Barcode', 'cDNA_position', 'CDS_position', 'Protein_position', 'HGVSc', 'HGVSp', 'HGVSp_Short', 'Transcript_ID', 'RefSeq']
    gene_mut_pos = maf_data[need_colunms]
    gene_id_all = gene_mut_pos[gene_mut_pos['Hugo_Symbol'] == gene_id]
    nonsense_mut = gene_id_all[gene_id_all['Variant_Classification'] == 'Nonsense_Mutation']
    cds_pos = int(nonsense_mut['CDS_position'].iloc[0].split('/')[0])
    cds_len = int(nonsense_mut['CDS_position'].iloc[0].split('/')[1])
    HGVSc = nonsense_mut['HGVSc'].iloc[0]
    pattern = r'c\.\d+([AGCT]>[AGCT])'
    mut_mode = re.search(pattern, HGVSc).group(1)
    ref_base = mut_mode.split('>')[0]
    mut_base = mut_mode.split('>')[1]
    return nonsense_mut, ref_base, mut_base, cds_pos, cds_len

def extract_refseq_id(nonsense_mut):
    try:
        refseq_id = nonsense_mut['RefSeq'].iloc[0]
    except Exception as e:
        print(f"Error in extract_refseq_id(): {e}")
    return refseq_id

def main(hg38cds, maf, gene_id, out_path): 
    logging.basicConfig(filename=f'{out_path}/{gene_id}_nonsense_mutation_aa.log', level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s: %(message)s')
    # 1. get refseq_id from cds.fa
    # get cds_seq_true
    nonsense_mut, ref_base, mut_base, cds_pos, cds_len = get_nonsense_mut(maf,gene_id)
    refseq_id = extract_refseq_id(nonsense_mut)
    if ',' in refseq_id:
        list_refid = list()
        list_refid = refseq_id.split(',')
        for i in list_refid:
            cds_seq = extract_gene_cds(i.split('.')[0],hg38cds)
            mut_cal = cds_seq[cds_pos-1]
            cds_cal = len(cds_seq)
            if mut_cal == ref_base and cds_cal == cds_len:
                cds_seq_true = cds_seq
    else:
        cds_seq = extract_gene_cds(refseq_id.split('.')[0],hg38cds)
        mut_cal = cds_seq[cds_pos-1]
        cds_cal = len(cds_seq)
        if mut_cal == ref_base and cds_cal == cds_len:
            cds_seq_true = cds_seq
    # 2. replce the mutate base
    cds_list = list(cds_seq_true)
    cds_list[cds_pos - 1] = mut_base
    mut_seq = ''.join(cds_list)
    # 3. translate the base to aa
    aa_seq = Seq(mut_seq).translate()
    if '*' in aa_seq:
        aa_seq = aa_seq[:-1]
    with open(f'{out_path}/{gene_id}_aa_seq.fa', 'w') as fasta:
        fasta.write(f'>{gene_id}\n')
        for i in range(0, len(aa_seq), 60):
            fasta.write(str(aa_seq)[i:i+60] + '\n')
    with open(f'{out_path}/{gene_id}_mut_cds_seq.fa', 'w') as fasta:
        fasta.write(f'>{gene_id}\n')
        for i in range(0, len(mut_seq), 60):
            fasta.write(str(mut_seq)[i:i+60] + '\n')
    with open(f'{out_path}/{gene_id}_cds_seq.fa', 'w') as fasta:
        fasta.write(f'>{gene_id}\n')
        for i in range(0, len(cds_seq_true), 60):
            fasta.write(str(cds_seq_true)[i:i+60] + '\n')
    # 4. write raw aa sequence
    raw_aa_seq = Seq(cds_seq_true).translate()
    if '*' in raw_aa_seq:
        raw_aa_seq = raw_aa_seq[:-1]
    with open(f'{out_path}/{gene_id}_raw_aa_seq.fa', 'w') as fasta:
        fasta.write(f'>{gene_id}\n')
        for i in range(0, len(raw_aa_seq), 60):
            fasta.write(str(raw_aa_seq)[i:i+60] + '\n')
if __name__ == "__main__":
    # hg38cds ='/home/ljz/year1/amlwes/mutation_aa/hg38_cds.fa'
    # maf ='/home/ljz/year1/amlwes/mutation_aa/merged_vep.maf'
    # gene_id = 'RUNX1'
    # out_path = './'
    hg38cds = sys.argv[1]
    maf = sys.argv[2]
    gene_id = sys.argv[3]
    out_path = sys.argv[4]
    main(hg38cds, maf, gene_id, out_path)


