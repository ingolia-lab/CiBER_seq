datadir = "/mnt/ingolialab/rmuller1/CiBER_seq_package/all_raw_fasta_gz/bc_validation/"

from Bio import SeqIO
class FASTQtoBCcounts():
  #BC1 ACTGC
  #BC2 CTGAG
  #BC3 GACTA
  #BC4 TGACT
  #BC5 GTCAC
  #BC6 AAGCT
  def __init__(self,fastq):
    bc1 = 0
    bc2 = 0
    bc3 = 0
    bc4 = 0
    bc5 = 0
    bc6 = 0
    n = 0
    self.fastq = fastq
    BC1 = []
    BC2 = []
    BC3 = []
    BC4 = []
    BC5 = []
    BC6 = []
    for entries in SeqIO.parse(self.fastq,'fastq'):
      n = n + 1
      seq = str(entries.seq)
      if seq.find("ACTGC") == 10:
        BC1.append(seq)
        bc1 = bc1 + 1
      elif seq.find("CTGAG") == 10:
        BC2.append(seq)
        bc2 = bc2 + 1
      elif seq.find("GACTA") == 10:
        BC3.append(seq)
        bc3 = bc3 + 1
      elif seq.find("TGACT") == 10:
        BC4.append(seq)
        bc4 = bc4 + 1
      elif seq.find("GTCAC") == 10:
        BC5.append(seq)
        bc5 = bc5 + 1
      elif seq.find("AAGCT") == 10:
        BC6.append(seq)
        bc6 = bc6 + 1
    print('BC1%: ',float(bc1)/float(n)*100.0)
    print('BC2%: ',float(bc2)/float(n)*100.0)
    print('BC3%: ',float(bc3)/float(n)*100.0)
    print('BC4%: ',float(bc4)/float(n)*100.0)
    print('BC5%: ',float(bc5)/float(n)*100.0)
    print('BC6%: ',float(bc6)/float(n)*100.0)
    self.BC1 = self.get_count_dict(BC1)
    self.BC2 = self.get_count_dict(BC2)
    self.BC3 = self.get_count_dict(BC3)
    self.BC4 = self.get_count_dict(BC4)
    self.BC5 = self.get_count_dict(BC5)
    self.BC6 = self.get_count_dict(BC6)
  def get_count_dict(self, arr):
    count_dict = {}
    for each in arr:
      if each in count_dict:
        new_dict_val = count_dict[each] + 1
        count_dict[each] = new_dict_val
      else:
        count_dict[each] = 1
    return count_dict

def compare_rep_dicts(dict1,dict2,outfile):
  with open(outfile,'w') as outtab:
    # first, find all seqs that dict1 and dict2 share, if they share then remove from dict2, 
    # if a dict1 key isn't in dict2 then output a 0 for dict2. 
    for each in dict1:
      d1val = dict1[each]
      if each in dict2:
        d2val = dict2.pop(each)
        # Will catch all reads which appear in both replicates
        outtab.write(str(each) +'\t'+ str(d1val) +'\t'+ str(d2val) +'\n')
      else:
        # If not in replicate two, output a 0 value for that replicate two
        outtab.write(str(each) +'\t'+ str(d1val) +'\t'+ str(0) +'\n')
    
    # Now, for what is left of dict2, go ahead and confirm it's not in dict1, and write out
    # the counts. Assign 0 for dict1 if truly not in dict1. 
    for each in dict2:
      d2val = dict2[each]
      if each in dict1:
        d1val = dict1.pop(each)
        # This is a sanity check, all things caught here shoul have been caught earlier. 
        outtab.write(str(each) +'\t'+ str(d1val) +'\t'+ str(d2val) +'\n')
      else:
        # If not in replicate one, output a 0 value for that replicate one
        outtab.write(str(each) +'\t'+ str(0) +'\t'+ str(d2val) +'\n')
  return
  
def main():
  DNA1 = FASTQtoBCcounts(datadir + "/DNA_rep1.trimmed.fastq")
  DNA2 = FASTQtoBCcounts(datadir + "/DNA_rep2.trimmed.fastq")
  len(DNA1.BC1)
  compare_rep_dicts(DNA1.BC1,DNA2.BC1,datadir + "/DNA.BC1.tsv")
  compare_rep_dicts(DNA1.BC2,DNA2.BC2,datadir + "/DNA.BC2.tsv")
  compare_rep_dicts(DNA1.BC3,DNA2.BC3,datadir + "/DNA.BC3.tsv")
  compare_rep_dicts(DNA1.BC4,DNA2.BC4,datadir + "/DNA.BC4.tsv")
  compare_rep_dicts(DNA1.BC5,DNA2.BC5,datadir + "/DNA.BC5.tsv")
  compare_rep_dicts(DNA1.BC6,DNA2.BC6,datadir + "/DNA.BC6.tsv")
  del(DNA1)
  del(DNA2)
  RNA1 = FASTQtoBCcounts(datadir + "/RNA_rep1.trimmed.fastq")
  RNA2 = FASTQtoBCcounts(datadir + "/RNA_rep2.trimmed.fastq")
  compare_rep_dicts(RNA1.BC1,RNA2.BC1,datadir + "/RNA.BC1.tsv")
  compare_rep_dicts(RNA1.BC2,RNA2.BC2,datadir + "/RNA.BC2.tsv")
  compare_rep_dicts(RNA1.BC3,RNA2.BC3,datadir + "/RNA.BC3.tsv")
  compare_rep_dicts(RNA1.BC4,RNA2.BC4,datadir + "/RNA.BC4.tsv")
  compare_rep_dicts(RNA1.BC5,RNA2.BC5,datadir + "/RNA.BC5.tsv")
  compare_rep_dicts(RNA1.BC6,RNA2.BC6,datadir + "/RNA.BC6.tsv")
  del(RNA1)
  del(RNA2)
main()

def reconcile_dna_rna(dna_file,rna_file,out_file):
  dna_lines = open(dna_file,'r').readlines()
  rna_lines = open(rna_file,'r').readlines()
  out_obj = open(out_file,'w')
  
  #make dictionaries to related 1:1 the idx of the line# and the read on said line. 
  dna_line_idx_dict = {}
  n = 0
  for each in dna_lines:
    tmp = each.split('\t')
    dna_line_idx_dict[tmp[0]] = n
    n = n + 1
  rna_line_idx_dict = {}
  n = 0
  for each in rna_lines:
    tmp = each.split('\t')
    rna_line_idx_dict[tmp[0]] = n
    n = n + 1
    
  # Match lines of DNA and RNA files based on if the read is in the DNA file and RNA file. 
  # Concatenate counts.
  for dna_read in dna_line_idx_dict:
    dna_tab = dna_lines[dna_line_idx_dict[dna_read]].strip('\n').split('\t')
    if dna_read in rna_line_idx_dict:
      rna_tab = rna_lines[rna_line_idx_dict.pop(dna_read)].strip('\n').split('\t')
      out_obj.write(dna_tab[0]+'\t'+dna_tab[1]+'\t'+dna_tab[2]+'\t'+rna_tab[1]+'\t'+rna_tab[2]+'\n')
    else:
      out_obj.write(dna_tab[0]+'\t'+dna_tab[1]+'\t'+dna_tab[2]+'\t'+str(0)+'\t'+str(0)+'\n')
  
  for rna_read in rna_line_idx_dict:
    rna_tab = rna_lines[rna_line_idx_dict[rna_read]].strip('\n').split('\t')
    if rna_read in dna_line_idx_dict:
      dna_tab = dna_lines[dna_line_idx_dict.pop(rna_read)].strip('\n').split('\t')
      out_obj.write(rna_tab[0]+'\t'+dna_tab[1]+'\t'+dna_tab[2]+'\t'+rna_tab[1]+'\t'+rna_tab[2]+'\n')
    else:
      out_obj.write(rna_tab[0]+'\t'+str(0)+'\t'+str(0)+'\t'+rna_tab[1]+'\t'+rna_tab[2]+'\n')
  return
  

reconcile_dna_rna(datadir + "/DNA.BC1.tsv",datadir + "/RNA.BC1.tsv",datadir + "/DNA_RNA.BC1.tsv")
reconcile_dna_rna(datadir + "/DNA.BC2.tsv",datadir + "/RNA.BC2.tsv",datadir + "/DNA_RNA.BC2.tsv")
reconcile_dna_rna(datadir + "/DNA.BC3.tsv",datadir + "/RNA.BC3.tsv",datadir + "/DNA_RNA.BC3.tsv")
reconcile_dna_rna(datadir + "/DNA.BC4.tsv",datadir + "/RNA.BC4.tsv",datadir + "/DNA_RNA.BC4.tsv")
reconcile_dna_rna(datadir + "/DNA.BC5.tsv",datadir + "/RNA.BC5.tsv",datadir + "/DNA_RNA.BC5.tsv")
reconcile_dna_rna(datadir + "/DNA.BC6.tsv",datadir + "/RNA.BC6.tsv",datadir + "/DNA_RNA.BC6.tsv")
