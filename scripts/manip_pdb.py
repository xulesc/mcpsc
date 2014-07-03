#!/usr/bin/python

from Bio.PDB import PDBParser, Polypeptide
from Bio import pairwise2
from os import listdir
from os.path import isfile, join

def get_structure(pdb_fname):
  parser = PDBParser()
  return parser.get_structure(pdb_fname, pdb_fname)

def get_sequence(structure, model_id = 0):
  return Polypeptide.Polypeptide(structure[model_id].get_residues()).get_sequence()

## extreemly slow for complex pairs
def get_best_pairwise_alignment_score(seq1, seq2):
  global_score = max([a[2] for a in pairwise2.align.globalxx(str(seq1), str(seq2))])
  #local_score = max([a[2] for a in pairwise2.align.localxx(str(seq1), str(seq2))])
  return global_score
  
def write_sequences(pdb_dir):
  pdbfiles = [f for f in listdir(pdb_dir) if isfile(join(pdb_dir,f)) and f.endswith('.pdb')]
  
  for file in pdbfiles:
    pdb_file_name = join(pdb_dir, file)
    print pdb_file_name
    out = open('%s' %pdb_file_name.replace('.pdb','.seq'), 'w')
    out.write('%s\n' %str(get_sequence(get_structure(pdb_file_name))))
    out.close
    
   
#########################################################################################
#print str(get_sequence(get_structure('/home/xule/workspace/pdb_chew_kedem/1aa9.pdb')))
#seq1 = get_sequence(get_structure('/home/xule/workspace/pdb_chew_kedem/1flp.pdb'))
#print str(seq1)
#seq2 = get_sequence(get_structure('/home/xule/workspace/pdb_chew_kedem/4enl.pdb'))
#print (seq2)
#print get_best_pairwise_alignment_score(seq1, seq2)
#write_sequences('../pdb_chew_kedem')
write_sequences('../pdb_rost_sander')