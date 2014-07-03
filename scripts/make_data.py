#!/usr/bin/python

import sys
import numpy as NP
from BeautifulSoup import BeautifulSoup
from Bio import ExPASy, SwissProt
from collections import defaultdict

tm_id = 0
ce_id = 1

class TimingData:
  def __init__(self, str):
    self.line = str
    (self.p1, self.p2, self.l1, self.l2, self.tmsec) = str.replace('\n','').split(' ')

## writes out p1 p2 l1 l2 ptime
## makes 2 files one for tmalign and one for ce
## name of files <algo_name>.infile.dat
def get_features_from_log(infile):
  ##
  tmfile = open('tm.%s.dat' %infile,'w')
  cefile = open('ce.%s.dat' %infile,'w')
  ##
  for line in open(infile, 'r'):
    ##
    if not line.startswith('Algo'):
      continue
    ##
    data = line.replace('\n','').replace(',','').split(' ')
    algo = int(data[1]); p1 = data[3]; p2 = data[5]; l1 = data[7]; l2 = data[9]; ptime = data[23]
    ##
    out_str = '%s %s %s %s %s\n' %(p1,p2,l1,l2,ptime)
    if algo == tm_id:
      tmfile.write('%s' %out_str)
    else:
      cefile.write('%s' %out_str)
  ##
  tmfile.close()
  cefile.close()
  
def get_quartiles_from_features(infile):
  ##
  lqfile = open('%slq.dat' %infile.replace('dat',''), 'w')
  uqfile = open('%suq.dat' %infile.replace('dat',''), 'w')
  ##
  data = [TimingData(line) for line in open(infile, 'r')]
  times = [int(d.tmsec) for d in data]
  lquart = NP.percentile(times,5)
  uquart = NP.percentile(times,95)
  for d in data:
    if d.p1 == d.p2:
      continue
    if int(d.tmsec) <= lquart:
      lqfile.write('%s' %d.line)
    elif int(d.tmsec) >= uquart:
      uqfile.write('%s' %d.line)

  ##
  lqfile.close()
  uqfile.close()
  
def get_uniprot_summary_for_pdb(infile, uniprot_dump_dir):
  outfile = open('%s.uniprot.dat' %infile, 'w')
  ##
  for line in open(infile, 'r'):
    pdb_dom = line.replace('\n','').split('.')[0]
    soup = BeautifulSoup(' '.join([line for line in open('%s/%s' %(uniprot_dump_dir,pdb_dom),'r')]))    
    result_row_htmls = soup.find(id='results').findAll('tr')
    result_row_htmls = result_row_htmls[1:] ## exclude the header row
    uniprot_ids_string = ''
    for result_row_html in result_row_htmls:
      uniprot_ids_string += result_row_html.findAll('td')[1].text + '|'
    uniprot_ids_string = uniprot_ids_string[:len(uniprot_ids_string) - 1]
    outfile.write('%s\t%d\t%s\n' %(pdb_dom, len(result_row_htmls), uniprot_ids_string))
  ##
  outfile.close()

def gen_uniprot_features_for_pdb(infile):
  for line in open(infile,'r'):
    (pdb_dom, count, uniprot_ids) = line.replace('\n','').split('\t')
    uniprot_ids = uniprot_ids.split('|')
    for uniprot_id in uniprot_ids:
      data = SwissProt.read(ExPASy.get_sprot_raw(uniprot_id)).__dict__  
      keep = False
      go = []; interpro = ''; evo_trace = ''
      for xref in data['cross_references']:
        if xref[0] == 'GO':
          go.append(xref[1])
        if xref[0] == 'InterPro':
          interpro = xref[1]
        if xref[0] == 'EvolutionaryTrace':
          evo_trace = xref[1]
        if xref[0] == 'PDB' and xref[1].lower() == pdb_dom.lower():
          keep = True
      if keep == False:
        continue
      organism = data['organism']
      loc = ''
      for comment in data['comments']:
        if comment.startswith('SUBCELLULAR LOCATION'):
          loc = comment
      print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(pdb_dom,uniprot_id,'|'.join(go),interpro,evo_trace,organism,loc)

def gen_all_features_file(pair_file, blast_file, uniprot_file):
  blast_data_list = [line.replace('\n','').split(' ') for line in open(blast_file, 'r') if line.startswith('..')]
  blast_data_map = defaultdict(dict)
  for blast_data in blast_data_list:
    dom1 = blast_data[0].split('/')[2].split('.')[0]
    dom2 = blast_data[1].split('/')[2].split('.')[0]
    bdata = []
    for d in blast_data[2:]:
      if d.endswith('%'):
        bdata.append(str(int(d[:len(d)-1])/float(100)))
      else:
        bdata.append(d)
    blast_data_map[dom1][dom2] = bdata[2:]
  
  uniprot_data_list = [line.replace('\n','').split('\t') for line in open(uniprot_file, 'r')]
  uniprot_data_map = {}
  for uniprot_data in uniprot_data_list:
    if uniprot_data_map.get(uniprot_data[0]) == None:
      uniprot_data_map[uniprot_data[0]] = uniprot_data[1:]
    else:
      v = uniprot_data_map[uniprot_data[0]]
      n = uniprot_data[1:]
      v = ['%s|%s' %(v[i],n[i]) for i in range(0,len(v))]
      uniprot_data_map[uniprot_data[0]] = v
  
  for line in open(pair_file, 'r'):
    (dom1, dom2, l1, l2, tmsec) = line.replace('\n','').split(' ')
    dom1 = dom1.split('.')[0]; dom2 = dom2.split('.')[0]
    blast_data = blast_data_map[dom1][dom2]
    if uniprot_data_map.get(dom1) == None:
      continue
    uniprot_data_dom1 = uniprot_data_map[dom1]
    if uniprot_data_map.get(dom2) == None:
      continue
    uniprot_data_dom2 = uniprot_data_map[dom2]
    diff_go_12 = list(set(uniprot_data_dom1[1].split('|')) - set(uniprot_data_dom2[2].split('|')))
    diff_go_21 = list(set(uniprot_data_dom2[1].split('|')) - set(uniprot_data_dom1[2].split('|')))
    uniprot_data_dom1[1] = str(len(uniprot_data_dom1[1].split('|')))
    uniprot_data_dom2[1] = str(len(uniprot_data_dom2[1].split('|')))
    print "'%s','%s','%s','%s','%s','%s','%s','%s','%d','%d'" %(dom1,dom2,tmsec,l1,l2,"','".join(blast_data),"','".join(uniprot_data_dom1),"','".join(uniprot_data_dom2),len(diff_go_12),len(diff_go_21))
    
    
###################################
#get_features_from_log(sys.argv[1])  
#get_quartiles_from_features(sys.argv[1])
#get_uniprot_summary_for_pdb(sys.argv[1], sys.argv[2])
gen_uniprot_features_for_pdb(sys.argv[1])
#gen_all_features_file(sys.argv[1], sys.argv[2], sys.argv[3])
##