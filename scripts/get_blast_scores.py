#!/usr/bin/python

import fileinput

##
score = 0; expect = 0; identity = 0; 
positives = 0; gaps = 0;
##
for line in fileinput.input():
  line = line.replace('\n','')
  if line.startswith('BLASTP'):
    score = -1; expect = -1; identity = -1; 
    positives = -1; gaps = -1;
    continue
  if line.startswith('Window'):
    print '%s %s %s %s %s' %(score, expect, identity, positives, gaps)
    continue
  # do something
  if line.find('Score') != -1:
    data = line.split(' ')
    score = data[3]; expect = data[9][:len(data[9])-1]
  if line.find('Identities') != -1:
    data = line.split(' ')
    identity = data[4][1:len(data[4])-2]
    positives = data[8][1:len(data[8])-2]
    gaps = data[12][1:len(data[12])-1]


##
##