#!/usr/bin/env python

import sys
import dxpy
import gzip
from Bio import SeqIO

expected_colnames = [
  '_name',
  '_loci',
  "Human",
  "Chimp",
  "Gorilla",
  "Orangutan",
  "Gibbon",
  "Rhesus",
  "Crab_eating_macaque",
  "Baboon",
  "Green_monkey",
  "Marmoset",
  "Squirrel_monkey",
  "Bushbaby",
  "Chinese_tree_shrew",
  "Squirrel",
  "Lesser_Egyptian_jerboa",
  "Prairie_vole",
  "Chinese_hamster",
  "Golden_hamster",
  "Mouse",
  "Rat",
  "Naked_mole_rat",
  "Guinea_pig",
  "Chinchilla",
  "Brush_tailed_rat",
  "Rabbit",
  "Pika",
  "Pig",
  "Alpaca",
  "Bactrian_camel",
  "Dolphin",
  "Killer_whale",
  "Tibetan_antelope",
  "Cow",
  "Sheep",
  "Domestic_goat",
  "Horse",
  "White_rhinoceros",
  "Cat",
  "Dog",
  "Ferret",
  "Panda",
  "Pacific_walrus",
  "Weddell_seal",
  "Black_flying_fox",
  "Megabat",
  "David's_myotis",
  "Microbat",
  "Big_brown_bat",
  "Hedgehog",
  "Shrew",
  "Star-nosed_mole",
  "Elephant",
  "Cape_elephant_shrew",
  "Manatee",
  "Cape_golden_mole",
  "Tenrec",
  "Aardvark",
  "Armadillo"
]

def load_ccds_seqs(fn):
  h = gzip.open(fn)
  ans = {}
  for record in SeqIO.parse(h, 'fasta'):
    ccds_id = record.id.split('|')[0]
    ans[ccds_id] = str(record.seq)
  return ans

def similarity(s1,s2):
  s1 = s1.upper()
  s2 = s2.upper()
  n = len(s1)
  assert n == len(s2)
  m = 0
  k = 0
  for i in xrange(n):
    c1 = s1[i]
    c2 = s2[i]
    if c1 != '-' and c2 != '-':
      m = m+1
      if c1 == c2:
        k = k+1
  return m, k

def main():
  job=dxpy.DXJob(dxid=sys.argv[1])
  ccds_seqs=load_ccds_seqs(sys.argv[2])
  alignments = dxpy.get_handler(job.describe()['output']['alignments'])
  alignments_desc = alignments.describe()

  # verify basic table properties
  assert 'CrossSpeciesAlignments' in alignments_desc['types']
  assert alignments_desc['length'] == 318

  # verify table schema
  cols = alignments.get_columns()
  assert len(cols) == len(expected_colnames)
  for i in xrange(len(expected_colnames)):
    assert cols[i]['type'] == 'string'
    assert cols[i]['name'] == expected_colnames[i]

  # verify table contents
  names = []

  bps = 0
  alns = [0 for i in xrange(3, len(expected_colnames))]
  ids = [0 for i in xrange(3, len(expected_colnames))]

  for row in alignments.iterate_rows():
    ccds_id = row[1]
    names.append(ccds_id)
    assert len(row) == 1+len(expected_colnames)
    # check length of each consensus sequence
    for col in xrange(4,len(row)):
      assert len(row[col]) == len(row[3])
    # check that human sequence matches CCDS
    hs_seq = str(row[3]).translate(None, '-').upper()
    ccds_seq = ccds_seqs[ccds_id].upper()
    if len(ccds_seq) == 3+len(hs_seq): # stop codon
      ccds_seq = ccds_seq[:len(ccds_seq)-3]
    if hs_seq != ccds_seq:
      print 'Sequence mismatch for {}'.format(ccds_id)
      print 'Sequence found in alignments table: {}'.format(hs_seq)
      print 'Expected sequence from CCDS: {}'.format(ccds_seq)
      assert False
    # collect alignment statistics
    bps = bps + len(hs_seq)
    for i in xrange(4, 1+len(expected_colnames)):
      a, s = similarity(row[3],row[i])
      alns[i-4] = alns[i-4] + a
      ids[i-4] = ids[i-4] + s

  assert len(set(names)) == 318

  print 'Alignment quality statistics (make sure these look reasonable):'
  print 'Informant\tAligned %\tIdentity %'
  for i in xrange(len(expected_colnames)-3):
    print '{}\t{}\t{}'.format(expected_colnames[i+3].rjust(12), alns[i]*100/bps, ids[i]*100/alns[i])

main()
