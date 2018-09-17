# Script to generate collapse found extended alleles in database files.
#
# Use:
# python collapse.py
#
# Author: Apostolos Papadopoulos
# Last update: Check with git log
# Last release: <beta version>

import os, sys, csv, itertools
from collections import defaultdict
from pprint import pprint

filenames = list()
if len(sys.argv) == 1:
  print "collapse.py requires at least one filename as an argument"
  exit(1)
else:
  try:
    if sys.argv[1]:
      filenames.append(sys.argv[1])
    if sys.argv[2]:
      filenames.append(sys.argv[2])
    if sys.argv[3]:
      filenames.append(sys.argv[3])
  except IndexError:
    pass

# gen_file = sys.argv[3]
# gen_path = os.path.join(gen_file)

collapsed_alleles_dict = defaultdict(list)
disagreements_dict = defaultdict(list)

def compare(alist):
  counter = 0
  zipped = zip(*alist)
  zlist = list(zipped)
  ret_str = ""
  ret_list = list()
  for char_pos in zlist:
    diff = all(x == char_pos[0] for x in char_pos)
    if diff == True:
      ret_str += char_pos[0]
      counter += 1
    else:
      ret_str += "#"
      temp = list()
      for char in char_pos:
        temp.append(char)
      temp.append(counter)
      ret_list.append(list(temp))
      counter += 1
  return ret_str, ret_list

def collapse(file):
  filtered_path = os.path.join(file)
  higher_res_alleles = list()
  extended_alleles_dict = defaultdict(list)

  with open(filtered_path) as file:
    reader = csv.reader(file, delimiter=" ")
    higher_res_alleles.extend(next(reader))

  line_no = 2 # start line for each extended allele
  for allele in higher_res_alleles:
    with open(filtered_path) as file:
      reader = csv.reader(file, delimiter=" ")
      for row in itertools.islice(reader, line_no, None,
        len(higher_res_alleles)+1):
        # tentative fix to `|` bug
        # delete them \m/ \m/ \m/ \m/
        list_val = filter(None, row[2:])
        for index, obj in enumerate(list_val):
          if obj == "|":
            if len(list_val[index-1]) < 10 and len(list_val[index+1]) < 10:
              del list_val[index]
              list_val[index-1:index+1] = [''.join(list_val[index-1:index+1])]
            elif len(list_val[index-1]) == 10 and len(list_val[index+1]) == 10:
              del list_val[index]
            elif len(list_val[index-1]) == 10 and len(list_val[index+1]) < 10:
              del list_val[index]
        extended_alleles_dict[allele].append(list_val)
      line_no += 1

  # k = key and list of list values as a dict item/pair
  # k[0] = allele only (NOT iterable)
  # k[1] = list of list only (all lines of this allele) (NOT iterable)
  # k[1][0] = a specific whole line (0 is iterable; substitute with 1 =
  #           second line, etc)
  # k[1][0][0] = specific block from whole line (second [0] is iterable,
  #           substitute with 1 = second block, etc)

  # print len(v): number of sequence lines
  # max number of second []: 41
  # max number of third []: 9
  collapsed_allele_resolution_key = ""
  same_snp = list()
  counter = 0

  for i in range(0, 43): # replace ints with len(max) calls
    for j in range(0, 10): # replace ints with len(max) calls
      templist = list()
      for k in extended_alleles_dict.iteritems():
        try:
          templist.append(k[1][i][j])
          if k[0][7] == ":":
            collapsed_allele_resolution_key = k[0][0:7]
          else:
            collapsed_allele_resolution_key = k[0][0:8]
          collapsed_alleles_dict[collapsed_allele_resolution_key]
        except IndexError:
          pass
      diff = all(x == templist[0] for x in templist)
      if diff == True:
        try:
          same_snp.append(templist[0])
        except IndexError:
          pass
      else:
        snp_diff, snp_disag = compare(list(set(templist)))
        same_snp.append(snp_diff)
        snp_disag.append(i)
        snp_disag.append(j)
        disagreements_dict[collapsed_allele_resolution_key].append(
          list(snp_disag)
        )
        counter += (len(snp_disag) - 2)

  # split values list every 10 list items (equiv to "blocks" from db file)
  composite_list = [same_snp[x:x+10] for x in range(0, len(same_snp),10)]

  for i in composite_list:
    collapsed_alleles_dict[collapsed_allele_resolution_key].append(i)

  disagreements_dict[collapsed_allele_resolution_key].append(counter)




coding_mapped = list()
# to codons preserving allele data via coords
def get_sequence_coords(first_file):
 
  print("\n\n###########################")
  print("## get sequence coords() ##")
  print("###########################")

  line_start = 1 # start line for each extended allele
  file_path = os.path.join(first_file)


  higher_res_alleles = list()
  nuc_list = list()
  coords = list()

  with open(file_path) as file:
    reader = csv.reader(file, delimiter=" ")
    higher_res_alleles.extend(next(reader))

  step = len(higher_res_alleles) + 1

  with open(file_path) as file:
    reader = csv.reader(file, delimiter=" ")
    cnt = 0

    for row in itertools.islice(reader, line_start, None, step):
      seqrow = filter(None, row[2:])

      for index, obj in enumerate(seqrow):
        if obj == "|":
          if len(seqrow[index-1]) < 10 and len(seqrow[index+1]) < 10:
            seqrow[index-1:index+1] = [''.join(seqrow[index-1:index+1])]
            seqrow[index-1:index+1] = [''.join(seqrow[index-1:index+1])]
            # del seqrow[index]
          elif len(seqrow[index-1]) == 10 and len(seqrow[index+1]) == 10:
            seqrow[index+1] = "|" + seqrow[index+1]
            del seqrow[index]

      nuc_list.append(seqrow)

      try:
        if "|" in "".join(seqrow):
          block_idx = 0
          block_pos = 0
          pipe_coords = list()
          for block_cnt, block in enumerate(seqrow):
            if "|" in block:
              block_idx = block_cnt
              block_pos = block.index("|")
              pipe_coords.append([block_idx, block_pos])
          # print "".join(seqrow), cnt, pipe_coords
          coords.append([cnt, pipe_coords])
        else:
          # print "".join(seqrow)
          pass
      except ValueError:
        pass
      cnt += 1

  nucseq = "".join(map("".join, nuc_list)).replace(".", "").split("|")
  coding_region = nucseq[1::2]

  flat_list = list()

  for item in coords:
    if len(item[1]) < 2:
      for i in item[1]:
        flat_list.append([item[0], i])
    else:
      for i in item[1]:
        flat_list.append([item[0], i])

  coords_array = list()
  for index, obj in enumerate(flat_list):
    if index % 2 == 0:
      coords_array.append([flat_list[index], flat_list[index+1]])

  global coding_mapped
  coding_mapped = zip(coding_region, coords_array)

  for i in coding_mapped:
    print list(i)

def map_coding_to_dicts():
  print("\n\n###########################")
  print("## Map coding to dicts() ##")
  print("###########################")
  
  temp_dict = defaultdict(list)
  exon_dict = defaultdict(list)

  keylist = list()
  for k in collapsed_alleles_dict.iterkeys():
    keylist.append(k)

  print " "
  alleleslist = list()
  for region in coding_mapped:
    for item in region:
      if type(item) == list:
        print item

        # NEW NEW EXPERIMENT --------------------------------- START
        for i in range(0, 43):
          if i >= item[0][0] and i <= item[1][0]:
            for k in keylist:
              if i == item[0][0]: # start
                for j in range(0, 10):
                  if j == item[0][1][0]: # first block
                    for idx, cr in enumerate(collapsed_alleles_dict[k][i][j]):
                      if idx >= j:
                        temp_dict[k].append(cr)
                  elif j > item[0][1][0]:
                    for idx, cr in enumerate(collapsed_alleles_dict[k][i][j]):
                      temp_dict[k].append(cr)
              elif i == item[1][0]: # end
                for j in range(0, 10):
                  if j == item[1][1][0]: # last block
                    for idx, cr in enumerate(collapsed_alleles_dict[k][i][j]):
                      if idx <= j:
                        temp_dict[k].append(cr)
                  elif j <= item[1][1][0]:
                    for idx, cr in enumerate(collapsed_alleles_dict[k][i][j]):
                      temp_dict[k].append(cr)
              else: # whole blocks
                for j in range(0, 10):
                  for idx, cr in enumerate(collapsed_alleles_dict[k][i][j]):
                    temp_dict[k].append(cr)

        # NEW NEW EXPERIMENT --------------------------------- END
        for k, v in temp_dict.iteritems():
          temp_dict[k].append("|")

  print " "
  for region in coding_mapped:
    for k2, v2 in temp_dict.iteritems():
      # print k2, v2
      exon_dict[region[0]].append(v2)

  # print_dict(exon_dict)

  print " "
  c = 0
  for k, v in exon_dict.iteritems():
    t_list = list("".join(v[0]).replace(".", "").split("|"))
    t_list2 = list("".join(v[1]).replace(".", "").split("|"))
    print k, len(k), "LEN"
    # print t_list
    # print t_list2
    for index, val in enumerate(t_list):
      if index == c:
        print val
    for index, val in enumerate(t_list2):
      if index == c:
        print val
    c += 1

spatial_counter = 0
def search_across():
  total_disagreements = 0
  total_disag_alleles = defaultdict(list)
  global spatial_counter
  for key in disagreements_dict.iterkeys():
    total_disag_alleles[key].append(0)

  for i in range(0, 43):
    for j in range(0, 10):
      collapsed_diffs_list = list()
      for k in collapsed_alleles_dict.iteritems():
        try:
          # print k[1][i][j]
          collapsed_diffs_list.append(k[1][i][j])
        except IndexError:
          pass

      a, b = compare(collapsed_diffs_list)

      # print collapsed characters if non-empty list
      if len(b) > 0:
        print "ORIG: ", b, i, j
        total_disagreements += len(b)
        coordinate(b, spatial_counter)

      # print the diff disagreements
      for k, v in disagreements_dict.iteritems():
        for item in v:
          try:
            if item[-2] == i and item[-1] == j:
              print "DIFF: ", item, k
              total_disag_alleles[k][0] = (total_disag_alleles[k][0] + 1)
          except TypeError:
            pass

      spatial_counter += 10

  print "Total character disagreements: ", total_disagreements
  print "Ranking of disagreeing alleles: ", total_disag_alleles

spatial_list = list()
def coordinate(char_list, spc):
  for item in char_list:
    try:
      if spc > 0:
        # print spatial_counter, item[-1]
        # print (int(spatial_counter) + int(item[-1]))
        spatial_list.append(int(spc) + int(item[-1]))
      else:
        spatial_list.append(int(spc) + int(item[-1]))
    except TypeError:
      pass

codons = {
  "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
  "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
  "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
  "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
  "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
  "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
  "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
  "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
  "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
  "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
  "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
  "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
  "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
  "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
  "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
  "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M"
}

def read_amino():
  nuc_list = list()
  with open(gen_path) as file:
    reader = csv.reader(file, delimiter=" ")
    for row in reader:
      nuc_list.append(filter(None, row[2:]))

  nuqsec = "".join(map("".join, nuc_list)).replace(".", "").split("|")

  coding_region = nuqsec[1::2]
  coding_region_str = "".join(coding_region)
  codonlist = map("".join, zip(*[iter(coding_region_str)]*3))

  amino_list = list()
  for codon in codonlist:
    for k, v in codons.iteritems():
      if codon == k:
        # print nuq, v
        amino_list.append(v)

  amino_str = "".join(amino_list)
  return amino_str



def print_dict(adict):
  for k, v in adict.iteritems():
    print k
    for i in v:
      print i





# call functions, print the results
for file in filenames:
  collapse(file)

print "PRINTING DISAG"
print_dict(collapsed_alleles_dict)
print_dict(disagreements_dict)
print "================================"
search_across()
print "Absolute position for each disagreement:", spatial_list
# print "Aminoacid sequence:", read_amino()
get_sequence_coords(filenames[0])
map_coding_to_dicts()
