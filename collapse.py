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

numRows = 43
if(sys.argv[1].startswith("a")):
	numRows=38

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


#this function finds the differences between input files
#collapses higher res alleles together -- notes where there are disagreements ("#")
def collapse(filename):
  filtered_path = os.path.join(filename)
  higher_res_alleles = list()
  extended_alleles_dict = defaultdict(list)

  with open(filtered_path) as filename:
    reader = csv.reader(filename, delimiter=" ")
    higher_res_alleles.extend(next(reader))

  line_no = 2 # start line for each extended allele
  for allele in higher_res_alleles:
    with open(filtered_path) as filename:
      reader = csv.reader(filename, delimiter=" ")
      for row in itertools.islice(reader, line_no, None, len(higher_res_alleles)+1):
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

  for i in range(0, numRows): # replace ints with len(max) calls
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




spatial_counter = 0
def search_across():
  total_disagreements = 0
  total_disag_alleles = defaultdict(list)
  global spatial_counter
  for key in disagreements_dict.iterkeys():
    total_disag_alleles[key].append(0)

  for i in range(0, numRows):
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




coding_mapped = list()
# to codons preserving allele data via coords
# gets position of coding sequence -- the bases within "|"
# this region will be the same for all a*.txt alleles, and all b*.txt alleles
def get_protein_sequence_coords(first_file):
 
  print("\n\n###########################################################")
  print("## get protein coding sequence coordinates of first file ##")
  print("###########################################################")

  line_start = 1 # start line for each extended allele
  file_path = os.path.join(first_file)
  
  higher_res_alleles = list()
  nuc_list = list()
  coords = list()

  with open(file_path) as filename:
    reader = csv.reader(filename, delimiter=" ")
    higher_res_alleles.extend(next(reader))

  step = len(higher_res_alleles) + 1

  with open(file_path) as filename:
    reader = csv.reader(filename, delimiter=" ")
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

#  nucseq = "".join(map("".join, nuc_list)).replace(".", "").split("|")
  nucseq = "".join(map("".join, nuc_list)).split("|")
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


  #the coding regions all seem to be correct, but coords_array is not eg. 
  #TGTGA
  #[[33, [1, 0]], [33, [2, 5]]] -- this is 15 not 5, but should be 5
  #quick hack
  for coord in coords_array:
    if(coord[0][0] == coord[1][0]):
      if(coord[1][1][0] == coord[0][1][0] + 1):
        coord[1][1][0] = coord[1][1][0] - 1

  global coding_mapped
  coding_mapped = zip(coding_region, coords_array)

  for i in coding_mapped:
    print list(i), "start-end, [row, column, block index]"




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



#using positions from get_protein_sequence_coords(), map the bases/codons within this region to amino acid
def map_coding_to_dicts():
  print("\n\n###########################")
  print("## Map coding to dicts() ##")
  print("###########################")
  
  temp_dict = {}
#  exon_dict = defaultdict(list)

  print "Each coding region: "
  for region in coding_mapped:
    print region[1] # this prints out each coding region 
    item=region[1]

    startrow=item[0][0]
    startblock=item[0][1][0]
    startind=item[0][1][1]

    endrow=item[1][0]
    endblock=item[1][1][0]
    endind=item[1][1][1]

    # NEW NEW EXPERIMENT --------------------------------- START
    for allele in collapsed_alleles_dict.iterkeys(): #for each allele
      if(allele not in temp_dict):
        temp_dict[allele] = ""
      for row in range(startrow, endrow + 1): #for each row
        for block in range(0, 10): #for each block
          for idx, ca in enumerate(collapsed_alleles_dict[allele][row][block]): #for each collapsed allele (enumerated)
            if( row == startrow == endrow and startblock <= block <= endblock): #if startrow == endrow
              if(block == startblock == endblock and startind <= idx < endind):
                temp_dict[allele] += ca
              elif(block == startblock != endblock and idx >= startind):
                temp_dict[allele] += ca
              elif(block == endblock != startblock and idx < endind): 
                temp_dict[allele] += ca
              elif (block != startblock and block != endblock):
                temp_dict[allele] += ca
            elif(row == startrow != endrow and startblock <= block): #if block == startblock and startblock != endblock
              if(startblock == block and idx >= startind):
                temp_dict[allele] += ca
              elif(startblock != block):
                temp_dict[allele] += ca
            elif(row == endrow != startrow and endblock >= block): #if block == endblock and endblock != startblock
              if(endblock == block and idx < endind):
                temp_dict[allele] += ca
              elif(endblock != block):
                temp_dict[allele] += ca
            elif(startrow < row < endrow): #between start and end
              temp_dict[allele] += ca

      temp_dict[allele] += "|"
      # NEW NEW EXPERIMENT --------------------------------- END
  print " "



  for i in range(len(coding_mapped)):
    ref_bases = coding_mapped[i][0]
    print "Reference : ", ref_bases, "LEN:", len(ref_bases)
    for allele in temp_dict:

      collapsed = temp_dict[allele].split("|")
      print allele, ": ", collapsed[i], "LEN:", len(collapsed[i])
    print



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





### call functions, print the results

#this finds the disagreements between input files
for file in filenames:
  collapse(file)

print "PRINTING DISAG"
print_dict(collapsed_alleles_dict)
print_dict(disagreements_dict)
print "================================"
search_across()
print "Absolute position for each disagreement:", spatial_list
# print "Aminoacid sequence:", read_amino()
get_protein_sequence_coords(filenames[0])
map_coding_to_dicts()
