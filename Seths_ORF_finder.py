### This program finds all open reading frames (ORFs) above a specified length.  Input is a fasta file of contigs, genomes or scaffolds
### and it returns a file of ORF coordinates and fasta file of ORF sequences.

import sys
from Bio import SeqIO

### make codon table 11
### any genetic code table from NCBI can be used here

table11 = '''FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
---M------**--*----M------------MMMM---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'''.split('\n')


### This creates python dictionaries of...
### codons = codons to amino acids
### starts = start codons
### stops = stop codons

codons,starts,stops = {},{},{}

for i in range(len(table11[0])):
	codons[table11[2][i]+table11[3][i]+table11[4][i]] = table11[0][i]
	if table11[1][i] == 'M':
		starts[table11[2][i]+table11[3][i]+table11[4][i]] = table11[1][i]
	if table11[1][i] == '*':
		stops[table11[2][i]+table11[3][i]+table11[4][i]] = table11[1][i]


def translate(seq,table11,scaffold_name,min_pro_len,strand,frame,dna):
	if strand == '+':
		orf,s = '',0
		for nt in range(0,len(seq)-2,3):
			if seq[nt:nt+3] in starts:
				if orf == '':
					orf += codons.get(seq[nt:nt+3])
					s = nt
				else:
					orf += codons.get(seq[nt:nt+3])
			elif seq[nt:nt+3] in stops:
				if orf == '':
					continue
				else:
					if len(orf) > min_pro_len:
						out_genes = open(sys.argv[1]+'.long_genes.fasta','a')
						out_coordinates = open(sys.argv[1]+".gene_coordinates.txt",'a')
						out_genes.write('>'+str(scaffold_name)+'\n'+str(orf)+'\n')
						out_coordinates.write(str(scaffold_name)+'\t'+strand+'\t'+str(s+frame)+'\t'+str(nt+frame)+'\n')
						out_genes.close()
						out_coordinates.close()
						orf = ''
						s = 0
                                        elif len(orf)>50:
                                                out_coordinates = open(sys.argv[1]+".gene_coordinates.txt",'a')
                                                out_coordinates.write(str(scaffold_name)+'\t'+strand+'\t'+str(len(dna)-nt-frame)+'\t'+str(len(dna)-s-frame)+'\n')
                                                out_coordinates.close()
                                                orf = ''
                                                s = 0
			else:
				if orf == '':
					continue
				else:
					orf += codons.get(seq[nt:nt+3],'X')
	elif strand == '-':
		orf,s = '',0
		for nt in range(0,len(seq)-2,3):
			if seq[nt:nt+3] in starts:
				if orf == '':
					orf += codons.get(seq[nt:nt+3])
					s = nt
				else:
					orf += codons.get(seq[nt:nt+3])
			elif seq[nt:nt+3] in stops:
				if orf == '':
					continue
				else:
					if len(orf) > min_pro_len:
						orfs.append(orf)
						out_genes = open(sys.argv[1]+'.long_genes.fasta','a')
						out_coordinates = open(sys.argv[1]+".gene_coordinates.txt",'a')
						out_genes.write('>'+str(scaffold_name)+'\n'+str(orf)+'\n')
						out_coordinates.write(str(scaffold_name)+'\t'+strand+'\t'+str(len(dna)-nt-frame)+'\t'+str(len(dna)-s-frame)+'\n')
						out_genes.close()
						out_coordinates.close()
						orf = ''
						s = 0
					elif len(orf)>50:
                                                out_coordinates = open(sys.argv[1]+".gene_coordinates.txt",'a')
                                                out_coordinates.write(str(scaffold_name)+'\t'+strand+'\t'+str(len(dna)-nt-frame)+'\t'+str(len(dna)-s-frame)+'\n')
                                                out_coordinates.close()
                                                orf = ''
                                                s = 0
			else:
				if orf == '':
					continue
				else:
					orf += codons.get(seq[nt:nt+3],'X')


min_pro_len = 500


for record in SeqIO.parse(sys.argv[1],'fasta'):
	dna = record.seq
	if len(record.seq) < 300:
		continue
	for strand, nuc in [('+', record.seq), ('-', record.seq.reverse_complement())]:
		for frame in range(3):
			length = 3 * ((len(record)-frame) // 3) #Multiple of three
			translate(nuc[frame:frame+length],table11,str(record.description),min_pro_len,strand,frame,dna)

