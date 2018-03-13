import sys
from Bio import SeqIO

kmers = {i.strip().split(',')[0]:0 for i in open('/lustre/projects/SethCommichaux/MitochondriaMysteries/Tribolium_castaneum_probes')}

seq_file = sys.argv[1]
results_file = seq_file[:seq_file.find('.')]+'_kmer_counts_dict.txt'

if seq_file.endswith('fasta'):
        type_parse = 'fasta'
else:
        type_parse = 'fastq'

for h,i in enumerate(SeqIO.parse(seq_file,type_parse)):
        dna = str(i.seq)
        rev_comp = str(i.seq.reverse_complement())
        for j in range(0,len(dna)+1-30):
                k = dna[j:j+30]
                l = rev_comp[j:j+30]
                if k in kmers:
                        kmers[k] += 1
                if l in kmers:
                        kmers[l] += 1

o = open(results_file,'w')

for k,v in sorted(kmers.items(),reverse=True,key=lambda x:x[1]):
        o.write(str(k)+'\t'+str(v)+'\n')

o.close()

