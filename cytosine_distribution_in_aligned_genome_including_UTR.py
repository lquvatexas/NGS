#Calculate cytosine distribution in aligned genome (conserved cytosines between diploid and tetraploid), including intron, exon, and different types of TE.
import sys
import numpy
import re
total_cytosine=0
def extract_cytosine_location(dir,file):
	global total_cytosine
	location={}
	IN_tmp=open(dir+'/'+file)
	for line_tmp in IN_tmp:
		total_cytosine+=1
		element=line_tmp.split('\t')
		location[int(element[0])]=1
	return location
IN=open('Gh_exon_intron_TE.txt')
IN_gene=open('Gh_gene_region.txt')
OUT=open('cytosine_distribution_in_aligned_genome_including_UTR.txt','w')
cytosine={}
chroms={}
for line in IN:
	tmp=line.strip().split('\t')
	if tmp[0] in chroms:
		chroms[tmp[0]].append(line.strip())
	else:
		chroms[tmp[0]]=[line.strip()]
IN.close()

for line in IN_gene:
	tmp=line.strip().split('\t')
	if tmp[4]=='+':
		chroms[tmp[1]].append(tmp[1]+'\t'+'upstream'+'\t'+str(int(tmp[2])-1000)+'\t'+str(int(tmp[2])-1))
		chroms[tmp[1]].append(tmp[1]+'\t'+'downstream'+'\t'+str(int(tmp[3])+1)+'\t'+str(int(tmp[3])+1000))
	else:
		chroms[tmp[1]].append(tmp[1]+'\t'+'upstream'+'\t'+str(int(tmp[3])+1)+'\t'+str(int(tmp[3])+1000))
		chroms[tmp[1]].append(tmp[1]+'\t'+'downstream'+'\t'+str(int(tmp[2])-1000)+'\t'+str(int(tmp[2])-1))
chrom_flag=''
for chrom_id in chroms:
	if chrom_flag != chrom_id:
		print(chrom_id)
		location=extract_cytosine_location('/scratch/02489/qs723/total_site_evolution/conserved_2000_diploid_to_tetraploid',chrom_id)
		chrom_flag=chrom_id
	for region in chroms[chrom_id]:
		tmp_region=region.split('\t')
		for i in range(int(tmp_region[2]),int(tmp_region[3])):
			if i in location:
				if tmp_region[1] in cytosine:
					cytosine[tmp_region[1]]+=1
				else:
					cytosine[tmp_region[1]]=1
OUT.write('total cytosines: '+str(total_cytosine)+'\n')
for i in cytosine:
	OUT.write(str(i)+'\t'+str(cytosine[i])+'\n')

		
