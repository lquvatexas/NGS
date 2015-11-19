import re
def calculate_Cnumber(genome_dir,region_file,gene_file, promoter_file):
	'''caculate methylated Cytosine number in gene body and gene promoter
	genome_dir: chromsome directory; region_file: gene region informitaion file; 
	gene_file: output methylated c number in gene body; promoter_file:output methylated c number in gene promoter.
	'''
	OUT_gene_file=open(gene_file,mode='w')
	OUT_promoter_file=open(promoter_file,mode='w')
	with open(region_file) as IN_file:
		chrom_flag=''
		#direction is used to determine the regin of promoter according to the gene direction
		direction={'+':[1,-1],'-':[2,1]}
		for line in IN_file:
			gene=line.strip().split('\t')
			number={'C':0,'C_cover':0,'CG':0,'CHG':0,'CHH':0}
			if re.search(r'Glyma(\d+)g\d+\.1',gene[0]):
				chrom_id='Gm'+re.search(r'Glyma(\d+)g\d+\.1',gene[0]).groups()[0]
			else:
				print (gene[0]+'is not right gene id')
			if chrom_id != chrom_flag:
				with open(genome_dir+'/'+chrom_id) as chrom_file:
					location={}
					for i in chrom_file:
						tmp_base=i.strip().split('\t')
						location[int(tmp_base[0])]=[tmp_base[2],(int(tmp_base[5])+int(tmp_base[6])),tmp_base[12]]
				chrom_flag=chrom_id
			#caculate C of all types number in gene body
			for base in range(int(gene[1]),int(gene[2])+1):
				if location.get(base):
					number['C']+=1
					if location[base][1] >1:
						number['C_cover']+=1
					if location[base][2] == 'Y':
						number[location[base][0]]+=1
			OUT_gene_file.write(gene[0]+'\t'+'\t'.join([str(number[stem]) for stem in ('C','C_cover','CG','CHG','CHH')])+'\n')
			#caculate C of all types number in gene promoter
			number={'C':0,'C_cover':0,'CG':0,'CHG':0,'CHH':0}
			for base in range((int(gene[direction[gene[3]][0]])+direction[gene[3]][1]),(int(gene[direction[gene[3]][0]])+direction[gene[3]][1]*1001)):
				if location.get(base):
					number['C']+=1
					if location[base][1] >1:
						number['C_cover']+=1
					if location[base][2] == 'Y':
						number[location[base][0]]+=1
			OUT_promoter_file.write(gene[0]+'\t'+'\t'.join([str(number[stem]) for stem in ('C','C_cover','CG','CHG','CHH')])+'\n')
	
def calculate_genebody_Cnumber(genome_dir,region_file,gene_file):
	'''caculate methylated Cytosine number in gene body
	genome_dir: chromsome directory; region_file: gene region informitaion file; 
	gene_file: output methylated c number in gene body.
	'''
	OUT_gene_file=open(gene_file,mode='w')
	with open(region_file) as IN_file:
		chrom_flag=''
		#direction is used to determine the regin of promoter according to the gene direction
		for line in IN_file:
			gene=line.strip().split('\t')
			number={'C':0,'C_cover':0,'CG':0,'CHG':0,'CHH':0}
			if re.search(r'Glyma(\d+)g\d+\.1',gene[0]):
				chrom_id='Gm'+re.search(r'Glyma(\d+)g\d+\.1',gene[0]).groups()[0]
			else:
				print (gene[0]+'is not right gene id')
			if chrom_id != chrom_flag:
				with open(genome_dir+'/'+chrom_id) as chrom_file:
					location={}
					for i in chrom_file:
						tmp_base=i.strip().split('\t')
						location[int(tmp_base[0])]=[tmp_base[2],(int(tmp_base[5])+int(tmp_base[6])),tmp_base[12]]
				chrom_flag=chrom_id
			#caculate C of all types number in gene body
			for base in range(int(gene[1]),int(gene[2])+1):
				if location.get(base):
					number['C']+=1
					if location[base][1] >1:
						number['C_cover']+=1
					if location[base][2] == 'Y':
						number[location[base][0]]+=1
			OUT_gene_file.write(gene[0]+'\t'+'\t'.join([str(number[stem]) for stem in ('C','C_cover','CG','CHG','CHH')])+'\n')

def calculate_promoter_Cnumber(genome_dir,region_file, promoter_file):
	'''caculate methylated Cytosine number in gene body and gene promoter
	genome_dir: chromsome directory; region_file: gene region informitaion file; 
	promoter_file:output methylated c number in gene promoter.
	'''
	OUT_promoter_file=open(promoter_file,mode='w')
	with open(region_file) as IN_file:
		chrom_flag=''
		#direction is used to determine the regin of promoter according to the gene direction
		direction={'+':[1,-1],'-':[2,1]}
		for line in IN_file:
			gene=line.strip().split('\t')
			number={'C':0,'C_cover':0,'CG':0,'CHG':0,'CHH':0}
			if re.search(r'Glyma(\d+)g\d+\.1',gene[0]):
				chrom_id='Gm'+re.search(r'Glyma(\d+)g\d+\.1',gene[0]).groups()[0]
			else:
				print (gene[0]+'is not right gene id')
			if chrom_id != chrom_flag:
				with open(genome_dir+'/'+chrom_id) as chrom_file:
					location={}
					for i in chrom_file:
						tmp_base=i.strip().split('\t')
						location[int(tmp_base[0])]=[tmp_base[2],(int(tmp_base[5])+int(tmp_base[6])),tmp_base[12]]
				chrom_flag=chrom_id
			#caculate C of all types number in gene promoter
			for base in range((int(gene[direction[gene[3]][0]])+direction[gene[3]][1]),(int(gene[direction[gene[3]][0]])+direction[gene[3]][1]*1001)):
				if location.get(base):
					number['C']+=1
					if location[base][1] >1:
						number['C_cover']+=1
					if location[base][2] == 'Y':
						number[location[base][0]]+=1
			OUT_promoter_file.write(gene[0]+'\t'+'\t'.join([str(number[stem]) for stem in ('C','C_cover','CG','CHG','CHH')])+'\n')

def compare_paralog_methylation(paralog_file,Cnumber_file,Outfile):
	OUT_file=open(Outfile,mode='w')
	with open(Cnumber_file) as IN_file1:
		number={}
		for line in IN_file1:
			gene=line.strip().split('\t')
			number[gene.pop(0)]=gene
	with open(paralog_file) as IN_file2:
		for line in IN_file2:
			gene=line.strip().split('_')
			for i in range(len(gene)):
				if i==0:
					OUT_file.write(gene[i]+'\t'+'\t'.join(number[gene[i]]))
				else:
					OUT_file.write('\t'+gene[i]+'\t'+'\t'.join(number[gene[i]]))
			OUT_file.write('\n')
	OUT_file.close()


if __name__ =="__main__":
	import sys
	calculate_Cnumber(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
