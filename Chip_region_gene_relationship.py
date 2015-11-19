import sys
import re
IN_FL1=open(sys.argv[1])
gene={}
chrom={}
for i in IN_FL1:
	element=i.strip().split('\t')
	gene[element[0]]=[element[1],element[2],element[3]]
	(tmp,)=re.search('^AT(\w)G',element[0]).groups()
	tmp='Chr'+tmp
	if tmp in chrom.keys():
		chrom[tmp].append(element[0])
	else:
		chrom[tmp]=[element[0]]
IN_FL1.close()

IN_FL2=open(sys.argv[2])
OUT_FL=open(sys.argv[3],mode='w')
OUT_FL.write('\t'.join(('chr','start','end','length','summit','tags','-10*log10(pvalue)','fold_enrichment','FDR(%)','RPM','promoter','gene_body')))
OUT_FL.write('\n')
for i in IN_FL2:
	promoter=[]
	gene_body=[]
	if not re.search('^Chr',i):
		if re.search('# tags after filtering in treatment:',i):
			(total_reads,)=re.search('# tags after filtering in treatment: (\d+)',i).groups()
			total_reads=float(total_reads)/1000000
		continue
	element=i.strip().split('\t')
	element[1]=int(element[1])
	element[2]=int(element[2])
	RPM=int(element[5])/total_reads
	for j in chrom[element[0]]:
		if gene[j][2]=='+' and (((int(gene[j][0])-2000) < element[1] < int(gene[j][1])) or ((int(gene[j][0])-2000) < element[2] < int(gene[j][1]))):
			if ((int(gene[j][0])-2000) < element[1] < int(gene[j][0])) or ((int(gene[j][0])-2000) < element[2] < int(gene[j][0])):
				promoter.append(j)
			if (int(gene[j][0]) < element[1] < int(gene[j][1])) or (int(gene[j][0]) < element[2] < int(gene[j][1])):
				gene_body.append(j)
		if gene[j][2]=='-' and ((int(gene[j][0]) < element[1] < (int(gene[j][1])+2000)) or (int(gene[j][0]) < element[2] < (int(gene[j][1])+2000))):
			if (int(gene[j][1]) < element[1] < (int(gene[j][1])+2000)) or (int(gene[j][1]) < element[2] < (int(gene[j][1])+2000)):
				promoter.append(j)
			if (int(gene[j][0]) < element[1] < int(gene[j][1])) or (int(gene[j][0]) < element[2] < int(gene[j][1])):
				gene_body.append(j)
	OUT_FL.write(i.strip()+'\t'+str(RPM)+'\t'+' '.join(sorted(promoter))+'\t'+' '.join(sorted(gene_body))+'\n')
