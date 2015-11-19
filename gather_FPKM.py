import sys
gene={}
for i in range(1,len(sys.argv)-1):
	IN=open(sys.argv[i])
	tmp=IN.readline()
	for line in IN:
		tmp=line.split('\t')
		if tmp[0] in gene:
			gene[tmp[0]].append(tmp[9])
		else:
			gene[tmp[0]]=[tmp[9]]
	IN.close()
OUT=open(sys.argv[-1],mode='w')
for i in sorted(gene.keys()):
	OUT.write(i+'\t'+'\t'.join(gene[i])+'\n')
