import re
import sys
IN=open(sys.argv[1])
OUT_map=open(sys.argv[2],mode='w')
OUT_unmap1=''
OUT_unmap2=''
if len(sys.argv)>3:
	OUT_unmap1=open(sys.argv[3],mode='w')
	OUT_unmap2=open(sys.argv[4],mode='w')
for i in IN:
	if re.search(r'^@',i):
		continue
	else:
		i2=IN.next()
		tmp=i.split('\t')
		tmp2=i2.split('\t')
		if tmp[1]=='83' or tmp[1]=='99':
			if re.search(r'XS:i:',tmp[12]) and re.search(r'XS:i:',tmp2[12]):
				AS1=int(re.search(r'AS:i:(\S+)',tmp[11]).groups()[0])
				XS1=int(re.search(r'XS:i:(\S+)',tmp[12]).groups()[0])
				AS2=int(re.search(r'AS:i:(\S+)',tmp2[11]).groups()[0])
				XS2=int(re.search(r'XS:i:(\S+)',tmp2[12]).groups()[0])
				if AS1>(XS1+4) and AS2>(XS2+4) and AS1>-12 and AS2>-12:
					OUT_map.write(i+i2)	
				continue
			else:
				AS1=int(re.search(r'AS:i:(\S+)',tmp[11]).groups()[0])
				AS2=int(re.search(r'AS:i:(\S+)',tmp2[11]).groups()[0])
				if AS1>-12 and AS2>-12:
					OUT_map.write(i+i2)
		else:
			if OUT_unmap1 and OUT_unmap2:
				OUT_unmap1.write('>'+tmp[0]+'\n'+tmp[9]+'\n')
				OUT_unmap2.write('>'+tmp2[0]+'\n'+tmp2[9]+'\n')
