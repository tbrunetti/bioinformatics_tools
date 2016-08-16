import os
import sys

def concatenateFiles():
	f=open('bamStats_RNAseq.txt', 'a')
	os.chdir(sys.argv[1])
	for files in os.listdir(sys.argv[1]):
		with open(files) as input:
			storage={}
			for line in input:
				line=line.rstrip('\n')
				line=line.split(':')
				if len(line)>1:
					storage[line[0].lstrip()]=line[1].lstrip()
			f.write(str(files)+'\t')
		for key in storage:	
			f.write(storage[key]+'\t')
		f.write('\n')

	print storage





if __name__=='__main__':
	concatenateFiles();