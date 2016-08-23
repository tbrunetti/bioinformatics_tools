import re
import os

def parse_flagstat(pathToFlagstats):
	#exact everything up to the first space
	extract_first_value = re.compile(r"([^\s]+)")
	extract_key = re.compile(r"^(\d+\s\W+\s\d+\s)+")  # regex to extract header info
	remove_parenthesis = re.compile(r"(\((.+)\))") # regex to header parenthesis specific to each sample
	os.chdir(pathToFlagstats)

	for files in os.listdir(pathToFlagstats):
		extract_bid = re.compile(r"^(\d+-\d+)+") # regex to extract BID
		bid=re.search(extract_bid, files).group(1)
		
		with open(files) as flagstat:
			headerStorage = ['BID']
			stats = [str(bid)]
			
			for line in flagstat:
				value = re.search(extract_first_value, line).group(1)
				header = re.sub(extract_key, "", line).strip('\n')
				newHeader = re.search(remove_parenthesis, header)
				
				# if None, means parenthesis not found and moves to else
				if newHeader != None:
					header = re.sub(remove_parenthesis,"", header).strip('\n') # subs regex for parenthesis with empty string
					headerStorage.append(header)
					stats.append(value)
				else:
					headerStorage.append(header)
					stats.append(value)
				
		# output all samples with flagstats in tab-delimited format
		with open(str(pathToOutput)+'ATAC-seq-unique-flagstats.txt', 'a+') as output:
			
			if os.stat(str(pathToOutput)+'ATAC-seq-unique-flagstats.txt').st_size == 0:
				output.write('\t'.join(headerStorage)+'\n')
				output.write('\t'.join(stats)+'\n')
			else:
				output.write('\t'.join(stats)+'\n')

if __name__=='__main__':
	pathToFlagstats = '/home/tonya/Desktop/fake_directory/'
	pathToOutput = '/home/tonya/Desktop/'
	parse_flagstat(pathToFlagstats);