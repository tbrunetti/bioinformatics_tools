import sys
import argparse
import os

def convert_multi_to_single(outDir):
	
	with open(sys.argv[1]) as multi_line_fasta:
		for line in multi_line_fasta:
			if line[0] == ">":
				new_file = open(str(line.rstrip('\n')[1:]) + '.fa', 'w')
			else:
				new_file.write(line.rstrip('\n'))


if __name__ == "__main__":
	parser=argparse.ArgumentParser(description='')
	parser.add_argument('-outDir', default=os.getcwd(), dest='outputDir', help='Full path to output directory, default is current working directory')
	args=parser.parse_args()

	convert_multi_to_single(outDir=args.outputDir)