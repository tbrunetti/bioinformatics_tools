import argparse
import sys
import os

def FPKM_to_TPM():
	pass;

def RPKM_to_TPM():
	pass;

def raw_to_RPKM():
	pass;

def raw_to_TPM():
	pass;

def raw_to_FPKM():
	pass;



if __name__=='__main__':
	parser=argparse.ArgumentParser(Description="Converts RNA seq quantification units")
	parser.add_argument('-input', required=True, dest='', help='list of full paths to all quantification files of samples that need to be calculated, one sample path per line')
	parser.add_argument('-method', required=True, dest='methodType', type=str, help='FPKM_TO_TPM, RPKM_TO_TPM', 'RAW_TO_RPKM', 'RAW_TO_FPKM', 'RAW_TO_TPM')
	parser.add_argument('-output', default=os.getcwd(), dest='outputResults', type=str, help='location of output directory for results.  Default=current working directory')
	args=parser.parse_args()

	if args.methodType=='FPKM_TO_TPM':
		FPKM_TO_TPM();
	elif args.methodType=='RPKM_TO_TPM':
		RPKM_TO_TPM();
