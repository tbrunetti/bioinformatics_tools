import pandas
import time
import os
import sys
import argparse
import logging
from memory_profiler import profile

@profile(precision=4)
def prepData(genoFile, genoFileType, inputTable, chrm, LDrange, outdir, config):
	logger = logging.getLogger('prepData')
	logging.debug('Running function prepData()')

	allAssocData = pandas.read_csv(inputTable, use_cols = ['CHROM', 'BEGIN', 'END', 'PVALUE'], skipInitialSpace = True, sep='\t', 
		dtype={'CHROM':'str', 'BEGIN':'int64', 'END':'int64', 'PVALUE':'float64'})

	tableSubset = allAssocData.loc[(allAssocData['CHROM'] == chrm) & (allAssocData['BEGIN'] >= int(LDrange.split('-').rstrip()[0]))
					& (allAssocData['END'] <= int(LDrange.split('-').rstrip()))]
	tableSubset['MARKER_ID'] = tableSubset['CHROM'].astype(str) + ':' +tableSubset['BEGIN'].astype(str)
	logging.info('Finshed subsetting summary statistics data...preparing to clear memory to free RAM')
	
	
	# START FREE MEMORY OF LARGE DATAFRAME
	allAssocData = pandas.DataFrame()
	freeMem = [allAssocData]
	del allAssocData
	del freeMem
	# END FREE MEMORY OF LARGE DATAFRAME

	
	if genoFileType == 'vcf':
		pass;
	elif genoFileType == 'plink':
		logging.debug('Subsetting plink file to match tableSubet')
		subprocess.call([config['plink'],
			'--bfile', genoFile[:-4],
			'--chr', chrm,
			'--from-bp', LDrange.split('-').rstrip()[0],
			'--to-bp', LDrange.split('-').rstrip()[1],
			'--make-bed',
			'--out', os.path.join(outdir, '{}_{}_{}'.format(genoFile[:-4], LDrange.split('-').rstrip()[0], LDrange.split('-').rstrip()[1]))
		])

	try:
		assert os.path.isfile(outdir, '{}_{}_{}.bed'.format(genoFile[:-4], LDrange.split('-').rstrip()[0], LDrange.split('-').rstrip()[1]))
	except AssertionError:
		logging.critical('Unable to generate plink file {}_{}_{}.bed.  System Exiting'.format(genoFile[:-4], LDrange.split('-').rstrip()[0], LDrange.split('-').rstrip()[1]))
		print('Critical Error Reached!  Unable to subset data from plink file.  System Exiting')
		sys.exit()

	logging.info('Writing results to {}'.format(os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile[:-4], LDrange.split('-').rstrip()[0], LDrange.split('-').rstrip()[1]))))
	tableSubset.to_csv(os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile[:-4], LDrange.split('-').rstrip()[0], LDrange.split('-').rstrip()[1])),
			index=False, sep='\t')


	findRefSnp(tableSubset=tableSubset, outdir=outdir)


def findRefSnp(tableSubset, outdir):
	logger = logging.getLogger('findRefSnp')
	logging.debug('Running function findRefSnp()')

	refSnp = tableSubset.loc[tableSubset['PVALUE'].idxmax()]['MARKER_ID']
	print('The reference SNP for LD in this subset is {}'.format(refSnp))


	calculateLD()

def calculateLD(chrm, LDrange, genoFileSubset, outdir, config):
	logger = logging.getLogger('calculateLD')
	logging.debug('Running function calculateLD()')

	pass;

def runLZ(outdir, config):
	logger = logging.getLogger('runLZ')
	logging.debug('Running function runLZ()')

	pass;


if __name__ == '__main__':
	import json

	parser = argparse.ArgumentParser(Description='Find reference snps and calculates LD for genomic regions-output \
		can be used as direct input into locuszoom')

	parser.add_argument('--config', required=True, type=str, help='Full path to generate_LD_config.json')
	parser.add_argument('--chrm', type=str, help='Chromosome to calculate -- must be present in dataset')
	parser.add_argument('--LDrange', type=str, help='base pair position range expressed as int-int inclusive. Example:\
		500-999 would equate to base pair 500bp through 999bp including both 500 and 999')
	parser.add_argument('--inputTable', type=str, help='Tab-delimited file with following headers: CHROM, BEGIN, END, PVALUE.  MARKER_ID is optional\
		but will be used if available, any additional column names will be ignored; CHROM, BEGIN, END should be integers and PVALUE a float')
	parser.add_argument('--outdir', type=str, default=os.getcwd(), help='Output directory to write results')

	genotypes = parser.add_mutually_exclusive_group(required=True)
	genotypes.add_argument('--vcfGeno', type=str, default=None, help="Full path to vcf file to use for genotype extraction")
	genotypes.add_argument('--plinkGeno', type=str, default=None, help="Full path to binary plink file (ends in .bed) for genotype extraction")

	parser.parse_args()
	

	logging.basicConfig(filename='generateLD.log', level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	logger = logging.getLogger('main_runtime')

	logging.info('Checking table format and headers')
	
	try:
		readTableHeader = open(args.inputTable, 'r')
		header = readTableHeader.readlines().split('\t')
		for names in ['CHROM', 'BEGIN', 'END', 'PVALUE']:
			assert names in header
		logging.info('confirmed, {} found in header of input table'.format(names))
	
	except AssertionError:
		logging.critical('FAIL, {} not found in header of input table; Program exiting.'.format(names))
		sys.exit(1)
	
	# read and store json config file as dictionary
	paths = json.loads(config)[0]


	try:
		assert os.path.exists(args.outdir)
	except AssertionError:
		logging.info('{} does not exist, creating new directory...'.format(args.outdir))
		os.mkdir(args.outdir)



	if args.vcfGeno != None:
		logging.info("Input genotype file is dectected as VCF")
		prepData(genoFile=args.vcfGeno, genoFileType='vcf', inputTable=args.inputTable, chrm=args.chrm, LDrange=args.LDrange, outdir=args.outdir, config=paths)
	elif args.plinkGeno != None:
		logging.info("Input genotype file is detected as plink format")
		prepData(genoFile=args.plinkGeno, genoFileType='plink', inputTable=args.inputTable, chrm=args.chrm, LDrange=args.LDrange, outdir=args.outdir, config=paths)
	else:
		pass




