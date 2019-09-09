import pandas
import time
import os
import sys
import argparse
import logging
from memory_profiler import profile
import subprocess

@profile(precision=4)
def prepData(genoFile, genoFileType, inputTable, chrm, LDrange, outdir, config):
	logger = logging.getLogger('prepData')
	logging.debug('Running function prepData()')

	allAssocData = pandas.read_csv(inputTable, usecols = ['CHROM', 'BEGIN', 'END', 'PVALUE'], skipinitialspace = True, sep='\t', 
		dtype={'CHROM':'str', 'BEGIN':'int64', 'END':'int64', 'PVALUE':'float64'})

	tableSubset = allAssocData.loc[(allAssocData['CHROM'] == chrm) & (allAssocData['BEGIN'] >= int(LDrange.split('-')[0].rstrip()))
					& (allAssocData['END'] <= int(LDrange.split('-')[1].rstrip()))]
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
			'--from-bp', LDrange.split('-')[0].rstrip(),
			'--to-bp', LDrange.split('-')[1].rstrip(),
			'--make-bed',
			'--out', os.path.join(outdir, '{}_{}_{}'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip()))
		])

	try:
		assert os.path.isfile(os.path.join(outdir, '{}_{}_{}.bed'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip())))
	except AssertionError:
		logging.critical('Unable to generate plink file {}_{}_{}.bed.  System Exiting'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip()))
		print('Critical Error Reached!  Unable to subset data from plink file.  System Exiting')
		sys.exit()

	logging.info('Writing results to {}'.format(os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip()))))
	tableSubset.to_csv(os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip())),
			index=False, sep='\t')


	findRefSnp(tableSubset=tableSubset, data=os.path.join(outdir, '{}_{}_{}'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip())), outdir=outdir)


@profile(precision=4)
def findRefSnp(tableSubset, data, outdir):
	logger = logging.getLogger('findRefSnp')
	logging.debug('Running function findRefSnp()')

	refSnp = tableSubset.loc[tableSubset['PVALUE'].idxmin()]['MARKER_ID']
	pval = tableSubset.loc[tableSubset['PVALUE'].idxmin()]['PVALUE']
	print('The reference SNP for LD in this subset is {} and has a p-value of {}'.format(refSnp, pval))
	logging.info('The reference SNP for LD in this subset is {} and has a p-value of {}'.format(refSnp, pval))

	# TO DO: if refSNP is more than one in same position take the first index (doesn't matter because they will both be in LD)


	#TODO: populate and trim arguments
	calculateLD(data=data, refSnp=refSnp, outdir=outdir, config=config)


@profile(precision=4)
def calculateLD(data, refSnp, outdir, config):
	logger = logging.getLogger('calculateLD')
	logging.debug('Running function calculateLD()')

	subprocess.call([
		config['plink'],
		'--bfile', data,
		'--r2', 'dprime',
		'--ld-snp', refSnp,
		'--ld-window-kb', str(1000),
		'--ld-window', str(99999),
		'--ld-window-r2', str(0),
		'--out', os.path.join(outdir, '{}_plinkCalculatedLD'.format(data))
		]
	)


	#TODO: change header names to match expectation of LZ
	plinkLD = pandas.read_csv(os.path.join(outdir, '{}_plinkCalculatedLD.ld'.format(data)), usecols=['SNP_B', 'SNP_A', 'DP', 'R2'], skipinitialspace=True,  delim_whitespace=True)

	lzLDformat = plinkLD.rename(columns={'SNP_B':'snp1', 'SNP_A':'snp2', 'DP':'dprime', 'R2':'rsquare'})

	lzLDformat.to_csv(os.path.join(outdir, '{}_plinkCalculatedLD_convertedForLZinput.ld'.format(data)), sep="\t", index=False)
	


@profile(precision=4)
def runLZ(outdir, config):
	logger = logging.getLogger('runLZ')
	logging.debug('Running function runLZ()')

	pass;


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Find reference snps and calculates LD for genomic regions-output \
		can be used as direct input into locuszoom')

	parser.add_argument('--config', required=True, type=str, help='Full path to configuration of software file (generate_LD_config.cfg)')
	parser.add_argument('--chrm', type=str, help='Chromosome to calculate -- must be present in dataset')
	parser.add_argument('--LDrange', type=str, help='base pair position range expressed as int-int inclusive. Example:\
		500-999 would equate to base pair 500bp through 999bp including both 500 and 999')
	parser.add_argument('--inputTable', type=str, help='Tab-delimited file with following headers: CHROM, BEGIN, END, PVALUE.  MARKER_ID is optional\
		but will be used if available, any additional column names will be ignored; CHROM, BEGIN, END should be integers and PVALUE a float')
	parser.add_argument('--outdir', type=str, default=os.getcwd(), help='Output directory to write results')

	genotypes = parser.add_mutually_exclusive_group(required=True)
	genotypes.add_argument('--vcfGeno', type=str, default=None, help="Full path to vcf file to use for genotype extraction (imputed if possible)")
	genotypes.add_argument('--plinkGeno', type=str, default=None, help="Full path to binary plink file (ends in .bed) for genotype extraction (imputed if possible)")

	args = parser.parse_args()
	

	logging.basicConfig(filename='generateLD.log', level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	logger = logging.getLogger('main_runtime')

	logging.debug('Checking table format and headers')
	
	try:
		readTableHeader = open(args.inputTable, 'r')
		header = readTableHeader.readline().rstrip().split('\t')
		for names in ['CHROM', 'BEGIN', 'END', 'PVALUE']:
			assert names in header
			logging.info('confirmed, {} found in header of input table'.format(names))
	
	except AssertionError:
		logging.critical('FAIL, {} not found in header of input table; Program exiting.'.format(names))
		sys.exit(42)
	
	# read and store config file as dictionary
	config = {}
	with open(args.config) as confFile:
		for line in confFile:
			(key, value) = line.split(':')
			config[key.strip()] = value.strip()

	try:
		assert os.path.exists(args.outdir)
	except AssertionError:
		logging.info('{} does not exist, creating new directory...'.format(args.outdir))
		os.mkdir(args.outdir)



	if args.vcfGeno != None:
		logging.info("Input genotype file is dectected as VCF")
		prepData(genoFile=args.vcfGeno, genoFileType='vcf', inputTable=args.inputTable, chrm=args.chrm, LDrange=args.LDrange, outdir=args.outdir, config=config)
	elif args.plinkGeno != None:
		logging.info("Input genotype file is detected as plink format")
		prepData(genoFile=args.plinkGeno, genoFileType='plink', inputTable=args.inputTable, chrm=args.chrm, LDrange=args.LDrange, outdir=args.outdir, config=config)
	else:
		pass




