import pandas
import time
import os
import sys
import argparse
import logging
from memory_profiler import profile
import subprocess
import datetime

@profile(precision=4)
def prepData(genoFile, genoFileType, inputTable, chrm, LDrange, outdir, config):
	logger = logging.getLogger('prepData')
	logger.setLevel(logging.DEBUG)

	logger.debug('Running function prepData()')

	allAssocData = pandas.read_csv(inputTable, usecols = ['CHROM', 'BEGIN', 'END', 'PVALUE'], skipinitialspace = True, sep='\t', 
		dtype={'CHROM':'str', 'BEGIN':'int64', 'END':'int64', 'PVALUE':'float64'})

	tableSubset = allAssocData.loc[(allAssocData['CHROM'] == chrm) & (allAssocData['BEGIN'] >= int(LDrange.split('-')[0].rstrip()))
					& (allAssocData['END'] <= int(LDrange.split('-')[1].rstrip()))]
	tableSubset['MARKER_ID'] = tableSubset['CHROM'].astype(str) + ':' +tableSubset['BEGIN'].astype(str)
	tableSubset.rename(columns={'CHROM':'#CHROM'}, inplace=True) # for epacts formatting requirement for locuszoom
	logger.debug('Finshed subsetting summary statistics data...preparing to clear memory to free RAM')
	
	
	# START FREE MEMORY OF LARGE DATAFRAME
	allAssocData = pandas.DataFrame()
	freeMem = [allAssocData]
	del allAssocData
	del freeMem
	# END FREE MEMORY OF LARGE DATAFRAME

	

	if genoFileType == 'vcf':
		pass;
	elif genoFileType == 'plink':
		logger.debug('Subsetting plink file to match tableSubet')
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
		logger.critical('Unable to generate plink file {}_{}_{}.bed.  System Exiting'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip()))
		print('Critical Error Reached!  Unable to subset data from plink file.  System Exiting')
		sys.exit()

	logger.info('Writing results to {}'.format(os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip()))))
	tableSubset.to_csv(os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip())),
			index=False, sep='\t')

	findRefSnp(tableSubset=tableSubset, pathToTableSubset=os.path.join(outdir, 'subsetData_{}_{}_{}.txt'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip())), 
		data=os.path.join(outdir, '{}_{}_{}'.format(genoFile.split('/')[-1][:-4], LDrange.split('-')[0].rstrip(), LDrange.split('-')[1].rstrip())), outdir=outdir)


@profile(precision=4)
def findRefSnp(tableSubset, pathToTableSubset, data, outdir):
	logger = logging.getLogger('findRefSnp')
	logger.setLevel(logging.DEBUG)
	logger.debug('Running function findRefSnp()')

	refSnp = tableSubset.loc[tableSubset['PVALUE'].idxmin()]['MARKER_ID']
	pval = tableSubset.loc[tableSubset['PVALUE'].idxmin()]['PVALUE']
	print('The reference SNP for LD in this subset is {} and has a p-value of {}'.format(refSnp, pval))
	logger.info('The reference SNP for LD in this subset is {} and has a p-value of {}'.format(refSnp, pval))

	# TO DO: if refSNP is more than one in same position take the first index (doesn't matter because they will both be in LD)


	#TODO: populate and trim arguments
	calculateLD(pathToTableSubset=pathToTableSubset, data=data, refSnp=refSnp, outdir=outdir, config=config)


@profile(precision=4)
def calculateLD(pathToTableSubset, data, refSnp, outdir, config):
	logger = logging.getLogger('calculateLD')
	logger.setLevel(logging.DEBUG)
	logger.debug('Running function calculateLD()')

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

	logger.debug("Finished generating LD file from plink")


	#change header names to match expectation of LZ
	logger.debug("Updating headers to match locuszoom expected input...")
	plinkLD = pandas.read_csv(os.path.join(outdir, '{}_plinkCalculatedLD.ld'.format(data)), usecols=['SNP_B', 'SNP_A', 'DP', 'R2'], skipinitialspace=True,  delim_whitespace=True)

	lzLDformat = plinkLD.rename(columns={'SNP_B':'snp1', 'SNP_A':'snp2', 'DP':'dprime', 'R2':'rsquare'})

	lzLDformat.to_csv(os.path.join(outdir, '{}_plinkCalculatedLD_convertedForLZinput.ld'.format(data)), sep="\t", index=False)
	logger.debug("Finsihed updating headers for locuszoom compatability")

	#TODO: call runLZ()
	calcLDflagStatus = True
	# user input refSnp takes precendence over inferred from calculation refSnp
	if args.refSnp is None:
		logger.info('Using reference snp {} calculated from calculateLD() function'.format(refSnp))
		runLZ(epacts=pathToTableSubset, refSnp=refSnp, ldFile=os.path.join(outdir, '{}_plinkCalculatedLD_convertedForLZinput.ld'.format(data)), 
			calcLDflag=calcLDflagStatus, build=args.build, flank=args.flank, refGene=args.refgene, specificRegion=args.specificRegion, 
			outdir=outdir, config=config)
	else:
		logger.warning('Using reference snp {} specified by user at runtime using --refSnp argument.  User runtime input takes precedence over inferred refSnp from calculateLD() function'.format(args.refSnp))
		runLZ(epacts=pathToTableSubset, refSnp=args.refSnp, ldFile=os.path.join(outdir, '{}_plinkCalculatedLD_convertedForLZinput.ld'.format(data)), 
			calcLDflag=calcLDflagStatus, build=args.build, flank=args.flank, refGene=args.refgene, specificRegion=args.specificRegion, 
			outdir=outdir, config=config)


@profile(precision=4)
def runLZ(epacts, refSnp, ldFile, calcLDflag, build, flank, refGene, specificRegion, outdir, config):
	logger = logging.getLogger('runLZ')
	logger.setLevel(logging.DEBUG)
	logger.debug('Running function runLZ()')


	if calcLDflag == True:
		# infer from calcuated LD file
		if (flank is None) and (specificRegion is None):
			logger.info("Using LD generated from program and inferring, chromosome, start, and end arguments based on LD calulation inputs at runtime")
			subprocess.check_call([
				config['locuszoom'],
				'--epacts', epacts,
				'--build', build,
				'--ld', ldFile,
				'--refsnp', refSnp,
				'--chr', args.chrm,
				'--start', args.LDrange.split('-')[0].rstrip(),
				'--end', args.LDrange.split('-')[1].rstrip(),
				'--prefix', 'chr{}_{}-{}_refSnp{}_{}_{}'.format(args.chrm, args.LDrange.split('-')[0].rstrip(), args.LDrange.split('-')[1].rstrip(), refSnp, build, args.cohort)
				]
			)

		# use specific region specificed at runtime
		elif (flank is None) and (specificRegion is not None):
			logger.info("Using LD generated from program and chromosome, start, and end arguments are determined by --specificRegion argument at runtime")
			subprocess.call([
				config['locuszoom'],
				'--epacts', epacts,
				'--build', build,
				'--ld', ldFile,
				'--refsnp', refSnp,
				'--chr', specificRegion.split(':')[0].rstrip(),
				'--start', specificRegion.split(':')[1].rstrip().split('-')[0].rstrip(),
				'--end', specificRegion.split(':')[1].rstrip().split('-')[1].rstrip(),
				'--prefix', 'chr{}_{}-{}_refSnp{}_{}_{}'.format(specificRegion.split(':')[0].rstrip(), specificRegion.split(':')[1].rstrip().split('-')[0].rstrip(),
				 args.LDrange.split('-')[1].rstrip(), specificRegion.split(':')[1].rstrip().split('-')[1].rstrip(), refSnp, build, args.cohort),
				'--plotonly'
				]
			)
		
		# use flank region specified at runtime
		elif (flank is not None) and (specificRegion is None):
			logger.info("Using LD generated from program and region is determined by --flank argument at runtime")
			subprocess.call([
				config['locuszoom'],
				'--epacts', epacts,
				'--build', build,
				'--ld', ldFile,
				'--refsnp', refSnp,
				'--flank', flank,
				'--prefix', 'chr{}_flank{}_refSnp{}_{}_{}'.format(args.chrm, flank, refSnp, build, args.cohort),
				'--plotonly'
				]
			)

		else:
			print("Not working...")

	#TODO
	elif calcLDflag == False:
		pass;


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Find reference snps and calculates LD for genomic regions-output \
		can be used as direct input into locuszoom')

	parser.add_argument('--config', required=True, type=str, help='Full path to configuration of software file (generate_LD_config.cfg)')
	parser.add_argument('--chrm', type=str, help='Chromosome to calculate -- must be present in dataset')
	parser.add_argument('--LDrange', type=str, help='base pair position range expressed as int-int inclusive. Example:\
		500-999 would equate to base pair 500bp through 999bp including both 500 and 999')
	parser.add_argument('--inputTable', type=str, help='Tab-delimited file with following headers: CHROM, BEGIN, END, PVALUE.  MARKER_ID is optional\
		but will be used if available, any additional column names will be ignored; CHROM, BEGIN, END should be integers and PVALUE a float, if using standalone flag --lzOnly,\
		also mark header beginning with # and no space following #')
	parser.add_argument('--outdir', type=str, default=os.getcwd(), help='Output directory to write results')
	parser.add_argument('--build', type=str, default='hg19', choices=['hg19', 'hg18', 'hg38'], help="[default: hg19] Human genome build positions are based")
	parser.add_argument('--cohort', type=str, default=datetime.datetime.now().strftime("%Y-%m-%d_%H:%M"), help='[default:date/time stamp] string defining a project name or population for naming output files')
	

	reference = parser.add_mutually_exclusive_group()
	reference.add_argument('--refSnp', type=str, default=None, required='--lzOnly' in sys.argv, help="Name of snp to use as a reference; input take prescedence over inferred refsnp when using full program to caluclate LD")
	reference.add_argument('--refgene', type=str, default=None, help="Plot all data that falls in the listed reference gene; example: using MEF2 here would plot all snps that fall with MEF2 gene")


	parser.add_argument('--lzOnly', action='store_true', help="This flag tells the program to only run locuszoom and not calculate a reference snp")	
	ldParser = parser.add_mutually_exclusive_group()
	ldParser.add_argument('--ldFile', type=str, default=None, required='--lzOnly' in sys.argv, help="User-supplied LD if you want locuszoom(LZ) only, must follow LZ ld file format")
	ldParser.add_argument('--referenceLD', type=str, default=None, required='--lzOnly' in sys.argv, 
		choices=['hg19,ASN,1000G_March2012', 'hg19,AFR,1000G_March2012', 'hg19,EUR,1000G_March2012', 'hg19,AMR,1000G_March2012', 
		'hg19,ASN,1000G_Nov2010', 'hg19,AFR,1000G_Nov2010', 'hg19,EUR,1000G_Nov2010', 'hg18,CEU,1000G_June2010', 'hg18,YRI,1000G_June2010',
		'hg18,JPT+CHB,1000G_June2010', 'hg18,CEU,1000G_Aug2009', 'hg18,YRI,1000G_Aug2009', 'hg18,JPT+CHB,1000G_Aug2009'],
		help="choose one of choices")

	

	genotypes = parser.add_mutually_exclusive_group(required=True)
	genotypes.add_argument('--vcfGeno', type=str, default=None, help="Full path to vcf file to use for genotype extraction (imputed if possible)")
	genotypes.add_argument('--plinkGeno', type=str, default=None, help="Full path to binary plink file (ends in .bed) for genotype extraction (imputed if possible)")


	locuszoomRegions = parser.add_mutually_exclusive_group()
	locuszoomRegions.add_argument('--specificRegion', type=str, default=None, help="Region to plot for locuszoom, example: 1:100-200 would plot chromome 1 starting at base pair 100 and ending at base pair 200")
	locuszoomRegions.add_argument('--flank', type=str, default=None, help="Total base pairs in kb to flank the reference snp; example: 500kb would  plot snps +/- 500kb from refsnp, no unit, MB and kb are accetable units")
	args = parser.parse_args()
	



	logging.basicConfig(filename='generateLD.log', level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	logger = logging.getLogger('main_runtime')
	logger.setLevel(logging.DEBUG)

	'''
	logger.debug("The current working directory is {} and user has specified output to be located at {}".format(os.getcwd(), args.outdir))
	try:
		os.chdir(args.outdir)
	except OSError:
		logger.critical("Error trying to change current working directory to {}.  Please check directory exists or that you have permission to access the directory".format(args.outdir))
		sys.exit(42)
	'''

	logger.debug('Checking table format and headers')
	
	try:
		readTableHeader = open(args.inputTable, 'r')
		header = readTableHeader.readline().rstrip().split('\t')
		for names in ['CHROM', 'BEGIN', 'END', 'PVALUE']:
			assert names in header
			logger.info('confirmed, {} found in header of input table'.format(names))
	
	except AssertionError:
		logger.critical('FAIL, {} not found in header of input table; Program exiting.'.format(names))
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
		logger.info('{} does not exist, creating new directory...'.format(args.outdir))
		os.mkdir(args.outdir)



	if args.vcfGeno != None and args.lzOnly == False:
		logger.info("Input genotype file is dectected as VCF")
		logger.warning("Human genome build is set to {}".format(args.build))
		prepData(genoFile=args.vcfGeno, genoFileType='vcf', inputTable=args.inputTable, chrm=args.chrm, LDrange=args.LDrange, outdir=args.outdir, config=config)
	elif args.plinkGeno != None and args.lzOnly == False:
		logger.info("Input genotype file is detected as plink format")
		logger.warning("Human genome build is set to {}".format(args.build))
		prepData(genoFile=args.plinkGeno, genoFileType='plink', inputTable=args.inputTable, chrm=args.chrm, LDrange=args.LDrange, outdir=args.outdir, config=config)
	elif args.plinkGeno != None and args.lzOnly:
		logger.info("Input genotype file is detected as plink format")
		logger.warning("Human genome build is set to {}".format(args.build))
		runLZ(genoFile=args.plinkGeno, refSnp=args.refsnp, ldFile=args.ldFile, calcLDflag=False, build=args.build, outdir=args.outdir, config=config)
	else:
		logger.critical("Doesn't match an option")
		sys.exit(42)





# Warning: The "test" version of EPACTS changed the format of the output. To make LZ work, you'll also need to add --epacts-beg-col BEG to your command line. 