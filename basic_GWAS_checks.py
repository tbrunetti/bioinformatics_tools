import pandas as pd
import numpy as np
import argparse
import subprocess

def check_AT_GC_snps(bimFile):
	out_file = open('at_gc_snps_to_remove.txt', 'w')
	
	bim = pd.read_table(bimFile, delim_whitespace=True, header=None, names=['chr', 'snp', 'gen_dist', 'pos', 'allele1', 'allele2'], dtype=str)
	at_snps = bim.loc[((bim['allele1'] == 'A') & (bim['allele2'] == 'T')) | ((bim['allele1'] == 'T') & (bim['allele2'] == 'A'))]
	gc_snps = bim.loc[((bim['allele1'] == 'G') & (bim['allele2'] == 'C')) | ((bim['allele1'] == 'C') & (bim['allele2'] == 'G'))]
	snps_to_remove = list(at_snps['snp']) + list(gc_snps['snp'])
	#print at_snps
	#print gc_snps
	
	for snps in snps_to_remove:
		out_file.write(snps + '\n')

	'''
	REMOVE SNPS???
	'''

	out_file.close()

def pre_qc_imputation_filtering(bimFile, plink, vcftools, bgzip):
	subprocess.call([plink, '--bfile', bimFile[:-4], '--maf', '0.0001', '--geno', '0.05', '--hwe', '0.0001', '--make-bed', '--out', bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC'])
	for chrm in range(1, 23):
		vcf_out = open(bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC_chr_'+str(chrm)+'_sorted.vcf.gz', 'w')
		subprocess.call([plink, '--bfile', bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC', '--chr', str(chrm), '--make-bed', '--out', bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC_chr_'+str(chrm)])
		subprocess.call([plink, '--bfile', bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC_chr_'+str(chrm), '--recode', 'vcf', '--out', bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC_chr_'+str(chrm)+'_vcf'])
		sort = subprocess.Popen((vcftools, bimFile[:-4]+'_prefiltered_for_input_into_imputation_QC_chr_'+str(chrm)+'_vcf.vcf'), stdout=subprocess.PIPE)
 		compress = subprocess.Popen((bgzip, '-c'), stdin=sort.stdout, stdout=vcf_out)
 		sort.wait()
 		compress.wait()
 		vcf_out.close()

def pre_imputation_filtering_post_QC(bimFile, plink, vcftools, bgzip, flips, removeSnps):
	subprocess.call([plink, '--bfile', bimFile[:-4], '--flip', flips, '--exclude', removeSnps, '--make-bed', '--out', bimFile[:-4]+'_postQC_strand_flipped'])
	for chrm in range(1, 23):
		vcf_out = open(bimFile[:-4]+'_postQC_strand_flipped_chr_'+str(chrm)+'_sorted.vcf.gz', 'w')
		subprocess.call([plink, '--bfile', bimFile[:-4]+'_postQC_strand_flipped', '--chr', str(chrm), '--make-bed', '--out', bimFile[:-4]+'_postQC_strand_flipped_chr_'+str(chrm)])
		subprocess.call([plink, '--bfile', bimFile[:-4]+'_postQC_strand_flipped_chr_'+str(chrm), '--recode', 'vcf', '--out', bimFile[:-4]+'_postQC_strand_flipped_chr_'+str(chrm)+'_vcf'])
		sort = subprocess.Popen((vcftools, bimFile[:-4]+'_postQC_strand_flipped_chr_'+str(chrm)+'_vcf.vcf'), stdout=subprocess.PIPE)
 		compress = subprocess.Popen((bgzip, '-c'), stdin=sort.stdout, stdout=vcf_out)
 		sort.wait()
 		compress.wait()
 		vcf_out.close()	





def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('--bimFile', dest='bim', type=str, help='Full path to PLINK bim file')
	parser.add_argument('--ATGCsnps', dest='remove_snps', action='store_true', help='This flag means to return a list of snp names that have AT or GC allele pairs')
	parser.add_argument('--preQCImputeFilter', dest='filter', action='store_true', help = 'This flag means to call PLINK and filter PLINK files based on pre-imputation filtering guidelines')
	parser.add_argument('--postQCFilter', dest='postQC', action='store_true', help = 'This flag means to call PLINK and flip strands that fail pre impute QC and then compress back into sorted VCF for input into imputation')
	parser.add_argument('--plink', dest='plink_exe', type=str, help = 'Full path to PLINK executable')
	parser.add_argument('--vcftools', dest='vcf_exe', type=str, help='Full path to vcftools sort executable')
	parser.add_argument('--bgzip', dest='bgzip_exe', type=str, help='Full path to bgzip executable')
	parser.add_argument('--flips', dest='flips', type=str, help='Full path to list of snps that need to be flipped')
	parser.add_argument('--remove', dest='remove', type=str, help='Full path to list of snps that should be removed post imputation QC')
	args = parser.parse_args()

	if args.remove_snps == True:
		check_AT_GC_snps(bimFile = args.bim)
	if args.filter == True:
		pre_qc_imputation_filtering(bimFile = args.bim, plink=args.plink_exe, vcftools=args.vcf_exe, bgzip=args.bgzip_exe)
	if args.postQC == True:
		pre_imputation_filtering_post_QC(bimFile = args.bim, plink=args.plink_exe, vcftools=args.vcf_exe, bgzip=args.bgzip_exe, flips=args.flips, removeSnps=args.removea)		


if __name__ == '__main__':
	main()


