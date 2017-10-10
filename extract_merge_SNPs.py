import argparse
import subprocess
import time
import os
import pandas
import collections

def extract_merge_vcf(vcftools, plink, reference, chr, snps, out):
	merge_filename = open(os.path.join(out, 'chr_' + str(chr) + '_merge_these_files.txt'), 'w')
	submitted = [] # submitted jobs
	job_pool = [] # actively running jobs
	
	with open(snps) as pos:
		for line in pos:
			filename='chr_' + str(chr) + '_' + str(line)
			job_pool.append(subprocess.Popen([vcftools, '-gzvcf', reference, '--chr', str(chr), '--from-bp', str(line), '--to-bp', str(line), '--out', os.path.join(out, filename), '--plink']))
			submitted.append(os.path.join(out, str(filename)+'.ped'))
			merge_filename.write(str(os.path.join(out, str(filename)+'.ped')) + '\t' + str(os.path.join(out, str(filename)+'.map')) + '\n')

	merge_filename.flush() # flush out buffer		
	
	# wait for all jobs to finish running	
	while True:
		if all(os.path.isfile(str(completed_file)) for completed_file in submitted):
			break;
		else:
			print("Waiting for all extractions to finish...")
			time.sleep(60)

	print('All jobs finished')		


def extract_plink(plink, reference, chr, snps, out):
	pos_not_found = open(os.path.join(out, 'chr_' + str(chr) + '_positions_not_found.txt'), 'w')
	merge_filename = open(os.path.join(out, 'chr_' + str(chr) + '_merge_these_files.txt'), 'w')
	available_positions = []
	with open(reference[:-4]+'.bim') as get_pos:
		for line in get_pos:
			line = line.rstrip().split() # remove any right trailing whitespace and then split by whitespace
			available_positions.append(str(line[3]))


	with open(snps) as pos:
		for line in pos:
			line = line.rstrip()
			filename='chr_' + str(chr) + '_' + str(line)
			if str(line) in available_positions:
				subprocess.call([plink, '--bfile', reference[:-4], '--chr', str(chr), '--from-bp', str(line), '--to-bp', str(line), '--make-bed', '--out', str(os.path.join(out, filename))])
				merge_filename.write(str(os.path.join(out, str(filename)+'.bed')) + '\t' + str(os.path.join(out, str(filename)+'.bim')) + '\t' + str(os.path.join(out, str(filename)+'.fam')) + '\n')
			else:
				pos_not_found.write(str(line) + '\n')


	merge_filename.flush() # flush out buffer
	pos_not_found.flush() # flush out buffer
	merge_filename.close()
	pos_not_found.close()

	print('All jobs finished')		
	
	# grep through all log file for the following: Error: All variants excluded.


def merge_plinks(plink, ref, mergelist):
	pass


def duplicate_snps(plink, dup_file, reference, out):
	total_times_duplicated = open(os.path.join(out, 'total_times_duplicated.txt'), 'w')
	subprocess.call(['uniq', '-c', dup_file], stdout=total_times_duplicated)
	total_times_duplicated.close()
	merge_these = open(os.path.join(out, 'all_duplicates_extracted_merged_updated_illumina_snp_names.txt'), 'w')
	merge_list_plink = open(os.path.join(out, 'merge_list_input_plink.txt'), 'w')
	with open(total_times_duplicated.name) as snp_to_dup:
		for line in snp_to_dup:
			line = line.rstrip().split()
			for i in range(int(line[0])): #repeat this the number of times that it needs to be duplicated
				snp_names = line[3].split(',')
				subprocess.call([plink, '--bfile', reference[:-4], '--chr', str(line[1]), '--from-bp', str(line[2]), '--to-bp', str(line[2]), '--recode tab', '--out', str(os.path.join(out, 'chr_'+str(line[1])+'_pos_'+str(line[2])+'_dup_'+str(i)))])
				with open(str(os.path.join(out, 'chr_'+str(line[1])+'_pos_'+str(line[2])+'_dup_'+str(i))+'.map')) as map_file:
					for info in map_file:
						info = info.rstrip().split()
						info[1] = snp_names[i]
				updated_map = open(str(os.path.join(out, 'chr_'+str(line[1])+'_pos_'+str(line[2])+'_dup_'+str(i))+'.map'), 'w')
				updated_map.write('\t'.join(info))
				updated_map.close()
				merge_these.write(str(os.path.join(out, 'chr_'+str(line[1])+'_pos_'+str(line[2])+'_dup_'+str(i)+'.ped')) + '\t' + str(os.path.join(out, 'chr_'+str(line[1])+'_pos_'+str(line[2])+'_dup_'+str(i)+'.map')) +'\n') # create a merge list for plink
	
	merge_these.close()

	base_file = subprocess.check_output(['sed', '1q', merge_these.name])
	base_file = base_file.split()[0][:-4]
	subprocess.call(['sed', '1d', merge_these.name], stdout=merge_list_plink)
	merge_list_plink.close() # flushes buffer and closes file

	subprocess.call([plink, '--file', base_file, '--merge-list', merge_list_plink.name, '--recode', 'tab', '--out', str(os.path.join(out, 'all_duplicates_extracted_merged_updated_illumina_snp_names'))])


def generate_duplicate_files(merge_bim, out):
	duplicate_file = open(os.path.join(out, 'snps_to_duplicate_formatted.txt'), 'w')
	merged_bim_files = pandas.read_table(merge_bim, dtype=str)
	duplicate_snps = []
	total_snp_dups = collections.Counter(merged_bim_files['SNP_name_to_update']) # the columns with duplicate snps names of set that needs duplicating
	for snp,total in total_snp_dups.iteritems():
		if total > 1:
			duplicate_snps.append(snp)
	for dup in duplicate_snps:
		subset_dataframe = merged_bim_files.loc[merged_bim_files['SNP_name_to_update'] == str(dup), ['chrm', 'pos', 'ref_SNP_name']]
		for repeat in range(len(list(subset_dataframe['ref_SNP_name']))):
			duplicate_file.write(str(list(subset_dataframe['chrm'])[0]) + '\t' + str(list(subset_dataframe['pos'])[0]) + '\t' + ','.join(list(subset_dataframe['ref_SNP_name'])) + '\n')
	duplicate_file.flush()
	duplicate_file.close()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='extracts and merges SNPs in VCF or PLINK formats')
	parser.add_argument('--method', required=True, dest='method', type=str, help='method to use.  Options:  vcf_extract, plink_extract, merge_plink, duplicate_snps, generate_dup_file')
	parser.add_argument('--merge_bim_file', dest='merge_dups_bim', type=str, help='merged bim file of the two merge files of interest; required headers = chrm, pos, ref_SNP_name, SNP_name_to_update')
	parser.add_argument('--file', dest='dup_file', type=str, help='file of duplicates (1 snp per line) format is tab-delimted chr and position and comma-separated list of snps names , i.e. if a snps is seen twice two lines of identical chr, position and identical comma-separated list of snp renames should be seen')
	parser.add_argument('--plink', dest='plink_exe', type=str, help='full path to PLINK executable')
	parser.add_argument('--vcftools', dest='vcf_exe', type=str, help='full path to VCF executable')
	parser.add_argument('--ref', dest='reference', type=str, help='reference file to extract given snp list postion from(vcf.gz or .bed)')
	parser.add_argument('--SNPs', dest='snp_file', type=str, help='text file of SNP bp positions; one per line')
	parser.add_argument('--chr', dest='chr', type=str, help='chromosome number to extract from')
	parser.add_argument('--outDir', dest='outDir', type=str, default=os.getcwd(), help='full path directory to store output')
	args = parser.parse_args()

	if args.method == 'vcf_extract':
		extract_merge_vcf(vcftools=args.vcf_exe, plink=args.plink_exe, reference=args.reference, chr=args.chr, snps=args.snp_file, out=args.outDir)
	elif args.method == 'plink_extract':
		extract_plink(plink=args.plink_exe, reference=args.reference, chr=args.chr, snps=args.snp_file, out=args.outDir)
	elif args.method == 'merge_plink':
		merge_plinks(plink, ref, mergelist)
	elif args.method == 'duplicate_snps':
		duplicate_snps(plink=args.plink_exe, dup_file=args.dup_file, reference=args.reference, out=args.outDir)
	elif args.method == 'generate_dup_file':
		generate_duplicate_files(merge_bim=args.merge_dups_bim, out=args.outDir)
	else:
		print("Invalid method.  Please choose one of the following: vcf_extract, plink_extract, merge_plink")