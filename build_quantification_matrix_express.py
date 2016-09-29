import os
import re
import operator

def build_matrix(express_dir_locations):
	# creates three files; each line transcript ID, each col is a BID with effcounts/tpm/fpkm data
	effCounts_data_matrix = open('express_effective_counts_matix', 'w')
	tpm_data_matrix = open('express_tpm_matrix', 'w')
	fpkm_data_matrix = open('express_fpkm_matrix', 'w')
	
	# regex to selectr proper express quant file
	quantFileName = re.compile("^e.*xprs")
	# regex to extract BID of sample
	bid = re.compile(r'(\d+-\d+)+')	
	# first key = transcriptID, 2nd key = BID, value = quant data
	transcript_quants = {}
	#stores BIDs of all samples with a quantifcation file
	all_sample_names = []
	os.chdir(express_dir_locations)
	for directories in os.listdir(express_dir_locations):
		os.chdir(directories)
		for files in os.listdir(express_dir_locations + directories):
			if quantFileName.match(files):
				print quantFileName.match(files)
				with open(files) as quantData:
					header = next(quantData)
					extract_bid = re.search(bid, directories).group(1)
					all_sample_names.append(extract_bid)
					
					for line in quantData:
						line = line.rstrip('\n').split('\t')
						transcript, eff_counts, eff_length, fpkm, tpm = line[1], float(line[7]), float(line[3]), float(line[10]), float(line[14])
						if transcript in transcript_quants:
							transcript_quants[transcript].update({extract_bid:[eff_counts, eff_length, fpkm, tpm]})
						else:
							transcript_quants[transcript]={extract_bid:[eff_counts, eff_length, fpkm, tpm]}
					os.chdir(express_dir_locations)
					break;
			
			else:
				continue;
	
	# adds headers to each file			
	effCounts_data_matrix.write('transcript'+'\t'+'\t'.join(all_sample_names)+'\n')
	tpm_data_matrix.write('transcript'+'\t'+'\t'.join(all_sample_names)+'\n')
	fpkm_data_matrix.write('transcript'+'\t'+'\t'.join(all_sample_names)+'\n')
	
	for key in transcript_quants:
		effCounts_data_matrix.write(str(key)+'\t')
		tpm_data_matrix.write(str(key)+'\t')
		fpkm_data_matrix.write(str(key)+'\t')

		for i in range(0, len(all_sample_names)):
			if transcript_quants[key][all_sample_names[i]] not in transcript_quants:
				effCounts_data_matrix.write('\t'.join(str(0)))
				tpm_data_matrix.write('\t'.join(str(0)))
				fpkm_data_matrix.write('\t'.join(str(0)))
			else:
				effCounts_data_matrix.write('\t'.join(transcript_quants[key][all_sample_names[i]][0]))
				tpm_data_matrix.write('\t'.join(transcript_quants[key][all_sample_names[i]][3]))
				fpkm_data_matrix.write('\t'.join(transcript_quants[key][all_sample_names[i]][2]))
		
		effCounts_data_matrix.write('\n')
		tpm_data_matrix.write('\n')
		fpkm_data_matrix.write('\n')
				

if __name__ == '__main__':
	express_dir_locations = '/lustre/beagle/tbrunetti/RNA_quantification/quantification_express/'
	build_matrix(express_dir_locations);