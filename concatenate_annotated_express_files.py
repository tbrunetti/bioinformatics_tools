import sys

def small_quant_matrix():
	
	final_file = open('VIP-concatenated-express-file.txt', 'w')
	gene_counts = {}

	for files in range(1, len(sys.argv)):
		with open(sys.argv[files]) as express_file:
			header = next(express_file)
			for line in express_file:
				line = line.split('\t')
				if line[0] in gene_counts:
					gene_counts[line[0]] = gene_counts[line[0]] + [str(int(float(line[9])))]
				else:
					gene_counts[line[0]] = [str(int(float(line[9])))]
		final_file.write('\t' + str(sys.argv[files]))

	for key in gene_counts:
		final_file.write('\n' + str(key) +'\t' + '\t'.join(gene_counts[key]))


if __name__ == '__main__':
	small_quant_matrix()