import pandas as pd
import numpy as np
import sys

def transcript_to_gene_level_counts():
	counts_data_matrix = pd.read_table('/home/tonya/express_effective_counts_matix_rounded', sep=' |\t', index_col = 0, header = 0)
	# create a new column for gene name corresponding to transcript ID
	counts_data_matrix['gene_name'] = np.nan

	# makes a dictionary of transcript IDs as keys and matching gene as values
	conversion_dictionary = {}
	with open(sys.argv[1]) as transcript_conversions:
		for line in transcript_conversions:
			enst, gene = line.rstrip('\n').split('\t')
			conversion_dictionary[enst] = gene

	for row_name in counts_data_matrix.index:
		if row_name in conversion_dictionary:
			counts_data_matrix.ix[row_name, 'gene_name'] = conversion_dictionary[row_name]
		else:
			counts_data_matrix.ix[row_name, 'gene_name'] = str(row_name)


	counts_data_matrix.sort_values('gene_name', axis=0, inplace=True)
	counts_data_matrix.to_csv('effective_counts_matrix_with_gene_names.csv')
    
    # gene dataframe
    gene_df = pd.DataFrame()
    # number of genes
    genes = list(set(counts_data_matrix['gene_name'].tolist()))
    for x in genes:
    	concatenate_rows = counts_data_matrix.loc[counts_data_matrix['gene_name'] == x]
    	concatenate_rows = concatenate_rows.sum(axis=0)
    	# convert the sum to a dataframe
    	concatenate_rows_transpose = pd.DataFrame(concatenate_rows).transpose()
    	concatenate_rows_transpose.ix[:, 'gene_name'] = x
    	gene_df = pd.concat([gene_df, concatenate_rows_transpose])
    
    gene_df.to_csv('express_effective_counts_matrix_rounded_gene_level.csv')


if __name__ == '__main__':
	transcript_to_gene_level_counts()