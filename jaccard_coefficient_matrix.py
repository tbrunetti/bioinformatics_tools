import sys
import numpy as np
from sklearn.metrics import jaccard_similarity_score
import pandas as pd

def calculate_jaccard():
	binaryArray = []
	samples = []
	with open(sys.argv[1]) as binary_matrix:
		for line in binary_matrix:
			line = line.strip('\n').split('\t')
			samples.append(line[0])
			line.pop(0)
			binaryArray.append(line)

	binaryArray = np.array(binaryArray)
	similarity_matrix = np.zeros((len(samples), len(samples)))
	
	for i in range(0, len(samples)):
		for j in range(0, len(samples)):
			similarity_matrix[i, j] = jaccard_similarity_score(binaryArray[i], binaryArray[j], normalize=True)		
	
	labeled_sim_matrix = pd.DataFrame(similarity_matrix, index = samples, columns = samples)
		
	labeled_sim_matrix.to_csv('/home/tonya/Desktop/jaccard_similarity_score_ATAC_seq_data_89.csv', doublequote=False)
if __name__ == '__main__':
	calculate_jaccard();