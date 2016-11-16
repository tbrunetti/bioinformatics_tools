import sys
import operator
import math
import matplotlib.pyplot as plt
import numpy as np

def parse_capture():
	regions = {}
	for x in range(1, len(sys.argv)):
		with open(sys.argv[x]) as input:
			for line in input:
				chrm, start, end, desc, normal, primary, met, adj_normal, adj_primary, adj_met, ratio_primary_normal, ratio_met_normal, ratio_met_primary = line.rstrip('\n').split('\t')
				if chrm in regions:
					regions[chrm] = regions[chrm] + [[int(start), int(end), desc, float(ratio_primary_normal),float(ratio_met_normal), float(ratio_met_primary)]]
				else:
					regions[chrm] = [[int(start), int(end), desc, float(ratio_primary_normal), float(ratio_met_normal), float(ratio_met_primary)]]


	for key in regions_sorted:
		chrm_regions = []
		ratios = []
		for i in range(0,len(regions_sorted[key])):
			expand_regions = list(range(regions_sorted[key][i][0], regions_sorted[key][i][1] +1))
			expand_ratios = [regions_sorted[key][i][3] for x in range(0, len(expand_regions))]
			chrm_regions = chrm_regions + expand_regions
			ratios = ratios + expand_ratios
		 
		all_ratios = all_ratios + ratios
       	log_ratiios = [math.log(ratios[x], 2) if ratios[x] !=0 else 0 for x in range(0, len(ratios))]
        #plt.scatter(chrm_regions, log_ratios)
       	#plt.xlabel(str(key), fontsize=18)
       	#plt.ylabel('log2(ratio(primary/normal))', fontsize=16)
       	#plt.axhline(linewidth=2, color='r')
        #plt.show()

        #plt.scatter(chrm_regions, ratios)
        #plt.xlabel(str(key), fontsize=18)
        #plt.ylabel('ratio(primary/normal)', fontsize=16)
        #plt.axhline(y=1, linewidth=2, color='r')
        #plt.show()

        print np.mean(all_ratios)
        print np.std(all_ratios)
        print min(all_ratios)
        print max(all_ratios)


		log_ratios = [math.log(ratios[x], 2) for x in range(0, len(ratios))]
		plt.scatter(chrm_regions, log_ratios)
		plt.show()



if __name__ == '__main__':
	parse_capture()