import argparse
import os
import typing

'''
function: get_top_n(inFile, outDir, prefixName, topN)
arguments:
    inFile: a tab-delimited file of all the super enhancers concatenated across samples and sorted by "stitchedPeakRank" column (column 9). Input data should not contain any headers;
            col 1 = REGION_ID
            col 2 = CHROM
            col 3 = START
            col 4 = STOP    
            col 5 = NUM_LOCI        
            col 6 = CONSTITUENT_SIZE        
            col 7 = ranking bam density
            col 8 = control bam density (optional)  -- if no control bam is used remaining columns move up
            col 9 = stitchedPeakRank  [col 8 if no control bam is used]      
            col 10 = isSuper [col 9 if no control bam is used]
    outDir: a path of where to write data 
    prefixName: a string of a prefix to prepend to the output files
    topN: an integer specifying how many of the top enhancers to include in the output
Output
    bed file of the topN super enhancers, or all super enhancers possible if < topN
'''

def get_top_n(inFile:str, outDir:str, prefixName:str, topN:int) -> None:
   
    import pandas

    topEnhancers = {}
    
    # assumes the input file is already sorted by column 8 (no control) or column 9 (with control) in ascending order
    with open(inFile) as rankedEnhancers:
        for line in rankedEnhancers:
            if len(topEnhancers.keys()) != topN:
                enhancer, chrom ,start, stop, *_ = line.split('\t') # split to only get bed components
                if enhancer in topEnhancers:
                    print('Duplicate enhancer found ... skipping ...')
                else:
                    topEnhancers[enhancer] = {'chrom':chrom, 'start':start, 'stop':stop}
            else:
                print('Top {} enhancers have been identified!'.format(str(topN)))
      
    # if file closes before top N is reached, this if statement is exectuted to let user know that topN > # unique enhancers present in dataset
    if len(topEnhancers.keys()) != topN:
        print('{} super enhancers identified. Your data had fewer than {} super enhancers'.format(str(len(topEnhancers.keys())), str(topN)))

    # convert to bed format from dict of dict
    topEnhancersDf = pandas.DataFrame.from_dict(topEnhancers, orient='index')
    topEnhancersDf = topEnhancersDf.rename_axis('enhancer').reset_index()
    topEnhancersDf = topEnhancersDf[['chrom', 'start', 'stop', 'enhancer']]
    topEnhancersDf.to_csv(os.path.join(outDir, '{}_top{}_super_enhancers.bed'.format(prefixName, str(topN))), sep = '\t', index=False, header=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='ROSE super enhancer functions',description='A series of functions to help manipulate files and output from ROSE super enhancer software')
    parser.add_argument('-i', '--input', type=str, help = "Full path to input file")
    parser.add_argument('-o', '--output', type = str, default = os.getcwd(), help = "Full path to output directory")
    parser.add_argument('-p', '--prefix', type = str, help = " string specifying a name/prefix to prepend to all output files")
    #parser.add_argument('-c', '--controlUsed', type = bool, action = 'store_true', help = "specify this flag at runtime if a control bam was used to remove background density")
    parser.add_argument('-n', '--topN', type = int, default = 1000, help = "An integer specifying how many top n unique super enhancers to select")
    args = parser.parse_args()

    get_top_n(inFile = args.input, outDir = args.output, prefixName = args.prefix, topN = args.topN)