from chunkypipes.components import *

class Pipeline(BasePipeline):
	def dependencies(self):
		return []

	def description(self):
		return {'Calculates reads in each gene region and for each intergenic region in a sample'}

	def configure(self):
		return {
			'samtools':{
				'path':'Full path to samtools executable'
			},
			'bedtools':{
				'path':'Full path to bedtools coverageBed executable'
			}
		}

	def add_pipeline_args(self, parser):
		parser.add_argument('-inputBAM', required=True, help="Full path to BAM file to be analyzed")
		parser.add_argument('-genic', required=True, help="Full path to genic regions in BED format")
		parser.add_argument('-intergenic', required=True, help="Full path to intergenic regions in BED format")
		parser.add_argument('-minOverlap', default='0.51', help="[FLOAT] between 0 and 1 that is the mininum percent a read must overlap region to be considered a hit")
	

	def run_pipeline(self, pipeline_args, pipeline_config):
		genicFile=open(pipeline_args['inputBAM'][:-4]+'_genic_overlapping_reads.txt', 'w')
		intergenicFile=open(pipeline_args['inputBAM'][:-4]+'_intergenic_overlapping_reads.txt', 'w')
		
		# software initiation
		samtools_sort=Software('samtools', pipeline_config['samtools']['path']+' sort')
		bedtools_cov=Software('bedtools', pipeline_config['bedtools']['path'])

		# sorts bam file
		samtools_sort.run(
			Parameter(pipeline_args['inputBAM']),
			Parameter(pipeline_args['inputBAM'][:-4]+'.sorted')
			)
		
		# calculates number of reads in genic region
		bedtools_cov.run(
					Parameter('-a', pipeline_args['genic']),
					Parameter('-b', pipeline_args['inputBAM'][:-4]+'.sorted.bam'),
					Parameter('-F', pipeline_args['minOverlap']),
					Redirect(stream=Redirect.STDOUT, dest=pipeline_args['inputBAM'][:-4]+'_genic_overlapping_reads.txt')
				)
	
		# calculates number of reads in intergenic regions	
		bedtools_cov.run(
					Parameter('-a', pipeline_args['intergenic']),
					Parameter('-b', pipeline_args['inputBAM'][:-4]+'.sorted.bam'),
					Parameter('-F', pipeline_args['minOverlap']),
					Redirect(stream=Redirect.STDOUT, dest=pipeline_args['inputBAM'][:-4]+'_intergenic_overlapping_reads.txt')
				)
		# removes sorted bam file to save space
		subprocess.call(['rm', '-rf', pipeline_args['inputBAM'][:-4]+'.sorted.bam'])
