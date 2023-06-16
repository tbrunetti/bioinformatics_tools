import argparse
import typing

def add_transcript(gtf:str, out:str) -> None:
    import re

    pattern = re.compile('.*(transcript_id ".*?");.*')
    exclusion_pattern = re.compile('(.*)exon_number \"[0-9]*\";.*')
    
    header = ''
    transcript_data = {}
    
    with open(gtf) as gtf_parse:
        for line in gtf_parse:
            if line.split('\t')[2] == 'exon':
                print(line.split('\t'))
                metadata = line.split('\t')[8]
                match = pattern.match(metadata)
                transcript_id = match.group(1).split("\"")[1]
                strand = line.split('\t')[6]
                chromosome = line.split('\t')[0]
                start_exon = int(line.split('\t')[3])
                end_exon = int(line.split('\t')[4])
                ref_db = line.split('\t')[1]
                if transcript_id in transcript_data:
                    if start_exon < transcript_data[transcript_id]['start_exon']:
                        transcript_data[transcript_id]['start_exon'] = start_exon
                    if end_exon > transcript_data[transcript_id]['end_exon']:
                        transcript_data[transcript_id]['end_exon'] = end_exon
                else:
                    exclude_exon_info = exclusion_pattern.match(metadata)
                    transcript_data[transcript_id] = {'strand':strand, 
                    'chromosome':chromosome,
                    'ref_db':ref_db,
                    'start_exon':start_exon,
                    'end_exon':end_exon,
                    'metadata':exclude_exon_info.group(1)}

            elif line.startswith('#'):
                header = header + line.strip() + '\n'
    
    with open(out, 'w') as gtf_write:
        for transcripts in transcript_data:
            gtf_write.write(transcript_data[transcripts]['chromosome'] +
            '\t' + transcript_data[transcripts]['ref_db'] + 
            '\ttranscript\t' + str(transcript_data[transcripts]['start_exon']) +
            '\t' + str(transcript_data[transcripts]['end_exon']) + '\t.\t' +
            transcript_data[transcripts]['strand'] + '\t.\t' + 
            transcript_data[transcripts]['metadata'] + '\n'
            )





if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='GTF modification',description='functions to modify gtf files')
    parser.add_argument('-g', '--gtf', type=str, help = "Full path to gtf file")
    parser.add_argument('-o', '--out', type = str, help = "Full path and name of final gtf to generate")
    args = parser.parse_args()

    add_transcript(gtf = args.gtf, out = args.out)