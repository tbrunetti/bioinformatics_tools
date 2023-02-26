import pandas
import re
import csv

def geneid_mapping_table(gtf):

    gene_map = {}

    with open(gtf) as annotationFile:
        for line in annotationFile:
            if line.startswith("#"):
                print('Skipping header...')
            elif (line.split("\t")[2] == 'exon') or (line.split('\t')[2] == "CDS"):
                geneId = re.match(r".*(gene_id.*?;).*", line.split("\t")[-1])
                #print(geneId.group(1))
                if geneId.group(1).split("\"")[1] in gene_map.keys():
                    print("Already seen {}".format(geneId.group(1).split("\"")[1]))
                else:
                    tmp = {}
                    tmp['feature'] = line.split("\t")[2]
                    tmp['strand'] = line.split("\t")[6]
                    gene = re.match(r".*(gene\s.*?;).*", line.split("\t")[-1])
                    transcriptId = re.match(r".*(transcript_id.*?;).*", line.split("\t")[-1])
                    note = re.match(r".*(note.*?;).*", line.split("\t")[-1])
                    product = re.match(r".*(product.*?;).*", line.split("\t")[-1])
                    try:
                        tmp['gene'] = gene.group(1).split("\"")[1]
                    except:
                        tmp['gene'] = ""

                    try: 
                        tmp['transcript_id'] = transcriptId.group(1).split("\"")[1]
                    except AttributeError:
                        tmp['transcript_id'] = ""

                    try:
                        tmp['note'] = note.group(1).split("\"")[1]
                    except:
                        tmp['note'] = ''
                        
                    try:
                        tmp['product'] = product.group(1).split("\"")[1]
                    except:
                        tmp['product'] = ''
                    gene_map[geneId.group(1).split("\"")[1]] = tmp
    
    gtf_mapping = pandas.DataFrame.from_dict(gene_map, orient = 'index')
    gtf_mapping.to_csv("cVibrio_cholerae_Genbank_GCA_003097695.1_ASM309769v1_geneId_mappings.tsv", index = True, sep = "\t")


def vibrio_cholerae_strain_A1552():
    mappingTable = pandas.read_table("/home/tonya/pluto_projects/clients/Pukatzki_lab/vibrio_cholerae_strain_A1552/ncbi-genomes-2023-02-25/genbank/GCA_003097695.1_ASM309769v1/Vibrio_cholerae_Genbank_GCA_003097695.1_ASM309769v1_geneId_mappings_cleaned_final.tsv", sep = "\t", index_col = None)
    gtf = pandas.read_table("/home/tonya/pluto_projects/clients/Pukatzki_lab/vibrio_cholerae_strain_A1552/ncbi-genomes-2023-02-25/genbank/GCA_003097695.1_ASM309769v1/GCA_003097695.1_ASM309769v1_genomic.gtf", skiprows=3, header=None, index_col=None)

    mappingTable['find_gene_id'] = "gene_id \"" + mappingTable['gene_id'] + "\""
    mappingTable['find_transcript_id'] = "transcript_id \"" + mappingTable['transcript_id'] + "\""
    mappingTable['replace_gene_id'] = "gene_id \"" + mappingTable['updated_gene_id'] + "\""
    mappingTable['replace_trascript_id'] = "transcript_id \"" + mappingTable['updated_transcript_id'] + "\""

    update_gene_id_pairings = dict(zip(mappingTable['find_gene_id'], mappingTable['replace_gene_id']))
    update_transcript_id_pairings = dict(zip(mappingTable['find_transcript_id'], mappingTable['replace_trascript_id']))

    gtf.replace({8:update_gene_id_pairings},regex = True, inplace=True)
    gtf.replace({8:update_transcript_id_pairings},regex = True, inplace=True)
    gtf.replace({2:{"CDS":"exon"}}, inplace = True)

    gtf.to_csv("/home/tonya/pluto_projects/clients/Pukatzki_lab/vibrio_cholerae_strain_A1552/ncbi-genomes-2023-02-25/genbank/GCA_003097695.1_ASM309769v1/Vibrio_cholerae_Genbank_GCA_003097695.1_ASM309769v1_updated.gtf", index = False, sep = "\t", header=False, quoting=csv.QUOTE_NONE, quotechar="")
