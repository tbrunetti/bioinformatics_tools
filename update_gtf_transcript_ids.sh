#!/bin/bash

new_gtf=test.gtf
orig_gtf=GCF_000006765.1_ASM676v1_genomic_CDS_replaced_with_exon_confirmed_no_overlap.gtf

while read line; do
	if [[ $(echo "${line}" | cut -f3) == "exon" ]] || [[ $(echo "${line}" | cut -f3) == "start_codon" ]] || [[ $(echo "${line}" | cut -f3) == "stop_codon" ]] ; then
		preDetails=$(echo "${line}" | awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}')
		details=$(echo "${line}" | awk -F"\t" '{print $9}')
		gene_id=$(echo "${details}" | sed -n 's/.*gene_id "\([^"]*\)".*/\1/p')
		gene_name=$(echo "${details}" | sed -n 's/.*gene "\([^"]*\)".*/\1/p')
		transcript_id=$(echo ${details} | sed -n 's/.*transcript_id "\([^"]*\)".*/\1/p')
		length=${#gene_name}
		if [ "${length}" -gt 0 ]; then
			new_tx_id=$(echo "${gene_id}_${gene_name}")
			updated_details=$(echo ${details} | sed "s/${transcript_id}/${new_tx_id}/") # double quotes in sed expands the variable unlike single quotes treated as literals
			echo -e "${preDetails}\t${updated_details}" >> ${new_gtf}
		else
			updated_details=$(echo ${details} | sed "s/${transcript_id}/${gene_id}/")
			echo -e "${preDetails}\t${updated_details}" >> ${new_gtf}
		fi
	else
		echo "${line}" >> ${new_gtf}
	fi
done < ${orig_gtf}

