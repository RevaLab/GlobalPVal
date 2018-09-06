import csv
import pickle
import sys

import numpy as np

from pathway import Pathway
from sample import Sample
from tumor import Tumor
from user_pathways import UserPathways

yes_or_no = {'y': True, 'n': False}

input_file = sys.argv[1]
db_name = sys.argv[2]
ascending = yes_or_no[sys.argv[3]]
base_dir = sys.argv[4]
custom_pw_file = sys.argv[5]
custom_pw_name = sys.argv[6]
custom_pw_exists = sys.argv[7]

# input_file = 'ucec.normal.txt'
# db_name = 'KEGG'
# ascending = 'n'
# base_dir = '/Users/anna/PycharmProjects/Pathway_Enrichment_Ranked_Expression'
# custom_pw_file = '0'
# custom_pw_name = '0'
# custom_pw_exists = 'n'

# Create custom pathway database
if custom_pw_file != '0':
    custom_pathways = UserPathways(custom_pw_file,
                                   db_name,
                                   base_dir,
                                   custom_pw_name,
                                   custom_pw_exists)
    custom_pw_db = custom_pathways.combined_user_pathways_and_db
else:
    custom_pw_db = False
    custom_pw_name = False

# Instantiate tumor
print("Instantiating tumor")
tumor = Tumor(
    input_file,
    db_name,
    ascending,
    base_dir,
    custom_pw_db,
    custom_pw_name
)
# Process samples in tumor
sample_count = 1

all_samples = len(tumor.gene_expression_table.columns)

for sample_id in tumor.gene_expression_table:
    sample = Sample(sample_id,
                    tumor.output_dir,
                    tumor.gene_expression_table[sample_id],
                    tumor.ascending)

    with open(sample.output_file, 'w') as tsvout:
        tsvout = csv.writer(tsvout, delimiter='\t')
        for pw in tumor.db['dict']:
            pw_data = tumor.db['dict'][pw]

            pathway = Pathway(
                pw_data['genes'],
                sample.genes_by_rank,
                tumor.all_genes,
                pw,
                sample_id
            )

            if np.isnan(pathway.geom_mean_p_vals):
                continue

            row = [pw,
                   pw_data['db'],
                   'Count:',
                   len(pw_data['genes']),
                   'PW_u_input:',
                   pathway.bg,
                   'Geom_mean_p_val:',
                   pathway.geom_mean_p_vals,
                   'Rank:',
                   pathway.rank
                   ]
            tsvout.writerow(row)
            tumor.add_pathway_to_final_summary(pw,
                                               sample_id,
                                               pathway.geom_mean_p_vals,
                                               len(pw_data['genes']),
                                               pw_data['db'])
            tumor.add_pathway_to_final_summary_rank(pw,
                                                    sample_id,
                                                    pathway.rank,
                                                    len(pw_data['genes']),
                                                    pw_data['db'])
            tumor.add_pathway_to_final_summary_enrichment(pw,
                                                          sample_id,
                                                          pathway.avg_enrichment,
                                                          len(pw_data['genes']),
                                                          pw_data['db'])
            tumor.add_pathway_to_final_summary_log(pw,
                                                   sample_id,
                                                   pathway.geom_mean_p_vals,
                                                   len(pw_data['genes']),
                                                   pw_data['db'])
    print("Finished with sample: {}/{}".format(sample_count, all_samples))
    sample_count += 1

tumor.write_final_summary_table()
tumor.write_final_summary_table_rank()
tumor.write_final_summary_table_enrichment()
tumor.write_final_summary_table_log()

tumor.find_unique_metabolites()
tumor.write_metabolite_summary()
tumor.record_metabolite_sample_summary()

# print(tumor.final_summary)