import os
import pickle

import pandas as pd


class Tumor:
    def __init__(self,
                 input_file,
                 db_name,
                 ascending,
                 base_dir,
                 custom_pw_db,
                 custom_pw_name
                 ):

        yes_or_no = {'y': True, 'n': False}
        """
            pathway: {percent: [samples]}
        """
        self.final_summary = {}

        """
            percent: [genes]
        """
        self.gene_summary = {}

        self.ascending = yes_or_no[ascending]
        # ascending = rank from bottom

        self.input_file = input_file
        self.db_name = db_name
        self.base_dir = base_dir

        if not custom_pw_db:
            self.db = pickle.load(open('{}/databases/{}.pkl'.format(self.base_dir, self.db_name), 'rb'))
        else:
            self.db = custom_pw_db
            self.db_name = custom_pw_name

        print("Creating output directory")
        self.output_dir = self.create_output_dir()

        print("Reading input file")
        self.gene_expression_table = self.read_file()

        self.all_genes = self.gene_expression_table.index

        self.final_summary = {
            'count': {},
            'db': {}
        }
        self.final_summary_rank = {
            'count': {},
            'db': {}
        }

        self.final_summary_enrichment = {
            'count': {},
            'db': {}
        }

    def create_output_dir(self):
        input_file_basename = os.path.basename(self.input_file)
        basename_wo_ext = '.'.join(input_file_basename.split(".")[:-1])
        if self.ascending:
            top_or_bottom = 'ascending'
        else:
            top_or_bottom = 'descending'

        output_dir = "{}/output_dir/{}_{}_{}".format(self.base_dir,
                                                                  basename_wo_ext,
                                                                  self.db_name,
                                                                  top_or_bottom)

        if not os.path.isdir(output_dir):
            os.mkdir(output_dir)

        return output_dir

    def read_file(self):
        gene_expression_table = pd.read_csv(self.input_file, header=0, index_col=0, sep="\t")
        return gene_expression_table.round(3)

    def add_pathway_to_final_summary(self, pw, sample, p_value, count, db):
        if sample not in self.final_summary:
            self.final_summary[sample] = {}

        self.final_summary[sample][pw] = p_value
        self.final_summary['count'][pw] = count
        self.final_summary['db'][pw] = db

    def add_pathway_to_final_summary_rank(self, pw, sample, rank, count, db):
        if sample not in self.final_summary_rank:
            self.final_summary_rank[sample] = {}

        self.final_summary_rank[sample][pw] = rank
        self.final_summary_rank['count'][pw] = count
        self.final_summary_rank['db'][pw] = db

    def add_pathway_to_final_summary_enrichment(self, pw, sample, enrichment, count, db):
        if sample not in self.final_summary_enrichment:
            self.final_summary_enrichment[sample] = {}

        self.final_summary_enrichment[sample][pw] = enrichment
        self.final_summary_enrichment['count'][pw] = count
        self.final_summary_enrichment['db'][pw] = db

    def write_final_summary_table(self):
        final_summary_df = pd.DataFrame(self.final_summary)
        output_summary = '{}/final_summary.txt'.format(self.output_dir)
        with open(output_summary, 'w') as tsvout:
            final_summary_df.to_csv(tsvout, sep="\t")

        print("finished writing summary: {}".format(output_summary))

    def write_final_summary_table_rank(self):
        final_summary_df = pd.DataFrame(self.final_summary_rank)
        output_summary = '{}/final_summary_rank.txt'.format(self.output_dir)
        with open(output_summary, 'w') as tsvout:
            final_summary_df.to_csv(tsvout, sep="\t")

        print("finished writing summary rank: {}".format(output_summary))

    def write_final_summary_table_enrichment(self):
        final_summary_df = pd.DataFrame(self.final_summary_enrichment)
        output_summary = '{}/final_summary_enrichment.txt'.format(self.output_dir)
        with open(output_summary, 'w') as tsvout:
            final_summary_df.to_csv(tsvout, sep="\t")

        print("finished writing summary enrichment: {}".format(output_summary))
