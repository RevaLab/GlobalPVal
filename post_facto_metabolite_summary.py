import csv
import os
import pickle

import pandas as pd


class AnalyzedTumor:
    def __init__(self,
                 input_file):

        self.input_file = input_file
        self.output_dir = os.path.dirname(self.input_file)
        if 'KEGG' in self.output_dir:
            self.db_name = 'KEGG'
        elif 'reactome' in self.output_dir:
            self.db_name = 'reactome'
        elif 'HMDB' in self.output_dir:
            self.db_name = 'HMDB_SMPDB'
        self.final_summary = self.read_file()
        self.metabolite_significant_pathways = {}
        self.find_unique_metabolites()
        self.write_metabolite_summary()

    def read_file(self):
        pathways_with_p_vals = pd.read_csv(self.input_file, header=0, index_col=0, sep="\t")
        pathways = pathways_with_p_vals.to_dict()
        return pathways

    def find_unique_metabolites(self):
        significant_pathways_samples = {}
        for sample in self.final_summary:
            if sample == 'count' or sample == 'db':
                continue
            for pathway in self.final_summary[sample]:
                if self.final_summary[sample][pathway] < 0.01:
                    if pathway not in significant_pathways_samples:
                        significant_pathways_samples[pathway] = set()
                    significant_pathways_samples[pathway].add(sample)

        significant_pathways = significant_pathways_samples.keys()
        all_metabolites = pickle.load(open('databases/{}_metabolites.pkl'.format(self.db_name), 'rb'))

        metabolite_significant_pathways = {}

        for metabolite in all_metabolites:
            significant_pathways_with_this_metabolite = all_metabolites[metabolite].intersection(significant_pathways)
            if len(significant_pathways_with_this_metabolite):
                if metabolite not in metabolite_significant_pathways:
                    metabolite_significant_pathways[metabolite] = {
                        'pathways': set(),
                        'samples': set()
                    }

                metabolite_significant_pathways[metabolite]['pathways'] = \
                    metabolite_significant_pathways[metabolite]['pathways'].union(
                        significant_pathways_with_this_metabolite)

                for pathway in significant_pathways_with_this_metabolite:
                    samples_where_this_pathway_is_significant = significant_pathways_samples[pathway]
                    metabolite_significant_pathways[metabolite]['samples'] = \
                        metabolite_significant_pathways[metabolite]['samples'].union(
                            samples_where_this_pathway_is_significant)
        self.metabolite_significant_pathways = metabolite_significant_pathways

    def write_metabolite_summary(self):
        output_summary = '{}/metabolite_summary.txt'.format(self.output_dir)
        with open(output_summary, 'w') as tsvout:
            tsvout = csv.writer(tsvout, delimiter='\t')
            for metabolite in self.metabolite_significant_pathways:
                row = [metabolite,
                       self.db_name,
                       'Samples:',
                       len(self.metabolite_significant_pathways[metabolite]['samples']),
                       'Pathways:',
                       len(self.metabolite_significant_pathways[metabolite]['pathways'])]
                row += self.metabolite_significant_pathways[metabolite]['pathways']
                tsvout.writerow(row)
        print("finished writing metabolite summary: {}".format(output_summary))


if __name__ == '__main__':
    for subdir, dirs, files in os.walk('output_dir'):
        for file in files:
            # print os.path.join(subdir, file)
            filepath = subdir + os.sep + file

            if filepath.endswith("final_summary.txt"):
                analyzed_tumor = AnalyzedTumor(filepath)
