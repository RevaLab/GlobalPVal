import scipy.stats as stats
import numpy as np


class Pathway:
    def __init__(self,
                 genes,
                 sample_genes_by_rank,
                 all_genes,
                 pw,
                 sample_id):

        self.pw = pw
        self.sample_id = sample_id

        self.genes = genes

        self.all_genes = set(all_genes)
        self.duplicates = len(all_genes) - len(self.all_genes)
        self.bg = len(set(all_genes).union(self.genes)) + self.duplicates
        self.sample_genes_by_rank = sample_genes_by_rank
        # self.pathway_count = len(self.genes)
        self.pathway_count = len(self.genes.intersection(set(self.sample_genes_by_rank)))
        self.pathway_gene_ranks = self.rank_pathway_genes()

        self.all_p_vals, \
        self.ranks_by_p_val, \
        self.p_val_reciprocals, \
        self.ranks_by_enrichment, \
        self.all_enrichments = self.calculate_all_p_vals()

        self.natural_logs = []

        # summary statistics
        self.geom_mean_p_vals = self.geo_mean_overflow(self.all_p_vals)
        self.rank = self.calculate_rank_given_p_vals()
        self.avg_enrichment = self.calculate_average_enrichment()

    def rank_pathway_genes(self):
        pathway_gene_ranks = []

        for gene in self.genes:
            try:
                pathway_gene_ranks.append(self.sample_genes_by_rank[gene]['rank'])
            except KeyError:
                pass

        pathway_gene_ranks.sort()

        return pathway_gene_ranks

    def calculate_all_p_vals(self):
        all_p_vals = []
        ranks_by_p_val = []  # (Rm * (1/pm))
        p_val_reciprocals = []
        ranks_by_enrichment = []
        enrichment_reciprocals = []
        all_enrichments = []
        for i, rank in enumerate(self.pathway_gene_ranks):

            overlap = i + 1
            rank_only = rank - overlap
            pathway_only = self.pathway_count - overlap
            bg_only = self.bg - rank - self.pathway_count + overlap
            p_value = self.calculate_p_val_for_each_gene_in_pw(overlap,
                                                               rank_only,
                                                               pathway_only,
                                                               bg_only)

            enrichment_coefficient = self.find_enrichment_coefficient(overlap,
                                                                      self.pathway_count,
                                                                      rank,
                                                                      self.bg)

            rank_by_p_val = rank / p_value
            ranks_by_p_val.append(rank_by_p_val)

            rank_by_enrichment = enrichment_coefficient / p_value
            ranks_by_enrichment.append(rank_by_enrichment)

            p_val_reciprocal = 1/p_value
            p_val_reciprocals.append(p_val_reciprocal)

            all_p_vals.append(p_value)
            all_enrichments.append(enrichment_coefficient)

        return all_p_vals, \
               ranks_by_p_val, \
               p_val_reciprocals, \
               ranks_by_enrichment, \
               all_enrichments

    def calculate_rank_given_p_vals(self):
        a = np.sum(self.ranks_by_p_val)
        b = np.sum(self.p_val_reciprocals)
        rank = a/b
        return rank

    def calculate_average_enrichment(self):
        a = np.sum(self.ranks_by_enrichment)
        b = np.sum(self.p_val_reciprocals)
        rank = a/b
        return rank

    # @staticmethod
    def calculate_p_val_for_each_gene_in_pw(self,
                                            overlap,
                                            rank_only,
                                            pathway_only,
                                            bg_only):
        try:
            odds_ratio, pvalue = stats.fisher_exact(
                [
                    [overlap, pathway_only],
                    [rank_only, bg_only]
                ],
                alternative='greater'
            )
            return pvalue
        except ValueError:
            print('Something is wrong with this 2x2 table: \n')
            print([overlap, pathway_only], [rank_only, bg_only])
            return 'NaN'

    @staticmethod
    def geo_mean_overflow(iterable):
        a = np.log(iterable)
        return np.exp(a.sum() / len(a))

    @staticmethod
    def find_enrichment_coefficient(overlap, pathway, rank, bg):
        enrichment_threshold = (pathway * rank) / bg

        coefficient = overlap / enrichment_threshold
        #
        # if coefficient < 1:
        #     return 1 / coefficient * -1
        return coefficient
