import pandas as pd


class Sample:
    def __init__(self,
                 sample_id,
                 output_dir,
                 gene_expressions,
                 ascending):

        self.output_dir = output_dir
        self.sample_id = sample_id
        self.ascending = ascending

        self.output_file = '{}/{}.txt'.format(self.output_dir, self.sample_id)

        self.gene_expressions = gene_expressions.sort_values(ascending=self.ascending)

        self.genes_by_rank = self.rank_genes()

    # {gene: rank}
    def rank_genes(self):
        expression_level_ranks = {}
        rank = 0
        for gene, expression in self.gene_expressions.iteritems():
            rank += 1
            expression_level_ranks[gene] = {
                'rank': rank,
                'expression': expression
            }
        return expression_level_ranks

