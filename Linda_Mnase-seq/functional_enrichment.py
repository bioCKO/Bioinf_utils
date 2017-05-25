### functional_enrichment.py
### Written by Yulia Mostovoy for 2012 CGRL RNA-Seq workshop
### 4/30/2012

import sys, argparse
import rpy2.robjects as robjects

def main():
    parser = argparse.ArgumentParser(description="Given a list of interesting genes, a list of all genes for which there was valid data in the study, and a file containing group definitions, tests whether the interesting genes are significantly over-represented in any of the groups.")
    parser.add_argument(metavar='Genes_of_interest', dest='genes_of_interest_filename', help='One gene name per line')
    parser.add_argument(metavar='List_of_all_genes', dest='all_genes_filename', help='One gene name per line')
    parser.add_argument(metavar='Groups_to_test', dest='group_filename', help='One group per line, format: Group name|group description [tab] list of genes (separated by tabs)')

    args = parser.parse_args()
    r = robjects.r

    ### open all input files ###
    genes_of_interest_file = openfile(args.genes_of_interest_filename, 'r')
    all_genes_file = openfile(args.all_genes_filename,'r')
    group_file = openfile(args.group_filename,'r')
    
    ### process list of all genes ###
    all_valid_genes = []
    for gene in all_genes_file:
        gene = gene.strip()
        all_valid_genes.append(gene)
    all_valid_genes = list(set(all_valid_genes)) # remove duplicates
    all_genes_file.close()

    ### process genes to test for enrichment ###
    genes_of_interest = []
    discarded = 0
    for gene in genes_of_interest_file:
        gene = gene.strip()
        if gene in all_valid_genes:
            genes_of_interest.append(gene)
        else:
            discarded += 1
    genes_of_interest = list(set(genes_of_interest)) # remove duplicates
    genes_of_interest_file.close()

    ### process group file ###
    groups = process_group_file(group_file, all_valid_genes)

    ### process and print enrichment ###
    enrichment(genes_of_interest, all_valid_genes, groups)

    sys.stderr.write("Genes in sample: %d\n" % len(genes_of_interest))
    sys.stderr.write("Genes in population: %d\n" % len(all_valid_genes))
    sys.stderr.write("Genes in sample rejected for not being present in population: %d\n" % discarded)
    sys.exit()

def enrichment(genes_of_interest, all_valid_genes, groups):
    r = robjects.r
    r("pvals = c()")
    group_list = []

    for group in groups:
        group_list.append(group)
        genes_of_interest_in_group = 0

        for gene_in_group in groups[group]["genes"]:
            if gene_in_group in genes_of_interest:
                genes_of_interest_in_group += 1

        groups[group]["found"] = genes_of_interest_in_group

        # calculate values for matrix
        genes_of_interest_not_in_group = len(genes_of_interest) - genes_of_interest_in_group
        genes_in_group_remainder = len(groups[group]["genes"]) - genes_of_interest_in_group
        genes_remainder = len(all_valid_genes) - len(groups[group]["genes"]) - genes_of_interest_not_in_group

        # get Fisher's exact test p-value for enrichment of this group, append to vector of p-values stored in R
        r("table = matrix(c(%d,%d,%d,%d), nrow=2, ncol=2)" %
           (genes_of_interest_in_group, genes_in_group_remainder, genes_of_interest_not_in_group, genes_remainder))
        r("pval = fisher.test(table, alternative='g')$p.value")
        r("pvals = c(pvals, pval)")

    # multiple testing correction of p-values
    pvals_raw = r("pvals")
    pvals_adj = r("p.adjust(pvals, method='bonferroni')")

    joined_data = zip(pvals_raw, pvals_adj, group_list)
    joined_data.sort()

    # print results
    print "Group\tFrequency\tpval_raw\tpval_adjusted\tDescription"
    
    for pval_raw, pval_adj, group in joined_data:
        print "%s\t%d/%d (%.2f%%)\t%e\t%e\t%s" % (group, groups[group]["found"], len(groups[group]["genes"]), 100.0*groups[group]["found"]/len(groups[group]["genes"]), 
               pval_raw, pval_adj, groups[group]["description"])
    return


def process_group_file(group_file, all_valid_genes):
    final_groups = {}
    for group in group_file:
        line = group
        group = group.strip().split('\t')
        if len(group)==1:
            sys.stderr.write("Error: Line in group file unable to be split - check format of group file!\n")
            sys.stderr.write("Problematic line: %s" % line)
            sys.exit()
        temp_group = {"description":'', "genes":[]}
        group_name = group.pop(0)
        if len(group_name.split('|')) > 1:
            temp_group['description'] = group_name.split('|')[1].strip()
            group_name = group_name.split('|')[0].strip()
        else:
            temp_group['description'] = ''
        genes = group
        for gene in genes:
            if gene in all_valid_genes:
                temp_group['genes'].append(gene)
        if len(temp_group['genes']) >= 5:
            final_groups[group_name] = temp_group
    return final_groups


def openfile(filename, mode):
    try:
        newfile = open(filename, mode)
    except IOError:
        sys.stderr.write("Unable to open file in mode %s: %s\n" % (mode, filename))
        sys.exit()

    return newfile



########################
if __name__=="__main__": main()
