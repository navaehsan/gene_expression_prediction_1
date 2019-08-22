# Prediction of gene expression
This script could be used to predict gene expression using allelic fold change (AFC). The method is based on Mohammadi, P., Castel, S. E., Brown, A. A., & Lappalainen, T. (2017). Quantifying the regulatory effect size of cis-acting genetic variation using allelic fold change. Genome research, 27(11), 1872-1884 


# aFC_Whole_Blood.txt
This is a tab delimited file providing allelic fold change for each eQTL. In this file for each variant in columnn "variant_id" we have the following information:

gene_id : gene associated to that variant

rank : rank of that variant in eQTL calling process

CHROM : the chromosome of the variant

POS : position of that variant

REF : The Reference allele

ALT : The Alternative allele

log2_aFC : The effect size of the variant

# gene_expression_estimation_function.ipynb

This notebooks provide two functions in R that will read the AFC file and predict the gene expression.

## AFC_gene_vector function:

Reading the AFC file "aFC_Whole_Blood.txt" and saving it in a dataframe named AFC_dt, This function will retrieve the AFC vector for a gene which will be used in the next function to predict the expression of a gene.

## gene_expression_estimation function:
This function predicts the gene expression using :

h1 : genotype of variants for the first haplotype

h2 : genotype of variants for the second haplotype

s : vector of effect size of variants

returns a vector : 

result[1] : log expression for haplotype 1 in log2 scale

result[2] : gene expression for haplotype 2 in log2 scale

result[3] : total expression in log2 scale

for using this function a phased vcf file is required to extract the genotypes of individuals.

In VCF file if genotype information is present, then the same types of data must be present for all samples. First a FORMAT   field is given specifying the data types and order (colon-separated alphanumeric String). This is followed by one field per sample, with the colon-separated data in this  field corresponding to the types specified in the format. The  first sub-field must always be the genotype (GT)

GT : genotype, encoded as allele values separated by either of / or |. " /" means genotype unphased and "|" means genotype phased. The allele values are 0 for the reference allele (what is in the REF  field), 1 for the  
allele listed in ALT. For diploid calls examples could be 0/1, 1|0. 

**The REF and ALT information should match the REF and ALT information in aFC_Whole_Blood.txt**

For example if a gene has three eQTLs naming chr6_143498271_ATTGAACAAAGTCC_A_b38, chr6_143501317_C_G_b38, chr6_143504621_C_T_b38 and vcf file has the following information for an individual:

ID  \t                                 individual_1

chr6_143498271_ATTGAACAAAGTCC_A_b38   0|1

chr6_143501317_C_G_b38                0|0

chr6_143504621_C_T_b38                1|0

a sample genotype for that gene would be : 

genotype_h1 = c(0,0,1)

genotype_h2 = c(1,0,0)

Let's call AFC_vector as the output of the AFC_gene_vector function then the result could be obtained by calling the function as follows:

result<-gene_expression_estimation(genotype_h1,genotype_h2,AFC_vector)

