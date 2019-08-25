# Prediction of gene expression
This script is designed to predict allele-specific expressian and total gene expression using allelic fold change (AFC). See the [manuscript](https://genome.cshlp.org/content/27/11/1872.short) for method description.

## Inputs

### AFC file
In the data folder there is a file named aFC_Whole_Blood.txt. This is a tab delimited file providing allelic fold change for each eQTL. In this file for each variant in columnn "variant_id" we have the following information:

- gene_id : gene associated to that variant

- rank : rank of that variant in eQTL calling process

- CHROM : the chromosome of the variant

- POS : position of that variant

- REF : The Reference allele

- ALT : The Alternative allele

- log2_aFC : The effect size of the variant

### VCF file

For using this script a phased vcf file is required to extract the genotypes of individuals. In the data folder there is a vcf_smaple.txt file. 

genotype, encoded as allele values separated by either of / or |. " /" means genotype unphased and "|" means genotype phased. The allele values are 0 for the reference allele (what is in the REF  field), 1 for the  
allele listed in ALT. For diploid calls examples could be 0/1, 1|0. 

**The REF and ALT information should match the REF and ALT information in aFC_Whole_Blood.txt**

In order to extract the subset of the vcf file you can use the following command which IDs.txt repesents the required variants a sample of this file is in data folder.

```Shell
for line in $(cat ~/data/IDs.txt)
    do
    id=`echo ${line} | sed 's/\"/\t/g'`
    chr=`echo ${id} | sed 's/_/\t/g' | cut -f 1`
    pos=`echo ${id} | sed 's/_/\t/g' | cut -f 2`
    echo ${id}
    echo ${chr}
    echo ${pos}
    tabix vcf_file ${chr}:${pos}-${pos} | grep -w $id >> vcf_file_subset.txt
    done


```

Here is an example to predict the expression of a specific gene for an individual using these inputs

```R
# AFC file
AFC_file = "data/aFC_Whole_Blood.txt"

#VCF file
vcf_file = "data/sample_vcf.txt"


# Name of gene that we want to predict its expression
gene_id = "gene_id"

# Name of individual that we want to predict its expression
individual_id = "individual_id"

#########################

expression_prediction_gene_individual<-function(AFC_file,vcf_file,gene_id,individual_id){
    # read AFC file
    AFC_dt=read.table(AFC_file, header=TRUE, sep="\t")

    #read vcf file 
    genotype_info= read.table(vcf_file, header=TRUE, sep = " ")
   

    # get the afc vector for a specific gene 
    # the function definition is available in R folder
    AFC_vector<-AFC_gene_vector(gene_id,AFC_dt)
   

    # get the afc vector for a specific gene 
    # the function definition is available in R folder
    variant_vector<-variant_gene_vector(gene_id,AFC_dt)


    #null genotypes
    genotype_h1 = c()
    genotype_h2 = c()

# assuming variant and individual_id exist in vcf file
    for (variant in variant_vector){
        
        genotype<-genotype_info[genotype_info$ID == variant,individual_id]      
        genotype_split_full <-unlist(strsplit(as.character(genotype),":"))[1]
        genotype_split<-unlist(strsplit(as.character(genotype),"|"))
        genotype_h1<-c(genotype_h1,genotype_split[1])
        genotype_h2<-c(genotype_h2,genotype_split[3])
    }

    
    #get the expression output vector for two haplotypes
    # the function definition is available in R folder
    #output[1] : log expression for haplotype 1 in log2 scale
    #output[2] : gene expression for haplotype 2 in log2 scale
    #output[3] : total expression in log2 scale
    result<-gene_expression_estimation(as.numeric(genotype_h1),as.numeric(genotype_h2),AFC_vector)
    
    return(result)
    }
    ```

# R resources

## gene_expression_estimation_functions.ipynb

This notebook provides functions in R used in the above example that will read the AFC file and predict the gene expression.


## gene_expression_lookupTable.R

This R script gets a sorted AFC file (**sorted based on gene_id**), counts the number of variants for each gene and produces expression values for all genotypes. To run the script use the following command:

```Shell
    Rscript gene_expression_lookupTable.R data\aFC_Whole_Blood.txt
```    



