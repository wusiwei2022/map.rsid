# Tutorial on how to convert genetic coordinate in GWAS summary statistics to RSIDs

## Background
Genetic variants can be presented in different formats in different GWAS summary statistics. The most common identifier for genetic variants is RSID. In some other GWAS, genetic variants may be expressed as chr:pos or other identifiable marker.However, some GWAS downstream analysis need to work on multiple GWAS summary statistics, where compatible identifiers are essential across all the GWAS summary statistics. RSID is the most common and stable identifier for genetic variants. Therefore, in this tutorial, we presented different approaches to label RSID to genetic variants using the coordinate and allele information of the variants.

## GWAS summary data
The example data that we will use in this tutorial is the publicly available GWAS summmary data by Genes & Health study (GRCh37) and the GWAS summary by UK Biobank Plasma Protein Project, UKBPPP (GRCh38). The GWAS summary data (GRCh37) of Genes & Health Study can be downloaded at https://www.genesandhealth.org/researchers/data/. The GWAS summary data (GRCh38) of UKBPPP can be downloaded at https://metabolomips.org/ukbbpgwas/. In this tutorial, we will use a subset of the GWAS summary data, which is available at /map.rsid/data/

## Approach 1: map with ensembl-vep
Ensembl-vep can be used to annotate genetic variants with RSIDs. You can download and install the ensembl-vep tool by following the official tutorial at https://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html. We recommended to install the VEP with Singularity if you are working on HPC.
The GWAS summary statistics need to be converted to ensembl-vep accepted format first. Please refer to https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html#input to understand the accepted format by ensembl-vep. We tend to format the GWAS summary statistics into vcf format because the required information is available in most GWAS summary and it is quite easy to handle. The necessary information are genetic coordinate (chromsome and position), reference allele, and alternative allele. It's also good to provide a unique identifier for each variant.<br>
Please note ensembl-vep can only navigate forwardly. Therefore, please do run ensembl-vep at the parent directory that accommodates the vep tool, the input, and ouput directory.

### Format the GWAS summary into vcf format
You may refer to /results/gwas_bmi_vep_subset.input and /results/gwas_ukbppp_A1BG_vep.input as the template for the input of ensembl-vep. There are 5 columns in the input file including chromsome, genetic position, identifier, reference allele, and alternative allele. The identifier is not mandatory. You may fill chromsome:position:reference allele:alternative allele into this column or simply fill in with ".".

### Map to RSID with ensembl-vep on Genes & Health BMI GWAS subset data (GRCh37)
```
singularity exec softwares/vep.sif\
 vep\
 --dir softwares/vep_data\
 --cache\
 --assembly GRCh37\
 -i tutorials/map.rsid/data/gwas_bmi_vep_subset.input\
 --check_existing\
 --af_1kg\
 --pick\
 --force_overwrite\
 --tab\
 --fields Uploaded_variation,Location,Allele,SAS_AF,Existing_variation\
 -o tutorials/map.rsid/results/gwas_bmi_vep/gwas_bmi_vep_subset.annotation
```

### Map to RSID with ensembl-vep on UKBPPP A1BG GWAS subset data (GRCh38)
```
singularity exec softwares/vep.sif\
 vep\
 --dir softwares/vep_data\
 --cache\
 --assembly GRCh38\
 -i tutorials/map.rsid/data/gwas_ukbppp_A1BG_vep_subset.input\
 --check_existing\
 --af_1kg\
 --pick\
 --force_overwrite\
 --tab\
 --fields Uploaded_variation,Location,Allele,SAS_AF,Existing_variation\
 -o tutorials/map.rsid/results/gwas_ukbppp_vep/gwas_ukbppp_A1BG_vep_subset.annotation
```

## Approach 2: map with bridge file (snp151common.txt)
We can also map RSID to genetic coordinate of the variants using a bridge file. The bridge file can be downloaded at https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ for reference genome version GRch37 (hg19). The full variants are available in snp151.txt.gz while the common variants are available in snp151common.txt.gz. The schema for the brdige file can be found at https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=snp151Common&hgta_table=snp151Common. The bridge file for reference genome version GRch38 can be downloaded at https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/. You can download the bridge file using wget.
```
# wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp151.txt.gz # complete variants (GRCh37)
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp151Common.txt.gz # common variants (GRCh37)
# wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/snp151.txt.gz # complete variants (GRCh38)
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/snp151Common.txt.gz # common variants (GRCh38)
```

### Process the snp151common bridging genetic coordinate and RSID
The bridge file snp151common includes 26 columns, the field name and descriptions are as follows. In this tutorial, we removed cDNA and unknown variants. We will only map genetic variants on chromosome 1 - 22, and  therefore removed genetic variants X chromosome, Y chromosome, and random sequence scaffold. We only map SNPs in this tutorial, and you may try to map all the genetic variants by not filtering by class == "single". Only colnames "chrom", "chromStart", "chromEnd", "name", and "alleles" will be used for mapping.
```{r}
snp151common = fread("snp151Common.txt.gz")
```

```{r}
colnames(snp151common) = c("bin", #  Indexing field to speed chromosome range queries.
                           "chrom", # Reference sequence chromosome or scaffold
                           "chromStart", # Start position in chrom
                           "chromEnd", # End position in chrom
                           "name", # dbSNP Reference SNP (rs) identifier
                           "score", # Not used
                           "strand", # Which DNA strand contains the observed alleles
                           "refNCBI", # Reference genomic sequence from dbSNP
                           "refUCSC", # Reference genomic sequence from UCSC lookup of chrom,chromStart,chromEnd
                           "observed", # The sequences of the observed alleles from rs-fasta files
                           "molType", # Sample type from exemplar submitted SNPs (ss)
                           "class", # Class of variant (single, in-del, named, mixed, etc.)
                           "valid", # Validation status of the SNP
                           "avHet", # Average heterozygosity from all observations. Note: may be computed on small number of samples.
                           "avHetSE", # Standard Error for the average heterozygosity
                           "func", # Functional category of the SNP (coding-synon, coding-nonsynon, intron, etc.)
                           "locType", # Type of mapping inferred from size on reference; may not agree with class
                           "weight", # The quality of the alignment: 1 = unique mapping, 2 = non-unique, 3 = many matches
                           "exceptions", # Unusual conditions noted by UCSC that may indicate a problem with the data
                           "submitterCount", # Number of distinct submitter handles for submitted SNPs for this ref SNP
                           "submitters", # List of submitter handles
                           "alleleFreqCount", # Number of observed alleles with frequency data
                           "alleles", # Observed alleles for which frequency data are available
                           "alleleNs", # Count of chromosomes (2N) on which each allele was observed. Note: this is extrapolated by dbSNP from submitted frequencies and total sample 2N, and is not always an integer.
                           "alleleFreqs", # Allele frequencies
                           "bitfields" # SNP attributes extracted from dbSNP's SNP_bitfield table
                           )
```

```{r}
table(snp151common$molType)
snp151common = snp151common %>% filter(molType == "genomic") # remove cDNA & unknown genetic variants (V11)
```

```{r}
# table(snp151common$chrom)
snp151common = snp151common %>% filter(chrom != "chrX" & chrom != "chrY") # remove genetic variants on chromosome X and Y
snp151common = snp151common[!grepl("_", snp151common$chrom), ] # remove the random sequence scaffold
```

```{r}
snp151common = snp151common %>% filter(class == "single") # Keep only SNPs
```

```{r}
snp151common = snp151common[, c("chrom", "chromStart", "chromEnd", "name", "alleles")]
snp151common = snp151common %>% mutate(chrom= gsub("chr", "", chrom))
```

```{r}
sum(duplicated(snp151common$rsid))  # check the uniqueness of rsids: no duplicates
```

### Format and map for Genes & Health BMI GWAS subset data
We only did the mapping for the genetic variants on autosome. Since column "rsid" have the same information as SNPID, we removed this column. The genetic position in GWAS summary data is usually the coordinate end of the genetic variant. Therefore, we use "chrom" and "chromEnd" to map rsid. There are 9,527,863 chromosomal variants carried for RSID mapping. 7,105,989 variants are successfully mapped while 6,509,534 SNPs have alleles matched with alleles in reference data. You can use the same pipeline if you need to map the fullset GWAS summary data. You can also do the conversion based on GRCh38 genetic coordinate by using the GRCh38 bridge file.
```{r}
bmi_gwas = fread("../data/gwas_bmi_snp151common.input")
```

```{r}
# table(bmi_gwas$CHR)
bmi_gwas = bmi_gwas %>% filter(!(CHR == "X"))
bmi_gwas = bmi_gwas %>% select(-c("rsid"))
dim(bmi_gwas)
```

```{r}
bmi_gwas = right_join(snp151common %>% select(-chromStart), bmi_gwas, by = c("chrom" = "CHR", "chromEnd" = "POS"))
# paste0(dim(bmi_gwas %>% filter(is.na(name)))[1], " genetic variants are unmappped")
bmi_gwas = bmi_gwas %>% filter(!is.na(name))
dim(bmi_gwas)[1]
```

```{r}
bmi_gwas = bmi_gwas[str_detect(bmi_gwas$alleles, bmi_gwas$Allele1) & str_detect(bmi_gwas$alleles, bmi_gwas$Allele2), ]
dim(bmi_gwas)[1]
```

```{r}
bmi_gwas = bmi_gwas %>% select(-alleles)
colnames(bmi_gwas) = c("CHR", "POS", "RSID", "SNPID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "imputationInfo", "N", "BETA", "SE", "Tstat", "p.value", "varT", "varTstar")
fwrite(bmi_gwas, "../results/gwas_bmi_snp151common/bmi_gwas.rsid")
```
