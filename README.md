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
