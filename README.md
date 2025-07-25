# Tutorial on how to convert genetic coordinate in GWAS summary statistics to RSIDs.

## Background
Genetic variants can be presented in different formats in different GWAS summary statistics. The most common identifier for genetic variants is RSID. In some other GWAS genetic variants may be expressed as chr:pos or other identifiable marker.However, some GWAS downstream analysis need to work on multiple GWAS summary statistics, where compatible identifiers are essential across all the GWAS summary statistics. RSID is the most common and stable identifier for genetic variants. Therefore, in this tutorial, we presented different approaches to label RSID to genetic variants using the coordinate and allele information of the variants.

## GWAS summary data
The example data that we will use in this tutorial is the publicly available GWAS summmary data by Genes & Health study (GRCh37) and the GWAS summary by UK Biobank Plasma Protein Project, UKBPPP (GRCh38).  

## Approach 1: map with ensembl-vep
The GWAS summary statistics need to be converted to ensembl-vep accepted format first. Please refer to https://asia.ensembl.org/info/docs/tools/vep/vep_formats.html#input to understand the accepted format by ensembl-vep. We tend to format the GWAS summary statistics into vcf format because the required information is available in most GWAS summary and it is quite easy to handle. The necessary information are genetic coordinate (chromsome and position), reference allele, and alternative allele. It's also good to provide a unique identifier for each variant.<br>
Please note ensembl-vep can only navigate forwardly. Therefore, please do run ensembl-vep at the parent directory that accommodates the vep tool, the input, and ouput directory.

### Format the GWAS summary into vcf format
You may refer to ../results/gwas_bmi_vep/gwas_bmi_vep_subset.input as the template for the input of ensembl-vep. There are 5 columns in the input file including chromsome, genetic position, identifier, reference allele, and alternative allele. The identifier is not mandatory. You may fill chromsome:position:reference allele:alternative allele into this column or simply fill in with ".".
```
zcat /rds/general/project/medbio-epiukb-archive-2018/live/pickup/sw5122/resources/GenesHealth/results_GWAS_public/2025_05_13_BMI.max_singlevariant51kGSA-TOPMEDr3-GWAS_BMI.max.regenie.gz |\
awk 'NR!=1 && $1!=23 {print $1" "$2" "$3" "$4" "$5}' > ../results/gwas_bmi_vep/gwas_bmi_vep.input
head -100 ../results/gwas_bmi_vep/gwas_bmi_vep.input > ../results/gwas_bmi_vep/gwas_bmi_vep_subset.input
```

### Map to RSID with ensembl-vep on subset of data
```
singularity exec ~/softwares/vep.sif\
  vep\
  --dir ~/softwares/vep_data\
  --cache\
  --assembly GRCh37\
  -i ../results/gwas_bmi_vep/gwas_bmi_vep_subset.input\
  --check_existing\
  --af_1kg\
  --pick\
  --force_overwrite\
  --tab\
  --fields Uploaded_variation,Location,Allele,SAS_AF,Existing_variation\
  -o ../results/gwas_bmi_vep/gwas_bmi_vep_subset.annotation
```

## Approach 2: map with bridge file (snp151common.txt)
The bridge file can be downloaded at https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ for reference genome version GRch37 (hg19). The full variants are available in snp151.txt.gz while the common variants are available in snp151common.txt.gz. The schema for the brdige file can be found at https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=snp151Common&hgta_table=snp151Common. You can download the bridge file using wget:
```
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/snp151.txt.gz
```
