# acc_regions_mammals_aves

<hr/>

## Roadmap

Analysis was performed following vignettes from https://github.com/CshlSiepelLab/RPHAST

<hr/>

### aves - sarcopterygii

result file from this analysis: `/aves/data/output/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv`

generated with `/aves/r_source/functional_regions/phastCons.analyze.R`

neutral model for conservation and acceleration `galGal6.phastCons77way.mod` downloaded from 
[https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/phastCons77way/galGal6.phastCons77way.mod](https://hgdownload.soe.ucsc.edu/goldenPath/galGal6/phastCons77way/galGal6.phastCons77way.mod)

conservation was computed with the code in `/aves/r_source/conservation`

conservation partial results are in `/aves/data/phastCons.zip`

acceleration was computed with the code in `/aves/r_source/acceleration`

acceleration results can be found in `/aves/data/acceleration/non_param_sim_phyloP/acc_elements_bed`

<hr/>

### mammals - sarcopterygii

result file from this analysis: `/mammals/data/output/join_filtered_elements_norm_cod_nonCod_oneacczerogap_functionalregions.csv`

generated with `/mammals/functional_regions/phastCons.analyze.R`

neutral model for conservation and acceleration in mammals - sarcopterygii

`hg38.phastCons100way.mod` downloaded from [https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.mod](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons100way/hg38.phastCons100way.mod)

conservation was computed with the code in `/mammals/conservation/scripts/conservation`

conservation partial results can be found in `/mammals/conservation/data/output`

acceleration was computed with the code in `/mammals/acceleration`

acceleration results can be found in `/mammals/acceleration/data/non_param_sim_phyloP/sarcopterygii_mammals/acc_elements_bed`

<hr/>

Code dependencies

* R depedencies are listed in `setup_2024.R`

* Python dependencies are:
    * boto3
    * botocore
    * shutil
    * biopython
    * pandas
    
* Binaries:
  
    * `mafSpeciesSubset` from [https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/](https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/)
    

<hr/>


If you require further information or have any questions about reusing this code for your analysis, feel free to contact me at paulati.ingebi (at) gmail.com

