
# ncontRanalyser
<!-- badges: start --> 
<!-- badges: end -->

The goal of ncontRanalyser is to take CCP4 NCONT output .txt files and generate dataframes which can be used for downstream analysis. This package contains functions to both carry out analysis of NCONT antibody and antigen data, generating AlphaFold-Predicted Epitopes (AFPE) and AlphaFold-Predicted Paratopes (AFPP, more below). Results in .txt file formats can be generated using command line CCP4 NCONT. For sample code on how to run this, scroll to the very bottom.

Note the functions in this package work for the output from CCP4 NCONT v8.0.019, and does not support other versions of CCP4 NCONT.

## Installation

You can install ```ncontRanalyser``` using the following code:

``` r
# easiest to do using the remotes package. Run the first line if you don't already have the remots package.
install.packages(remotes)
library(remotes)
# then use
remotes::install_github("sproc2104/ncontRanalyser")

```

## Reference-based analysis

When a reference structure is provided and analysed using the ```reference_analysis()``` command, it is used to calculate antibody and antigen n/d (ratio of number of contacts to avergae distance) quartiles for the NCONT contacts. These quartiles can then be used to detect changes in the predicted epitope and paratope residues by using the ```mutant_analysis()``` command.

Note that the ```reference_analysis()``` command **must be assigned to a variable named ```ref```**. If it is not, then the ```mutant_analysis()``` command will not recognise the quartiles generated from the reference structure, and will fail.

The ```reference_analysis()``` and ```mutant_analysis()``` commands take as input the .txt files containing NCONT ouptputs for the reference or mutant structures, respectively.
``` r
library('ncontRanalyser')

ref<- reference_analysis("ref_heavychain.txt", "ref_lightchain.txt")
ref ##to see output

mut<- mutant_analysis("mutant_heavychain.txt", "mutant_lightchain.txt")
mut
```

The significance of the quartiles is arbitrary, but in previous analysis, I have selected those in the Q1 and above as composing the AFPE or AFPP.


## Reference-free analysis (own quartile input)

In the case where you want to choose the quartiles or not run a reference structure first, you can provide a set of values to use as quartiles.

Quartiles must be provided in the vector form, where:
``` r
ab_quartiles<- c(max_val, Q1, Q2, Q3)
ag_quartiles<- c(max_val, Q1, Q2, Q3)
```
These can then be passed to the ```mutant_analysis()``` function, which will assign antigen and antibody residue quartiles based on these values.
``` r
mut<- mutant_analysis("mutant_heavychain.txt", "mutant_lightchain.txt", ab_quartiles, ag_quartiles)
```

## Saving files

The outputs from above commands can be saved using the ```save_txt()``` function, which generate two separate .csv files for antibody and antigen residues.

```r
save_csv(ref, "abX_agY")

## The command will generate the files:
## "abX_agY_afpp_antibody.csv"
## "abX_agY_afpe_antigen.csv"
```


## Example code for running command line CCP4 NCONT

If you do not already have CCP4 NCONT, follow the installation instructions [here] (https://www.ccp4.ac.uk/download/index.php#os=macos). 
The two commands below will generate two separate .txt results files, as is needed for the ncontRanalyser package.

```unix
## where A is the heavy chain, B is the light chain, C is the antigen, and 4 is the distance threshold for contacts in A.
ncont XYZIN structure_1.pdb <<EOF >heavy_output.txt
source A
target C 
maxdist 4
EOF

ncont XYZIN structure_1.pdb <<EOF >light_output.txt
source B
target C 
maxdist 4
EOF
```
The ```heavy_output.txt``` and ```light_output.txt``` are now ready to be used as input for the functions above.
