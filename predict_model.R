###############################################################################

# Set up libraries and options/files

###############################################################################
library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tidyr)
library(ggplot2)
library(tools)
library(lme4)
library(tibble)
library(pROC)
library(caret)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column names of identifier column", action = 'store'),
  make_option('--ms', type = 'character', help = 'File path to metabolic score file (made using MetS_calc.R)'),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno # AD exposure (phenotype of cohort)
MetS_fp=opt$ms # AD MetS (predictor) 
outdir <- opt$outdir # File path of output directory

# sinking all output to a log file 

sink(paste0(outdir, cohort, "_MetS_AD_assoc.log"))

###############################################################################

# Read in covariates, phenotype and MRS files 

###############################################################################

MetS <- read.table(MetS_fp, header = T) # File with ID and MetS
print(colnames(MRS))

# check that there is a MetS column in the file 

if('AD_MetS' %in% colnames(MetS) == FALSE){
  stop('No AD_MetS column in the MetS file')
} else {
  print('AD_MetS column in the MetS file')
}

# support phenotype .rds file

if (endsWith(pheno_fp, '.rds')){
  ad_pheno <- readRDS(pheno_fp)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .rds file')
}

# check that there is an antidep column in the file 

if('antidep' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep column in the phenotype file')
} else {
  print('antidep column in the phenotype file')
}

# remove missing values (if any?)
ad_pheno <- ad_pheno %>% filter(!is.na(antidep)) 

# logging phenotype characteristics 
print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' : Number of cases: ',
             nrow(ad_pheno %>% 
                    filter(antidep==1)), 
             'Number of controls: ',
             nrow(ad_pheno%>% 
                    filter(antidep==0))))

#merge the phenotype and MetS file together 

MetS_pheno <- merge(MetS, ad_pheno, by = id_col)

# logging phenotype characteristics after merging 

print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' after merging with MetS and pheno: Number of cases: ',
             nrow(MetS_pheno %>% 
                    filter(antidep==1)), 
             ' \n Number of controls: ',
             nrow(MetS_pheno %>% 
                    filter(antidep==0))))

###############################################################################
  
  # Generalised linear model (GLM)
  
###############################################################################

# Fit a logistic model 
# logit link function and binomial family 
# Outcome - Antidepressant exposure phenotype 
# Predictor - Antidepressant MetS (from MetS_calc.R)

MetS_pheno$antidep <- as.factor(MetS_pheno$antidep)

assoc_mod <- glm(antidep ~ scale(AD_MetS), family=binomial (link=logit), data = MetS_pheno)

# Extract the effect estimates, standard errors and p-value 
warnings()
print(assoc_mod)
print(summary(assoc_mod)$coefficients %>% as.data.frame())

assoc_coefs <- summary(assoc_mod)$coefficients %>% as.data.frame()
assoc_coefs <- rownames_to_column(assoc_coefs, var = "Coefficient")

# save the coefficients 
outfile <- file.path(outdir, paste0(cohort, "_MetS_AD_coefficients.rds"))
print(paste0('Saving the model coefficients to ', outfile))
saveRDS(assoc_coefs, outfile)


