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
  make_option('--covs', type = 'character', help = 'File path to covariate file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno # AD exposure (phenotype of cohort)
MetS_fp=opt$ms # AD MetS (predictor) 
covs_fp=opt$covs 
outdir <- opt$outdir # File path of output directory

# sinking all output to a log file 

sink(paste0(outdir, cohort, "_MetS_AD_assoc.log"))

###############################################################################

# Read in covariates, phenotype and MetS files 

###############################################################################

MetS <- read.table(MetS_fp, header = T) # File with ID and MetS
print(colnames(MetS))

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

all_covs <- readRDS(covs_fp)

print(paste0("Covariates read in ", paste(colnames(all_covs %>% dplyr::select(-all_of(id_col))), collapse = ", ")))

#merge the phenotype and MetS file together 

MetS_pheno <- merge(MetS, ad_pheno, by = id_col)
MetS_pheno_covs <- merge(MetS_pheno, all_covs, by = id_col)


# logging phenotype characteristics after merging 

print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ' after merging with MetS and pheno: Number of cases: ',
             nrow(MetS_pheno_covs %>% 
                    filter(antidep==1)), 
             ' \n Number of controls: ',
             nrow(MetS_pheno_covs %>% 
                    filter(antidep==0))))

###############################################################################
  
  # Generalised linear model (GLM)
  
###############################################################################

# Fit a logistic model 
# logit link function and binomial family 
# Outcome - Antidepressant exposure phenotype 
# Predictor - Antidepressant MetS (from MetS_calc.R)

assoc_mod <- glm(as.factor(antidep) ~ scale(AD_MetS) + as.integer(age) + as.factor(assessment_centre) +
                 as.factor(spectrometer) + as.factor(sex) + as.factor(smoking_status) + as.factor(education) + 
                 as.factor(ethnicity) + as.numeric(bmi) + as.factor(mdd) + as.factor(alcohol_drinking), 
                 family=binomial (link=logit), data = MetS_pheno_covs)

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

###############################################################################

# Calculating the AUC and ROC graph

###############################################################################

predicted_probs <- predict(assoc_mod, type = 'response')

# Take the true outcomes (MetS_pheno_covs$antidep)
# and the predicted probabilities for the 1 ('case') class
# returns false positive and true positive rates for different classification thresholds

roc_curve <- roc(MetS_pheno_covs$antidep, predicted_probs)
auc_value <- auc(roc_curve)

# save ROC curve object for plotting all cohorts together
print(paste0('Saving the ROC curve object for plotting all cohorts together to rds object: ', outdir, cohort, '_roc_curve.rds'))
saveRDS(roc_curve, paste0(outdir, cohort, '_roc_curve.rds'))

# ROC Graph 
print(paste0('Saving the ROC curve for the cohort alone to ', outdir, cohort, '_assoc_ROC_curve.pdf'))
cairo_pdf(file = paste0(outdir, cohort, '_assoc_ROC_curve.pdf'), width = 8, height = 6)
plot.roc(roc_curve, col = "blue", lwd =2, main = paste0('ROC Curve: ', cohort))
dev.off()

###############################################################################

# Precision Recall Curve 

###############################################################################

pr_curve <- pr.curve(MetS_pheno_covs$antidep, predicted_probs, curve = T)

cairo_pdf(file = paste0(outdir, cohort, '_assoc_precision_recall.pdf'), width = 8, height = 6)
plot(pr_curve, main= paste0(cohort, ' : Precision Recall Curve'), col = 'red')
dev.off()

sink()
