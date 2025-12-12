###############################################################################

# Demographic information about sample from the cohort 
# age, sex, smoking_status, education, bmi, ethnicity, mdd, alcohol_drinking

###############################################################################

# Set up libraries and options/files

###############################################################################

library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tidyr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort, ideally no spaces (for graphs and documentation)", action='store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column names of identifier column in phenotype and covariate files", action = 'store'),
  make_option('--ms', type = 'character', help = 'File path to metabolic score file (made using MetS_calc.R)'),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file, column names (ID and antidep_expo),
  make_option('--demo', type = 'character', help = 'File path to demographics file, column names,
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
id_col <- opt$id_column # Vector of identifier columns 
pheno_fp=opt$pheno # AD exposure (phenotype of cohort)
MetS_fp=opt$ms # AD MetS (predictor) 
demo_fp=opt$demo # Demographic variables
outdir <- opt$outdir # File path of output directory

# sinking all output to a log file 

sink(paste0(outdir, cohort, "_cohort_demographics.log"))


###############################################################################

# Read in files

###############################################################################

# read in file with additional demographic information for manuscript 

demographics <- readRDS(demo_fp)
AD_MetS <- readRDS(MetS_fp)

if (endsWith(pheno_fp, '.rds')){
  ad_pheno <- readRDS(pheno_fp)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .rds')
}


###############################################################################

# File checks 

###############################################################################

if('MetS' %in% colnames(AD_MetS) == FALSE){
  stop('No MetS column in the AD_MetS file')
} else {
  print('MetS column in the AD_MetS file')
}


# check that there is an antidep column in the file 

if('antidep_expo' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep_expo column in the phenotype file, please change name')
} else {
  print('antidep_expo column in the phenotype file')
}

# filter out any missing phenotype values (if any?)
print('Missing values in the antidep_expo phenotype? ')
table(is.na(ad_pheno$antidep_expo))

ad_pheno <- ad_pheno %>% filter(!is.na(antidep_expo))

# Demographics file 
# If applicable: age, sex, smoking_status, educatio, bmi, ethnicity, mdd, alcohol_drinking
# load in vector of required demographic variables 

req_demo_vars <- c('age', 'sex', 'smoking_status', 'education', 'bmi', 'ethnicity', 'mdd', 'alcohol_drinking')
if(all(req_demo_vars %in% colnames(demographics))){
  print('All demographic variables loaded in and named correctly')
} else {
    print('Not all demographic variables present / not loaded correctly')
    cols_missing = req_demo_vars[!(req_demo_vars %in% colnames(demographics))]
    print(paste0('Incorrect/missing columns: ', cols_missing))
}

# merge with the AD_MetS file 

demographics <- merge(demographics, AD_MetS, by= id_col, all = TRUE)

###############################################################################

# Summarising demographics function 

###############################################################################

# create a summary of this phenotype 
# summarise function with the continuous covs 

# summary function for cases and controls 
# phenotype is the phenotype column name in the phenotype file 
# file is the dataframe with the phenotype and covariate information 


phenotype_summary <- function(phenotype, file) {
  file <- file[!is.na(file[,phenotype]),]
  print(colnames(file))
  summary <- file %>% 
    group_by(across(all_of(phenotype)))  %>% 
    summarise(age = paste0(signif(mean(age, na.rm = T),3), 
                           ' (', signif(sd(age, na.rm = T),3), ')'), 
              bmi = paste0(signif(mean(bmi, na.rm = T),3), ' (', 
                           signif(sd(bmi, na.rm = T),3),')'), 
              MetS = paste0(signif(mean(MetS, na.rm = T),3),' (', 
                                signif(sd(MetS, na.rm = T),3), ')'), 
              
              n = n()) %>%
           
           Female = c(file %>% 
                        group_by(across(all_of(c(phenotype, 'sex')))) %>% 
                        summarise(count = n()) %>% 
                        mutate(percentage = (count/sum(count))*100) %>%
                        mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                        filter(!!sym(phenotype)==0 & sex == 'Female') %>% 
                        pull(val), 
                      file %>% group_by(across(all_of(c(phenotype, 'sex')))) %>% 
                        summarise(count = n()) %>% 
                        mutate(percentage = (count/sum(count))*100) %>%
                        mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                        filter(!!sym(phenotype)==1 & sex == 'Female') %>% 
                        pull(val)), 
           Male = c(file %>% 
                      group_by(across(all_of(c(phenotype, 'sex')))) %>% 
                      summarise(count = n()) %>% 
                      mutate(percentage = (count/sum(count))*100) %>%
                      mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                      filter(!!sym(phenotype)==0 & sex == 'Male') %>% 
                      pull(val), 
                    file %>% group_by(across(all_of(c(phenotype, 'sex')))) %>% 
                      summarise(count = n()) %>% 
                      mutate(percentage = (count/sum(count))*100) %>%
                      mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                      filter(!!sym(phenotype)==1 & sex == 'Male') %>% 
                      pull(val))) %>% 
    as.data.frame()
  
  summary <- summary %>% select(phenotype, n, age ,bmi,  Female, Male, MetS)
  colnames(summary) <- c(phenotype, 'N', 'Age (SD)','BMI (SD)', 'Sex Female (%)', 'Sex Male (%)', 'AD MetS (%)')
  return(summary)
}

###############################################################################

# Merging with phenotype and summarising

###############################################################################

demographics_pheno <- merge(ad_pheno, demographics, by = id_col)
cleaned_demographics_pheno <- demographics_pheno[complete.cases(demographics_pheno$antidep_expo), ]
demo_summary <- phenotype_summary(phenotype = 'antidep_expo', file = cleaned_demographics_pheno)

# add MDD summary if that info is available 

if('mdd' %in% colnames(demographics_pheno)){
  print('MDD data available for the phenotype')
  # summarise the MDD data 
  # print the table of mdd variable to log (to sanity check the conversion)
  print('The mdd variable:')
  table(demographics_pheno$mdd)
  demographics_pheno <- demographics_pheno %>%
    mutate(mdd_pheno = case_when(mdd == 0 ~ 'Control', mdd==1 ~ 'Case', TRUE ~ NA))
  
  print('The mdd named variable:')
  table(demographics_pheno$mdd_pheno)
  MDD_summary <- data.table(antidep_expo = c(0, 1),
    Cases = c(demographics_pheno %>%
                             group_by(across(all_of(c('antidep_expo', 'mdd_pheno')))) %>% 
                             summarise(count = n()) %>% 
                             mutate(percentage = (count/sum(count))*100) %>%
                             mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                             filter(antidep_expo==0 & mdd_pheno == 'Case') %>% 
                             pull(val), 
              demographics_pheno  %>%
                             group_by(across(all_of(c('antidep_expo', 'mdd_pheno')))) %>% 
                             summarise(count = n()) %>% 
                             mutate(percentage = (count/sum(count))*100) %>%
                             mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                             filter(antidep_expo==1 & mdd_pheno == 'Case') %>% 
                             pull(val)),
  Controls = c(demographics_pheno  %>%
                 group_by(across(all_of(c('antidep_expo', 'mdd_pheno')))) %>% 
                 summarise(count = n()) %>% 
                 mutate(percentage = (count/sum(count))*100) %>%
                 mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                 filter(antidep_expo==0 & mdd_pheno == 'Control') %>% 
                 pull(val), 
               demographics_pheno  %>%
                 group_by(across(all_of(c('antidep_expo', 'mdd_pheno')))) %>% 
                 summarise(count = n()) %>% 
                 mutate(percentage = (count/sum(count))*100) %>%
                 mutate(val = paste0(count, ' (', signif(percentage,2), '%)')) %>%
                 filter(antidep_expo==1 & mdd_pheno == 'Control') %>% 
                 pull(val)
  )) %>% as.data.frame()
  colnames(MDD_summary) <- c('antidep_expo','MDD cases', 'MDD controls')

  demo_summary <- merge(demo_summary, MDD_summary, by = 'antidep_expo')

} else{
  print('No MDD phenotype data available')
}

###############################################################################

# Save the demographic summary 

###############################################################################

write.table(demo_summary, paste0(outdir, cohort, '_demo_summary.tsv'), sep = '\t',quote = F, row.names = F)
print(paste0('Saved the demographic summary to ', cohort, '_demo_summary.txt'))

sink()

