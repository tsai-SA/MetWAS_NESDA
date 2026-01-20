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
library(ggpubr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
# ukb_weights contain two columns: Metabolites and weights
# pheno should have two columns: ID and antidep_expo

option_list <- list(
  make_option('--cohort', type='character', help="Cohort", action='store'),
  make_option('--std_met', type='character', help="The filepath for standardised metabolite file", action='store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column names of participant identifier column", action = 'store'),
  make_option('--ukb_weights', type = 'character', default = 'UKB_AD_met_weights.rds', help = "Weights file provided for calculating MS"),
  make_option('--pheno', type = 'character', help = 'File path to antidepressant exposure phenotype file'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)
cohort <- opt$cohort
std_met_filepath=opt$std_met
ukb_weights_filepath=opt$ukb_weights # Weights file from GS
id_col <- opt$id_column # Vector of identifier columns 
pheno_filepath=opt$pheno
outdir <- opt$outdir

# Create a log file
sink(paste0(outdir, cohort, "_AD_MS.log"))
print(paste0('Calculating metabolic scores for ', cohort))
print(paste0('Read in the processed metabolite file from: ', std_met_filepath))
print(paste0('Read in the weights file from: ', ukb_weights_filepath))
print(paste0('Read in the pheno file from: ', pheno_filepath))
print(paste0('Output to be saved in: ', outdir))

###############################################################################

# Read the files

###############################################################################

std_met <- readRDS(std_met_filepath) 
print(paste0('The metabolite file has data for ', nrow(std_met), ' participants and ', std_met %>% select(-ID) %>% ncol(), " Metabolites"))  # Remove ID col
ukb_weights <- readRDS(ukb_weights_filepath) 
print(paste0('The weights file has weights for ', nrow(weights), ' Metabolites'))

###############################################################################

# Assess missingness of metabolites

###############################################################################

if(std_met %>% select(-ID) %>% ncol() != nrow(ukb_weights)){
  print(paste0('Number of non-zero coef metabolites read in for metabolite file (', std_met %>% select(-ID) %>% ncol(),
               ') does not match the number of weights provided (', nrow(ukb_weights), ')'))
  # If not all weights included, save the metabolites included in the MS for the cohort 
  readRDS(std_met %>% select(-ID) %>% colnames(), paste0(outdir,cohort, '_metinMS.txt'))
  print(paste0('Metabolites used in the MS are saved in: ', outdir, cohort, '_metinMS.txt'))
} else {
  print(paste0('Number of probes read in for std_met matches the number of weights provided: n = ', nrow(ukb_weights)))
}

missing_percentage <- std_met %>% select(-ID) %>%
  summarise_all(~ mean(is.na(.)) * 100) %>%
  gather(Metabolites, MissingPercentage) %>% arrange(desc(MissingPercentage))

print(paste0('The metabolite with the highest level of missingness is: ', missing_percentage$Metabolites[1]))
print(paste0('There are ', nrow(missing_percentage %>% filter(MissingPercentage > 5)), ' Metabolites with more than 5% missingness and ',
             nrow(missing_percentage %>% filter(MissingPercentage > 50)), ' with more than 50% missingness'))

print(missing_percentage)

missing_hist <- ggplot(missing_percentage, aes(x = MissingPercentage, fill = MissingPercentage < 50)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = '% of Missingness', y = 'Number of metabolite')+
  geom_vline(xintercept = 50, color = "red", linetype = 'dashed')+
  ggtitle(cohort)

missing_plot <- ggplot(missing_percentage %>% filter(MissingPercentage > 50),
                           aes(x = reorder(Metabolites, -MissingPercentage), y = MissingPercentage)) +
    geom_bar(stat='identity', fill = 'skyblue', color = 'black') + 
    labs(title = paste0(cohort, ': Metabolite missingness in MS'), x = "Metabolite", y = "% Missing") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# Plot the % of missingness against the absolute value of the MS coefficient 

missingness_weights <- merge(missing_percentage, ukb_weights, by = 'Metabolites')

missing_weights_plt <- ggplot(missingness_weights, aes(x = MissingPercentage, y = abs(weights))) + 
  geom_point() + 
  theme_minimal() + 
  labs(x = '% of Missingness', y ='MS weight', 
       title = paste0(cohort, ': % missing vs MS weights for metabolites > 50% missingness'))

all_missing_plots <- ggarrange(missing_hist, missing_weights_plt, missing_plot, nrow = 3, ncol = 1)

ggsave(paste0(outdir, cohort, "_MS_met_missingness.png"), missing_plot, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_MS_met_missingness_hist.png"), missing_hist, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_MS_met_missingness_weights.png"), missing_weights_plt, width = 8, height = 6, device='png', dpi=300)
ggsave(paste0(outdir, cohort, "_MS_met_missingness_allplots.png"), all_missing_plots, width = 8, height = 6, device='png', dpi=300)

write.table(missing_percentage, paste0(outdir, cohort, '_met_missingness.txt'), row.names = F, quote = F)

###############################################################################

# Calculate the metabolic scores 

###############################################################################
print('Change std_met colnames into UKB ID format)
rename_vec <- setNames(ukb_weights$ukb_id, ukb_weights$nesda_abbre)
std_met <- std_met %>%
  rename(any_of(rename_vec))  

print('Converting std_met to long format')
long_std_met <- std_met %>% 
  pivot_longer(-c(all_of(id_col)), 
               names_to = "Metabolites", 
               values_to= "Mval")
print('Merging with the weights')
long_std_met <- merge(long_std_met, ukb_weights, by = 'Metabolites')

print('Calculating the metabolic scores')

# Group by ID
# and then calculate a weighted sum of all the metabolites per ID

MetS <- long_std_met %>% group_by(!!sym(id_col)) %>%
  summarise(weighted_sum = sum(weights*Mval, na.rm = T)) %>% as.data.frame()

print(paste0('Metabolic scores calculated for ', nrow(MetS), ' participants in the ', cohort, ' cohort'))

###############################################################################

# Distribution of the number of metabolites within the MetS (non-missing)

###############################################################################

# Showing as similar metric to the above missing plots, but another way of looking at it 

num_met <- long_std_met %>%
  group_by(!!sym(id_col)) %>%
  summarise(met_ms = sum(!is.na(Mval)))

met_ms_hist <- ggplot(num_met, aes(x = met_ms)) + 
  geom_histogram() +
  theme_minimal() + 
  labs(x = 'Number of metabolites included in the MetS', y = 'Frequency') + 
  ggtitle(paste0(cohort, ': Histogram of number of metabolites included in the MetS'))

ggsave(paste0(outdir, cohort, "_indiv_nummet_MS.png"), met_ms_hist, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Metabolic Score Distributions

###############################################################################

print(paste0('Plotting distributions across the whole ', cohort, ' sample'))

# Distribution of the metabolic score across the cohort 

MS_dist <- ggplot(MetS, aes(x = weighted_sum)) + 
  geom_histogram() + 
  theme_minimal() + 
  labs(x = 'Metabolic Profile Score', y = 'Count')+
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_AD_MS_overalldist.png"), MS_dist, width = 8, height = 6, device='png', dpi=300)

# look at Distribution in AD exposed and AD not exposed (violin plots)

if (endsWith(pheno_filepath, '.rds')){
ad_pheno <- readRDS(pheno_filepath)
} else {
  stop('Unsupported phenotype file, please provide the phenotype as a .rds file')
}

if('antidep_expo' %in% colnames(ad_pheno) == FALSE){
  stop('No antidep column in the phenotype file')
} else {
  print('antidep column in the phenotype file')
}

ad_pheno <- ad_pheno %>% filter(!is.na(antidep_expo)) # remove missing values if there are any 

print(paste0('Read in the Antidepressant exposure phenotype for ', cohort, ': Number of cases: ',
             nrow(ad_pheno %>% 
                    filter(antidep_expo==1)), 
             ' Number of controls: ',
             nrow(ad_pheno%>% 
                    filter(antidep_expo==0))))

ad_pheno_MS <- merge(ad_pheno, MetS, by = id_col)

print('Plotting MetS distributions for AD exposure cases and controls ')
MS_pheno_dists <- ggplot(ad_pheno_MS, aes(x = weighted_sum, fill = as.factor(antidep_expo))) + 
  geom_histogram(alpha = 0.8) + 
  theme_minimal() + 
  labs(x = 'Metabolic Profile Score', y = 'Count', fill = 'AD use') +
  ggtitle(cohort)

ggsave(paste0(outdir, cohort, "_AD_MS_phenodist.png"), MS_pheno_dists, width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Saving the metabolic score 

###############################################################################
outfile <- file.path(outdir, paste0(cohort, "_AD_MetS.rds"))
print(paste0('Saving the metabolic score to ', outfile))

colnames(MetS)[2] <- 'AD_MetS'
saveRDS(MetS, outfile)
sink()

