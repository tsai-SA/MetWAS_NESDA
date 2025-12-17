
###### Set up libraries and options/files

library(data.table)
library(dplyr)
library(optparse)
library(readr)
library(tools)
library(ggplot2)
library(tidyr)

parse <- OptionParser()

# setting up options for the filepaths to the correct files
option_list <- list(
  make_option('--cohort', type='character', help="Cohort", action='store'),
  make_option('--metabolites', type='character', help="The filepath for metabolites file", action='store'), 
  make_option('--probe', type = 'character', help= "The filepath for the list of probe metabolites", action = 'store'),
  make_option('--id_column', type = 'character', default="ID", help = "Column name for identifier column", action = 'store'),
  make_option('--outdir', type = 'character', help = 'The filepath for output directory', action = 'store')
)

args = commandArgs(trailingOnly=TRUE)
opt <- parse_args(OptionParser(option_list=option_list), args=args)


# setting up arguments from the options 
print('Setting up the options')
cohort <- opt$cohort
met_filepath=opt$metabolites # metabolite file
probe_filepath=opt$probe # List of metabolite names
id_col <- opt$id_column # Vector of identifier column
out_dir <- opt$outdir

sink(paste0(out_dir, cohort, "_", "_std_metabolites.log"))  #log file

print(paste0('Metabolite file from : ', met_filepath))
print(paste0('List of probe metabolites from : ', probe_filepath))
print(paste0('ID column : ', id_col))
print(paste0('Output to be saved in : ', out_dir))

###############################################################################

# Read the files

###############################################################################

if (endsWith(met_filepath, ".rds")){
  metabolites <- readRDS(met_filepath)
} else {
  stop("Unsupported file format. Please provide a file with .rds format")
}

met_list <- readRDS(probe_filepath)
print('Read in files')

###############################################################################

# Filter the metabolites to just the metabolites used in the LASSO training 

###############################################################################
print('Filter to just the LASSO metabolites')
ms_met <- metabolites %>% select(c(all_of(id_col), intersect(names(metabolites), met_list)))

print(paste0('Filtered to ', ms_met %>% select(-id_col) %>% ncol())) 
rm(metabolites) # remove large metabolite object 

###############################################################################

# Scale the metabolite levels 

###############################################################################
print('Scaling the metabolite columns')

met_std <- ms_met %>% mutate(across(-c(all_of(id_col)), scale))

# Plotting the distribution of unstandardised and standardised values
# for 3 randomly selected metabolites

met_std_compare <- sample(setdiff(names(met_std), id_col), 3)
met_both <- rbind(
  ms_met %>% 
  select(c(all_of(id_col), all_of(met_std_compare))) %>%
    pivot_longer(cols = -c(all_of(id_col)), 
    names_to = "Metabolites", 
  values_to = "Mval") %>% 
  as.data.frame() %>% 
  mutate(Values = 'Unstandardised'),
met_std %>% 
  select(c(all_of(id_col), all_of(met_std_compare))) %>%
  pivot_longer(cols = -c(all_of(id_col)), 
               names_to = "Metabolites", 
               values_to = "Mval") %>%
  mutate(Values = 'Standardised')
)

met_dists <- ggplot(met_both, aes(x = Mval, fill = Metabolites)) +
  geom_histogram() + 
  facet_grid(Values~Metabolites) +
  ggtitle(paste0(cohort, ': Random sample of metabolites - standardisation'))

ggsave(filename=paste0(out_dir, cohort, "_", "_std_unstd_metabolites.png"),met_dists, 
       width = 8, height = 6, device='png', dpi=300)

###############################################################################

# Write out the filtered and standardised metabolite data 

###############################################################################

outfile <- paste0(out_dir, cohort, "_", "_std_metabolites.rds")
saveRDS(met_std, outfile)
sink()


