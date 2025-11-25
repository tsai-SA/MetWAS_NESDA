# MetWAS-NESDA
External validation of metabolome-wide association study of antidepressant exposure

Shared files required
UKB_AD_met_weights.rds: Metabolites and their LASSO weights trained using UKB cohort
Met_list.rds: List of metabolite with non-zero coefficient (LASSO) used in calculating the metabolic scores

## Metabolite preprocessing: Z-score standardisation
The metabolic scores were trained using standardised metabolic levels using z-score standardisation. Therefore we would like the metabolic scores to be calculated also using standardised metabolic levels.

The R script preprocessing_metabolites.R will read the dataframe (as .rds format) containing the metabolite levels, and will filter it to the metabolites that with non-zero coefficients from our LASSO training model. The metabolite dataframe should have rows as participant ID and columns as metabolite names.

Arguments:
--cohort : Cohort name, e.g 'UKB' or 'NESDA'
--metabolites : The file path to the metabolite file (rds format)
--list : The file path for the list of metabolites from LASSO(non-zero coefficients) provided by us 
--id_column : The name of the participant identifier column in the data 
--analysis: Either 'sig' (plot unstandardised vs standardised distribution of metabolites) or 'mrs' (just standardisation)
--outdir : The directory where the outputs will be saved

Example:
Rscript preprocessing_metabolites.R \
  --cohort "UKB" \
  --metabolites "/Users/Desktop/test_met_data.rds" \
  --list "/Users/Desktop/probes_list.rds" \
  --id_column "ID" \
  --analysis "mrs" \
  --outdir "/Users/Desktop/"

## Calculating metabolic scores
The R script MetS_calc.R: will calculate metabolic scores for participants.
Arguments:
--cohort : Cohort name
--std_met : The file path for the preprocessed DNAm file from process_DNAm_MRS.R, e.g /Users/data/DNAm/AD_MRS/GS_mrs_DNAm_preproc.txt
--id_column : The column name of the identifier column (default == IID)
--ukb_weights : The file path for the MRS weights files provided by GS e.g /Users/data/DNAm/GS_AD_MRS_weights.txt
--pheno : The file path to the antidepressant exposure phenotype file for your cohort. Script can read in either a .csv, or .txt file and should follow a format of FID(Optional), IID (Required), and antidep, with antidep being coded as 0 (no exposure) and 1 (exposure). e.g /Users/data/phenos/AD_pheno.csv

--outdir : The directory where the results and graphs will be saved e.g /Users/data/DNAm/AD_MRS/

