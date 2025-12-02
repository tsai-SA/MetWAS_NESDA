# MetWAS-NESDA
External validation of metabolome-wide association study of antidepressant exposure

Shared files required
UKB_AD_met_weights.rds: Metabolites and their LASSO weights trained using UKB cohort
Met_list.rds: List of metabolite with non-zero coefficient (LASSO) used in calculating the metabolic scores

## Metabolite preprocessing: Z-score standardisation
The metabolic scores were trained using standardised metabolic levels using z-score standardisation. Therefore we would like the metabolic scores to be calculated also using standardised metabolic levels.

The R script preprocessing_metabolites.R will read the dataframe (as .rds format) containing the metabolite levels, and will filter it to the metabolites that with non-zero coefficients from our LASSO training model. The metabolite dataframe should have rows as participant ID and columns as metabolite names.

Arguments:
--cohort : Cohort name, e.g 'UKB' or 'NESDA' \
--metabolites : The file path to the metabolite file (rds format) \
--list : The file path for the list of metabolites from LASSO(non-zero coefficients) provided by us \
--id_column : The column name of the identifier column (default == ID) \
--analysis: Either 'sig' (plot unstandardised vs standardised distribution of metabolites) or 'mrs' (just standardisation) \
--outdir : The directory where the outputs will be saved

Example:
Rscript preprocessing_metabolites.R \
  --cohort "UKB" \
  --metabolites "/Users/Desktop/test_met_data.rds" \
  --list "/Users/Desktop/probes_list.rds" \
  --id_column "ID" \
  --analysis "mrs" \
  --outdir "/Users/Desktop/"

## Calculate metabolic scores
The R script MetS_calc.R will read the rds file from preprocessing_metabolites.R, the LASSO weights provided by us, and the antidepressant exposure phenotype in your cohort (0 = no exposure, 1 = exposure). Then, it will calculate metabolic score for each participant.

Arguments:
--cohort : Cohort name, e.g 'UKB' or 'NESDA' \
--std_met : The file path for the standardised metabolite file from preprocessing_metabolites.R \
--id_column : The column name of the identifier column (default == ID) \
--ukb_weights : The file path for the metabolite weights file provided by us \
--pheno : The file path to the antidepressant exposure phenotype file for your cohort (rds format). The file should have two columns: ID and antidep_expo with antidep_expo coded as 0 (no exposure) and 1 (exposure) \
--outdir : The directory where the results and graphs will be saved 

Example:
Rscript MetS_calc.R \
--cohort "UKB" \
--std_met "/Users/Desktop/UKB_mrs_preproc_metabolites.rds" \
--id_column "ID" \
--ukb_weights "/Users/Desktop/lasso_coef.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--outdir "/Users/Desktop/"

## Predictive models

Key elements of the model:

Phenotype : Antidepressant exposure phenotype (0/1)
Predictor : Metabolic scores output from MetS_calc.R
Covariates: Should have rows as participant ID and columns as covariates.
- Technical covariates = assessment centre (factor variable) and spectrometer (factor variable).
- Demographic variables = age (numeric), sex (factor), smoking status (factor), socioeconomic status (factor), education , BMI, ethnicity, diagnosis of major depressive diisorder, alcohol drinking status)




General Output

We would like the coefficients of the model, alongside the standard errors, t values, P values. We will assess the model's performance using AUC, including ROC and PR curves.


