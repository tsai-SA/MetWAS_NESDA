# MetWAS-NESDA
External validation of metabolome-wide association study of antidepressant exposure

Shared files required
UKB_AD_met_weights.txt: Text file of the metabolites and their LASSO weights
Met_list.txt: Text file of the list of metabolite names and groups used in calculating the metabolic scores

## Metabolite preprocessing
The metabolic scores were trained using standardised metabolic levels using z-score standardisation. Therefore we would like the metabolic scores to be calculated also using standardised metabolic levels.

R:

preprocessing_metabolites.R: which will read the dataframe (as .rds format) containing the metabolite levels, and will filter it to the metabolites that with non-zero coefficients from our LASSO training model. The metabolite dataframe should have rows as participant ID and columns as metabolite names.
