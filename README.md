# MetWAS-NESDA
External validation of metabolome-wide association study of antidepressant exposure

Shared files required:
1) UKB_AD_met_weights.rds: Metabolites and their LASSO weights trained using UKB cohort
2) Met_list.rds: List of metabolite with non-zero coefficient (LASSO) used in calculating the metabolic scores

## Overall summary:
1) Run preprocessing_metabolites.R to standardise the metabolites
2) Run MetS_calc.R to apply the weights from our training model to your cohort and produce a metabolic score for each participant in your cohort
3) Run predict_model.R to use the metabolic scores created in the previous step to predict antidepressant exposure status in your cohort
Please refer to the following information for further details and please let us know if you have any questions, thank you so much for your help!

## Metabolite preprocessing: Z-score standardisation
The metabolic scores were trained using standardised metabolic levels using z-score standardisation. Therefore we would like the metabolic scores to be calculated also using standardised metabolic levels.

The R script preprocessing_metabolites.R will read :
1) The dataframe (as .rds format) containing your cohort's metabolite levels, it should have rows as participant ID and columns as metabolite names.
2) The Met_list.rds file which contains a list of probe metabolites provided by us
Then the R script will filter your cohort's metabolites to the probe metabolites from our training model and scale them

Arguments:
--cohort : Cohort name, e.g 'UKB' or 'NESDA' \
--metabolites : The file path to the metabolite file (rds format) \
--list : The file path for the list of metabolites from LASSO(non-zero coefficients) provided by us \
--id_column : The column name of the identifier column (default == ID) \
--analysis: Either 'sig' (plot unstandardised vs standardised distribution of metabolites) or 'mrs' (just standardisation) \
--outdir : The directory where the outputs will be saved

Example:
```bash
Rscript preprocessing_metabolites.R \
  --cohort "UKB" \
  --metabolites "/Users/Desktop/test_met_data.rds" \
  --list "/Users/Desktop/probes_list.rds" \
  --id_column "ID" \
  --analysis "mrs" \
  --outdir "/Users/Desktop/"
```

## Calculate metabolic scores
The R script MetS_calc.R will read: 
1) The _mrs_preproc_metabolites.rds file from preprocessing_metabolites.R,
2) The LASSO weights provided by us
3) The antidepressant exposure phenotype in your cohort (0 = no exposure, 1 = exposure), this file should be .rds format with two columns, ID and antidep_expo
Then, it will calculate metabolic score for each participant.

Arguments:
--cohort : Cohort name, e.g 'UKB' or 'NESDA' \
--std_met : The file path for the standardised metabolite file from preprocessing_metabolites.R \
--id_column : The column name of the identifier column (default == ID) \
--ukb_weights : The file path for the metabolite weights file provided by us \
--pheno : The file path to the antidepressant exposure phenotype file for your cohort (rds format). The file should have two columns: ID and antidep_expo with antidep_expo coded as 0 (no exposure) and 1 (exposure) \
--outdir : The directory where the results and graphs will be saved 

Example:
```bash
Rscript MetS_calc.R \
--cohort "UKB" \
--std_met "/Users/Desktop/UKB_mrs_preproc_metabolites.rds" \
--id_column "ID" \
--ukb_weights "/Users/Desktop/lasso_coef.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 1 : basic model in all participants
The R script basic_model_all.R will read:
1) The _AD_MetS.rds file from MetS_calc.R output
2) The antidepressant exposure phenotype in your cohort (0 = no exposure, 1 = exposure), this file should be .rds format with two columns, ID and antidep_expo
3) A .rds file with ID and the **basic covariates** including sex, age, techinical covariates such as assessment_centre and/or spectrometer (if applicable) as columns.

Then the metabolic scores created in the previous step will be used to predict antidepressant exposure status in your cohort. 

Key elements of the model: Phenotype (antidepressant exposure phenotype), predictor (metabolic scores output from MetS_calc.R), basic covariates (should have rows as participant ID and columns as covariates)

General Output: We would like the coefficients of the model, alongside the standard errors, t values, P values. We will assess the model's performance using AUC, including ROC and PR curves.

Arguments:
--cohort : Cohort name, e.g 'UKB' or 'NESDA' \
--id_column : The column name of the identifier column (default == ID) \
--ms : file path to the _AD_MetS.rds file from MetS_calc.R output
--pheno : The file path to the antidepressant exposure phenotype file for your cohort (rds format). The file should have two columns: ID and antidep_expo with antidep_expo coded as 0 (no exposure) and 1 (exposure) \
--basic_covs : The file path to the basic covariate file for your cohort (rds format). It should include ID and each covariate as column \
--outdir : The directory where the results and graphs will be saved 


```bash
Rscript basic_model_all.R \
--cohort "UKB" \
--id_column "ID" \
--ms "/Users/Desktop/UKB_AD_MetS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--basic_covs "/Users/Desktop/basic_covs.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 2 : basic model in MDD cases only
The basic_model_mdd.R script is similar to model 1, but with additional mdd column (MDD diagnosis) and it will filter out participants with no MDD diagnosis
1) The _AD_MetS.rds file from MetS_calc.R output
2) The antidepressant exposure phenotype in your cohort (0 = no exposure, 1 = exposure), this file should be .rds format with two columns, ID and antidep_expo
3) A .rds file with ID and the **basic covariates** including sex, age, mdd, techinical covariates such as assessment_centre and/or spectrometer (if applicable) as columns. The mdd column should be coded into 0/1 (1 = mdd, 0 = controls), the basic_model_mdd.R script will filter out the controls.
   
```bash
Rscript basic_model_mdd.R \
--cohort "UKB" \
--id_column "ID" \
--ms "/Users/Desktop/UKB_AD_MetS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--basic_covs "/Users/Desktop/basic_covs.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 3 : complex model in all participants
This is also similar to model 1 but with additional lifestyle-related covariates such as smoking status, socioeconomic status, educational level, MDD diagnosis, ethnicity, depressive symptom scores, alcohol drinking (only when applicable)

The R script complex_model_all.R will read:
1) The _AD_MetS.rds file from MetS_calc.R output
2) The antidepressant exposure phenotype in your cohort (0 = no exposure, 1 = exposure), this file should be .rds format with two columns, ID and antidep_expo
3) A .rds file with ID and the **complex covariates** including sex, age, techinical covariates such as assessment_centre and/or spectrometer, and lifestyle covariates (if applicable) as columns.

```bash
Rscript complex_model_all.R \
--cohort "UKB" \
--id_column "ID" \
--ms "/Users/Desktop/UKB_AD_MetS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--complex_covs "/Users/Desktop/complex_covs.rds" \
--outdir "/Users/Desktop/"
```

## Predictive model 4 : complex model in MDD cases only
This is also similar to model 3 (with lifestyle-related covariates) but in MDD cases only

The R script complex_model_mdd.R will read:
1) The _AD_MetS.rds file from MetS_calc.R output
2) The antidepressant exposure phenotype in your cohort (0 = no exposure, 1 = exposure), this file should be .rds format with two columns, ID and antidep_expo
3) A .rds file with ID and the **complex covariates** including sex, age, techinical covariates such as assessment_centre and/or spectrometer, and lifestyle covariates (if applicable) as columns. The mdd column should be coded into 0/1 (1 = mdd, 0 = controls), the complex_model_mdd.R script will filter out the controls.

```bash
Rscript complex_model_mdd.R \
--cohort "UKB" \
--id_column "ID" \
--ms "/Users/Desktop/UKB_AD_MetS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--complex_covs "/Users/Desktop/complex_covs.rds" \
--outdir "/Users/Desktop/"
```

## Demographic information about your cohort
The file cohort_demograph_info.R will generate a table of demographic information of your cohort (useful for the manuscript and interpretating our results). It formats information on age, sex, bmi, smoking, MDD diagnosis, and AD MetS (generated from MetS_calc.R).

--cohort: Cohort name, e.g 'GS' or 'NESDA'
--id_column: The column name of the identifier column (default == ID)
--ms: The filepath to the _AD_MetS.rds file from MetS_calc.R output (colnames: ID, AD_MetS)
--pheno: The filepath to the AD phenotype file (colnames: ID, antidep_expo(0 = no exposure/1 = antidepressant exposure))
--demo: Filepath to the file containing demographic variables: age(numeric), sex(factor: Male/Female), bmi(numeric), mdd(factor: 0 = controls/1 = cases), AD_MetS)

```bash
Rscript cohort_demograph_info.R \
--cohort "UKB" \
--id_column "ID" \
--ms "/Users/Desktop/UKB_AD_MetS.rds" \
--pheno "/Users/Desktop/pheno.rds" \
--demo "/Users/Desktop/covs.rds" \
--outdir "/Users/Desktop/"
```
