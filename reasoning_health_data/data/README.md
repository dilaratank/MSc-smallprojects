# Data dictionary 

## BMI dataset (ProcessedGlobalMeanBMI.csv):
Field name: | Data type: | Description:
--- | --- | --- |
region | Character | Full country name, like “Afghanistan”
ISO | Character | 3166-1 alpha-3, a three-letter alphabetic code used to uniquely name and identify countries
mean_BMI | Numeric | Mean BMI for a country in the year 2016, averaged over men and women 


## AFib dataset (ProcessedGlobalAFibPrevalence.csv) :
Field name: | Data type: | Description:
--- | --- | --- |
Country | Character | Full country name, like “Afghanistan”
Standardized Prevalence | Character | 2017 Age-standardized prevalence rate of atrial fibrillation per 100.000 people
Prevalence | Character | Prevalence group. “Low” if the prevalence rate was 400 or less. “Medium” if the prevalence was between 400 and 600. “High” if the prevalence was 600 or above.  
iso_code | Character | 3166-1 alpha-3, a three-letter alphabetic code used to uniquely name and identify countries

## Combined dataset (Combined_ds.csv):
Field name: | Data type: | Description:
--- | --- | --- |
49436004 | Numeric | SNOMED-CT code for atrial fibrillation which covers the 2017 Age-standardized prevalence rate of atrial fibrillation per 100.000 people
prevalence_group | Character | Prevalence group. “Low” if the prevalence rate was 400 or less. “Medium” if the prevalence was between 400 and 600. “High” if the prevalence was 600 or above. 
ISO | Character | ISO 3166-1 alpha-3, a three-letter alphabetic code used to uniquely name and identify countries
region | Character | Full country name, like “Afghanistan”
60621009 | Numeric | SNOMED-CT code for BMI which covers the mean BMI for a country in the year 2016, averaged over men and women 

