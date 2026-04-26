import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

from scipy.stats import ttest_ind

csv_path = '/data/rudko/vgrouza/invivomouse/pddata/stats'
df_all_data = pd.read_csv(os.path.join(csv_path, 'pd_mouse_data_emd_only.csv'))

# Convert categorical variables to categories
df_all_data['Genotype'] = df_all_data['Genotype'].astype('category')
df_all_data['Infected'] = df_all_data['Infected'].astype('category')
df_all_data['Timepoint'] = df_all_data['Timepoint'].astype('category')

# Separate data by timepoints
df_timepoint1 = df_all_data[df_all_data['Timepoint'] == 1]
df_timepoint2 = df_all_data[df_all_data['Timepoint'] == 2]

# Define ROIs
rois = ['EMDstriatum', 'EMDthalamus', 'EMDpssc', 'EMDdg']

results = {}
posthoc_results = {}

# Perform two-way ANOVA for each ROI at each timepoint
for roi in rois:
    # Timepoint 1
    model1 = ols(f'{roi} ~ C(Genotype) * C(Infected)', data=df_timepoint1).fit()
    anova_table1 = sm.stats.anova_lm(model1, typ=2)
    results[f'Timepoint 1 - {roi}'] = anova_table1

    # Timepoint 2
    model2 = ols(f'{roi} ~ C(Genotype) * C(Infected)', data=df_timepoint2).fit()
    anova_table2 = sm.stats.anova_lm(model2, typ=2)
    results[f'Timepoint 2 - {roi}'] = anova_table2

    # Perform post-hoc tests for significant effects at Timepoint 1
    if anova_table1['PR(>F)']['C(Genotype)'] < 0.05:
        print(f"\nPost-hoc analysis for {roi} at Timepoint 1 - Genotype")
        tukey_genotype1 = pairwise_tukeyhsd(df_timepoint1[roi], df_timepoint1['Genotype'])
        print(tukey_genotype1)
        posthoc_results[f'Timepoint 1 - {roi} - Genotype'] = tukey_genotype1

    if anova_table1['PR(>F)']['C(Infected)'] < 0.05:
        print(f"\nPost-hoc analysis for {roi} at Timepoint 1 - Infected")
        tukey_infected1 = pairwise_tukeyhsd(df_timepoint1[roi], df_timepoint1['Infected'])
        print(tukey_infected1)
        posthoc_results[f'Timepoint 1 - {roi} - Infected'] = tukey_infected1

    if anova_table1['PR(>F)']['C(Genotype):C(Infected)'] < 0.05:
        print(f"\nPost-hoc analysis for {roi} at Timepoint 1 - Interaction")
        tukey_interaction1 = pairwise_tukeyhsd(df_timepoint1[roi],
                                               df_timepoint1['Genotype'].astype(str) + ':' + df_timepoint1[
                                                   'Infected'].astype(str))
        print(tukey_interaction1)
        posthoc_results[f'Timepoint 1 - {roi} - Interaction'] = tukey_interaction1

    # Perform post-hoc tests for significant effects at Timepoint 2
    if anova_table2['PR(>F)']['C(Genotype)'] < 0.05:
        print(f"\nPost-hoc analysis for {roi} at Timepoint 2 - Genotype")
        tukey_genotype2 = pairwise_tukeyhsd(df_timepoint2[roi], df_timepoint2['Genotype'])
        print(tukey_genotype2)
        posthoc_results[f'Timepoint 2 - {roi} - Genotype'] = tukey_genotype2

    if anova_table2['PR(>F)']['C(Infected)'] < 0.05:
        print(f"\nPost-hoc analysis for {roi} at Timepoint 2 - Infected")
        tukey_infected2 = pairwise_tukeyhsd(df_timepoint2[roi], df_timepoint2['Infected'])
        print(tukey_infected2)
        posthoc_results[f'Timepoint 2 - {roi} - Infected'] = tukey_infected2

    if anova_table2['PR(>F)']['C(Genotype):C(Infected)'] < 0.05:
        print(f"\nPost-hoc analysis for {roi} at Timepoint 2 - Interaction")
        tukey_interaction2 = pairwise_tukeyhsd(df_timepoint2[roi],
                                               df_timepoint2['Genotype'].astype(str) + ':' + df_timepoint2[
                                                   'Infected'].astype(str))
        print(tukey_interaction2)
        posthoc_results[f'Timepoint 2 - {roi} - Interaction'] = tukey_interaction2

# Print results for each timepoint and ROI
for key, value in results.items():
    print(f"\nANOVA results for {key}:")
    print(value)

# Applying Bonferroni correction to p-values from post-hoc tests
all_pvals = []
for key, tukey_result in posthoc_results.items():
    all_pvals.extend(tukey_result.pvalues)

_, corrected_pvals, _, _ = multipletests(all_pvals, alpha=0.05, method='bonferroni')

# Printing the Bonferroni-corrected p-values
print("\nBonferroni-corrected p-values for post-hoc tests:")
for i, key in enumerate(posthoc_results.keys()):
    print(
        f"{key}: Corrected p-values - {corrected_pvals[i * len(posthoc_results[key].pvalues):(i + 1) * len(posthoc_results[key].pvalues)]}")


## Tests for Timepoint 1 vs Timepoint 2

# Store results
comparison_results = {}

# Perform t-test for each ROI, comparing Timepoint 1 and Timepoint 2
for roi in rois:
    # Perform independent t-test
    t_stat, p_val = ttest_ind(df_timepoint1[roi].dropna(), df_timepoint2[roi].dropna(), equal_var=False)

    # Store results
    comparison_results[roi] = {'t-statistic': t_stat, 'p-value': p_val}

# Print the comparison results
print("Comparison of Timepoint 1 vs Timepoint 2 for each ROI (pooled across Genotype and Infected):")
for roi, result in comparison_results.items():
    print(f"{roi}: t-statistic = {result['t-statistic']:.4f}, p-value = {result['p-value']:.4f}")


# Filter data by genotype
df_KO = df_all_data[df_all_data['Genotype'] == 'KO']
df_WT = df_all_data[df_all_data['Genotype'] == 'WT']

# Store results
interaction_results = {}

# Perform two-way ANOVA with interaction for each ROI and genotype
for roi in rois:
    # KO group
    model_KO = ols(f'{roi} ~ C(Infected) * C(Timepoint)', data=df_KO).fit()
    anova_table_KO = sm.stats.anova_lm(model_KO, typ=2)
    interaction_results[f'KO - {roi}'] = anova_table_KO

    # WT group
    model_WT = ols(f'{roi} ~ C(Infected) * C(Timepoint)', data=df_WT).fit()
    anova_table_WT = sm.stats.anova_lm(model_WT, typ=2)
    interaction_results[f'WT - {roi}'] = anova_table_WT

# Print results for each genotype and ROI
print("Testing for difference in effect of infection between Timepoint 1 and Timepoint 2 for KO and WT groups separately:")
for key, value in interaction_results.items():
    print(f"\nANOVA results for {key}:")
    print(value)
