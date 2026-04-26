import os, glob
import pandas as pd
import numpy as np

masterDir = os.path.dirname(__file__)

all_df = []

for root, dirsList, filesList in os.walk(masterDir):
    for dir in dirsList:
        if dir == 'Iba1 soma':
            dirPath = os.path.join(root, dir) 
            for fileNamePath in glob.glob(dirPath + '/*.csv'): 
                with open(fileNamePath, 'r') as csvfile:
                    df = pd.read_csv(fileNamePath, sep=',', index_col=False) 
                    df.convert_dtypes() 

                    df['Brain #'] = fileNamePath.split('\\')[-4].split('.')[1]
                    df['Brain #'] = df['Brain #'].str.replace(' ', '') 
                    df['Brain region'] = fileNamePath.split('\\')[-1].split('.oif.csv')[0].split('-')[-1][0] 
                    
                    Iba1CountDorsal = df.loc[df['Brain region'] == 'D', 'Area'].count()/(0.21197*0.21197)
                    Iba1CountVentral = df.loc[df['Brain region'] == 'V', 'Area'].count() /(0.21197*0.21197)

                    Iba1areaDorsal = df.loc[df['Brain region'] == 'D', 'Area'].mean() 
                    Iba1areaVentral = df.loc[df['Brain region'] == 'V', 'Area'].mean() 

                    Iba1MeanDorsal = df.loc[df['Brain region'] == 'D', 'Mean'].mean()
                    Iba1MeanVentral = df.loc[df['Brain region'] == 'V', 'Mean'].mean()

                    Iba1PerimDorsal = df.loc[df['Brain region'] == 'D', 'Perim.'].mean()
                    Iba1PerimVentral = df.loc[df['Brain region'] == 'V', 'Perim.'].mean()

                    Iba1CircDorsal = df.loc[df['Brain region'] == 'D', 'Circ.'].mean()
                    Iba1CircVentral = df.loc[df['Brain region'] == 'V', 'Circ.'].mean()

                    Iba1RoundDorsal = df.loc[df['Brain region'] == 'D', 'Round'].mean()
                    Iba1RoundVentral = df.loc[df['Brain region'] == 'V', 'Round'].mean()

                    Iba1SolidityDorsal = df.loc[df['Brain region'] == 'D', 'Solidity'].mean()
                    Iba1SolidityVentral = df.loc[df['Brain region'] == 'V', 'Solidity'].mean()

                    df['Count (dorsal striatum)'] = Iba1CountDorsal 
                    df['Count (ventral striatum)'] = Iba1CountVentral
                    df['Area (dorsal striatum)'] = Iba1areaDorsal 
                    df['Area (ventral striatum)'] = Iba1areaVentral
                    df['Intensity (dorsal striatum)'] = Iba1MeanDorsal 
                    df['Intensity (ventral striatum)'] = Iba1MeanVentral
                    df['Perimeter (dorsal striatum)'] = Iba1PerimDorsal 
                    df['Perimeter (ventral striatum)'] = Iba1PerimVentral
                    df['Circularity (dorsal striatum)'] = Iba1CircDorsal 
                    df['Circularity (ventral striatum)'] = Iba1CircVentral
                    df['Roundness (dorsal striatum)'] = Iba1RoundDorsal 
                    df['Roundness (ventral striatum)'] = Iba1RoundVentral
                    df['Solidity (dorsal striatum)'] = Iba1SolidityDorsal 
                    df['Solidity (ventral striatum)'] = Iba1SolidityVentral

                    df = df.drop_duplicates(subset=['Brain #', 'Area (dorsal striatum)', 'Area (ventral striatum)', 'Intensity (dorsal striatum)', 'Intensity (ventral striatum)'], keep='first')

                    all_df.append(df)

df_merged = pd.concat(all_df, ignore_index=False, sort=True) 
columns = ['Brain #', 'Brain region', 'Count (dorsal striatum)', 'Count (ventral striatum)', 'Area (dorsal striatum)', 'Area (ventral striatum)', 'Intensity (dorsal striatum)', 'Intensity (ventral striatum)', 'Perimeter (dorsal striatum)', 'Perimeter (ventral striatum)', 'Circularity (dorsal striatum)', 'Circularity (ventral striatum)', 'Roundness (dorsal striatum)', 'Roundness (ventral striatum)', 'Solidity (dorsal striatum)', 'Solidity (ventral striatum)']
df_merged = df_merged.reindex(columns=columns)

df_merged = df_merged.replace(0, np.nan).dropna(axis=0, how='all')

df_merged = df_merged.groupby(['Brain #', 'Brain region']).mean().reset_index()

df_merged = df_merged.groupby(['Brain #']).mean().reset_index()

df_merged.to_excel(masterDir + "/compilation-Iba1-soma.xlsx", index=False) 