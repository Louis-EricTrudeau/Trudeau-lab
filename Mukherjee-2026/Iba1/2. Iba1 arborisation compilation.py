import os, glob
import pandas as pd
import numpy as np

masterDir = os.path.dirname(__file__)

all_df = []

for root, dirsList, filesList in os.walk(masterDir):
    for dir in dirsList:
        if dir == 'Iba1 arborisation':
            dirPath = os.path.join(root, dir) 
            for fileNamePath in glob.glob(dirPath + '/*.csv'): 
                with open(fileNamePath, 'r') as csvfile:
                    df = pd.read_csv(fileNamePath, sep=',', index_col=False) 
                    df.convert_dtypes() 
                    
                    df['Brain #'] = fileNamePath.split('\\')[-4].split('.')[1].replace(' ', '') 
                    df["Brain region"] = df.apply(lambda row: row.Image.split('-')[-1][0], axis = 1)
                    
                    Iba1AreaDorsal = df.loc[df['Brain region'] == 'D', 'Area'].mean() 
                    Iba1AreaVentral = df.loc[df['Brain region'] == 'V', 'Area'].mean() 

                    Iba1MeanDorsal = df.loc[df['Brain region'] == 'D', 'Mean'].mean()
                    Iba1MeanVentral = df.loc[df['Brain region'] == 'V', 'Mean'].mean()

                    df['Area (dorsal striatum)'] = Iba1AreaDorsal
                    df['Area (ventral striatum)'] = Iba1AreaVentral
                    df['Intensity (dorsal striatum)'] = Iba1MeanDorsal
                    df['Intensity (ventral striatum)'] = Iba1MeanVentral

                    df = df.drop_duplicates(subset=['Brain #', 'Area (dorsal striatum)', 'Area (ventral striatum)', 'Intensity (dorsal striatum)', 'Intensity (ventral striatum)'], keep='first')

                    all_df.append(df)

df_merged = pd.concat(all_df, ignore_index=False, sort=True)
columns = ['Brain #', 'Brain region', 'Area (dorsal striatum)', 'Area (ventral striatum)', 'Intensity (dorsal striatum)', 'Intensity (ventral striatum)']
df_merged = df_merged.reindex(columns=columns) 

df_merged = df_merged.replace(0, np.nan).dropna(axis=0, how='all')

df_merged = df_merged.groupby(['Brain #', 'Brain region']).mean().reset_index()

df_merged = df_merged.groupby(['Brain #']).mean().reset_index()

df_merged.to_excel(masterDir + "/compilation-Iba1-arborisation.xlsx", index=False) 