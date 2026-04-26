import os, glob
import pandas as pd
import numpy as np

masterDir = os.path.dirname(__file__)

all_df = []

for root, dirs_list, files_list in os.walk(masterDir):
    for file_name in files_list:
        if os.path.splitext(file_name)[-2] + os.path.splitext(file_name)[-1] == "intensity.csv":
            file_name_path = os.path.join(root, file_name)
            with open(file_name_path, 'r') as csvfile:
                df = pd.read_csv(file_name_path, sep=',', index_col=False) 
                df.convert_dtypes() 

                df['Brain #'] = file_name_path.split('\\')[-2].split('.')[1] 
                df["Brain region"] = df.apply(lambda row: row.Slice.split('-')[-1][0], axis = 1)
                        
                THDorsalStriatum = df.loc[df['Brain region'] == 'D', 'Mean'].mean() 
                THVentralStriatum = df.loc[df['Brain region'] == 'V', 'Mean'].mean() 

                tdTomatoDorsalStriatum = df.loc[df['Brain region'] == 'D', 'tdTomato'].mean() 
                tdTomatoVentralStriatum = df.loc[df['Brain region'] == 'V', 'tdTomato'].mean() 

                DATDorsalStriatum = df.loc[df['Brain region'] == 'D', 'DAT'].mean() 
                DATVentralStriatum = df.loc[df['Brain region'] == 'V', 'DAT'].mean() 

                THareaDorsalStriatum = df.loc[df['Brain region'] == 'D', '%Area'].mean()
                THareaVentralStriatum = df.loc[df['Brain region'] == 'V', '%Area'].mean()
                tdTomatoareaDorsalStriatum = df.loc[df['Brain region'] == 'D', 'tdTomato area'].mean()
                tdTomatoareaVentralStriatum = df.loc[df['Brain region'] == 'V', 'tdTomato area'].mean()
                DATareaDorsalStriatum = df.loc[df['Brain region'] == 'D', 'DAT area'].mean()
                DATareaVentralStriatum = df.loc[df['Brain region'] == 'V', 'DAT area'].mean()


                df['TH intensity (dorsal striatum)'] = THDorsalStriatum
                df['TH intensity (ventral striatum)'] = THVentralStriatum
                df['tdTomato intensity (dorsal striatum)'] = tdTomatoDorsalStriatum
                df['tdTomato intensity (ventral striatum)'] = tdTomatoVentralStriatum
                df['DAT intensity (dorsal striatum)'] = DATDorsalStriatum
                df['DAT intensity (ventral striatum)'] = DATVentralStriatum
                df['TH %Area (dorsal striatum)'] = THareaDorsalStriatum
                df['TH %Area (ventral striatum)'] = THareaVentralStriatum
                df['tdTomato %Area (dorsal striatum)'] = tdTomatoareaDorsalStriatum
                df['tdTomato %Area (ventral striatum)'] = tdTomatoareaVentralStriatum
                df['DAT %Area (dorsal striatum)'] = DATareaDorsalStriatum
                df['DAT %Area (ventral striatum)'] = DATareaVentralStriatum

                all_df.append(df)

df_merged = pd.concat(all_df, ignore_index=False, sort=True) 
df_merged = df_merged.reindex(columns=['Brain #', 'TH intensity (dorsal striatum)', 'TH intensity (ventral striatum)', 'tdTomato intensity (dorsal striatum)', 'tdTomato intensity (ventral striatum)', 'DAT intensity (dorsal striatum)', 'DAT intensity (ventral striatum)', 'TH %Area (dorsal striatum)', 'TH %Area (ventral striatum)', 'tdTomato %Area (dorsal striatum)', 'tdTomato %Area (ventral striatum)', 'DAT %Area (dorsal striatum)', 'DAT %Area (ventral striatum)']) # Reorder the columnns
df_merged = df_merged.drop_duplicates() 
df_merged.to_excel(masterDir + "/compilation-ec-gaussian-threshold.xlsx", index=False) 