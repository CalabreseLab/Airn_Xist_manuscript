# parse MEME restults .txt file to extract all motifs and their stats

import pandas as pd
import re
import os

# get all foldernames in this directory D:\Cygwin\home\Think\airnupdatemotif
foldernames = os.listdir('D:\\Cygwin\\home\\Think\\XAK_classic_allw')
#foldernames = os.listdir('/Users/sl/XAKmotif')

# get the foldernames starts with XAK
x_allw = [x for x in foldernames if x.startswith('mXist')]
a_allw = [x for x in foldernames if x.startswith('airn')]
k_allw = [x for x in foldernames if x.startswith('kcnq1ot1')]

# define a function to parse the meme.txt file in each folder
def parse_meme_txt(foldername):
    print(foldername)
    # get the path of the meme.txt file
    meme_txt_path = os.path.join('D:\\Cygwin\\home\\Think\\XAK_classic_allw', foldername, 'meme.txt')
    #meme_txt_path = os.path.join('/Users/sl/XAKmotif', foldername, 'meme.txt')

    # check if the file exists
    if os.path.exists(meme_txt_path):
        # Open the file and read the lines
        with open(meme_txt_path, 'r') as f:
            lines = f.readlines()

        # Regular expression pattern to match lines starting with "MOTIF"
        pattern_motif = r'^MOTIF .+ MEME-\d+\twidth'
        pattern_following = r'^Motif .+ MEME-\d+ regular expression'

        # List to store matched lines
        data = []
        last_motif_index = -1

        for i in range(len(lines)-1):
            line = lines[i].strip()
            # If the line matches the motif pattern
            if re.match(pattern_motif, line):
                parts = re.split(r'[ \t]+', line)
                item = {
                    'MOTIF': parts[1],
                    'MEME': parts[2],
                    'Width': parts[5],
                    'Sites': parts[8],
                    'LLR': parts[11],
                    'E-value': parts[14],
                    'PatternFull': None,  # Set 'PatternFull' to None by default
                    'block': foldername[:-18]  # Set 'block' to the foldername but leaves out the -meme-classic-allw part
                }
                data.append(item)
                last_motif_index = len(data) - 1
            # If the line matches the following line pattern
            elif last_motif_index != -1 and re.match(pattern_following, line):
                # Get the second line following the matched line
                second_line = lines[i+2].strip() if i+2 < len(lines) else None
                # Append the second line to the last motif item as 'PatternFull'
                if second_line is not None:
                    data[last_motif_index]['PatternFull'] = second_line

        # Convert the data to a DataFrame
        df = pd.DataFrame(data)

        # convert df['P-value'] to numeric
        df['E-value'] = pd.to_numeric(df['E-value'])

        # keep only the motifs with E-value < 0.05 and save as a new df
        df = df[df['E-value'] < 0.05]

        # sort the df_fl by E-value, ascending
        df_sort = df.sort_values(by='E-value', ascending=True)

        return df_sort

# loop through all foldernames_de and parse the meme.txt file in each folder
# concatename all dfs into one df

x_motif = pd.concat([parse_meme_txt(x) for x in x_allw])
a_motif = pd.concat([parse_meme_txt(x) for x in a_allw])
k_motif = pd.concat([parse_meme_txt(x) for x in k_allw])

# dave df_de as a csv file
x_motif.to_csv('xist_classic_allw_combined.csv', index=False)
a_motif.to_csv('airn_classic_allw_combined.csv', index=False)
k_motif.to_csv('kot1_classic_allw_combined.csv', index=False)

