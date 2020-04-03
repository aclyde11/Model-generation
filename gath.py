import glob
import pandas as pd
import sys

files = glob.glob(sys.argv[1] + "/*/")

dfs = []
for file in files:
    try:
        df = pd.read_csv(file + "metrics.csv")
        df['file'] = file
        dfs.append(df)
        print(df)
    except:
        pass

pd.concat(dfs).to_csv("combined_adpr_2.csv")

