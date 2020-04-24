import pandas as pd
import numpy as np
import sys
df = pd.read_csv(sys.argv[1], error_bad_lines=False, header=None, skiprows=2)
print(len(np.unique(df.iloc[:, 1].astype(str))))