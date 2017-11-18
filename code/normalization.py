import pandas as pd
import numpy as np

vec = pd.read_table('G:\project2\\NPM201507\\code\\1024term_vector.txt', index_col=0)
print vec.head(5)
df =(vec - vec.min()) / (vec.max() - vec.min())

df.to_csv('1024vector.txt', sep='\t', index=True)