import pandas as pd
a = pd.DataFrame([[1,2,3],
                  [3,1,3],
                  [1,2,1]])
result = a.apply(pd.value_counts)
print result



