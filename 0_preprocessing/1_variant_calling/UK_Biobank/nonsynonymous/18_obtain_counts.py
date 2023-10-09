import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

# Get variant burden counts for each sample
df=pd.read_csv('tables/17_gene_ids.csv')

out=df.Sample.value_counts()

# Save to file
out.to_csv('tables/18_burden_table.csv', header=False)

# Also make a histogram of burden
plot_df=pd.DataFrame(out)
sns.histplot(plot_df, x='Sample')
plt.xlabel('SNV burden')
plt.savefig('Figures/18_burden_histogram.pdf')
plt.close()

print(plot_df.Sample.value_counts())
