import sys
import pandas as pd

fname = sys.argv[1]

df=pd.read_table(fname)

df['mismatch_prop'] = df.eval('mismatches/reads_all')

df['percent_A'] = df.eval('(A/reads_all) * 100')

df['percent_C'] = df.eval('(C/reads_all) * 100')

df['percent_G'] = df.eval('(G/reads_all) * 100')

df['percent_T'] = df.eval('(T/reads_all) * 100')

df2 = df[df.mismatch_prop >= 0.3]

df3 = df[df.reads_all == 0]

df2.to_csv(fname + '.degenerates.csv')

df3.to_csv(fname + '.gaps.csv')

print('data filtered')
