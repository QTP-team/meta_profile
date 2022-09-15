import argparse
import pandas as pd

parser = argparse.ArgumentParser(
        usage = 'python rules/profile_merge.py -i */*_profile.tsv -o merge.tsv -s .tsv')
parser.add_argument("-i", "--infiles", nargs = "+", help = "One or more profile files")
parser.add_argument("-o", "--output")
parser.add_argument("-s", "--suffix", default = ".tsv", help = "suffix of profile files")

args = parser.parse_args()

df = pd.DataFrame(columns=["Genome_ID","rel_ab"])
for file in args.infiles:
    data = pd.read_csv(file, sep = "\t", header = 0)
    data = data[['Genome_ID', 'rel_ab']]
    data.columns = ['Genome_ID',file.split('/')[1].split('_profile')[0]]
    df = pd.merge(data,df,on='Genome_ID',how='outer')
del df['rel_ab']
df = df.fillna(0)
df.rename(columns={'Genome_ID':'Sample'}, inplace = True)
df = df.T
df.to_csv(args.output,header = 0)
