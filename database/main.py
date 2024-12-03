import pandas as pd

class Database:
    def __init__(self, path, sep='\t'):
        self.df = pd.read_csv(path, sep=sep)


df = pd.read_csv("biogrid-sample.txt", sep='\t')
# df = pd.read_csv("BIOGRID-ALL-4.4.239.tab3.txt", sep='\t')
to_list = ["Synonyms Interactor A", "Synonyms Interactor B" ]
print(df[to_list])

for i in to_list:
    df[i] = df[i].apply(lambda x: x[0:-1].split('|'))


find = "YAL001C"

#print(df.columns)

#interest = ['Systematic Name Interactor A',
#       'Systematic Name Interactor B', 'Official Symbol Interactor A',
#       'Official Symbol Interactor B', 'Synonyms Interactor A',
#       'Synonyms Interactor B',]
#
#for i in interest:
#    print(df[i].head())

print(df[to_list])

#print(df.loc[find in df[to_list[0]]])
#print(df.loc[find in df[to_list[1]]])

# for col in df.columns:
# 	print(f"{col = }")
# 	print(df[col].head(3))
# 	print("\n\n\n")
