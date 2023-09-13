
#start with a dataframe of the VMH identified metabolites and their logged fold change value 
## I have called this the "sample" dataframe, so you can also name yours sample 


import pandas as pd


#use the path wherever you save the useful names file i sent you:
useful_names = {}
with open('/Users/theogoodchild-michelman/Desktop/research/summer_research/metaboliteVMH_readable.txt') as f:
    next(f) 
    for key, *values in csv.reader(f, delimiter='\t'):
        key_met = key.replace("EX_", "_")
        key_met = key_met.replace("[e]", "[u]")

        useful_names[key_met] = values

sample_realnames = []
counter = 0

sample["real names"] = sample_realnames
sample.set_index("real names")


for i in sample.index:
    met = i.replace("_", "", 1)
    met = met.replace("[u]", "")
    met = met.replace("[", "")
    met = met.replace("]", "")
    if met in useful_names:
        sample_realnames.append(useful_names[met])

#saving the logged fold change data to look at
sample.sort_values("abs", ascending=False).to_csv("logged_foldchange.csv")

#import your fold change data
logged_df = pd.read_csv("logged_foldchange.csv")
logged_df = logged_df.sort_values("abs", ascending=False)


#you can change 45 to however many metabolites you want to show
logged_df = logged_df.head(45)
logged_df = logged_df.set_index("real names")

logged_df['positive'] = logged_df['median logged values'] > 0
logged_df = logged_df.sort_values('median logged values')

#making the graph titles look nice
logged_df.index = [str(i) for i in logged_df.index]
logged_df.index = [i.replace("[","") for i in logged_df.index]
logged_df.index = [i.replace("]","") for i in logged_df.index]
logged_df.index = [i.replace("'","") for i in logged_df.index]


#actually plotting the graph now"
import matplotlib.pyplot as plt


plt.rcParams["figure.figsize"] = (27,13)
plt.rcParams.update({'font.size': 18})

fig, ax = plt.subplots()

logged_df['median logged values'].plot(kind='barh',
                             color=logged_df.positive.map({True: 'r', False: 'b'}))

plt.savefig('foo.png', facecolor='w',bbox_inches="tight")

