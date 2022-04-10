import streamlit as st
from threading import activeCount
import matplotlib.pyplot as plt
#import umap
import io
import numpy as np
import pandas as pd

dfcurrent=pd.read_csv("G25_Current_DNA.csv")
Xcurrent=dfcurrent.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])

dfancient=pd.read_csv("G25_Ancient_DNA.csv")
Xancient=dfancient.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])
#Xcurrent

c=pd.read_csv('clustergmm15.csv')
c=c.drop(columns=['Unnamed: 0'])

frames = [dfcurrent, dfancient]
dfcombined = pd.concat(frames)
Xcombined=dfcombined.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])
#dfcombined

dfcurrentgroup=dfcurrent.groupby(['DNA sample ethnicity']).mean().reset_index()
Xcurrentgroup=dfcurrentgroup.drop(columns=['DNA sample ethnicity'])

dfadnalineages=pd.read_csv("adnalineages.csv")
dfancientpcadna = pd.merge(dfadnalineages,dfancient)
dfancienthpg=dfancientpcadna.groupby(['Assigned Mutation']).mean().reset_index()

def euclidean_distance(p, q):
    n = 25 # dimensions
    return sum([(p[i] - q[i]) ** 2 for i in range(n)]) ** 0.5

p=Xcombined.iloc[735:740]

dfdistances=dfcombined
distances=[]
#Induvidual

for i in range(len(Xcombined)):
  distances.append(euclidean_distance(Xcombined.iloc[i],p.iloc[3]))
dfdistances['distances']=distances

dfdistances=dfdistances.sort_values(by=['distances'])

#p=pd.read_csv('')
p=Xcombined.iloc[735:740]

import numpy as np
p1 = np.zeros((len(p),len(c)))
p2 = np.zeros((len(p),len(c)))

for i in range(len(p)):
  for j in range(len(c)):
    p1[i][j]=euclidean_distance(c.iloc[j],p.iloc[i])
distmat=pd.DataFrame(p1)

for q in range(len(p)):
  tot=distmat.iloc[q].sum()

  for w in range(len(c)):
    p2[q,w]=1-((tot-distmat.iloc[q,w])/tot)

ancestry=pd.DataFrame(p2)


#dfdistances
#print(dfdistances['DNA sample ethnicity and id'].iloc[:3])

Tools = st.selectbox("Choose your Tool", ["Distance Tool", "ML Ancestry Tool", "Ancient DNA Lineage Tool"]) 

if page == "Distance Tool":
  st.dataframe(dfdistances)

elif page == "ML Ancestry Tool":
  st.dataframe(ancestry)

elif page == "Ancient DNA Lineage Tool":

  fig2, ax2 = plt.subplots(figsize=(30, 24))
  ax2.scatter(dfancienthpg['1'], dfancienthpg['2'],s = 1)
  #ax.scatter(point['1'],point['2'],s=500)
  for i in range(len(dfancienthpg)):
    ax2.annotate(dfancienthpg['Assigned Mutation'][i], (dfancienthpg['1'][i], dfancienthpg['2'][i]))
  st.pyplot()

elif page == "Page 4":
