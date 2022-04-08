# -*- coding: utf-8 -*-
"""st1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1J6XoM5A3K10-5VhqbcsJWk5Ro9t_bXXD
"""

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

frames = [dfcurrent, dfancient]
dfcombined = pd.concat(frames)
Xcombined=dfcombined.drop(columns=['DNA sample ethnicity and id','DNA sample ethnicity','sample id'])
#dfcombined

dfcurrentgroup=dfcurrent.groupby(['DNA sample ethnicity']).mean().reset_index()
Xcurrentgroup=dfcurrentgroup.drop(columns=['DNA sample ethnicity'])


st.title('ML App')

st.dataframe(dfcurrentgroup)  

def euclidean_distance(p, q):
    n = 25 # dimensions
    return sum([(p[i] - q[i]) ** 2 for i in range(n)]) ** 0.5

#p=pd.read_csv('')
p=Xcombined.iloc[735:740]

dfdistances=dfcombined
distances=[]
#Induvidual

for i in range(len(Xcombined)):
  distances.append(euclidean_distance(Xcombined.iloc[i],p.iloc[3]))
dfdistances['distances']=distances

dfdistances=dfdistances.sort_values(by=['distances'])

dfdistances
print(dfdistances['DNA sample ethnicity and id'].iloc[:3])



c=pd.read_csv('clustergmm15.csv')
c=c.drop(columns=['Unnamed: 0'])

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

st.dataframe(ancestry)


uploaded_file = st.file_uploader("Choose a file")
if uploaded_file is not None:
     # To read file as bytes:
     bytes_data = uploaded_file.getvalue()
     st.write(bytes_data)

     # To convert to a string based IO:
     stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
     st.write(stringio)

     # To read file as string:
     string_data = stringio.read()
     st.write(string_data)

     # Can be used wherever a "file-like" object is accepted:
     dataframe = pd.read_csv(uploaded_file)
     st.write(dataframe)
    
    
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(30, 24))
ax.scatter(dfcurrentgroup['1'], dfcurrentgroup['2'],s = 1)
#ax.scatter(point['1'],point['2'],s=500)

for i in range(len(dfcurrentgroup)):
  ax.annotate(dfcurrentgroup['DNA sample ethnicity'][i], (dfcurrentgroup['1'][i], dfcurrentgroup['2'][i]))
