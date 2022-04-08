# -*- coding: utf-8 -*-
"""st1.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1J6XoM5A3K10-5VhqbcsJWk5Ro9t_bXXD
"""

pandas
matplotlib
pydeck
streamlit
flexidate

import streamlit as st
from threading import activeCount
import streamlit.components.v1 as components
import plotly.figure_factory as ff
#import umap
import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


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

def PCAplot():
  fig, ax = plt.subplots()
  ax.scatter(dfcurrentgroup['1'], dfcurrentgroup['2'],s = 1)
  #ax.scatter(point['1'],point['2'],s=500)

  for i in range(len(dfcurrentgroup)):
    ax.annotate(dfcurrentgroup['DNA sample ethnicity'][i], (dfcurrentgroup['1'][i], dfcurrentgroup['2'][i]))

  st.pyplot
PCAplot()
