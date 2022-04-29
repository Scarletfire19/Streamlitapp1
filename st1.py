import streamlit as st
from threading import activeCount
import matplotlib.pyplot as plt
import pickle
import umap
import io
import numpy as np
import pandas as pd
import streamlit_authenticator as stauth
import os
import plotly.express as px
from scipy.spatial import distance
import plotly.graph_objects as go



st.set_option('deprecation.showPyplotGlobalUse', False)

#https://wallpapercave.com/wp/wp4594030.jpg
#https://wallpapercave.com/wp/wp2170040.jpg
#https://blog.23andme.com/wp-content/uploads/2008/09/novembreblogpostfig.jpg
#https://media.istockphoto.com/photos/background-picture-id185314478?b=1&k=20&m=185314478&s=612x612&w=0&h=XTWvVxe3IdiV6B5bJDi2DAICfqmyaZxWnWyXY_Lw5kk=

page_bg_img = '''
<style>
body {
background-image: url("https://images.unsplash.com/photo-1542281286-9e0a16bb7366")
background-size: cover;
}
</style>
'''

st.markdown(page_bg_img, unsafe_allow_html=True)

st.title("Machine Learning Genomics App")


# streamlit_app.py


#names = ['Blue Bear','Cream Rabbit']
#usernames = ['bhalu','khwabekhargosh']
#passwords = ['bear','rabbit']
#hashed_passwords = stauth.Hasher(passwords).generate()
#authenticator = stauth.Authenticate(names,usernames,hashed_passwords,'some_cookie_name','some_signature_key',cookie_expiry_days=30)
#name, authentication_status, username = authenticator.login('Login','main')

#if authentication_status:
#    st.write('Welcome *%s*' % (name))
#    st.title('Some content')
#elif authentication_status == False:
#    st.error('Username/password is incorrect')
#elif authentication_status == None:
#    st.warning('Please enter your username and password')

#if authentication_status:
#    authenticator.logout('Logout', 'main')
#    st.write('Welcome *%s*' % (name))
#    st.title('Some content')
#elif authentication_status == False:
#    st.error('Username/password is incorrect')

dfadnalatlonglineage=pd.read_csv('adnalatlonglineage.csv')
dfadnalatlonglineage['Long'] = pd.to_numeric(dfadnalatlonglineage['Longitude'],errors='coerce')
dfadnalatlonglineage['Lat']=pd.to_numeric(dfadnalatlonglineage['Latitude'],errors='coerce')

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

dfumap=pd.read_csv('umap.csv')
dftsne=pd.read_csv('tSNE.csv')


uploaded_file = st.file_uploader("Enter G25 co-ordinates")
if uploaded_file is not None:
     input = pd.read_csv(uploaded_file)
     st.write(input)
      
          
inputhaplogroup = st.text_input("Enter Y Haplogroup")

try:
  if len(input)>10:
     st.write("Proceed with ML Genomics Tools")
except:
     st.write("Please enter input")

def euclidean_distance(w, q):
    n = 25 
    return sum([(w[i] - q[i]) ** 2 for i in range(n)]) ** 0.5

p=Xcombined.iloc[735:740]
dfdistances=dfcombined
distances=[]



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
k=0


#dfdistances
#print(dfdistances['DNA sample ethnicity and id'].iloc[:3])

Tools = st.selectbox("Choose your Tool", ["Genetic World Map"]) 

if Tools == "Genetic World Map":
     st.title("Genetic World Map")
     figworldmap = px.scatter_mapbox(dfadnalatlonglineage, lat="Lat", lon="Long", hover_name="Sample ID/Group/Community", hover_data=["Location","Y chromosome haplogroup", "mtDNA haplogroup","Average Confidence Interval Date","Dating","Gender","Published Study"], color_discrete_sequence=["maroon"], zoom=3, height=1200)
     figworldmap.update_layout(mapbox_style="stamen-terrain")
     st.plotly_chart(figworldmap, use_container_width=True)
