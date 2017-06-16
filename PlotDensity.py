import numpy as np
import pandas as pd
import seaborn as sns
import numpy as np
import scipy.stats as st
import plotly.plotly as py
from plotly.graph_objs import *   

import matplotlib.pyplot as plt

def kde_scipy( vals1, vals2, (a,b), (c,d), N ):
    
    #vals1, vals2 are the values of two variables (columns)
    #(a,b) interval for vals1; usually larger than (np.min(vals1), np.max(vals1))
    #(c,d) -"-          vals2 
    
    x=np.linspace(a,b,N)
    y=np.linspace(c,d,N)
    X,Y=np.meshgrid(x,y)
    positions = np.vstack([Y.ravel(), X.ravel()])

    values = np.vstack([vals1, vals2])
    kernel = st.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    
    return [x, y, Z]

    
def make_kdeplot(varX, varY, (a,b), (c,d), N, colorsc, title):
    #varX, varY are lists, 1d numpy.array(s), or dataframe columns, storing the values of two variables
   
    x, y, Z = kde_scipy(varY, varX, (a,b), (c,d), N )
    
    data = Data([
       Contour(
           z=Z, 
           x=x,
           y=y,
           colorscale=colorsc,
           #reversescale=True,
           opacity=0.9,    
           contours=Contours(
               showlines=False)      
        ),        
     ])

    layout = Layout(
        title= title,  
        font= Font(family='Georgia, serif',  color='#635F5D'),
        showlegend=False,
        autosize=False,
        width=650,
        height=650,
        xaxis=XAxis(
            range=[a,b],
            showgrid=False,
            nticks=7
        ),
        yaxis=YAxis(
            range=[c,d],
            showgrid=False,
            nticks=7
        ),
        margin=Margin(
            l=40,
            r=40,
            b=85,
            t=100,
        ),
    )
     
    return Figure( data=data, layout=layout )

                
N=1000
A = np.loadtxt('Mu_0_NoStrip_s0.dat')
c1 = A[:,7][400:]
c2 = A[:,1][400:]

fig=make_kdeplot(c1, c2, (min(c1),max(c1)), (min(c2),max(c2)), 
                 N, 'terrain_r','kde plot of two sets of data' )

py.sign_in('empet', 'my_api_key')
py.iplot(fig, filename='kde-2D-CSCE')                               
                                                
                                                                                