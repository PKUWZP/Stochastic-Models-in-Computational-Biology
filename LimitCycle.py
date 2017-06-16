#!/usr/bin/python

from numpy import *
from pylab import *
import pandas as pd
from scipy import stats, integrate
import seaborn as sns

A = loadtxt('Mu_3_std6_alpha0250.dat')   ## Fast Decoys
c1 = array(A[:,7][8400:]) ## mRNA
c2 = array(A[:,1][8400:]) ## Nn
c3 = array(A[:,2][8400:]) ## In 
c4 = array(A[:,10][8400:]+A[:,11][8400:]+A[:,12][8400:]+A[:,13][8400:]+A[:,14][8400:]+A[:,15][8400:]+A[:,16][8400:]+A[:,17][8400:]+A[:,18][8400:]+A[:,19][8400:]+A[:,20][8400:]) ## Db

#B = loadtxt('Mu_0_NoStrip_s0.dat')   ## Resonance Decoys
#b1 = array(B[:,7][400:50000]) ## mRNA
#b2 = array(B[:,1][400:50000]) ## Nn
#b3 = array(B[:,2][400:50000]) ## In 
#b4 = array(B[:,11][400:50000]) ## Db

#C = loadtxt('Mu_-4_NoStrip_s0.dat')
#d1 = array(C[:,7][400:50000]) ## mRNA
#d2 = array(C[:,1][400:50000]) ## Nn
#d3 = array(C[:,2][400:50000]) ## In 
#d4 = array(C[:,11][400:50000]) ## Db


#C1,C2 = [],[]
#for i in range(400,16401,200):
#    C1.append(sort(c1[i:i+200]))
#    C2.append(sort(c2[i:i+200]))

#C1 = array(C1)
#C2 = array(C2)

#r = []
#rho = []
#r_temp = []
#r_temp2 = []
#r_temp3 = []
#r_temp4 = []

#for theta in arange(0,2.05*pi,0.01*pi):
#    for i in range(0,len(c1)):
#        if c1[i]<=200000 and c2[i]>=10000:
#           theta1 = arctan(abs(c2[i]-10000)/abs(c1[i]-200000)) 
#           if theta1>=theta and theta1<=theta+0.01*pi:
#              #rho.append(theta)
#              r_temp.append(((c1[i]-200000)**2 + (c2[i]-10000)**2)**0.5)

 #       if c1[i]>=200000 and c2[i]>=10000:
 #          theta2 = arctan(abs(c2[i]-10000)/abs(c1[i]-200000))
 #          theta2 = pi - theta2
 #          if theta2>=theta and theta2<=theta+0.01*pi:
 #             #rho.append(theta)
 #             r_temp2.append(((c1[i]-200000)**2 + (c2[i]-10000)**2)**0.5)

 #       if c1[i]>=200000 and c2[i]<=10000:
  #         theta3 = arctan(abs(c2[i]-10000)/abs(c1[i]-200000))
  #         theta3 = pi + theta3
  #         if theta3>=theta and theta3<=theta+0.01*pi:
              #rho.append(theta)
  #            r_temp3.append(((c1[i]-200000)**2 + (c2[i]-10000)**2)**0.5)
           
  #      if c1[i]<=200000 and c2[i]<=10000:
  #         theta4 = arctan(abs(c2[i]-10000)/abs(c1[i]-200000))
   #        theta4 = 2*pi - theta4
    #       if theta4>=theta and theta4<=theta+0.01*pi:
     #         #rho.append(theta)
      #        r_temp4.append(((c1[i]-200000)**2 + (c2[i]-10000)**2)**0.5)

  #  if len(r_temp4)>0:
   #    r.append(mean(r_temp4))
   #    rho.append(theta)
   #    r_temp4 = []
    
   # if len(r_temp)>0:
   #    r.append(mean(r_temp))
   #    rho.append(theta)
   #    r_temp = []

 #   if len(r_temp2)>0:
 #      r.append(mean(r_temp2))
 #      rho.append(theta)
 #      r_temp2 = []
    
 #   if len(r_temp3)>0:
 #      r.append(mean(r_temp3))
 #      rho.append(theta)
 #      r_temp3 = []
       

def av(X):
    av_X=0
    for i in range(len(X)):
        av_X += X[i]/len(X)
    return av_X

#print len(rho)
#print len(r)
#c1m = av(C1)
#c2m = av(C2)
data1 = column_stack((c1,c4)) #1. mRNA, 2. Nn 3. In 4. Db. (Db-mRNA) 
data2 = column_stack((c2,c4))  ## Db-Nn
data3 = column_stack((c3,c4))  ## Db-In
data4 = column_stack((c1,c2))  ## mRNA-Nn
data5 = column_stack((c1,c3))  ## mRNA-In
data6 = column_stack((c2,c3))  ## Nn-In

df = pd.DataFrame(data1, columns=["x", "y"])
#print c1m
#sns.jointplot(x="x", y="y", data=df, kind="kde");

#cmap = sns.cubehelix_palette(as_cmap=True,'PuBu',rot=-.1,dark=0.3,light=1,reverse=True) #dark=0, light=1, reverse=True)
sns.kdeplot(df.x, df.y, cmap='terrain_r', n_levels=60, shade=True)
ylabel(r"$NF\kappa B-DNA$",fontsize=20,labelpad=0)
xlabel(r"$N_n$",fontsize=20,labelpad=0)
#xlabel(r"$mRNA$",fontsize=20,labelpad=0)
title(r"$k_{doff} = 0.018 min^{-1}, k_{s} = 10 \mu M^{-1} min^{-1}$",fontsize=20)
ax = gca()
ax.set_axis_bgcolor('white')
ax.tick_params(axis='both', which='major', labelsize=15)
#ax=subplot(111,projection='polar')
#ax.plot(rho, r, color='r', linewidth=3)
#plot(r*cos(rho)+200000,r*sin(rho)+10000,linewidth=3)
#plot(c1[450:500],c2[450:500])
#plot(c1,c2)
#xlim([130000,240000])
#xlim(-4000,30000)
#ylim([-3000,23000])
show()
