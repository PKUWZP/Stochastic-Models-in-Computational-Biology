#!/usr/bin/python
#!/opt/apps/software/MPI/GCC/4.4.7/OpenMPI/1.8.8/Python/2.7.9/bin/python


from scipy.integrate import odeint
from pylab import *
import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D

# Nn 0  In 1  N  2  I  3  (NI)n  4  (NI)  5  Im  6  BAD  7  UAD 8  BD  9  UD  10  OFF 11  BD_2-10 12-20 UD_2-10 21-29 
## BAD = bound artificial decoys
## UAD = unbound artificial decoys
## BD = bound natural decoys
## UD = unbound natural decoys

def deriv (y,t,k1,k2,k5,k6,k7,k8,k9,k10,k11,k12,k13,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,k49,k50,k51,k52,k53):
	
	#i = 51
	#file = open('param.in','r')
	#data = file.readlines()
	#for line in data:
	 #   words = line.split()
	  #  globals()["k"+str(i)] = float(words[1])
	  #  i+=1
	    
      #	for cv in range(1,25):
	#    globals()["k"+str(cv+26)] = 0.0001
	    ##globals()["k"+str(cv+50)] = 1
	 #   globals()["k"+str(cv+74)] = 0.0000
	    

    # all of above parameters come from Supplemental Materials in Krishna et al, PNAS 2006
        #for jc in range(1,25):
            
         #   sum1 =0
          #  sum1 = sum1 + globals()["k"+str(jc+50)]*y[jc+11]  ## unbinding reactions of decoys
           # sum2 =0
           # sum2 = sum2 + globals()["k"+str(jc+26)]*y[0]*y[jc+35] ## binding reactions of decoys
            #sum3 =0
           # sum3 =sum3 + globals()["k"+str(jc+74)]*y[1]*y[jc+11]  ## stripping reactions of decoys
            
            ##bounded decoys 
         #   globals()["dy"+str(jc+11)] = -globals()["k"+str(jc+50)]*y[jc+11] + globals()["k"+str(jc+26)]*y[jc+35]*y[0] - globals()["k"+str(jc+74)]*y[1]*y[jc+11]
            ## unbounded decoys
         #   globals()["dy"+str(jc+35)] = globals()["k"+str(jc+50)]*y[jc+11] - globals()["k"+str(jc+26)]*y[jc+35]*y[0] + globals()["k"+str(jc+74)]*y[1]*y[jc+11]
            

	dy0 = k1*y[2]-k2*y[0]*y[1]+k13*y[4]+k15*y[7]-k16*y[8]*y[0]+k17*y[9]-k18*y[0]*y[10]-k21*y[11]*y[0]+k22*(1-y[11])+k24*y[7]+k36*y[12]+k37*y[13]+k38*y[14]+k39*y[15]+k40*y[16]+k41*y[17]+k42*y[18]+k43*y[19]+k44*y[20]-k27*y[0]*y[21]-k28*y[0]*y[22]-k29*y[0]*y[23]-k30*y[0]*y[24]-k31*y[0]*y[25]-k32*y[0]*y[26]-k33*y[0]*y[27]-k34*y[0]*y[28]-k35*y[0]*y[29]

	dy1 = k9*y[3]-k10*y[1]-k2*y[0]*y[1]+k13*y[4]-k19*y[1]*y[7]-k20*y[1]*y[9]-k26*y[1]*(1-y[11])-k45*y[1]*y[12]-k46*y[1]*y[13]-k47*y[1]*y[14]-k48*y[1]*y[15]-k49*y[1]*y[16]-k50*y[1]*y[17]-k51*y[1]*y[18]-k52*y[1]*y[19]-k53*y[1]*y[20]
	
	dy2 = -k7*y[2]*y[3]+(k8+k11)*y[5]-k1*y[2]

	dy3 = k6*y[6]-k7*y[2]*y[3]+k8*y[5]-k9*y[3]+k10*y[1]

	dy4 = k2*y[0]*y[1]-(k13+k12)*y[4]+k19*y[1]*y[7]+k20*y[1]*y[9]+k26*y[1]*(1-y[11])+k45*y[1]*y[12]+k46*y[1]*y[13]+k47*y[1]*y[14]+k48*y[1]*y[15]+k49*y[1]*y[16]+k50*y[1]*y[17]+k51*y[1]*y[18]+k52*y[1]*y[19]+k53*y[1]*y[20]
	
	dy5 = k7*y[2]*y[3]-(k8+k11)*y[5]+k12*y[4]

	dy6 = k23*(1-y[11])-k5*y[6]

	dy7 = -k15*y[7]+k16*y[8]*y[0]-k19*y[1]*y[7]-k24*y[7]

	dy8 = k15*y[7]-k16*y[8]*y[0]+k19*y[1]*y[7]-k25*y[8]

	dy9 = -k17*y[9]+k18*y[0]*y[10]-k20*y[1]*y[9]

	dy10 = k17*y[9]-k18*y[0]*y[10]+k20*y[1]*y[9]

	dy11 = -k21*y[11]*y[0]+k22*(1-y[11])+k26*y[1]*(1-y[11])
	
	dy12 = -k36*y[12]+k27*y[21]*y[0]-k45*y[1]*y[12]
	
	dy13 = -k37*y[13]+k28*y[22]*y[0]-k46*y[1]*y[13]
	
	dy14 = -k38*y[14]+k29*y[23]*y[0]-k47*y[1]*y[14]
	
	dy15 = -k39*y[15]+k30*y[24]*y[0]-k48*y[1]*y[15]
	
	dy16 = -k40*y[16]+k31*y[25]*y[0]-k49*y[1]*y[16]
	
	dy17 = -k41*y[17]+k32*y[26]*y[0]-k50*y[1]*y[17]
	
	dy18 = -k42*y[18]+k33*y[27]*y[0]-k51*y[1]*y[18]
	
	dy19 = -k43*y[19]+k34*y[28]*y[0]-k52*y[1]*y[19]
	
	dy20 = -k44*y[20]+k35*y[29]*y[0]-k53*y[1]*y[20]
	
	dy21 = k36*y[12]-k27*y[21]*y[0]+k45*y[1]*y[12]
	
	dy22 = k37*y[13]-k28*y[22]*y[0]+k46*y[1]*y[13]
	
	dy23 = k38*y[14]-k29*y[23]*y[0]+k47*y[1]*y[14]
	
	dy24 = k39*y[15]-k30*y[24]*y[0]+k48*y[1]*y[15]
	
	dy25 = k40*y[16]-k31*y[25]*y[0]+k49*y[1]*y[16]
	
	dy26 = k41*y[17]-k32*y[26]*y[0]+k50*y[1]*y[17]
	
	dy27 = k42*y[18]-k33*y[27]*y[0]+k51*y[1]*y[18]
	
	dy28 = k43*y[19]-k34*y[28]*y[0]+k52*y[1]*y[19]
	
	dy29 = k44*y[20]-k35*y[29]*y[0]+k53*y[1]*y[20]
	
	
	return [dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16,dy17,dy18,dy19,dy20,dy21,dy22,dy23,dy24,dy25,dy26,dy27,dy28,dy29]

dt = 0.5   # time step
k1 = 5.4     # N -> Nn   kNin
k2 = 0.006    # Nn + In -> (NI)n   kf_n
#k3 = 0.03    #kbn 
#k4 = 1.03   # Nn -> Nn + mRNA 
k5 = 0.017   # gamma_m   mRNA -> 0
k6 = 0.24    # kt   mRNA->mRNA + Ic
k7 = 0.006    # kf    N + I -> (NI)
k8 = 0.03    # kb     (NI)-> N + I
k9 = 0.018   # kIin  Ic -> In
k10 = 0.012  # kIout   In -> Ic
k11 = 0.25   # alpha    (NI)-> N
k12 = 0.83   # kNIout    NIn -> NIc
k13 = 0.03   # kbn       (NI)n -> Nn + In 
#k14 = 1.030 # BA -> BA + Im 
############################################################
k15 = 0     # BAD -> UAD + Nn   unbinding of artificial decoys     k_adoff
k16 = 0.000  # UAD + Nn -> BAD       binding of artificial decoys   k_adon
############################################################
#k17 = 6.737947e-03     # BD -> UD + Nn     Decoy unbinding Kdoff
k17 = 1.234098e-04    # BD-> UD + Nn  Decoy unbinding Kdoff
k18 = 0.002  # Nn + UD -> BD   Decoy binding Kdon
############################################################
k19 = 0.000 # In + BAD -> UAD + INn    stripping of artificial decoys
k20 = 0.000  # In + BD -> INn + UD   stripping of decoys
############################################################
k21 = 0.002  # OFF + Nn -> ON    Promoter Binding
k22 = 1    # ON -> OFF + Nn      Promoter unbinding
############################################################
k23 = 5000  # ON -> ON + mRNA      transcription  
#k24 = 0.0  # BA + OFF -> ON + UD      antenna transfer
#k25 = 0  # ON + UA -> BA + OFF        antenna transfer
k24 = 0.00 # BAD -> Nn  degradation of bound-aritificial decoys
k25 = 0.00 # UAD -> 0   degradation of unbound-artificial decoys
k26 = 0.000  # In + ON -> OFF + INn      stripping on promoter
############################################################# binding rates  k_don
k27 = 0.002
k28 = 0.002
k29 = 0.002
k30 = 0.002
k31 = 0.002
k32 = 0.002
k33 = 0.002
k34 = 0.002
k35 = 0.002
############################################################# unbinding rates k_doff
#k36 = 4.978707e-02
k36 = 9.118820e-04
#k37 = 3.678794e-01
k37 = 6.737947e-03
#k38 = 2.718282e+00
k38 = 4.978707e-02
#k39 = 2.008554e+01
k39 = 3.678794e-01
#k40 = 1.484132e+02
k40 = 2.718282e+00
#k41 = 1.096633e+03
k41 = 2.008554e+01
#k42 = 8.103084e+03
k42 = 1.484132e+02
#k43 = 5.987414e+04
k43 = 1.096633e+03
#k44 = 4.424134e+05
k44 = 8.103084e+03
		
############################################################## stripping rates k_s
k45 = 0.000
k46 = 0.000
k47 = 0.000
k48 = 0.000
k49 = 0.000
k50 = 0.000
k51 = 0.000
k52 = 0.000
k53 = 0.000
	

def FindPeaks(list):
	result = []
	for i in range(1,len(list)-1):
		if list[i]>list[i-1] and list[i]>list[i+1]:
			result.append(list[i])
	return result
	
def CalculatePeriod(list,time):
    result = []
    for i in range(1,len(list)-1):
        if list[i]>list[i-1] and list[i]>list[i+1]:
           result.append(i)
    Period = time[result[2]]-time[result[1]]
    return Period

def DecayRate(list):
    result = []
    for i in range(1,len(list)-1):
        if list[i]>list[i-1] and list[i]>list[i+1]:
           result.append(list[i])
    DecayRate = -(result[len(result)-1]-result[1])/result[1]
    return DecayRate

def AveragePeakValue(list):
    result =[]
    for i in range(1,len(list)-1):
        if list[i]>list[i-1] and list[i]>list[i+1]:
           result.append(list[i])
    MeanPeak = mean(result[1:])
    return MeanPeak                            

def ListAllPeriod(list,time):
	periods = []
	times = []
	for i in range(1,len(list)-1):
		if list[i]>list[i-1] and list[i]>list[i+1]:
			times.append(time[i])
	for j in range(1,len(times)-1):
		periods.append((times[j+1]-times[j]))
	return periods	

## Read copy numbers of species from the input file
f1 = open('input400.in','r')
data1 = f1.readlines()
gh = 0
for line in data1:
    words = line.split()
    globals()["z"+str(gh)] = float(words[1])
    gh+=1

#time = np.linspace(0.0, 10000.0, 10000)
#time1 = np.linspace(0.0,2000,2000)
X=0  ## artifical decoys
#yinit = array([0,0,0,0,0,100000,0,0,X,0,z0,1,z1,0,0,0,0,0,0,0,0,0,z2,z3,z4,z5,z6,z7,z8,z9])
yinit = np.array([0,0,0,0,0,80000,0,0,X,z0,0,1,z1,z2,z3,z4,z5,z6,z7,z8,z9,0,0,0,0,0,0,0,0,0])
#time = np.linspace(0.0, 6000.0, 6000/dt)
time = np.linspace(0.0,10000,10001)
## Set up the grid size
#GS=101
#GI = 1001
#mamp = np.zeros((GS,GS)) ## Oscillation Amplitude
#Aped = np.zeros((GS,GS)) ## Oscillation Period

#S1 = np.linspace(0.001,1,GS) ## k_off
#S2 = np.linspace(0.0,100,GI) ## k_stripping
#S3 = np.linspace(0.001,1,GS) ## k_adoff 
#S4 = np.linspace(0,100000,GI) ## Decoy Number
#S5 = np.linspace(0,0.02,GS) ## gamma degradation rate
#cc = 0  ## Progress representation
#perc = np.array([0])

#for sc in range(GS):
#    for sc1 in range(GS):
#        k15 = S3[sc]
#        k24 = S5[sc1]
#        k25 = S5[sc1]
y = odeint (deriv, yinit, time, args = (k1,k2,k5,k6,k7,k8,k9,k10,k11,k12,k13,k15,k16,k17,k18,k19,k20,k21,k22,k23,k24,k25,k26,
        k27,k28,k29,k30,k31,k32,k33,k34,k35,k36,k37,k38,k39,k40,k41,k42,k43,k44,k45,k46,k47,k48,k49,k50,k51,k52,k53), mxstep=500000000)
#        data = y[:,0][1000:]
#        stime = time [1000:]
#        mamp[sc][sc1] = max(data) - min(data)
#        cc+=1
#        perc[0] = (cc*101.0/(GS*GS))
#        savetxt('progress_bar.txt',perc)

#savetxt('amp_out.txt',mamp)
#savetxt('period_out.txt',Aped)
#savetxt('scan1.out',S3)
#savetxt('scan2.out',S5)
                
#Z = (y[:,0] + y[:,4])/(y[:,2] + y[:,5])
#Z1 = (y[:,0])/(y[:,2])
#Z2 = y[:,0] + y[:,2] + y[:,4] + y[:,5]
#np.savetxt('NnNc_Ratio.txt',)
#y0mean = mean (y[:,0][1000:6000])
#y0max = max (y[:,0][1000:6000])
#y0min = min (y[:,0][1000:6000])
#Z = (y0max - y0min)/y0mean
#count = 0
#for i in range (0, len(y[:,0])):
#	if y[:,0][i]>y0mean:
#		count=count+1	
#print count

#print "Spikeness is:", Z

#print "Average Peak Value is:",AveragePeakValue(y[:,0]/1000)
#print "Decay Rate is:", DecayRate(y[:,0]/1000)
#print "Period is (min):", CalculatePeriod(y[:,0]/1000,time)
#print "Peak Time:", ListAllPeriod(y[:,0]/1000,time)
#print FindPeaks(y[:,0]/1000)

#print Z
#print y0min

#figure()
#plot (time,y[:,0]/100000,linewidth=3.0) #label = r"IkB-induced stripping rate: $40\mu M^{-1} min^{-1}$")
#plot (time1,y[:,3]/10000, label = r'$I\kappa B$')
#plot(y[:,0]/100000, y[:,1]/100000)
#xlim([0,0.05])
#rc('text',usetex=True)
#rc('font',family='serif')
#plot (time, y[:,3]/1000)
#xlim([-1,1])
#ylim([-1,4])
#xlabel (r'\textbf{time} (min)',fontsize = 15)
#ylabel (r"$NF\kappa B_n (\mu M)$",fontsize = 16)
#ylabel (r"$I\kappa B (\mu M)$", fontsize = 16)
#ylabel (r"$NF\kappa B_n (\mu M);I\kappa B (10\mu M)$", fontsize=16)
#legend (loc = 'upper right')
#legend (loc = 'upper right')
#savefig('Fig3_6_1000.eps',format='eps',dpi=1000)
#savefig('PNASFIG4.eps')
#show()
savetxt('mRNA_std6.txt',y[:,6])
savetxt('Db_std6.txt',y[:,9]+y[:,12]+y[:,13]+y[:,14]+y[:,15]+y[:,16]+y[:,17]+y[:,18]+y[:,19]+y[:,20])