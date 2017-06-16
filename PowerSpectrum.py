#!/usr/bin/python

from numpy import *
from pylab import *
from scipy import signal

style.use('classic')
rc("axes", linewidth=3.0)
rc("lines", markeredgewidth=4.0)
rc('xtick', labelsize=28)
rc('ytick', labelsize=28)
rcParams['lines.linewidth'] = 5
rcParams['figure.figsize'] = 11, 10

sim_t = 1

def pow_spec(X):
    from numpy import *
    Z=fft.fft(X)
    S=abs(Z)**2
    Omega = fft.fftfreq(len(S),sim_t)
    #f,Pxx_spec = signal.periodogram(X, len(X))
    #return (f,Pxx_spec)
    return (Omega, S)

#X1 = loadtxt('Autocorr_-4_NoStrip_s0.txt')
#X2 = loadtxt('Autocorr_0_NoStrip_s0.txt')
#X3 = loadtxt('Autocorr_3_NoStrip_s0.txt')
#X1 = loadtxt('Mu_0_NoStrip_s0.dat')
X1 = loadtxt('Autocorr_0_NoStrip_s6.txt')
#X1 = loadtxt('Autocorr_3_std6_alpha025.txt')
X2 = loadtxt('Autocorr_0_stp10_s6.txt')
#X2 = loadtxt('Autocorr_3_std6_stp10_alpha025.txt')

f1 = pow_spec(X1)
f2 = pow_spec(X2)
#f3 = pow_spec(X3)
print len(f1[0])
print len(f2[1])

#plot(f1[0],f1[1],linewidth=3,label=r"$k_{doff} = 0.018 min^{-1}$")
#plot(f2[0],f2[1],linewidth=3,label=r"$k_{doff} = 1 min^{-1}$")
#plot(f3[0],f3[1],linewidth=3,label=r"$k_{doff} = 20 min^{-1}$")
#plt.semilogy(f1[0], f1[1], basey=10,color='darkred', linewidth = 3)
plot(f1[0]/2.2,f1[1],linewidth=4,color='b',label=r"$\sigma^2 = 6, k_{s} = 0 \mu M^{-1} min^{-1}$")
plot(f2[0]/2.2,f2[1],linewidth=4,color='r',label=r"$\sigma^2 = 6, k_{s} = 10 \mu M^{-1} min^{-1}$")

#xlim([0.015,0.026])
xlim([0.007,0.015])
ylim([0,1.6e5])
#title(r"$k_{doff} = 1 min^{-1}$",fontsize=30, family="Times New Roman",y=1.04)
xlabel(r"$\Omega [min^{-1}]$",fontsize=36,family="Times New Roman")
ylabel(r"$S(\Omega)$",fontsize=40,family="Times New Roman",labelpad=0)
ticklabel_format(style='sci', axis='both', scilimits=(0,0))
tight_layout()
#ax = gca()
#ax.set_yscale('log')
lg=legend(loc = "upper right",prop={'family':'Times New Roman','size':25}, fancybox = True)
lg.set_frame_on(True)
savefig('NewFIG3_D.pdf',format='pdf')


show()
