from pylab import *

def Psi(rho):
    #return 1.0-exp(-rho)
    return exp(-1.0/rho)

X = range(10)

P1 = range(10)
P2 = range(10)
P3 = range(10)
P4 = range(10)
P5 = range(10)

ps = range(10)
G = -15
for i in range(len(X)):
    X[i] = X[i]*0.1+0.1
    r = X[i]
    P1[i] = r/3.0 + -1.0/6.0*Psi(r)**2
    P2[i] = r/3.0 + -2.0/6.0*Psi(r)**2
    P3[i] = r/3.0 + -4.0/6.0*Psi(r)**2
    P4[i] = r/3.0 + -8.0/6.0*Psi(r)**2
    P5[i] = r/3.0 + -16.0/6.0*Psi(r)**2
    ps[i]= Psi(r)

subplot(3,3,1)
plot(X,ps)
title('aaa')
subplot(3,3,4)
plot(X,P1)
subplot(3,3,5)
plot(X,P2)
subplot(3,3,6)
plot(X,P3)
subplot(3,3,7)
plot(X,P4)
show()



