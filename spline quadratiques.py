print('interpolation polymiale par morceaux:methode des splines quadratique')
N=int(input('veuillez entre le nombre de points d interpolation '))
A=float(input('veuillez entre la borne inferieure de votre intervalle'))
B=float(input('veuillez entre la borne superieure de votre intervalle'))
def fonction(t,n,a,b):
    import numpy as np
    def m(r):
        m=np.log(r)+np.cos(r)
        return m 
    def idct(x,y,r):
        if x<=r and r<y :
            idct=1
        else:
            idct=0
        return idct
    h=(b-a)/(n-1)
    T=np.zeros(n)
    R=np.zeros(n)
    S=np.zeros(n-1)
    G=np.zeros(n)
    for i in range(n):
        T[i]=a+i*h

    R[0]=2*((m(T[1])-m(T[0]))/h)-(1/a-np.sin(a))
    for i in range(1,n-1):
        R[i]=2*((m(T[i])-m(T[i-1]))/h)

    G[0]=R[0]
    for i in range(1,n-1):
        G[i]=R[i]-G[i-1]

    R[n-1]=2*((m(T[n-1])-m(T[n-2]))/h)-G[n-2]
    G[n-1]=R[n-1]
    for i in range(n-1):
        S[i]=(G[i+1]*(t-T[i])**2)/(2*h)-(G[i]*(t-T[i+1])**2)/(2*h)+m(T[i])+G[i]*h/2

    z=0
    for i in range(n-1):
        z+=idct(T[i],T[i+1],t)*S[i]

    return z
import numpy as np
import matplotlib.pyplot as plt 
k=np.linspace(A,B,300)
w=np.zeros(300)
for i in range(299):
    w[i]=fonction(k[i],N,A,B)

w[299]=np.log(B)+np.cos(B)
x=np.linspace(A,B,300)
y=np.log(x)+np.cos(x)
plt.figure()
plt.plot(k,w,c='red')
plt.plot(x,y)
plt.show()
