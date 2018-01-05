from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import odeint
from numpy.fft import irfft, rfft


def wanda01(x,t,k):
    rx=0.02
    wx=3
    ry=0.03
    wy=5
    nx0 = x[1]
    nx1 = (rx - x[0]**2) * x[1] - wx*x[0] + k*(x[2] - x[0])
    ny0 = x[3]
    ny1 = (ry - x[2]**2) * x[3] - wy*x[2]
    res = array([nx0, nx1, ny0, ny1])
    return res
    
    
tt=linspace(0, 1000, 100000)
k=[0, 0.2, 0.4, 0.6, 0.8, 1]
coef=zeros(len(k))

for i in range(len(k)):    
    res=odeint(wanda01,[-0.2, 0.1, 0.4, 0.2],tt,args=(k[i],))
    
    sopres=irfft(1j*rfft(res))
        
    # phase analisys 
    fi_x=zeros(len(res[:, 0]))
    ksi_x=arctan(sopres[:, 0]/res[:, 0])
    dob_x=0
        
    for el in range(len(res[:, 0])-1):
        if ksi_x[el+1]<ksi_x[el]:
            dob_x+=1
        fi_x[el]=ksi_x[el]+dob_x*pi
            
        
    fi_y=zeros(len(res[:, 2]))
    ksi_y=arctan(sopres[:, 2]/res[:, 2])
    dob_y=0
        
    for el in range(len(res[:, 2])-1):
        if ksi_y[el+1]<ksi_y[el]:
            dob_y+=1
        fi_y[el]=ksi_y[el]+dob_y*pi
        
        
    #coef phase synhro
    defi=fi_x - fi_y
    
    coef[i]=abs(mean(exp(1j*defi)))


xlabel('k')
ylabel('coef')
plot(k, coef, 'r', lw = 2)
grid()
show()