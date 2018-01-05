from numpy import *
from matplotlib.pyplot import *
from scipy.integrate import odeint
from numpy.fft import irfft, rfft
from mpl_toolkits.mplot3d import axes3d


#global parameters
a1=0.22
b1=0.15
c1=3
a2=0.2
b2=0.15
c2=3.5

cl = 10000


    
def ressler(r,t,k):
    global a1,b1,c1,a2,b2,c2
    x01 = r[0]
    y01 = r[1]
    z01 = r[2]
    x02 = r[3]
    y02 = r[4]
    z02 = r[5]
    
    dx1 = -y01-z01
    dx2 = -y02-z02 + k*(dx1-y02)
    dy1 = x01 + a1*y01
    dy2 = x02 + a2*y02
    dz1 = b1 - z01*(c1 - x01)
    dz2 = b2 - z02*(c2 - x02)
    res = array([dx1,dy1,dz1,dx2,dy2,dz2])
    return res
        
tt=linspace(0,500,100000)
k=linspace(0, 0.1, 10)
coef=zeros(len(k))
    
for i in range(len(k)):
    att=odeint(ressler,[0.1, 0.1, 0.1, 0.1, 0.1, 0.1],tt,args=(k[i],))
    sopres=irfft(1j*rfft(att))
    # phase analysis
    fi_x=zeros(len(att[cl:,0]))
    ksi_x=arctan(sopres[cl:,0]/att[cl:,0])
    dob_x=0
        
    for el in range(len(att[cl:,0])-1):
        if ksi_x[el+1]<ksi_x[el]:
            dob_x+=1
        fi_x[el]=ksi_x[el]+dob_x*pi
            
        
    fi_y=zeros(len(att[cl:,3]))
    ksi_y=arctan(sopres[cl:,3]/att[cl:,3])
    dob_y=0
        
    for el in range(len(att[cl:,3])-1):
        if ksi_y[el+1]<ksi_y[el]:
            dob_y+=1
        fi_y[el]=ksi_y[el]+dob_y*pi
        
        
    #coef phase synhro
    defi=fi_x - fi_y
    
    coef[i]=abs(mean(exp(1j*defi)))



# fig = figure(str(coef))
# ax = fig.gca(projection='3d')
# ax.plot(att[cl:,0],att[cl:,1],att[cl:,2],color='black')
# ax.plot(att[cl:,3],att[cl:,4],att[cl:,5],color='red')

figure('coef phase synhro')
plot(k, coef,  'm', lw = 2)
grid()
xlabel('k')
ylabel('coef')
show()