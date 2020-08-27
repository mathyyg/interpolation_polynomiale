import numpy as np
import numpy.polynomial.polynomial as nppol
import matplotlib.pyplot as plt

"""def g(x):
    return (2*x(**3)-5*x(**2)+3*x-4)
"""

def f(x):
    return (5*x*np.sin(x)/(1+x*x))

def eval_cubic_spline(x,cs) :
    for borne,poly in cs :
        if x <= borne :
            return nppol.polyval(x,poly)
    return 0

def get_system_to_solve(px,py) :
    t3,t2 = len(px),len(py)
    ##t1=len(px)
    rt=(4*t3-4,4*t3-4)
    matriceA=np.zeros(rt)
    ##interpolation
    matriceA[0,0]=px[0]
    matriceA[0,1]=px[0]
    matriceA[0,2]=px[0]
    matriceA[0,3]=1
    k=1
    #k++
    for i in range(1,t3-1):
        ##partie1
        matriceA[k,(i-1)*4+0]=px[i]**3
        matriceA[k,(i-1)*4+1]=px[i]**2
        matriceA[k,(i-1)*4+2]=px[i]
        matriceA[k,(i-1)*4+3]=1        
        ##partie2
        matriceA[k+1,(i-1)*4+0+4]=px[i]**3
        matriceA[k+1,(i-1)*4+1+4]=px[i]**2
        matriceA[k+1,(i-1)*4+2+4]=px[i]
        matriceA[k+1,(i-1)*4+3+4]=1
        k+=1
        k+=1
        #k++
        ##partie3
    matriceA[k,i*4+0]=px[i+1]**3
    matriceA[k,i*4+1]=px[i+1]**2
    matriceA[k,i*4+2]=px[i+1]
    matriceA[k,i*4+3]=1
    k+=1
    #k++
##derivabilite
    for i in range(0,t3-2):
        ##partie1
        matriceA[k,i*4+0]=3*px[i+1]**2
        matriceA[k,i*4+1]=2*px[i+1]
        matriceA[k,i*4+2]=1
        matriceA[k,i*4+3]=0
        ##partie2
        matriceA[k,i*4+4]=-3*px[i+1]**2
        matriceA[k,i*4+5]=-2*px[i+1]
        matriceA[k,i*4+6]=-1
        matriceA[k,i*4+7]=0
        k+=1
    """derivabilite 
    seconde
"""
    for i in range(0,t3-2):
        ##partie3
        matriceA[k,i*4+0]=6*px[i+1]
        matriceA[k,i*4+1]=2
        matriceA[k,i*4+4]=-6*px[i+1]
        matriceA[k,i*4+5]=-2
        k+=1
        #k++
##cond.arbitraire
##partie1
    t4=len(matriceA)
    matriceA[k,0]=3*px[0]**2
    matriceA[k,1]=2*px[0]
    matriceA[k,2]=1
    matriceA[k,t4-4]=-3*px[t3-1]**2
    matriceA[k,t4-3]=-2*px[t3-1]
    matriceA[k,t4-2]=-1
##partie2
    matriceA[k+1,0]=6*px[0]
    matriceA[k+1,1]=2
    matriceA[k+1,t4-4]=-6*px[t3-1]
    matriceA[k+1,t4-3]=-2
##matriceB
    matriceB=np.zeros(4*t2-4)
    matriceB[0]=py[0]
    h=1
    for i in range(1,t2-1):
        matriceB[h] = py[i]
        matriceB[h+1] = py[i]
        h+=1
        h+=1
    matriceB[h]=py[i+1]
    return matriceA,matriceB

def get_cubic_splines(px,py) :
    getSplines = []
    rt = len(px)
    mA,mB = get_system_to_solve(px,py)
    inverseA = np.linalg.inv(mA) 
    res = inverseA@mB
    for i in range(0,rt-1) :
        getSplines.append((px[i+1], [ res[3+i*4],res[2+i*4],res[1+i*4],res[0+i*4]]))
    return getSplines

"""f"""
x_1=np.linspace(0,10,100)
y_1=np.array(list(map(f,x_1)))
c_f1,=plt.plot(x_1,y_1,'r',linewidth=1,color="red")
c_f1.set_label("fonction de base")
"""px et py"""
px_alt=np.linspace(0,10,6)
px=list(px_alt)
rt2=len(px)
py=[0 for i in range(rt2)]
for i in range(rt2):
    py[i] = f(px[i])

"""vérification points"""
lv = plt.scatter(px,py,c="green")
lv.set_label("points de vérification")
"""fenetre"""
plt.gcf().canvas.set_window_title('Mathys Gagner - Mateo Esteves')

"""affichage fin"""
cs = get_cubic_splines(px,py)
x_2 = np.linspace(0,10,100)
y_2=[eval_cubic_spline(e,cs) for e in x_2]
c,=plt.plot(x_2,y_2)
c.set_label("approximation")
plt.grid(color='black',linestyle='--')
plt.legend()
plt.show()
