import numpy as np
import numpy.polynomial.polynomial as nppol
import matplotlib.pyplot as plt

def f(x) :
    return 5*x*np.sin(x)/(1+x*x)

def Q(racines,i) :
    temp=racines.copy()
    x_t=temp.pop(i)
    coef=1/nppol.polyvalfromroots(x_t,temp)
    return (coef*nppol.polyfromroots(temp))

def base_lagrange(px) :
    res_tab = []
    rt = len(px)
    for k in range(rt):
        res_tab.append(Q(px,k))
    return res_tab

def poly_lagrange(px,py) :
    tab_q=base_lagrange(px)
    rt = len(px)
    res=nppol.polyzero
    for h in range(rt):
        res=nppol.polyadd(res,py[h] * tab_q[h])
    return res

##affichage
x=np.linspace(0, 10, 100)
test_poly = poly_lagrange([0, 2, 4, 6, 8, 10], [f(0), f(2), f(4), f(6), f(8), f(10)])
y=f(x)
y1=nppol.polyval(x,test_poly)
legende1,=plt.plot(x,y)
legende2,=plt.plot(x,y1)
##legende
legende1.set_label("fonction de base")
legende2.set_label("approximation")
plt.legend()
plt.grid(color='black',linestyle='--')
plt.gcf().canvas.set_window_title('Mathys Gagner - Mateo Esteves')

##initialisation
plt.show()