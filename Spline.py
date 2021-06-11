
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.linalg import solve


def runge_func(x):
    return 1/(1+pow(x,2))

def L(k,xi,x):
    prod = 1
    for j in range(len(xi)):
        if j == k:
            continue
        prod*= (x-xi[j])/(xi[k]-xi[j])
    return prod

def pLagr(xi,fi,x):#wertet zu Datensatz (xi,fi) das Interpolationspolynom p(x) aus
    sum = 0
    for j in range(len(xi)):
        sum += fi[j] * L(j, xi, x)
    return sum

def getVector_l(xi):
    l_vector = [0]*(len(xi)-2)
    for i in range(len(xi)-2):
        l_vector[i] = xi[i+2]-xi[i+1]
    return l_vector
def getVector_r(xi):
    r_vector = [0]*(len(xi)-2)
    for i in range(len(xi)-2):
        r_vector[i] = xi[i+1]-xi[i]
    return r_vector
def getVector_d(xi):
    d_vector = [0]*(len(xi)-2)
    for i in range(len(xi)-2):
        d_vector[i] = 2*((xi[i+2]-xi[i+1])+(xi[i+1]-xi[i]))
    return d_vector
def getVector_R(xi,fi,K1,K2):#s1''(a) = K1 und s2''(b) = K2, die Randbedingungen!
    n = len(xi)
    R_vector = [0]*n
    R_vector[0] = 3*((fi[1]-fi[0])/(xi[1]-xi[0]))-K1*(xi[1]-xi[0]) #erste Element von K1 abhängig
    R_vector[n-1] = 3*((fi[n-1]-fi[n-2])/(xi[n-1]-xi[n-2]))-K2*(xi[n-1]-xi[n-2]) #letzte Element von K2 abhängig
    for i in range(1,n-1):
        R_vector[i] = 3*((fi[i+1]-fi[i])*(xi[i]-xi[i-1])/(xi[i+1]-xi[i])+(fi[i]-fi[i-1])*(xi[i+1]-xi[i])/(xi[i]-xi[i-1]))
    return R_vector


def solve_S(xi,fi,K1,K2):
    n = len(xi)
    M = np.zeros((n, n))
    #natürlichen kubischen Spline
    M[0][0] = 2
    M[0][1] = 1
    M[n-1][n-1] = 2
    M[n-1][n-2] = 1
    l = getVector_l(xi)
    d = getVector_d(xi)
    r = getVector_r(xi)
    for i in range(1,n-1):#Zeilen
        M[i][i-1] = l[i-1]
        M[i][i-1+1] = d[i-1]
        M[i][i-1+2] = r[i-1]
    R = getVector_R(xi,fi,K1,K2)
    S = solve(M, R)
    return S

def s_i(xi,fi,Si,Si_plus_1,untere_Grenze_i,obere_Grenze_i,x): #Pf,3 beschränkt auf das Intervall [xi,xi+1] wertet x aus
    c1_i = fi[untere_Grenze_i]
    c2_i = Si
    c3_i = (3*fi[obere_Grenze_i]-3*fi[untere_Grenze_i]-2*Si*(xi[obere_Grenze_i]-xi[untere_Grenze_i])-Si_plus_1*(xi[obere_Grenze_i]-xi[untere_Grenze_i]))/pow(xi[obere_Grenze_i]-xi[untere_Grenze_i],2)
    c4_i = (2 * fi[untere_Grenze_i] - 2 * fi[obere_Grenze_i] + Si*(xi[obere_Grenze_i]-xi[untere_Grenze_i])+Si_plus_1*(xi[obere_Grenze_i]-xi[untere_Grenze_i])) / pow(xi[obere_Grenze_i] - xi[untere_Grenze_i], 3)
    y = c1_i + c2_i*(x-xi[untere_Grenze_i])+c3_i*pow(x-xi[untere_Grenze_i],2)+c4_i*pow(x-xi[untere_Grenze_i],3)
    return y
def aquid_x(a,b,n):# Intervall [a,b] mit n Stützstellen
    xi = [0]*n
    for i in range(n):
        xi[i] = -5 + (10*((i+1)-1))/(n-1)
    return xi
def spline(xi,fi,S,x):
    #In welchem Intervall liegt x?
    #Fall 1: Linke Grenze
    #Fall 2: Rechte Grenze
    #Fall 3: ansonsten von inneren Knoten umrandet
    n = len(xi)
    obere_Grenze_i = 0
    untere_Grenze_i = 0
    Si = 0 #von S, welches mit solve_S gelöst wurde, sind nur Si und Si+1 relevant für die Koef.
    Si_plus_1 = 0
    for i in range(1,n):
        if x <= xi[i]:
            obere_Grenze_i = i
            untere_Grenze_i = i-1
            Si_plus_1 = S[i]
            Si = S[i-1]
            break
    #Grenzen gefunden, jetzt soll s_i ausgewertet werden,
    #wobei Si, Si+1 benötigt werden für die Koef.
    y = s_i(xi,fi,Si,Si_plus_1,untere_Grenze_i,obere_Grenze_i,x)
    return y


#Spline-Funktion
n = 15 #Anzahl der Stützstellen größer 1
xi = aquid_x(-5,5,n) #Erstellen von N äquidistanten Stützstellen
fi = [] #Stützstellenwerte berechnet
x = np.linspace(-5,5, 100) #all x values to be evaluated
for i in xi:
    fi.append(runge_func(i))
S = solve_S(xi,fi,0,0)
spline_y = [0]*len(x)

for j in range(len(x)):
    spline_y[j] = spline(xi,fi,S,x[j])
plt.plot(x,spline_y, label="Spline Interpolation", linewidth=2)



#Runge-Funktion
plt.plot(x,runge_func(x), label="Runge Funktion", linewidth=2)


plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.show()

