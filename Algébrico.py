import numpy as np
import  pandas
import math as c
import matplotlib.pyplot as plt

ax = plt.axes(projection='3d')

#número de pontos
pt = 30
#alpha da função, velocidade que o calor se espalha
alpha = 0.01
#tempo
t = 10
#limita de série
series=50

#pontos do período de 0 a pi em x a y
x = np.linspace(0, c.pi, pt)
y = np.linspace(0, c.pi, pt)

#é preciso do mesmo número de pontos psoteriormente, por isso crio 2 ndarray vazias
x2 = np.linspace(0, 0, pt)
y2 = np.linspace(0, 0, pt)

#criação das ndarray
X, Y = np.meshgrid(x, y)
X2, Y2 = np.meshgrid(x2, y2)

#A00 definido por Fourier
A00 = (1/(c.pi**2))*(c.e**c.pi -1)**2

#Série com Am0 Fourier
Xsum = []

for p in range(0, pt):
    Xsoma = 0
    for m in range(1, series):
        Am0 = (2/(c.pi**2))*(((c.e**(c.pi) - 1)*(c.e**(c.pi)*(-1)**(m)-1))/(m**(2)+1))*c.cos(X[0][p]*m)*c.e**(-alpha*m**(2)*t)
        Xsoma = Xsoma + Am0
    Xsum.append(Xsoma)

#Série com A0n Fourier
Ysum = []
for q in range(0, pt):
    Ysoma = 0
    for n in range(1, series):
        A0n = (2/(c.pi**2))*(((c.e**(c.pi) - 1)*(c.e**(c.pi)*(-1)**(n)-1))/(n**(2)+1))*c.cos(Y[q][0]*n)*c.e**(-alpha*n**(2)*t)
        Ysoma = Ysoma + A0n
    Ysum.append(Ysoma)

#Série com Amn Fourier. Aqui que preciso da ndarray vazio,
#para conseguir criar uma tabela com pt x pt e preenchê-la com seus respectivos valores de z
XYsum = X2
for r in range(0, pt):
    for s in range(0,pt):
        XYsoma = 0
        for n in range(1, series):
            for m in range(1, series):
                Amn = (4/(c.pi**2))*(((c.e**(c.pi)*(-1)**(m)-1)*((c.e**(c.pi)*(-1)**(n)-1)))/((n**(2)+1)*(m**(2)+1)))*c.cos(X[0][s]*m)*c.cos(Y[r][0]*n)*c.e**(-alpha*(m**2+n**2)*t)
                XYsoma = XYsoma + Amn
        XYsum[s][r] = XYsoma

#soma de todos os valores
SOMAFINAL = Y2
for t in range (0, pt):
    for u in range (0, pt):
        SOMAFINAL[u][t] = A00 + Xsum[u] + Ysum[t] +XYsum[u][t]

ax.plot_surface(X, Y, SOMAFINAL, cmap='hot')
ax.set_zlim(0, 500)
ax.set_xlim(0, c.pi)
ax.set_ylim(0, c.pi)
plt.show()
