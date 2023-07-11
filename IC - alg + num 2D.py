import numpy as np
import  pandas
import math as c
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#número de pontos
pt = 100
#alpha da função, velocidade que o calor se espalha
alpha = 0.01
#tempo
t = 10
#limita de série
series = 100
#Delta T
Dt = 0.005
#max_iter_time tem que ser verificado para os valores de Dt e T
max_iter_time = int(t/Dt)+1

plate_length = pt

h = c.pi/pt

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
Algébrico = Y2
for tn in range (0, pt):
    for un in range (0, pt):
        Algébrico[un][tn] = A00 + Xsum[un] + Ysum[tn] +XYsum[un][tn]


if (1 - (4*alpha*Dt)/h**2) <= 0:
    print("Valores errados para o método explícito")

else:
    u = np.empty((max_iter_time, plate_length, plate_length))
    #valores da temperatura na condição inicial
    u[0] = c.e**(X+Y)

    def calculate(u):
        #pontos internos
        for l in range (0, max_iter_time-1, 1):
            for i in range (0, pt-1):
                for j in range (0, pt-1):
                    u[l+1][i][j] = u[l][i][j]*(1 - (4*alpha*Dt)/h**2) + alpha * Dt * ((u[l][i+1][j]+u[l][i][j+1]+u[l][i-1][j]+u[l][i][j-1])/h**2)

            #bordas sem cantos
            for k in range (0, pt-1):
                u[l+1][0][k] = u[l][0][k] * (1 - (4 * alpha * Dt) / h ** 2) + alpha * Dt * ((2*u[l][1][k] + u[l][0][k+1] + u[l][0][k - 1]) / h ** 2)
                u[l+1][pt-1][k] = u[l][pt-1][k] * (1 - (4 * alpha * Dt) / h ** 2) + alpha * Dt * ((2 * u[l][pt-2][k] + u[l][pt-1][k + 1] + u[l][pt-1][k - 1]) / h ** 2)
                u[l+1][k][0] = u[l][k][0] * (1 - (4 * alpha * Dt) / h ** 2) + alpha * Dt * ((u[l][k+1][0] + 2* u [l][k][1] + u[l][k-1][0]) / h ** 2)
                u[l+1][k][pt-1] = u[l][k][pt-1] * (1 - (4 * alpha * Dt) / h ** 2) + alpha * Dt * ((u[l][k + 1][pt-1] + 2 * u[l][k][pt-2] + u[l][k - 1][pt-1]) / h ** 2)

            #4 cantos
            u[l+1][0][0]= u[l][0][0]*(1 - (4*alpha*Dt)/h**2) + alpha * Dt * ((2*u[l][1][0]+2*u[l][0][1])/h**2)
            u[l+1][pt-1][pt-1]= u[l][pt-1][pt-1]*(1 - (4*alpha*Dt)/h**2) + alpha * Dt * ((2*u[l][pt-2][pt-1]+2*u[l][pt-1][pt-2])/h**2)
            u[l+1][pt-1][0]= u[l][pt-1][0]*(1 - (4*alpha*Dt)/h**2) + alpha * Dt * ((2*u[l][pt-2][0]+2*u[l][pt-1][1])/h**2)
            u[l+1][0][pt-1]= u[l][0][pt-1]*(1 - (4*alpha*Dt)/h**2) + alpha * Dt * ((2*u[l][0][pt-2]+2*u[l][1][pt-1])/h**2)

        return u

    Numérico = calculate(u)

relativeerror = ((np.absolute(Algébrico - Numérico[int(t/Dt)]))/Numérico[int(t/Dt)])*100

def plotheatmap(relativeerror_k, k):
        # Clear the current plot figure
        plt.clf()
        plt.title(f"Erro relativo para t = {t} unidade de tempo")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.xlim(0, c.pi)
        plt.ylim(0, c.pi)
        # This is to plot u_k (u at time-step k)
        plt.imshow(relativeerror_k, cmap='hot', vmin=0, vmax=1,
                   extent=[0, c.pi, c.pi, 0])  # Definindo os limites dos eixos x e y
        plt.colorbar()
        return plt

def animate(k):
    plotheatmap(relativeerror, k)

anim = animation.FuncAnimation(plt.figure(), animate, interval=10, frames=max_iter_time, repeat=False)
plt.show()
print("Done!")
