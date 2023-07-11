import numpy as np
import  pandas
import math as c
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

#tamanha da placa, alpha, número de pontos, Dt, etc
max_iter_time = 1001
pt = 100
plate_length = pt
alpha = 0.01
Dt = 0.01
h = c.pi/pt

if (1 - (4*alpha*Dt)/h**2) <= 0:
    print("Valores errados para o método explícito")

else:
    #pontos do período de 0 a pi em x a y
    x = np.linspace(0, c.pi, pt)
    y = np.linspace(0, c.pi, pt)
    X, Y = np.meshgrid(x, y)

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

    u = calculate(u)




