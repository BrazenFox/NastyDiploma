import math
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go

C = 0.12
NI = 36.0
CO = 0.0
CU = 2.0
W = 3.1
MN = 1.5
SI = 0.5
ALL = 0.0
CR = 15.0
VA = 0.0
TI = 1.3
MO = 0.0
NB = 0.0
ZR = 0.0

UD = 24.0  # Вольт
J = 150.0  # Ампер

PSI = 0.7

D = 4.0  # толщина металла с см
V = 0.3  # скорость сварки см/с

EE = math.e
PP = math.pi

G = (1 / 4) * ((1 / 3) * (C + (1 / 5) * (NI + CO + CU + W / 3)) + MN / 14 + (SI + ALL) / 7 + (CR + V + TI) / 13 + (
        MO + NB + ZR) / 24)
L = (1 / 20) * (11 - 7 * math.pow(G, 1 / 4))
CG = (1 / 10) * (55 - 9 * math.pow(G, 1 / 4))
A = L / CG

TL1 = 5 * (4 * (77 - 3 * C) - MN) - 0.5 * (13 * math.sqrt(CR) + 25 * math.sqrt(NI))
TL2 = -2 * (4 * SI + MO + VA + CO + ALL + 3 * CU + 7 * TI + 9 * NB) - 13 * ZR
TL = TL1 + TL2
AL = (1 / 37 + TL / 28500) ** 2
#######################################################################
Q = PSI * J * UD
##################################################################
PE = V * D / A
E2 = Q / (L * D * TL)
E1 = E2 / PE
E3 = E2 * PE
E4 = E3 * PE
#########################################################################
C = 0.000001


def rasschet_T_by_X_Y_Z(X, Y, Z):
    A3 = Q / (2 * PP * L)
    B = V / (2 * A)
    R = math.sqrt(X ** 2 + Y ** 2 + Z ** 2)
    T = A3 * math.exp(-B * (R + X)) / R
    return T


def rasschet_Tmax(Z, T):
    return rasschet_Xm_Ym_by_Z_Tmax(Z, T)[3]


def rasschet_Xm_Ym_by_Z_Tmax(Z, T):
    A3 = Q / (2 * PP * L)
    B = V / (2 * A)
    D = A3 * B
    PO = D / T
    P = PO
    P3 = PO * math.exp(-P / (1 + P))
    R = P3 - P
    while R ** 2 > C:
        P = P3
        P3 = PO * math.exp(-P / (1 + P))
        R = P3 - P

    P1 = P3 * math.sqrt(1 + 2 * P3) / (1 + P3)
    YM = P1 / B
    XM = -(P3 ** 2) / (B * (1 + P3))
    TM = -XM / V
    R = math.sqrt(XM ** 2 + YM ** 2 + Z ** 2)
    T = A3 * math.exp(-B * (R + XM)) / R
    return XM, YM, TM, T



fig, ax = plt.subplots()
Z=0
count_of_graph = 5
T = rasschet_Tmax(Z, TL)
XM = [0] * count_of_graph
YM = [0] * count_of_graph
TT = [((i+1) * T/count_of_graph) for i in range(count_of_graph)]
x = np.arange(-5, 1, 0.1)
t_array = [[0] * len(x) for i in range(count_of_graph)]
for t in range(count_of_graph):
    xm, ym, tm, tt = rasschet_Xm_Ym_by_Z_Tmax(Z, TT[t])
    XM[t] = xm
    YM[t] = ym
    print(YM[t])
    for i in range(len(x)):
        t_array[t][i] = rasschet_T_by_X_Y_Z(x[i], YM[t], Z)
    ax.plot(x, t_array[t], label = u'Y '+str(t))

"""XM_plot = [0] * 100
YM_plot = [0] * 100
TT_plot = [((i+1) * T/100) for i in range(100)]
for t in range(100):
    xm, ym, tm, tt = rasschet_Xm_Ym_by_Z_Tmax(Z, TT_plot[t])
    XM_plot[t] = xm
    YM_plot[t] = ym
ax.plot(XM_plot, TT_plot)"""


plt.title(u'График распределения температур по оси X для различных значений Y по схеме ДТ-ПТ')
plt.xlabel(u'X, см')
plt.ylabel(u'Температура, T\N{DEGREE SIGN}C')
ax.plot(XM, TT, 'r--',label = u'Максимумы')
plt.legend()
plt.grid(True)
plt.show()



fig, ax = plt.subplots()

for i in range(count_of_graph):
    ax.plot(-x/V, t_array[i], label = u'Y '+str(t))
XM = np.array(XM)
ax.plot(-XM/V, TT, 'r--',label = u'Максимумы')

plt.title(u'График термического цикла по схеме ДТ-ПТ')
plt.xlabel(u't,c')
plt.ylabel(u'Температура, T\N{DEGREE SIGN}C')
plt.legend()
plt.grid(True)
plt.show()


x = np.arange(-10,2,0.01)
y = np.arange(-3,3,0.01)
X, Y = np.meshgrid(x, y)
t_array = [[0] * len(x) for i in range(len(y))]
for i in range(len(y)):
    for j in range(len(x)):
        t_array[i][j] = rasschet_T_by_X_Y_Z(x[j], y[i], 0.5)

fig, ax = plt.subplots()

ax.contour(X,Y,t_array, levels = 5)

max_index=[0]*len(t_array)
for i in range(len(t_array)):
    max_index[i] = 0
    for j in range(len(t_array[i])):
        if t_array[i][j]>t_array[i][max_index[i]]:
            max_index[i] = j
print(len(t_array), len(t_array[0]))
print(X[0][0], Y[0][0])
Xmax=[0]*len(max_index)
Ymax=[0]*len(max_index)
for i in range(len(max_index)):
    Xmax[i] = X[i][max_index[i]]
    Ymax[i] = Y[i][max_index[i]]
ax.plot(Xmax, Ymax, 'r--',label = u'Максимумы')
plt.xlabel(u'X,см')
plt.ylabel(u'Y,см')
plt.title("Температурное поле предельного состояния по схеме ДТ-ПТ")
plt.legend()
plt.grid(True)

plt.show()

x_count = 200
y_count = 100
z_count = 100

X, Y, Z = np.mgrid[-20:2:200j, -5:5:100j, -5:5:100j]
values=[[[0 for i in range(z_count)] for j in range(y_count)] for k in range(x_count)]
for i in range(x_count):
    for j in range(y_count):
        for k in range(z_count):
            values[i][j][k]=rasschet_T_by_X_Y_Z(X[i][j][k], Y[i][j][k], Z[i][j][k])

values = np.array(values)
fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=values.flatten(),
    isomin=0,
    isomax=1500,
    opacity=0.1, # needs to be small to see through all surfaces
    surface_count=21, # needs to be a large number for good volume rendering
    ))
fig.show()




