import matplotlib.pyplot as plt
import plotly.graph_objects as go
from DT_PT import DT_PT

print("Введите напряжение в Вольтах, например: 10.3)")
UD = float(input())
print("Введите ток в Амперах, например 140.0")
J = float(input())
print("Введите скорость сварки в см/с, например 0.3")
V = float(input())
dt_pt_main = DT_PT(UD, J, V)
print("Теплопроводность", dt_pt_main.L)
print("Объемная теплоемкость", dt_pt_main.CG)
print("Температуропроводность", dt_pt_main.A)
print("Расчетное значение температуры ликвидус", dt_pt_main.TL)
print("Коэффициен полной поверхности теплоотдачи", dt_pt_main.AL)
print("Расчетное значение эффективной тепловой мощности", dt_pt_main.Q)
#####################################################################################
fig, ax = plt.subplots()
Z=0
count_of_graph = 5
x, t_array, XM, TT = dt_pt_main.first_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i))
ax.plot(XM, TT, 'r--',label = u'Максимумы')
plt.title(u'График распределения температур по оси X для различных значений Y по схеме ДТ-ПТ при силе тока I='+str(J)+'А напряжении U=' + str(UD) + 'В и скорости v=' + str(V) + 'см/с')
plt.xlabel(u'X, см')
plt.ylabel(u'Температура, T\N{DEGREE SIGN}C')
plt.legend()
plt.grid(True)
plt.show()
#####################################################################################
fig, ax = plt.subplots()
Z=0
count_of_graph = 5
dt_pt = DT_PT(10.3, 140.0, 0.3)
x, t_array, XM, TT = dt_pt.first_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i), color='red')
ax.plot(XM, TT, 'r--', label=u'Максимумы', color='red')
dt_pt = DT_PT(10.4, 150.0, 0.4)
x, t_array, XM, TT = dt_pt.first_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i), color='green')
ax.plot(XM, TT, 'r--', label=u'Максимумы', color='green')
dt_pt = DT_PT(10.5, 160.0, 0.5)
x, t_array, XM, TT = dt_pt.first_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i), color='blue')
ax.plot(XM, TT, 'r--',label = u'Максимумы', color='blue')
plt.title(u'График распределения температур по оси X для различных значений Y по схеме ДТ-ПТ при различных силах тока силе тока (140, 150, 160) напряжениях U (10.3, 10.4, 10.5) и скоростях v (0.3, 0.4, 0.5) соответственно')
plt.xlabel(u'X, см')
plt.ylabel(u'Температура, T\N{DEGREE SIGN}C')
plt.legend()
plt.grid(True)
plt.show()
#####################################################################################
fig, ax = plt.subplots()
Z=0
count_of_graph = 5
x, t_array, XM, TT = dt_pt_main.second_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i))
ax.plot(XM, TT, 'r--',label = u'Максимумы')
plt.title(u'График термического цикла по схеме ДТ-ПТ при силе тока I='+str(J)+'А напряжении U=' + str(UD) + 'В и скорости v=' + str(V) + 'см/с')
plt.xlabel(u't,c')
plt.ylabel(u'Температура, T\N{DEGREE SIGN}C')
plt.legend()
plt.grid(True)
plt.show()
#####################################################################################
fig, ax = plt.subplots()
Z=0
count_of_graph = 5
dt_pt = DT_PT(10.3, 140.0, 0.3)
x, t_array, XM, TT = dt_pt.second_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i), color='red')
ax.plot(XM, TT, 'r--',label = u'Максимумы',color='red')

dt_pt = DT_PT(10.4, 150.0, 0.4)
x, t_array, XM, TT = dt_pt.second_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i), color='green')
ax.plot(XM, TT, 'r--',label = u'Максимумы', color='green')

dt_pt = DT_PT(10.5, 160.0, 0.5)
x, t_array, XM, TT = dt_pt.second_graph(Z, count_of_graph)
for i in range(count_of_graph):
    ax.plot(x, t_array[i], label=u'Y ' + str(i), color='blue')
ax.plot(XM, TT, 'r--',label = u'Максимумы',color='blue')

plt.title(u'График термического цикла по схеме ДТ-ПТ при различных силах тока силе тока (140, 150, 160) напряжениях U (10.3, 10.4, 10.5) и скоростях v (0.3, 0.4, 0.5) соответственно')
plt.xlabel(u't,c')
plt.ylabel(u'Температура, T\N{DEGREE SIGN}C')
plt.legend()
plt.grid(True)
plt.show()
#####################################################################################
fig, ax = plt.subplots()
X, Y, t_array, Xmax, Ymax = dt_pt_main.contour()
ax.contour(X,Y,t_array, levels = 5)
ax.plot(Xmax, Ymax, 'r--',label = u'Максимумы')
plt.xlabel(u'X,см')
plt.ylabel(u'Y,см')
plt.title('Температурное поле предельного состояния по схеме ДТ-ПТ при силе тока I='+str(J)+'А напряжении U=' + str(UD) + 'В и скорости v=' + str(V) + 'см/с')
plt.legend()
plt.grid(True)
plt.show()
#####################################################################################
fig, ax = plt.subplots()
dt_pt = DT_PT(10.3, 140.0, 0.3)
X, Y, t_array, Xmax, Ymax = dt_pt.contour()
ax.contour(X,Y,t_array, levels = 5, colors = ['red'])
ax.plot(Xmax, Ymax, 'r--',label = u'Максимумы', color = 'red')

dt_pt = DT_PT(10.4, 150.0, 0.4)
X, Y, t_array, Xmax, Ymax = dt_pt.contour()
ax.contour(X,Y,t_array, levels = 5, colors = ['green'])
ax.plot(Xmax, Ymax, 'r--',label = u'Максимумы', color = 'green')

dt_pt = DT_PT(10.5, 160.0, 0.5)
X, Y, t_array, Xmax, Ymax = dt_pt.contour()
ax.contour(X,Y,t_array, levels = 5, colors = ['blue'])
ax.plot(Xmax, Ymax, 'r--',label = u'Максимумы', color = 'blue')
plt.xlabel(u'X,см')
plt.ylabel(u'Y,см')
plt.title('Температурное поле предельного состояния по схеме ДТ-ПТ при различных силах тока силе тока (140, 150, 160) напряжениях U (10.3, 10.4, 10.5) и скоростях v (0.3, 0.4, 0.5) соответственно')
plt.legend()
plt.grid(True)
plt.show()
#####################################################################################
X,Y,Z,values = dt_pt_main.DDD()
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
fig.update_layout(
    title='Распределение температур движущегося точечного источника на поверхности полубесконечного тела\n\r'
          'при силе тока I='+str(J)+'А напряжении U=' + str(UD) + 'В и скорости v=' + str(V) + 'см/с',
    )
fig.show()