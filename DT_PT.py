import math
import numpy as np


class DT_PT:

    def __init__(self, UD, J, V, Cap=0.000001, C=0.12, NI=36.0, CO=0.0, CU=2.0, W=3.1, MN=1.5, SI=0.5, ALL=0.0, CR=15.0,
                 VA=0.0, TI=1.3, MO=0.0, NB=0.0, ZR=0.0, PSI=0.7, D=4.0):
        self.C = C
        self.NI = NI
        self.CO = CO
        self.CU = CU
        self.W = W
        self.MN = MN
        self.SI = SI
        self.ALL = ALL
        self.CR = CR
        self.VA = VA
        self.TI = TI
        self.MO = MO
        self.NB = NB
        self.ZR = ZR
        self.UD = UD
        self.J = J
        self.V = V
        self.PSI = PSI
        self.D = D
        self.Cap = Cap
        self.G = (1 / 4) * (
                (1 / 3) * (C + (1 / 5) * (NI + CO + CU + W / 3)) + MN / 14 + (SI + ALL) / 7 + (CR + V + TI) / 13 + (
                MO + NB + ZR) / 24)
        self.L = (1 / 20) * (11 - 7 * math.pow(self.G, 1 / 4))
        self.CG = (1 / 10) * (55 - 9 * math.pow(self.G, 1 / 4))
        self.A = self.L / self.CG
        self.TL1 = 5 * (4 * (77 - 3 * C) - MN) - 0.5 * (13 * math.sqrt(CR) + 25 * math.sqrt(NI))
        self.TL2 = -2 * (4 * SI + MO + VA + CO + ALL + 3 * CU + 7 * TI + 9 * NB) - 13 * ZR
        self.TL = self.TL1 + self.TL2
        self.AL = (1 / 37 + self.TL / 28500) ** 2
        self.Q = PSI * J * UD
        self.PE = V * D / self.A
        self.E2 = self.Q / (self.L * D * self.TL)
        self.E1 = self.E2 / self.PE
        self.E3 = self.E2 * self.PE
        self.E4 = self.E3 * self.PE

    def rasschet_T_by_X_Y_Z(self, X, Y, Z):
        A3 = self.Q / (2 * math.pi * self.L)
        B = self.V / (2 * self.A)
        R = math.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        T = A3 * math.exp(-B * (R + X)) / R
        return T

    def rasschet_Tmax(self, Z, T):
        return self.rasschet_Xm_Ym_by_Z_Tmax(Z, T)[3]

    def rasschet_Xm_Ym_by_Z_Tmax(self, Z, T):
        A3 = self.Q / (2 * math.pi * self.L)
        B = self.V / (2 * self.A)
        D = A3 * B
        PO = D / T
        P = PO
        P3 = PO * math.exp(-P / (1 + P))
        R = P3 - P
        while R ** 2 > self.Cap:
            P = P3
            P3 = PO * math.exp(-P / (1 + P))
            R = P3 - P

        P1 = P3 * math.sqrt(1 + 2 * P3) / (1 + P3)
        YM = P1 / B
        XM = -(P3 ** 2) / (B * (1 + P3))
        TM = -XM / self.V
        R = math.sqrt(XM ** 2 + YM ** 2 + Z ** 2)
        T = A3 * math.exp(-B * (R + XM)) / R
        return XM, YM, TM, T

    def first_graph(self, Z, count_of_graph):
        T = self.rasschet_Tmax(Z, self.TL)
        XM = [0] * count_of_graph
        YM = [0] * count_of_graph
        TT = [((i + 1) * T / count_of_graph) for i in range(count_of_graph)]
        x = np.arange(-5, 1, 0.01)
        t_array = [[0] * len(x) for i in range(count_of_graph)]
        for t in range(count_of_graph):
            xm, ym, tm, tt = self.rasschet_Xm_Ym_by_Z_Tmax(Z, TT[t])
            XM[t] = xm
            YM[t] = ym
            for i in range(len(x)):
                t_array[t][i] = self.rasschet_T_by_X_Y_Z(x[i], YM[t], Z)
        return x, t_array, XM, TT

    def second_graph(self, Z, count_of_graph):
        x, t_array, XM, TT = self.first_graph(Z, count_of_graph)
        x = -x/self.V
        XM = -np.array(XM)/self.V
        return x, t_array, XM, TT

    def contour(self, x = np.arange(-10, 2, 0.01), y= np.arange(-3, 3, 0.01)):
        X, Y = np.meshgrid(x, y)
        t_array = [[0] * len(x) for i in range(len(y))]
        for i in range(len(y)):
            for j in range(len(x)):
                t_array[i][j] = self.rasschet_T_by_X_Y_Z(x[j], y[i], 0.5)

        max_index = [0] * len(t_array)
        for i in range(len(t_array)):
            max_index[i] = 0
            for j in range(len(t_array[i])):
                if t_array[i][j] > t_array[i][max_index[i]]:
                    max_index[i] = j
        Xmax = [0] * len(max_index)
        Ymax = [0] * len(max_index)
        for i in range(len(max_index)):
            Xmax[i] = X[i][max_index[i]]
            Ymax[i] = Y[i][max_index[i]]
        return X, Y, t_array, Xmax, Ymax

    def DDD(self, x_count=200, y_count=100, z_count=100):
        X, Y, Z = np.mgrid[-20:2:x_count+0j, -5:5:y_count+0j, -5:5:z_count+0j]
        values = [[[0 for i in range(z_count)] for j in range(y_count)] for k in range(x_count)]
        for i in range(x_count):
            for j in range(y_count):
                for k in range(z_count):
                    values[i][j][k] = self.rasschet_T_by_X_Y_Z(X[i][j][k], Y[i][j][k], Z[i][j][k])
        values = np.array(values)
        return X,Y,Z,values