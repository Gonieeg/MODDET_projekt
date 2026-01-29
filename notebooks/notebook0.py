import numpy as np
alfa=0.000019
ro=1013
r=287.05
c=1005
lambda_wall=0.17
lambda_window=0.96
lambda_air=0.0262
P=1267

S = [7, 16, 20, 24, 28]
temp_outside = 10

ht = 1
T = 1
t = int(T/ht)

n = 1
od = 0
do = 4
Nx = 20 * n + 1
Ny = 20 * n + 1
N = Nx * Ny

x = np.linspace(od, do, Nx)
y = np.linspace(od, do, Ny)
X, Y = np.meshgrid(x, y)

hx = x[1] - x[0]
hy = y[1] - y[0]
Xf = X.flatten()
Yf = Y.flatten()

# indeksy okna pokoju 1 - od 1 do 3 na top ścianie
ind_window = np.where((Yf == y[-1]) & (Xf >= 1) & (Xf <= 3))[0]
#ideksy ścian pokoju 1 - bez okna
ind_wall = np.where((Xf == od) | (Xf == do) | (Yf == od) | (Yf == do) & ~((Xf >= 1) & (Xf <= 3)))[0]


rs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9]

