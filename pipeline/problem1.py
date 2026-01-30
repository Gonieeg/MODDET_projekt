import numpy as np
import funkcje

# macierz dyskretyzująca drugie pochodne
def D2(N):
  D2 = -2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)
  return D2

# macierz dyskretyzująca pierwsze pochodne backward (prawo/góra)
def D1B(N):
  D1 = -np.eye(N, k=-1) + np.eye(N)
  return D1

# macierz dyskretyzująca pierwsze pochodne forward (lewo/dół)
def D1F(N): # MOZE BYC ZLE
  D1 = np.eye(N, k=+1) - np.eye(N)
  return D1


######
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
######


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


Id = np.eye(Nx*Ny)
id_Nx, id_Ny = np.eye(Nx), np.eye(Ny)

D2x, D2y = D2(Nx), D2(Ny)



laplacian = np.kron(id_Ny, D2x) / hx**2 + np.kron(D2y, id_Nx) / hy**2 # sprawdzic w poprzednihc czy taka kolejnosc arguemntow kron

factor = 1 * ht #alfa tu

A = Id - factor * laplacian
#b = u^n + ht * f(x, u^n)
#b = u_current + ht * f(x, u_current) # f - farelka


# indeksy okna pokoju 1 - od 1 do 3 na top ścianie
ind_okno = np.where((Yf == y[-1]) & (Xf >= 1) & (Xf <= 3))[0]

#ideksy ścian pokoju 1 - bez okna
#ind_wall = np.where((Xf == od) | (Xf == do) | (Yf == od) | (Yf == do) & ~((Xf >= 1) & (Xf <= 3)))[0]
ind_sciany_N = np.where((Xf == do) & ~((Xf >= 1) & (Xf <= 3)))[0] # bez okna
ind_sciany_S = np.where((Yf == od))[0]
ind_sciany_E = np.where((Xf == do))[0]
ind_sciany_W = np.where((Xf == od))[0]
# komplet
ind_scian = np.concatenate([ind_sciany_N, ind_sciany_E, ind_sciany_S, ind_sciany_W])


# warunki brzegowe

# macierze zadające warunki brzegowe (każda ściana ma swój)(okno ma swój)
### kolejnosc kron?
BsW = np.kron(id_Ny, D1F(Nx)) + lambda_wall/lambda_air * Id
BsS = np.kron(D1F(Ny), id_Nx) + lambda_wall/lambda_air * Id
BsE = np.kron(id_Ny, D1B(Nx)) + lambda_wall/lambda_air * Id
BsN = np.kron(D1B(Ny), id_Nx) + lambda_wall/lambda_air * Id
BoknoB = np.kron(D1B(Ny), id_Nx) + lambda_window/lambda_air * Id #okno

A[ind_sciany_N, :] = BsN[ind_sciany_N, :]
A[ind_sciany_E, :] = BsN[ind_sciany_E, :]
A[ind_sciany_S, :] = BsN[ind_sciany_S, :]
A[ind_sciany_W, :] = BsN[ind_sciany_W, :]
A[ind_okno, :] = BoknoB[ind_okno, :]


#? jak i gdzie:
#A u^{n+1} = b^n # "na wierszach roznych od ind_brzegowe"

# przed petla liczaca
# grzejnik 0.1
ind_grzejnik = np.where((Yf == 4 - r) & (Xf >= 1.5) & (Xf <= 2.5))[0]
# odleglosci grzejnika od okna w metrach
rs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9]

# ustawienia grzałki ZINENIC I TO NA KELWINY CHYBA STARTOWA TEMPERATURE TEZ!!
S = {0: 7, 1: 16, 2: 20, 3: 24, 4: 28}


def f(grzejnik, u, ust_grzalki=0, tc=0):  # albo room zamiast x,y? indeksy pokoju
  A = 4*4
  #Si = S[ust_grzalki]

  u_n = np.zeros_like(u)
  nu = P * r / (ro * A * c)

  if np.mean(u) <= 20: #Si:
    # tylko tam gdzie grzejnik robic to dzialanie
    u_n[grzejnik] = u[grzejnik] * nu
  return u_n


u_0 = np.full_like(X, 20.0).flatten()
u_current = u_0.copy()

for _ in range(t):
  # równanie du/dt
  u_current += ht * f(ind_grzejnik, u_current) # zmodyfikowac f zeby bralo koordy grzejnika albo
  #u_current[ind_grzejnik] += ht * f(x, u_current)
  # chyba bardziej przyszlosciowo jest w funkcji? dac tam u[grzejnik]


  # warunki brzegowe
  u_current[ind_scian] = lambda_wall / lambda_air * temp_outside
  u_current[ind_okno] = lambda_wall / lambda_air * temp_outside

  # krok symulacji
  u_current = np.linalg.solve(A, u_current)


