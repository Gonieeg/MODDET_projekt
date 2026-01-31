import numpy as np
import matplotlib.pyplot as plt
from funkcje import * #Kel, Cel, D2, D1B, D1F
# import all - zmienne tez powinny, ale spr czy cos

######
alfa=0.000019
p=1013
r=287.05
c=1005
lambda_wall =  0.17 / 0.25#0.3/0.25
lambda_window = 0.96 / 0.005 #1.1/0.005  # nie moge mnozyc przez grubosc bo wtedy wieksze i bardziej mrozace jest wall
lambda_air = 0.0262
P=1267

#S = [7, 16, 20, 24, 28]
temp_outside = Kel(0)
######


ht = 0.1
T = 10#000 # [s]
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
#alfa = 0.019
factor = alfa * ht #alfa realistycznie jest mala i prawie sie nie zmienia powietrze (ale przez to rysunek sie praktycznie nie zmienia!!)

A = Id - factor * laplacian


# indeksy okna pokoju 1 - od 1 do 3 na top ścianie
ind_okno = np.where((Yf == y[-1]) & (Xf >= 1) & (Xf <= 3))[0]

#ideksy ścian pokoju 1 - bez okna
#ind_wall = np.where((Xf == od) | (Xf == do) | (Yf == od) | (Yf == do) & ~((Xf >= 1) & (Xf <= 3)))[0]
ind_sciany_N = np.where((Yf == do) & ~((Xf >= 1) & (Xf <= 3)))[0] # bez okna
ind_sciany_S = np.where((Yf == od))[0]
ind_sciany_E = np.where((Xf == do))[0]
ind_sciany_W = np.where((Xf == od))[0]
# komplet
ind_scian = np.concatenate([ind_sciany_N, ind_sciany_E, ind_sciany_S, ind_sciany_W])


# warunki brzegowe

# macierze zadające warunki brzegowe (każda ściana ma swój)(okno ma swój)
### kolejnosc kron?
BsW = np.kron(id_Ny, D1F(Nx)) / hx + lambda_wall/lambda_air * Id
BsS = np.kron(D1F(Ny), id_Nx) / hy + lambda_wall/lambda_air * Id
BsE = np.kron(id_Ny, D1B(Nx)) / hx + lambda_wall/lambda_air * Id
BsN = np.kron(D1B(Ny), id_Nx) / hy + lambda_wall/lambda_air * Id
BoknoB = np.kron(D1B(Ny), id_Nx) / hy + lambda_window/lambda_air * Id #okno

A[ind_sciany_N, :] = BsN[ind_sciany_N, :]
A[ind_sciany_E, :] = BsE[ind_sciany_E, :]
A[ind_sciany_S, :] = BsS[ind_sciany_S, :]
A[ind_sciany_W, :] = BsW[ind_sciany_W, :]
A[ind_okno, :] = BoknoB[ind_okno, :]



# przed petla liczaca
# grzejnik
#odleglosc = 0.2
# odleglosc nie dziala poniewaz dokladnosc floatow. stawiamy grzjnik w j-tym rzedzie a r polizcymt pozniej
j = Ny - 2
ind_grzejnik = np.where((Yf == y[j]) & (Xf >= 1.5) & (Xf <= 2.5))[0]
# odleglosci grzejnika od okna w metrach
rs = [0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 3.9] #PARZYSTE Z TAKIM HX HY

# ustawienia grzałki ZINENIC
S = {0: 7, 1: 16, 2: 20, 3: 24, 4: 28} # w *C
for k in S: # na kelwiny
    S[k] = Kel(S[k])

# na razie nie jest wprowadzone
def f(grzejnik, u, ust_grzalki):  # albo room zamiast x,y? indeksy pokoju
  """Funkcja implementująca działanie grzejnika"""
  Area = 4*4 # zrobic funkcje do tego w kolejnych problemach
  Si = S[ust_grzalki]

  u_n = np.zeros_like(u)
  nu = P * r / (p * Area * c)

  if np.mean(u) <= Si:
    # tylko tam gdzie grzejnik robic to dzialanie
    u_n[grzejnik] = u[grzejnik] * nu
  return u_n

#print("temp termostatu graniczna:", S[4])

u_0 = np.full_like(X, Kel(20.0)).flatten()
u_current = u_0.copy()

u_0 += ht * f(ind_grzejnik, u_0, 4) # zmodyfikowac f zeby bralo koordy grzejnika albo
  #u_current[ind_grzejnik] += ht * f(x, u_current)
#print(np.mean(u_0) <= S[4])

#print(f(ind_grzejnik, u_0, 4))#.mean())
#print(u_0)#.mean())


  # warunki brzegowe
u_0[ind_scian] = lambda_wall / lambda_air * temp_outside
u_0[ind_okno] = lambda_window / lambda_air * temp_outside

#print(lambda_wall / lambda_air * temp_outside)
#print(lambda_window / lambda_air * temp_outside)

print("mean 0: ", u_current.mean())
for _ in range(t):
  # równanie du/dt
  u_current += ht * f(ind_grzejnik, u_current, 4) # zmodyfikowac f zeby bralo koordy grzejnika albo
  #u_current[ind_grzejnik] += ht * f(x, u_current)

  # warunki brzegowe
  u_current[ind_scian] = lambda_wall / lambda_air * temp_outside
  u_current[ind_okno] = lambda_window / lambda_air * temp_outside # to cos sprawia ze okno grzeje? robi skrajnosc

  # krok symulacji
  u_current = np.linalg.solve(A, u_current)
  print(f"mean {_+1}: ", u_current.mean())


#u_0.mean()
#u_current.max()
u_0 = Cel(u_0)
u_current = Cel(u_current)

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

levels0 = np.linspace(u_0.min(), u_0.max(), 50)
'''
if Cel(temp_outside) > u_current.min():
  print(u_current < Cel(temp_outside))
  '''
z_min = u_current.min()
z_max = u_current.max()
levels = np.linspace(z_min, z_max, 50)

im1 = axs[0].contourf(X, Y, u_0.reshape(Nx, Ny), levels=levels0, cmap='viridis')
axs[0].set_title('Warunek początkowy')
fig.colorbar(im1, ax=axs[0], ticks=np.linspace(u_0.min(), u_0.max(), 6))

im2 = axs[1].contourf(X, Y, u_current.reshape(Nx, Ny), levels=levels, cmap='viridis')
axs[1].set_title('Wynik końcowy')
fig.colorbar(im2, ax=axs[1], ticks=np.linspace(z_min, z_max, 6))

plt.tight_layout()
plt.show()

# Czemu okno tak mało chłodzi - mimo że ma 1000x skrajniejsze wartości?
# Czemu grzejnik grzeje tylko tam gdzie jest i to do milionów?
# - do milionów bo czas symulacji i wyłącza się gdy jest osiagnięta średnia - a rosnie tylko w tym miejscu
# - rosnie tylko na grzejniku bo współczynnik dyfuzji jest malutki??