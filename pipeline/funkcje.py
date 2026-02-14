import numpy as np

# współczynnik dyfuzji
alfa = 0.000019

# ciśnienie Pa
p = 101300

# ind. stała gazowa powietrza
r = 287.05

# ciepło właściwe
c = 1005

# współczynnik przenikania ciepła przez materiał (podzielony przez jego grubość w m)
# breeze blocks - pustaki
lambda_wall = 0.1#/0.3 m
# szyba
lambda_window = 0.96#/0.01 m
# powietrze
lambda_air = 0.0262

# moc grzejnika
P = 1267


def D2(N):
  """Macierz dyskretyzująca drugie pochodne"""
  D2 = -2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)
  # cos wyzerowac na rogach?
  return D2


def D1B(N):
  """Macierz dyskretyzująca pierwsze pochodne 'backward' (do kierunków północ N, wschód E)"""
  D1 = -np.eye(N, k=-1) + np.eye(N)
  return D1


def D1F(N): # MOZE BYC ZLE
  """Macierz dyskretyzująca pierwsze pochodne 'forward' (kierunki południowy S, zachodni W)"""
  D1 = np.eye(N, k=+1) - np.eye(N)
  return D1


def Kel(Cel):
  """Przekształca stopnie Celsjusza na stopnie Kelvina"""
  return Cel + 273.15

def Cel(Kel):
  """Przekształca stopnie Kelvina na stopnie Celsjusza"""
  return Kel - 273.15

# ustawienia grzałki
S = {0: 7, 1: 16, 2: 20, 3: 24, 4: 28} # w *C
for k in S: # na kelwiny
    S[k] = Kel(S[k])


def f(grzejnik, u, wnetrze, ust_grzalki):  # 1. gdzie dodawac cieplo 2. macierz pokoju 3. przy jakiej sredniej koniec grzania
    """Funkcja implementująca działanie grzejnika"""
    Area = 4*4*2
    Si = S[ust_grzalki]

    u_n = np.zeros_like(u)
    nu = P * r / (p * Area * c)

    if np.mean(u[wnetrze]) <= Si:
        # tylko tam gdzie grzejnik robic to dzialanie
        u_n[grzejnik] = u[grzejnik] * nu
    return u_n

