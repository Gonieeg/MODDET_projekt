import numpy as np
alfa = 0.000019
p = 101300
r = 287.05
c = 1005
lambda_wall = 0.17 #bez podzielenia przez grubosc okno spada do 2 a nie 0
lambda_window = 0.96 # bez / spada mniej
lambda_air = 0.0262
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

# ustawienia grzałki ZINENIC
S = {0: 7, 1: 16, 2: 20, 3: 24, 4: 28} # w *C
for k in S: # na kelwiny
    S[k] = Kel(S[k])


def f(grzejnik, u, wnetrze, ust_grzalki):  # 1. gdzie dodawac cieplo 2. macierz pokoju 3. przy jakiej sredniej koniec grzania
    """Funkcja implementująca działanie grzejnika"""
    Area = 4*4  # zrobic funkcje do tego w kolejnych problemach
    Si = S[ust_grzalki]

    u_n = np.zeros_like(u)
    nu = P * r / (p * Area * c)

    if np.mean(u[wnetrze]) <= Si:  # ind_wnetrza dodac do funkcji?
        # tylko tam gdzie grzejnik robic to dzialanie
        u_n[grzejnik] = u[grzejnik] * nu
    return u_n

