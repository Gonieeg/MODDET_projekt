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

# t jako argument?
# funkcja grzejnika

def room_area(ind_room):
    ## A=b*c tylko jak to wziac po indeksach pokoju
    return (x[-1]-x[0]) * (y[-1]-y[0])


def f(room, u, ust_grzalki, tc): # albo room zamiast x,y? indeksy pokoju
    A = room_area(room)
    Si = S[ust_grzalki]
    if  np.mean(u[room]) <= Si:
        return u * P * r / (ro * A * c)


def przesuniecie_grzejnika(r=0.1): # d - odleglosc od okna in [0.1, 0.2, ..., 3.9]
    # przy zalozeniu ze hx, hy sa co 0.1?
    # metr dlugosci
    # powinien byc 0.1m width
    # 3.9, 3.8, ...., 0.2, 0.1
    ind_heater = np.where((Yf == 4 - r) & (Xf >= 1.5) & (Xf <= 2.5))[0]
    return ind_heater


# ściany zewnętrzne
def ho(wall, u): # wall - indeksy ściany której temperature liczymy
    # u - temperatura na ścianie
    return -lambda_wall * 0.25 / lambda_air * (u[wall] - temp_outside) # razy grubość ściany


# tylko co jesli to jest sciana dzialowa???
def hi(wall, u):
    return -lambda_wall * 0.12 / lambda_air * (u[wall] - !!!!!)  # wziąć indeksy po drugiej stronie ściany ??

# okna
def g(window, u):
    return -lambda_window * 0.005 / lambda_air * (u[window] - temp_outside)  # wziąć indeksy po drugiej stronie ściany ??



def start_condition(x, y):
    ... # jakies 20 stopni wszedzie
