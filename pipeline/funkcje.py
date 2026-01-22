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

def f(room, u, ust_grzalki, tc): # albo room zamiast x,y? indeksy pokoju
    A = powierzchnia_pokoju(room)
    Si = S[ust_grzalki]
    if  np.mean(u[room]) <= Si:
        return u * P * r / (ro * A * c)


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