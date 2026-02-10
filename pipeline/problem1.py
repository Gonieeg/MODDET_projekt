import numpy as np
import matplotlib.pyplot as plt
from funkcje import * #Kel, Cel, D2, D1B, D1F
# import all - zmienne tez powinny, ale spr czy cos

class Problem1:

  #def __init__(self):
  # importuj z dane.csv


  def symuluj(self, T, ht, n, temp_outside, j, strategia):

    # dyskretyzacja przestrzeni
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

    # indeksy ścian, okna, grzejnika
    #   - indeksy okna pokoju 1 - od 1 do 3 na północnej ścianie
    ind_okno = np.where((Yf == y[-1]) & (Xf >= 1) & (Xf <= 3))[0]

    #   - indeksy ścian pokoju - bez okna
    ind_sciany_N = np.where((Yf == y[-1]) & ~((Xf >= 1) & (Xf <= 3)))[0] # bez okna
    ind_sciany_S = np.where((Yf == y[0]))[0]
    ind_sciany_E = np.where((Xf == x[-1]))[0]
    ind_sciany_W = np.where((Xf == x[0]))[0]
    # komplet
    ind_scian = np.concatenate([ind_sciany_N, ind_sciany_E, ind_sciany_S, ind_sciany_W])

    # indeksy grzejnika
    ind_grzejnik = np.where((Yf == y[j]) & (Xf >= 1.5) & (Xf <= 2.5))[0]

    # indeksy wnętrza pokoju, do obliczania sredniej temperatury pomieszczenia - bez scian (tam gdzie startowo jest 20C)
    ind_wnetrza = np.setdiff1d(np.arange(Nx*Ny), np.concatenate([ind_scian, ind_okno]))


    # macierz ewolucji
    Id = np.eye(N)
    id_Nx, id_Ny = np.eye(Nx), np.eye(Ny)
    D2x, D2y = D2(Nx), D2(Ny)

    laplacian = np.kron(id_Ny, D2x) / hx**2 + np.kron(D2y, id_Nx) / hy**2
    factor = alfa * ht
    A = Id - factor * laplacian


    # implementacja warunków brzegowych

    # macierze zadające warunki brzegowe
    BsW = - (np.kron(id_Ny, D1F(Nx)) / hx) + lambda_wall/lambda_air * Id
    BsS = - (np.kron(D1F(Ny), id_Nx) / hy) + lambda_wall/lambda_air * Id
    BsE = np.kron(id_Ny, D1B(Nx)) / hx + lambda_wall/lambda_air * Id
    BsN = np.kron(D1B(Ny), id_Nx) / hy + lambda_wall/lambda_air * Id
    #okno
    BoknoB = np.kron(D1B(Ny), id_Nx) / hy + lambda_window/lambda_air * Id

    # warunki w macierzy ewolucji
    A[ind_sciany_N, :] = BsN[ind_sciany_N, :]
    A[ind_sciany_E, :] = BsE[ind_sciany_E, :]
    A[ind_sciany_S, :] = BsS[ind_sciany_S, :]
    A[ind_sciany_W, :] = BsW[ind_sciany_W, :]
    A[ind_okno, :] = BoknoB[ind_okno, :]


    # warunek początkowy
    temp_outside = Kel(temp_outside)
    u_0 = np.full_like(X, Kel(20.0)).flatten()
    u_current = u_0.copy()


    # warunki brzegowe dla wyświetlania początkowego
    u_0[ind_scian] = lambda_wall / lambda_air * temp_outside
    u_0[ind_okno] = lambda_window / lambda_air * temp_outside


    # pętla symulacji
    t = int(T/ht)
    zmianysredniej = []
    #print("mean 0: ", u_current[ind_wnetrza].mean())
    for _ in range(t):
      # równanie du/dt
      u_current += ht * f(ind_grzejnik, u_current, ind_wnetrza, strategia)

      # warunki brzegowe
      u_current[ind_scian] = lambda_wall / lambda_air * temp_outside
      u_current[ind_okno] = lambda_window / lambda_air * temp_outside # zimno z okna praktycznie nie robi dyfuzji do wewnątrz

      # krok symulacji
      u_current = np.linalg.solve(A, u_current)
      zmianysredniej.append(u_current[ind_wnetrza].mean())

    print(f"mean {_+1}: ", u_current[ind_wnetrza].mean())

    odch = np.std(u_current[ind_wnetrza])
    print(f"standard deviation:{odch} ")


    u_0 = Cel(u_0)
    u_current = Cel(u_current).reshape(Ny, Nx)

    # wykresy
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    levels0 = np.linspace(0, 30, 50)
    levels = np.linspace(u_current.min(), u_current.max(), 50)

    im1 = axs[0, 0].contourf(X, Y, u_0.reshape(Ny, Nx), levels=levels0, cmap='viridis')
    axs[0, 0].set_title('Warunek początkowy')
    fig.colorbar(im1, ax=axs[0, 0])

    im2 = axs[0, 1].contourf(X, Y, u_current, levels=levels, cmap='viridis')
    axs[0, 1].set_title(f'Wynik po {t} krokach')
    fig.colorbar(im2, ax=axs[0, 1])

    im3 = axs[1, 0].contourf(
      X, Y, u_current,
      levels=20, vmin=0, vmax=40, cmap='jet')
    axs[1, 0].set_title(f'Wynik po {t} krokach')
    fig.colorbar(im3, ax=axs[1, 0])


    zmianysredniej = [Cel(x) for x in zmianysredniej]
    axs[1, 1].plot(range(t), zmianysredniej)
    axs[1, 1].set_title('Średnia temperatura w czasie')
    axs[1, 1].set_xlabel('krok czasu')
    axs[1, 1].set_ylabel('Temperatura [°C]')

    plt.tight_layout()
    plt.show()

    return u_current, zmianysredniej, odch


sim1 = Problem1()
a=sim1.symuluj(1000, 1, 2, 0, 2, 0)



#### dodac parametr temperatury wewnetrznej?
#### lepsze stawianie grzejnika?
# finalnie mean i std ze scianami wlacznie powinno byc ale tbh bez scian (i moze grzejnika? bo on zawyza std tym ze bardziej sie nagrzewa zdala od okna) mowi wiecej
