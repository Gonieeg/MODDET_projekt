import numpy as np
import matplotlib.pyplot as plt
from funkcje import *

# do dokończenia bo źle symuluje :(

class Problem2:
  """Klasa do symulacji przestrzeni 3 pomieszczeń w Problemie 2."""
  def symuluj3(self, T, ht, n, temp_outside, strategiaP, strategia1, strategia2):

    # każdy pokój ma osobną macierz A do dyskretyzacji, z uwzględnionymi jego własnymi warunkami brzegowymi
    # główny pokój (środkowy) odpowiada za implementację warunków brzegowych na ścianach działowych

    # dyskretyzacja przestrzeni
    # -- pokój środkowy
    od, do = 0, 4 # pasożyt P
    Nx, Ny = 20 * n + 1, 20 * n + 1
    N = Nx * Ny

    x,y = np.linspace(od, do, Nx), np.linspace(od, do, Ny)
    X, Y = np.meshgrid(x, y)
    # sąsiad 1 (zachodni)
    od1, do1 = -4, 0
    x1 = np.linspace(od1, do1, Nx)
    X1, Y1 = np.meshgrid(x1, y)
    # sąsiad 2 (wschodni)
    od2, do2 = 4, 8
    x2 = np.linspace(od2, do2, Nx)
    X2, Y2 = np.meshgrid(x2, y)

    # macierz do wykresów na końcu
    Nx_total = Nx * 3 - 2  # bo łączymy bez powtarzania krawędzi
    X_total = np.zeros((Ny, Nx_total))
    Y_total = np.zeros((Ny, Nx_total))
    X_total[:, 0:Nx] = X1
    X_total[:, (Nx - 1):(2*Nx - 1)] = X
    X_total[:, (2 * Nx - 2):] = X2
    Y_total[:, 0:Nx_total] = np.tile(y.reshape(-1, 1), (1, Nx_total))

    hx, hy = x[1] - x[0], y[1] - y[0]
    Xf, Yf = X.flatten(), Y.flatten()

    # indeksy ścian, okna, grzejnika
    #   - indeksy okna pokoju: wycentrowane, zajmuje pół ściany
    ind_okno = np.where((Yf == y[-1]) & (Xf >= (x[0] + 0.25*(x[-1]-x[0]))) & (Xf <= (x[0] + 0.75*(x[-1]-x[0]))))[0]

    #   - indeksy ścian pokojów - bez okna
    ind_sciany_N = np.where((Yf == y[-1]) & ~((Xf >= (x[0] + 0.25*(x[-1]-x[0]))) & (Xf <= (x[0] + 0.75*(x[-1]-x[0])))))[0] # bez okna
    ind_sciany_S = np.where((Yf == y[0]))[0]
    ind_sciany_E = np.where((Xf == x[-1]))[0]
    ind_sciany_W = np.where((Xf == x[0]))[0]
    # komplet
    ind_scian = np.concatenate([ind_sciany_N, ind_sciany_E, ind_sciany_S, ind_sciany_W])

    # indeksy grzejników w pokoju
    ind_grzejnik = np.where((Yf == y[-2]) & (Xf >= (x[0] + 0.375*(x[-1]-x[0]))) & (Xf <= (x[0] + 0.875*(x[-1]-x[0]))))[0] # dodac drugi grzejnik???

    # indeksy wnętrza pokoju, do obliczania sredniej temperatury pomieszczenia - bez scian (tam gdzie startowo jest 20C)
    ind_wnetrza = np.setdiff1d(np.arange(Nx*Ny), np.concatenate([ind_scian, ind_okno]))
    #ind_wnetrzaP = np.where(np.isclose(u_0, Kel(20.0)))[0]

    # indeksy przy ścianach, z których będziemy brać temperaturę "zewnętrzną" czyli z wnętrza sąsiedniego pokoju
    ind_wnetrza_E = np.where((Xf == x[-2]))[0]
    ind_wnetrza_W = np.where((Xf == x[1]))[0]

    # macierz ewolucji
    Id, id_Nx, id_Ny = np.eye(N), np.eye(Nx), np.eye(Ny)
    D2x, D2y = D2(Nx), D2(Ny)

    laplacian = np.kron(id_Ny, D2x) / hx**2 + np.kron(D2y, id_Nx) / hy**2
    A = Id - alfa * ht * laplacian


    # Implementacja warunków brzegowych

    # macierze zadające warunki brzegowe - wszystkie pomieszczenia mają takie same (ściany i okna w takich samych miejscach)
    BsW = - (np.kron(id_Ny, D1F(Nx)) / hx) + lambda_wall/lambda_air * Id
    BsS = - (np.kron(D1F(Ny), id_Nx) / hy) + lambda_wall/lambda_air * Id
    BsE = np.kron(id_Ny, D1B(Nx)) / hx + lambda_wall/lambda_air * Id
    BsN = np.kron(D1B(Ny), id_Nx) / hy + lambda_wall/lambda_air * Id
    # okno
    BoknoB = np.kron(D1B(Ny), id_Nx) / hy + lambda_window/lambda_air * Id

    # warunki w macierzy ewolucji
    A[ind_sciany_N, :] = BsN[ind_sciany_N, :]
    A[ind_sciany_E, :] = BsE[ind_sciany_E, :]
    A[ind_sciany_S, :] = BsS[ind_sciany_S, :]
    A[ind_sciany_W, :] = BsW[ind_sciany_W, :]
    A[ind_okno, :] = BoknoB[ind_okno, :]

    #AP, A1, A2 = A.copy(), A.copy(), A.copy()  # !

    # warunek początkowy
    temp_outside = Kel(temp_outside)
    u_0 = np.full_like(X, Kel(20.0)).flatten()
    u_P, u_1, u_2 = u_0.copy(), u_0.copy(), u_0.copy()  # !


    # warunki brzegowe dla wyświetlania początkowego
    u_0[ind_scian] = lambda_wall / lambda_air * temp_outside
    print('pierdylion', lambda_wall / lambda_air * temp_outside)
    u_0[ind_okno] = lambda_window / lambda_air * temp_outside
    print('dwa pierdyliardy', lambda_window / lambda_air * temp_outside)


    # pętla symulacji
    t = int(T/ht)
    zmianysredniej = []

    warunek_zewn = lambda_wall / lambda_air * temp_outside

    for _ in range(t):
      # równanie du/dt
      u_P += ht * f(ind_grzejnik, u_P, ind_wnetrza, strategiaP)
      u_1 += ht * f(ind_grzejnik, u_1, ind_wnetrza, strategia1)
      u_2 += ht * f(ind_grzejnik, u_2, ind_wnetrza, strategia2)

      # warunki brzegowe
      # - ściany zewnętrzne
      # -- południowe i północne i okna
      for u in [u_1, u_P, u_2]:
        u[ind_sciany_N] = warunek_zewn
        u[ind_okno] = warunek_zewn
        u[ind_sciany_S] = warunek_zewn
      # -- zachodnia i wschodnia pokojów 1 i 2
      u_1[ind_sciany_W] = warunek_zewn
      u_2[ind_sciany_E] = warunek_zewn

      # - ściany wewnętrzne
      #lambda_wall / lambda_air * temp_outside
      #jakos wziac temperature z kazdego punkta z wewnatrz pokojow 1,2 i uzyc jej jako zewn u Pasożyta
      # potem ustawic takie same sciany u sasiadow?
      # czy to zadziala czy nie wezmie pod uwage chlodznie z pokoju w srodku? powinno sie robic w dyfuzji i think
      u_P[ind_sciany_W] = lambda_wall / lambda_air * u_1[ind_sciany_E-1]
      u_P[ind_sciany_E] = lambda_wall / lambda_air * u_2[ind_sciany_W+1]
      u_2[ind_sciany_W] = lambda_wall / lambda_air * u_P[ind_sciany_E-1]
      u_1[ind_sciany_E] = lambda_wall / lambda_air * u_P[ind_sciany_W+1]
      # krok symulacji
      u_P = np.linalg.solve(A, u_P)
      zmianysredniej.append(u_P[ind_wnetrza].mean())
      # moja propozycja: dopiero tutaj aktualizujemy warunki brzegowe u sasiadow, zeby u_P mial wplyw na nie ale tez zeby sciany byly podobne?
       #= u_P[ind_sciany_E]#lambda_wall / lambda_air * u_1[ind_sciany_E - 1]
       #= u_P[ind_sciany_W]#lambda_wall / lambda_air * u_2[ind_sciany_W + 1]
      u_1 = np.linalg.solve(A, u_1)
      u_2 = np.linalg.solve(A, u_2)

    #print(f"mean {_+1}: ", u_P[ind_wnetrza].mean())
    for u in [u_1, u_P, u_2]:
      print(u.min(),u.max(), u.mean())
    odch = np.std(u_P[ind_wnetrza])
    print(f"standard deviation:{odch} ")


    u_0 = Cel(u_0)
    u_P = Cel(u_P).reshape(Ny, Nx)
    u_1 = Cel(u_1).reshape(Ny, Nx)
    u_2 = Cel(u_2).reshape(Ny, Nx)

    fig, axs = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    vmin = min(u_1.min(), u_P.min(), u_2.min())
    vmax = max(u_1.max(), u_P.max(), u_2.max())

    levels = np.linspace(vmin, vmax, 30)  # gęstość konturów
    cmap = "jet"

    cs1 = axs[0].contourf(X, Y, u_1, levels=levels, cmap=cmap)
    axs[0].set_title("Pokój 1")
    axs[0].set_xlabel("x")
    axs[0].set_ylabel("y")
    fig.colorbar(cs1, ax=axs[0])

    csP = axs[1].contourf(X, Y, u_P, levels=levels, cmap=cmap)
    axs[1].set_title("Pokój P (środkowy)")
    axs[1].set_xlabel("x")
    fig.colorbar(csP, ax=axs[1])

    cs2 = axs[2].contourf(X, Y, u_2, levels=levels, cmap=cmap)
    axs[2].set_title("Pokój 2")
    axs[2].set_xlabel("x")
    fig.colorbar(cs2, ax=axs[2])


    plt.show()

    return u_P, zmianysredniej, odch


sim1 = Problem2()
a=sim1.symuluj3(2, 1, 1, 0, 0, 0, 0)
print(a[0])

