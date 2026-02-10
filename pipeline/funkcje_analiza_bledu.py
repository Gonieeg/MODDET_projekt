import numpy as np
import matplotlib.pyplot as plt
from funkcje import * #Kel, Cel, D2, D1B, D1F
#import tqdm

def test2kroku(T, kroki, n, temp_outside=0, j=2, strategia=0):

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


    # implementacja warunków brzegowych

    # macierze zadające warunki brzegowe
    BsW = - (np.kron(id_Ny, D1F(Nx)) / hx) + lambda_wall/lambda_air * Id
    BsS = - (np.kron(D1F(Ny), id_Nx) / hy) + lambda_wall/lambda_air * Id
    BsE = np.kron(id_Ny, D1B(Nx)) / hx + lambda_wall/lambda_air * Id
    BsN = np.kron(D1B(Ny), id_Nx) / hy + lambda_wall/lambda_air * Id
    #okno
    BoknoB = np.kron(D1B(Ny), id_Nx) / hy + lambda_window/lambda_air * Id

    # warunek początkowy
    temp_outside = Kel(temp_outside)
    u_0 = np.full_like(X, Kel(20.0)).flatten()

    Dwall = lambda_wall / lambda_air * temp_outside
    Dwindow = lambda_window / lambda_air * temp_outside

    bledy =[]
    for ht in kroki: #tqdm.tqdm(kroki):
        u_1 = u_0.copy()
        u_2 = u_0.copy()

        A1 = Id - alfa * ht * laplacian
        A2 = Id - alfa * ht/2 * laplacian

        for A in [A1,A2]:
        # warunki w macierzy ewolucji
            A[ind_sciany_N, :] = BsN[ind_sciany_N, :]
            A[ind_sciany_E, :] = BsE[ind_sciany_E, :]
            A[ind_sciany_S, :] = BsS[ind_sciany_S, :]
            A[ind_sciany_W, :] = BsW[ind_sciany_W, :]
            A[ind_okno, :] = BoknoB[ind_okno, :]

        t = int(T / ht)

        blad_ht=0

        for _ in range(t):
            # równanie du/dt
            u_1 += ht * f(ind_grzejnik, u_1, ind_wnetrza, strategia)

            # warunki brzegowe
            u_1[ind_scian] = Dwall
            u_1[ind_okno] = Dwindow

            # krok symulacji
            u_1 = np.linalg.solve(A1, u_1)

            for i in range(2): # dla ht/2 dwa obroty pętli
                # równanie du/dt
                u_2 += ht/2 * f(ind_grzejnik, u_2, ind_wnetrza, strategia)

                # warunki brzegowe
                u_2[ind_scian] = Dwall
                u_2[ind_okno] = Dwindow

                # krok symulacji
                u_2 = np.linalg.solve(A2, u_2)

            blad_ht += np.mean(np.abs(u_1 - u_2))

        bledy.append(blad_ht / t)


    plt.plot(nasze_kroki, np.array(bledy))

    plt.title("Wykres błędu podwojonego kroku")
    plt.xlabel("h_t")
    plt.ylabel("log(MAE)")  # lub MAE lub RMS

    plt.yscale('log')
    plt.show()
    return bledy

nasze_kroki = np.arange(0.1,10.1,0.1)
test2kroku(20,nasze_kroki,1,0,2,0)