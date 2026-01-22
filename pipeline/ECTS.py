import numpy as np
import matplotlib.pyplot as plt

#### Próba #1 dla przetestowania i zrozumienia struktury projektu - klasa estymująca utworzona z poprzedniej listy

class ETCS:
    def start(self, T, n, od, do, alfa=1):
        self.T = T
        self.n = n
        self.diff_coeff = alfa

        self.od = od
        self.do = do
        self.Nx = 20 * n + 1
        self.Ny = 20 * n + 1
        self.N = self.Nx * self.Ny

    @staticmethod
    def D2(N):
        return -2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)

    @staticmethod
    def D1B(N):
        return -np.eye(N, k=-1) + np.eye(N)

    @staticmethod
    def f(x, y):
        """Warunek początkowy"""
        return np.exp(- (x + 1) ** 2 - (y + 1) ** 2)

    def symuluj_dla(self, ht):
        """Symulacja dla danego ht, t obrotów pętli i zapisanie do self.wyniki"""

        self.ht = ht

        self.t = int(self.T/ht)

        self.x = np.linspace(self.od, self.do, self.Nx)
        self.y = np.linspace(self.od, self.do, self.Ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)

        hx = self.x[1] - self.x[0]
        hy = self.y[1] - self.y[0]
        Xf = self.X.flatten()
        Yf = self.Y.flatten()
        self.Xf = Xf
        self.Yf = Yf

        id_Nx, id_Ny = np.eye(self.Nx), np.eye(self.Ny)
        D2x, D2y = self.D2(self.Nx), self.D2(self.Ny)
        laplacian = np.kron(id_Ny, D2x) / hx ** 2 + np.kron(D2y, id_Nx) / hy ** 2

        factor = self.diff_coeff * ht
        A = np.eye(self.Nx * self.Ny) - factor * laplacian

        # indeksy brzegowe Neumanna
        ind_brzeg_x_right = np.where((Xf == self.x[-1]) & (Yf <= 0.5))[0]
        ind_brzeg_y_top = np.where((Yf == self.y[-1]) & (Xf <= 0))[0]
        ind_brzeg_x_pion = np.where((Xf == 0) & (Yf >= 0.5))[0]

        # indeksy brzegowe Dirichleta
        ind_brzeg_x_left = np.where(Xf == self.x[0])[0]
        ind_brzeg_y_bot = np.where(Yf == self.y[0])[0]
        ind_brzeg_y_poziom = np.where((Yf == 0.5) & (Xf > 0))[0]
        ind_dirichlet = np.concatenate([ind_brzeg_x_left, ind_brzeg_y_bot, ind_brzeg_y_poziom])
        self.ind_dirichlet = ind_dirichlet

        # Neumann
        Cx_right = np.kron(id_Ny, self.D1B(self.Nx) / hx ** 2)
        Cy_top = np.kron(self.D1B(self.Ny) / hy ** 2, id_Nx)
        A[ind_brzeg_x_right, :] = Cx_right[ind_brzeg_x_right, :]
        A[ind_brzeg_y_top, :] = Cy_top[ind_brzeg_y_top, :]
        A[ind_brzeg_x_pion, :] = Cx_right[ind_brzeg_x_pion, :]

        # Dirichleta
        N = self.Nx * self.Ny
        Id = np.eye(N)
        A[ind_dirichlet, :] = Id[ind_dirichlet, :]

        # ograniczamy macierz do obszaru L
        indeksyL = np.where(
            ((Xf <= 0) & (Yf >= -1) & (Yf <= 1)) |
            ((Xf >= 0) & (Yf >= -1) & (Yf <= 0.5))
        )[0]
        self.indeksyL = indeksyL
        self.A_L = A[np.ix_(indeksyL, indeksyL)]

        self.u_0 = self.f(self.X, self.Y).flatten()

        ### SYMULACJA

        self.wyniki = {}
        epsilony = [0, 0.01, 0.1]
        mask_eps = (np.isclose(self.Xf[self.indeksyL], 0.1)) & (np.isclose(self.Yf[self.indeksyL], 0.1))
        u_0 = self.u_0

        for e in epsilony:
            u_current = u_0.copy()
            u_current[self.ind_dirichlet] = 0
            u_current = u_current[self.indeksyL]

            for it in range(self.t):
                u_current[mask_eps] += e * self.ht
                u_current = np.linalg.solve(self.A_L, u_current)

            Ufull = np.full(self.N, np.nan)
            Ufull[self.indeksyL] = u_current
            self.wyniki[e] = Ufull.reshape(self.Ny, self.Nx)

    def wykresy(self):
        """Tworzy wykresy na podstawie wyników symulacji"""

        fig, axs = plt.subplots(2, 2, figsize=(12, 10))

        # --- wykres początkowy ---

        Ufull = np.full(self.N, np.nan)
        Ufull[self.indeksyL] = self.u_0[self.indeksyL]
        self.wyniki["przed"] = Ufull.reshape(self.Ny, self.Nx)

        # dane do wykresów
        dane = [self.wyniki["przed"], self.wyniki[0], self.wyniki[0.01], self.wyniki[0.1]]
        tytuly = ["Obszar L przed symulacją",
                  "Obszar L po symulacji dla ε = 0",
                  "Obszar L po symulacji dla ε = 0.01",
                  "Obszar L po symulacji dla ε = 0.1"]

        # wykresy
        for ax, Z, title in zip(axs.flat, dane, tytuly):
            im = ax.imshow(Z, origin='lower', extent=[self.x.min(), self.x.max(), self.y.min(), self.y.max()],
                           cmap='viridis', interpolation='nearest', aspect='equal')
            ax.set_title(title)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            fig.colorbar(im, ax=ax)

        plt.tight_layout()
        plt.show()


sim = ETCS()
sim.start(T=10, n=1, od=-1, do=1)
sim.symuluj_dla(ht=1)
sim.wykresy()
