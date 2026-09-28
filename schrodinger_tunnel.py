"""
Équation de Schrödinger 1D : effet tunnel et comparaison de schémas numériques.

Un paquet d'ondes gaussien est envoyé sur une barrière de potentiel plus
haute que son énergie moyenne. On propage la fonction d'onde avec quatre
schémas en temps (Euler explicite, Euler implicite, Crank-Nicolson,
Leapfrog) et une méthode spectrale (diagonalisation de l'hamiltonien),
puis on suit la conservation de la norme et la probabilité de transmission.
Une seconde figure compare la stabilité de Leapfrog et de Crank-Nicolson
pour deux pas de temps (dt = 0.001 et dt = 0.003).

Unités : hbar = m = 1. Bords : conditions de Dirichlet (psi = 0).

Projet M1 Physique, CY Cergy Paris Université (mars 2026).

Utilisation :
    python schrodinger_tunnel.py            # animation, puis figure de stabilité
    python schrodinger_tunnel.py --sauver   # enregistre le GIF et la figure dans figures/
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

SAUVER = "--sauver" in sys.argv

# ---------------------------------------------------------------------------
# Paramètres
# ---------------------------------------------------------------------------

# Domaine spatial : assez large pour que l'onde n'atteigne pas les bords avant t ~ 2
L = 10.0
N = 250
x = np.linspace(-L, L, N)
dx = x[1] - x[0]

# Pas de temps : limité par la stabilité de Leapfrog (dt < 1 / E_max ~ 0.00297,
# un peu moins que la borne dx^2 / 2 de la particule libre à cause de la
# barrière). Avec dt = 0.003, la norme explose.
dt = 0.001

# Barrière de potentiel de hauteur V0 entre x = 1 et x = 2
V0 = 30.0
V = np.zeros(N)
V[(x > 1.0) & (x < 2.0)] = V0

# État initial : paquet d'ondes gaussien qui se dirige vers la barrière
# (énergie centrale k0^2 / 2 = 24.5 < V0, mais une partie du paquet a une
# énergie supérieure à V0 à cause de sa dispersion en impulsion)
x0 = -4.0
sigma = 0.8
k0 = 7.0
psi_initial = np.exp(-0.5 * ((x - x0) / sigma) ** 2) * np.exp(1j * k0 * x)

# Normalisation : l'intégrale de |psi|^2 vaut 1
norme = np.sqrt(np.sum(np.abs(psi_initial) ** 2) * dx)
psi_initial /= norme

# ---------------------------------------------------------------------------
# Hamiltonien discrétisé (différences finies, matrice tridiagonale)
# H = -1/2 d²/dx² + V
# ---------------------------------------------------------------------------

H = np.zeros((N, N), dtype=complex)
for i in range(N):
    H[i, i] = 1.0 / dx ** 2 + V[i]
    if i > 0:
        H[i, i - 1] = -0.5 / dx ** 2
    if i < N - 1:
        H[i, i + 1] = -0.5 / dx ** 2

I = np.eye(N, dtype=complex)

# Pré-calcul des matrices (les inverser dans la boucle serait trop lent)
M_gauche = I + 1j * (dt / 2) * H
M_droite = I - 1j * (dt / 2) * H
inv_imp = np.linalg.inv(I + 1j * dt * H)      # Euler implicite
inv_M_gauche = np.linalg.inv(M_gauche)        # Crank-Nicolson

# Méthode spectrale (référence) : psi(t) = somme c_n exp(-i E_n t) phi_n
energies, modes_propres = np.linalg.eigh(H)
coefficients = np.dot(np.conj(modes_propres.T), psi_initial)

# États initiaux des différents schémas
psi_exp = np.copy(psi_initial)
psi_imp = np.copy(psi_initial)
psi_cn = np.copy(psi_initial)

# Leapfrog a besoin de deux pas initiaux : on fait le premier avec Crank-Nicolson
psi_lf_prev = np.copy(psi_initial)
psi_lf_curr = np.dot(inv_M_gauche, np.dot(M_droite, psi_initial))

temps_liste = [0.0]
prob_exp, prob_imp, prob_cn, prob_lf = [1.0], [1.0], [1.0], [1.0]

# ---------------------------------------------------------------------------
# Graphiques
# ---------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Gauche : densité de probabilité et barrière
ax1.set_xlim(-L, L)
ax1.set_ylim(0, 1.6)
ax1.set_title("Observation de l'effet tunnel")
ax1.set_ylabel(r"Densité de probabilité $|\psi|^2$")
ax1.plot(x, V / V0, color="gray", alpha=0.3, label="Barrière")

line_cn, = ax1.plot(x, np.abs(psi_cn) ** 2, color="blue", lw=2, alpha=0.8, label="Crank-Nicolson")
line_lf, = ax1.plot(x, np.abs(psi_lf_curr) ** 2, color="cyan", ls="-.", alpha=0.8, label="Leapfrog")
line_sp, = ax1.plot(x, np.abs(psi_initial) ** 2, color="magenta", ls=":", alpha=0.8, label="Spectrale (référence)")
ax1.legend(loc="upper left")

# Droite : conservation de la norme
ax2.set_xlim(0, 1.0)
ax2.set_ylim(0, 3.0)
ax2.set_title("Conservation de la probabilité")
ax2.set_xlabel("Temps t")

line_p_exp, = ax2.plot(temps_liste, prob_exp, color="red", label="Euler explicite (diverge)")
line_p_imp, = ax2.plot(temps_liste, prob_imp, color="orange", label="Euler implicite (dissipe)")
line_p_cn, = ax2.plot(temps_liste, prob_cn, color="green", label="Crank-Nicolson (unitaire)")
line_p_lf, = ax2.plot(temps_liste, prob_lf, color="cyan", label="Leapfrog", ls="--")
ax2.legend(loc="upper left")

t_global = 0.0


def update(frame):
    global psi_exp, psi_imp, psi_cn, psi_lf_prev, psi_lf_curr, t_global

    # 25 pas de calcul par image pour accélérer l'animation
    for _ in range(25):
        t_global += dt

        # Euler explicite (on arrête le calcul quand la norme dépasse 10)
        if prob_exp[-1] < 10.0:
            psi_exp = np.dot((I - 1j * dt * H), psi_exp)

        # Euler implicite
        psi_imp = np.dot(inv_imp, psi_imp)

        # Crank-Nicolson
        psi_cn = np.dot(inv_M_gauche, np.dot(M_droite, psi_cn))

        # Leapfrog : psi(t+dt) = psi(t-dt) - 2i dt H psi(t)
        if prob_lf[-1] < 10.0:
            psi_lf_next = psi_lf_prev - 2j * dt * np.dot(H, psi_lf_curr)
            psi_lf_prev = np.copy(psi_lf_curr)
            psi_lf_curr = np.copy(psi_lf_next)

        # Normes
        temps_liste.append(t_global)
        prob_exp.append(np.sum(np.abs(psi_exp) ** 2) * dx if prob_exp[-1] < 10.0 else prob_exp[-1])
        prob_imp.append(np.sum(np.abs(psi_imp) ** 2) * dx)
        prob_cn.append(np.sum(np.abs(psi_cn) ** 2) * dx)
        prob_lf.append(np.sum(np.abs(psi_lf_curr) ** 2) * dx if prob_lf[-1] < 10.0 else prob_lf[-1])

    # Méthode spectrale : solution exacte du problème discrétisé
    psi_sp = np.dot(modes_propres, coefficients * np.exp(-1j * energies * t_global))

    line_cn.set_ydata(np.abs(psi_cn) ** 2)
    line_lf.set_ydata(np.abs(psi_lf_curr) ** 2)
    line_sp.set_ydata(np.abs(psi_sp) ** 2)

    line_p_exp.set_data(temps_liste, prob_exp)
    line_p_imp.set_data(temps_liste, prob_imp)
    line_p_cn.set_data(temps_liste, prob_cn)
    line_p_lf.set_data(temps_liste, prob_lf)

    # Défilement de l'axe du temps
    if t_global > ax2.get_xlim()[1]:
        ax2.set_xlim(0, t_global + 0.5)

    # Probabilité de présence au-delà de la barrière (x > 2)
    transmission = np.sum(np.abs(psi_cn[x > 2.0]) ** 2) * dx
    ax1.set_xlabel(f"Position x   |   t = {t_global:.2f}   |   au-delà de la barrière : {transmission * 100:.2f} %")

    return line_cn, line_lf, line_sp, line_p_exp, line_p_imp, line_p_cn, line_p_lf


plt.tight_layout()

if SAUVER:
    os.makedirs("figures", exist_ok=True)
    # 400 images x 25 pas : jusqu'à t = 10. On voit la séparation du paquet
    # (t ~ 1), puis les réflexions sur les bords et les interférences.
    ani = animation.FuncAnimation(fig, update, frames=400, blit=False)
    ani.save("figures/effet_tunnel.gif", writer=animation.PillowWriter(fps=20), dpi=55)
    print("Animation enregistrée : figures/effet_tunnel.gif")
    # Transmission à t = 1.6, avant que les ondes n'atteignent les bords
    # (calculée avec la méthode spectrale, exacte pour H discrétisé)
    psi_16 = np.dot(modes_propres, coefficients * np.exp(-1j * energies * 1.6))
    print(f"Transmission à t = 1.6 : {np.sum(np.abs(psi_16[x > 2.0]) ** 2) * dx * 100:.1f} %")
    print(f"t = {t_global:.2f} : norme Crank-Nicolson = {prob_cn[-1]:.6f}, "
          f"Euler implicite = {prob_imp[-1]:.3f}")
else:
    ani = animation.FuncAnimation(fig, update, frames=200, interval=30, blit=False)
    plt.show()


# ---------------------------------------------------------------------------
# Stabilité : Leapfrog et Crank-Nicolson pour deux pas de temps
# ---------------------------------------------------------------------------

def normes(schema, pas, t_max=1.5):
    """Norme de psi au cours du temps pour un schéma et un pas de temps donnés."""
    n_pas = int(round(t_max / pas))
    temps = np.arange(n_pas + 1) * pas
    Mg = I + 1j * (pas / 2) * H
    Md = I - 1j * (pas / 2) * H
    cn = np.linalg.solve(Mg, Md)  # matrice d'un pas de Crank-Nicolson
    psi_prev = psi_initial.copy()
    psi = cn @ psi_initial        # premier pas de Leapfrog fait avec Crank-Nicolson
    resultat = [1.0, np.sum(np.abs(psi) ** 2) * dx]
    for _ in range(n_pas - 1):
        if schema == "Crank-Nicolson":
            psi = cn @ psi
        else:
            psi_prev, psi = psi, psi_prev - 2j * pas * (H @ psi)
        resultat.append(np.sum(np.abs(psi) ** 2) * dx)
        if resultat[-1] > 1e3:  # on arrête une fois l'explosion bien visible
            break
    return temps[:len(resultat)], np.array(resultat)


E_max = np.max(energies)
fig2, ax = plt.subplots(figsize=(9, 5))
for schema, pas, style in [("Crank-Nicolson", 0.003, "g-"),
                           ("Leapfrog", 0.001, "c--"),
                           ("Leapfrog", 0.003, "r-")]:
    t, n = normes(schema, pas)
    ax.plot(t, n, style, lw=2, label=f"{schema}, dt = {pas}  (dt·E_max = {pas * E_max:.2f})")
ax.set_yscale("log")
ax.set_ylim(0.5, 1e3)
ax.set_xlabel("Temps t")
ax.set_ylabel(r"Norme $\int |\psi|^2 dx$ (échelle log)")
ax.set_title(f"Stabilité : Leapfrog exige dt < 1/E_max ≈ {1 / E_max:.5f}, Crank-Nicolson est toujours stable")
ax.grid(True, which="both", alpha=0.3)
ax.legend(loc="upper left")
fig2.tight_layout()

if SAUVER:
    fig2.savefig("figures/stabilite.png", dpi=80)
    print("Figure enregistrée : figures/stabilite.png")
else:
    plt.show()
