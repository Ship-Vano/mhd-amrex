#!/usr/bin/env python3
"""
plot_verification.py — графики верификационных расчётов (CSV автономного
драйвера mhd2d_verify ИЛИ срезы из plotfile'ов AMReX, сконвертированные в CSV).

Использование:
    python3 scripts/plot_verification.py briowu  out_briowu.csv
    python3 scripts/plot_verification.py ot      out_ot.csv
    python3 scripts/plot_verification.py rotor   out_rotor.csv
    python3 scripts/plot_verification.py alfven  out_alfven_16.csv out_alfven_32.csv ...
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(fname):
    d = np.genfromtxt(fname, delimiter=",", names=True)
    return d


def plot_briowu(files):
    d = load(files[0])
    # срез вдоль y: берём первый ряд ячеек
    ys = np.unique(d["y"]); m = d["y"] == ys[0]
    x = d["x"][m]
    fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    panels = [("rho", r"$\rho$"), ("u", r"$u$"), ("By", r"$B_y$"), ("p", r"$p$")]
    for a, (k, lbl) in zip(ax.flat, panels):
        a.plot(x, d[k][m], "k.-", ms=2, lw=0.7)
        a.set_ylabel(lbl); a.grid(alpha=0.3)
    for a in ax[1]: a.set_xlabel("x")
    fig.suptitle("Brio–Wu, t = 0.1 (HLLD + CT, MUSCL-MC, RK2, CFL=0.4)")
    fig.tight_layout()
    fig.savefig("briowu.png", dpi=150)
    print("-> briowu.png")


def plot_ot(files):
    d = load(files[0])
    n = int(round(np.sqrt(len(d["x"]))))
    X = d["x"].reshape(n, n); Y = d["y"].reshape(n, n)
    R = d["rho"].reshape(n, n); P = d["p"].reshape(n, n)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for a, F, t in [(ax[0], R, r"$\rho$"), (ax[1], P, r"$p$")]:
        im = a.pcolormesh(X, Y, F, cmap="jet", shading="auto")
        a.set_aspect(1); a.set_title(t + ", t = 0.5")
        plt.colorbar(im, ax=a, shrink=0.85)
    fig.tight_layout()
    fig.savefig("ot_maps.png", dpi=150)
    print("-> ot_maps.png")

    # давление вдоль y = 0.3125 — сравнение с рис. 8 статьи Авдеевой–Лукина
    ys = np.unique(d["y"])
    j = np.argmin(np.abs(ys - 0.3125))
    m = d["y"] == ys[j]
    fig, a = plt.subplots(figsize=(8, 4))
    a.plot(d["x"][m], d["p"][m], "k.-", ms=3, lw=0.8)
    a.set_xlabel("x"); a.set_ylabel("p")
    a.set_title(f"Orszag–Tang: давление вдоль y = {ys[j]:.4f}, t = 0.5")
    a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("ot_pressure_slice.png", dpi=150)
    print("-> ot_pressure_slice.png")

    print(f"max |divB| в выходном файле: {np.max(np.abs(d['divB'])):.3e}")


def plot_rotor(files):
    d = load(files[0])
    n = int(round(np.sqrt(len(d["x"]))))
    X = d["x"].reshape(n, n); Y = d["y"].reshape(n, n)
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for a, k, t in [(ax[0], "rho", r"$\rho$"), (ax[1], "p", r"$p$")]:
        im = a.pcolormesh(X, Y, d[k].reshape(n, n), cmap="jet", shading="auto")
        a.set_aspect(1); a.set_title(t + ", t = 0.15")
        plt.colorbar(im, ax=a, shrink=0.85)
    fig.tight_layout(); fig.savefig("rotor_maps.png", dpi=150)
    print("-> rotor_maps.png")


def plot_alfven(files):
    """Профили B⊥ и таблица сходимости по серии сеток."""
    alpha = np.pi / 6
    ca, sa = np.cos(alpha), np.sin(alpha)
    fig, a = plt.subplots(figsize=(9, 4.5))
    errs = []
    for f in files:
        d = load(f)
        N = int(f.rstrip(".csv").split("_")[-1])
        bperp = d["By"] * ca - d["Bx"] * sa
        xi = d["x"] * ca + d["y"] * sa          # координата вдоль волны
        bex = 0.1 * np.sin(2 * np.pi * xi)      # точное решение при t = 1
        order = np.argsort(xi)
        a.plot(xi[order], bperp[order], ".", ms=2, label=f"N={N}")
        errs.append((N,
                     np.mean(np.abs(bperp - bex)),
                     np.sqrt(np.mean((bperp - bex) ** 2))))
    xi_s = np.linspace(0, 2, 400)
    a.plot(xi_s, 0.1 * np.sin(2 * np.pi * xi_s), "k-", lw=1, label="точное")
    a.set_xlabel(r"$\xi = x\cos\alpha + y\sin\alpha$")
    a.set_ylabel(r"$B_\perp$")
    a.set_title(r"Альфвеновская волна, $\alpha=30°$, t = 1 (один период)")
    a.legend(); a.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig("alfven_profiles.png", dpi=150)
    print("-> alfven_profiles.png")

    print("\n N    L1(Bperp)      L2(Bperp)     порядок(L1)")
    prev = None
    for N, l1, l2 in errs:
        r = f"{np.log2(prev / l1):.2f}" if prev else "  – "
        print(f"{N:4d}  {l1:.6e}  {l2:.6e}   {r}")
        prev = l1


if __name__ == "__main__":
    mode, files = sys.argv[1], sys.argv[2:]
    {"briowu": plot_briowu, "ot": plot_ot,
     "rotor": plot_rotor, "alfven": plot_alfven}[mode](files)
