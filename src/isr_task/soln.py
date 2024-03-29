import time

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numba as nb
import numpy as np
import scipy.constants as const

import isr_task.config as cf
import isr_task.isr2 as isr
import isr_task.py_plot as pylt


def prob3():
    t0 = time.perf_counter()
    params = {
        "N_f": int(1e4),
        "f_min": -2e6,
        "f_max": 2e6,
        "f": 430e6,
        "n_e": 2e10,
        "B": 3.5e-5,
        "M": 16,
        "T_e": 200,
        "T_i": 200,
        "aspect": 135,
    }
    x, y = isr.isr(params)
    y = 10 * np.log10(y)
    t1 = time.perf_counter()
    print(f"Took {t1-t0:.2f} seconds")

    plt.figure(figsize=(10, 6))
    plt.plot(x, y)
    plt.xlabel("Frequency $f$")
    plt.ylabel("Echo power [dB]")
    plt.grid(alpha=0.2)
    # plt.savefig('docs/gyrolines.pdf')
    plt.show()


def prob4():
    t0 = time.perf_counter()
    params = {
        "N_f": int(1e3),
        "f_min": 3.5e6,
        "f_max": 7.5e6,
        "f": 933e6,
        "n_e": 2e11,
        "B": 5e-5,
        "M": 16,
        "T_e": 2000,
        "T_i": 2000,
        "aspect": 180,
    }

    T_E = [8000, 6000, 4000, 2000]
    data = []
    for T_e in T_E:
        params["T_e"] = T_e
        x, y = isr.isr(params)
        data.append((x, y))
    t1 = time.perf_counter()
    print(f"Took {t1-t0:.2f} seconds")

    lab = [r"$T_\mathrm{e}=\,$" + f"{t_e}" + r"$\,\mathrm{K}$" for t_e in T_E]
    plt.figure("temps", figsize=(10, 6))
    pylt.ridge_plot(
        data,
        "squeeze",
        "grid",
        xlabel="Frequency $f$",
        ylabel="Echo power",
        labels=lab,
        figname="temps",
    )
    # for d, l in zip(data, lab):
    #     plt.plot(d[0], d[1], label=l)
    # plt.legend()
    # plt.savefig('docs/temps.pdf', dpi=200)
    plt.show()


def prob5():
    t0 = time.perf_counter()
    params = {
        "N_f": int(1e3),
        "f_min": -5e3,
        "f_max": 5e3,
        "f": 430e6,
        "n_e": 2e10,
        "B": 5e-5,
        "M": 16,
        "T_e": 300,
        "T_i": 200,
        "aspect": 180,
    }

    T_I = [400, 300, 200, 100]
    data = []
    for T_i in T_I:
        params["T_i"] = T_i
        params["T_e"] = T_i * 1.5
        x, y = isr.isr(params)
        data.append((x, y))
    t1 = time.perf_counter()
    print(f"Took {t1-t0:.2f} seconds")

    lab = [r"$T_\mathrm{e}=\,$" + f"{t_e}" + r"$\,\mathrm{K}$" for t_e in T_I]
    plt.figure("temps", figsize=(10, 6))
    c = ["r", "g", "b", "magenta"]
    plt.grid(True, which="major", ls="-", alpha=0.2)
    for (x, y), col, l in zip(data[::-1], c, lab[::-1]):
        plt.plot(x, y, f"{col}", label=l)
    plt.legend()
    # plt.tick_params(axis='both', which='both', labelbottom=False, labelleft=False)
    plt.xlabel(r"Frequency $f$")
    plt.ylabel("Echo power")
    # pylt.ridge_plot(data, 'squeeze', 'grid', xlabel='Frequency $f$', \
    #         ylabel='Echo power [dB]', labels=lab, figname='temps', ylim=(- 7e4, 6e5))
    # plt.savefig('docs/ionline_soln.pdf')
    plt.show()


def extra_mov(save=False):
    # w, x, y = extra_data()
    with np.load("docs/assets/video_data_ions_vectorized.npz", mmap_mode="r") as f:
        w = f["w"]
        x = f["x_ax"]
        y = f["y_ax"]
    if x.shape[0] != y.shape[0]:
        y = y.T

    fig, ax = plt.subplots(figsize=(12, 6))
    # Electrons
    # ax.set(xlim=(-1e-6, 6e-6), ylim=(-1.1, 1.1))
    # ax.set(xlim=(- 1e-5, 1.6e-4), ylim=(-1.1, 1.1))
    # Ions
    ax.set(xlim=(-2e-6, 1e-5), ylim=(-1.1, 1.1))
    # ax.set(xlim=(- 1e-5, 1.6e-4), ylim=(-1.1, 1.1))
    inityr = np.real(y[:, 0])
    inityi = np.imag(y[:, 0])
    line1 = ax.plot(x, inityr, "r--", lw=0.7, label="Real")[0]
    line2 = ax.plot(x, inityi, "b--", lw=0.7, label="Imag")[0]
    ax.legend()

    def draw(frame, add_colorbar):
        y_ = y[:, frame]
        y_r = np.real(y_)
        y_i = np.imag(y_)
        line1.set_ydata(y_r)
        line2.set_ydata(y_i)
        f = w[frame] / (2 * np.pi)
        title = f"f = {f:.2e}"
        ax.set_title(title)
        return ax

    def init():
        return draw(0, add_colorbar=True)

    def animate(frame):
        return draw(frame, add_colorbar=False)

    frames = len(w)

    ani = animation.FuncAnimation(
        fig, animate, frames, interval=100, init_func=None, repeat=False
    )
    plt.show()
    if save:
        ani.save("video_ions_zoom.mp4", writer=animation.FFMpegWriter(fps=60), dpi=400)
    plt.close(fig)


@nb.njit(parallel=True, cache=True)
def extra_run() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    M = const.m_e
    # M = 16 * (const.m_p + const.m_n) / 2
    B = 3.5e-5
    T = cf.T_e
    aspect = cf.aspect * np.pi / 180
    k = -4 * np.pi * cf.f / const.c
    w_c = const.e * B / M  # isr.gyro("i", cf.B, M)
    # w_c = const.e * B / M  # isr.gyro("i", cf.B, M)
    w = 2 * np.pi * np.linspace(0, int(1e7) ** (1 / cf.ORDER), int(1e4)) ** cf.ORDER
    y = np.linspace(0, 1.5e-4 ** (1 / cf.ORDER), cf.N_Y) ** cf.ORDER
    # y = np.linspace(0, 1.5e-2 ** (1 / cf.ORDER), cf.N_Y) ** cf.ORDER
    G = isr.maxwellian_integrand(y, 0, k, aspect, T, w_c, const.m_e)
    # G = isr.maxwellian_integrand(y, 0, k, aspect, T, w_c, M)
    ## Original implementation // Took 261 seconds
    # data = np.ndarray
    # for w_ in w:
    #     val = np.exp(1j * w_ * y) * G
    #     if isinstance(data, type):
    #         data = val
    #     else:
    #         data = np.c_[data, val]
    ## Vectorized // 5.367 seconds (more load on mem, 1e4 freqs was killed...)
    # data = np.exp(np.einsum('i,j->ij', w, 1j * y)) * G
    ## Numba implementation // Took 3 sek, 15 with 1e4 freqs (~8 sec with cache on successive)
    data = np.zeros((len(w), len(y)), dtype=np.complex64)
    for i in nb.prange(len(w)):
        val = np.exp(1j * w[i] * y) * G
        data[i, :] = val
    return y, np.asarray(data), w


def extra_data():
    t0 = time.perf_counter()
    y, data, w = extra_run()
    t1 = time.perf_counter()
    print(f"Took {t1-t0:.3f} seconds")
    # x_ax = y
    # y_ax = data
    # np.savez("video_data_ions_vectorized.npz", w=w, x_ax=x_ax, y_ax=y_ax)


def main():
    # extra_mov()
    # extra_data()
    prob3()
    prob4()
    prob5()

if __name__ == '__main__':
    main()
