import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import scipy.integrate as si

import isr_task.config as cf

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "DejaVu Sans",
        "axes.unicode_minus": False,
        "font.size": 15,
    }
)


def isr(params: dict = None) -> tuple[np.ndarray, np.ndarray]:
    """Calculate an incoherent scatter spectrum.

    Parameters
    ----------
    params : dict
        Give the physical parameters of the IS experiment.
        'N_f' -- number of data points along the frequency axis
        'f_min' -- lower bound of frequency axis
        'f_max' -- upper bound of frequency axis
        'f' -- radar frequency
        'T_e' -- electron temperature
        'T_i' -- ion temperature
        'n_e' -- electron number density
        'B' -- magnetic field strength
        'aspect' -- aspect angle
        'M' -- ion mass in amu

    Returns
    -------
    f_ax : np.ndarray
        Frequency axis
    IS : np.ndarray
        Echo power at each frequency
    """
    # Set physical parameters
    if params is None:
        N_f = cf.N_F
        f_min = -2e6
        f_max = 2e6
        f = cf.f  # 1/s - radar frequency
        T_e = cf.T_e  # K - electron temperature
        T_i = cf.T_i  # K - ion temperature
        n_e = cf.n_e  # 1/m^3 - electron number density
        B = cf.B  # T - magnetic field strength (towards Earth)
        aspect = cf.aspect  # degree - radar pointing direction to magnetic field line
        aspect = np.pi / 180 * aspect
        M_amu = cf.M  # amu - ion mass
        M = M_amu * (const.m_p + const.m_n) / 2  # Convert to kg
    else:
        N_f = params["N_f"]
        f_min = params["f_min"]
        f_max = params["f_max"]
        f = params["f"]
        T_e = params["T_e"]
        T_i = params["T_i"]
        n_e = params["n_e"]
        B = params["B"]
        aspect = params["aspect"]
        aspect = np.pi / 180 * aspect
        M_amu = params["M"]
        M = M_amu * (const.m_p + const.m_n) / 2  # Convert to kg
    nu = 0  # 1/s - collision frequency

    # Calculate constants
    k = -4 * np.pi * f / const.c
    l_D = debye(T_e, n_e)
    w_c = gyro("e", B)
    W_c = gyro("i", B, M_amu)

    # Susceptibility
    f_ax = np.linspace(f_min, f_max, N_f)  # Frequency axis
    # Integration variable of Gordeyev
    y_e = np.linspace(0, 1.5e-4 ** (1 / cf.ORDER), cf.N_Y) ** cf.ORDER
    y_i = np.linspace(0, 1.5e-2 ** (1 / cf.ORDER), cf.N_Y) ** cf.ORDER
    G_e = maxwellian_integrand(y_e, nu, k, aspect, T_e, w_c, const.m_e)
    G_i = maxwellian_integrand(y_i, nu, k, aspect, T_i, W_c, M)
    Fe = F(f_ax, y_e, nu, G_e)
    Fi = F(f_ax, y_i, nu, G_i)

    Xp = (1 / (2 * l_D**2 * k**2)) ** (1 / 2)
    chi_e = 2 * Xp**2 * Fe
    chi_i = 2 * Xp**2 * Fi

    # Calculate the IS spectrum (and ignore "divide" and "invalid" errors)
    with np.errstate(divide="ignore", invalid="ignore"):
        IS = (
            n_e
            / (np.pi * 2 * np.pi * f_ax)
            * (np.imag(Fe) * np.abs(1 + chi_i) ** 2 + np.imag(Fi) * np.abs(chi_e) ** 2)
            / (np.abs(1 + chi_e + chi_i) ** 2)
        )

    return f_ax, IS


def F(f_ax: np.ndarray, y: np.ndarray, nu: float, G: np.ndarray) -> np.ndarray:
    """Calculate the helper function 'F'.

    Parameters
    ----------
    f_ax : np.ndarray
        The frequency axis
    y : np.ndarray
        Axis of integration in the Gordeyev integral
    nu : float or int
        The collision frequency
    G : np.ndarray
        The integrand in the Gordeyev integral

    Returns
    -------
    func : np.ndarray
        The 'F' function
    """
    # Calculate the F functions that include susceptibility
    a = np.array([])
    for f in f_ax:
        w = 2 * np.pi * f
        sint = my_integration_method(w, y, G)
        a = np.r_[a, sint]

    func = 1 + (1j * 2 * np.pi * f_ax + nu) * a
    return func


def maxwellian_integrand(
    y: np.ndarray, nu: float, k: float, aspect: float, T: float, w_c: float, m: float
) -> np.ndarray:
    """Calculate a Maxwellian integrand for a Gordeyev integral.

    Parameters
    ----------
    y : np.ndarray
        Axis of integration in the Gordeyev integral
    nu : float or int
        The collision frequency
    k : float
        The radar wave number
    aspect : float
        The aspect anngle in radians
    T : float
        Temperature
    w_c : float
        Gyro frequency
    m : float
        Mass in kg

    Returns
    -------
    G : np.ndarray
        The Maxwellian integrand to be used in the Gordeyev integral
    """
    G = np.exp(
        -y * nu
        - k**2
        * np.sin(aspect) ** 2
        * T
        * const.k
        / (m * w_c**2)
        * (1 - np.cos(w_c * y))
        - 0.5 * (k * np.cos(aspect) * y) ** 2 * T * const.k / m
    )

    return G


def my_integration_method(w: float, y: np.ndarray, G: np.ndarray) -> np.ndarray:
    """A simple wrapper for integrating of the Gordeyev integral.

    Parameters
    ----------
    w : float
        The angular/signed frequency to be evaluated
    y : np.ndarray
        Axis of integration in the Gordeyev integral
    G : np.ndarray
        The integrand in the Gordeyev integral

    Returns
    -------
    sint : np.ndarray
        The value of the Gordeyev integral at each frequency data points
    """
    val = np.exp(1j * w * y) * G
    sint = si.simps(val, y)
    return sint


def debye(T: float, n: float) -> float:
    """Calculate the Debye length for elecrons.

    Parameters
    ----------
    T : float
        Temperature
    n : float
        Electron number density

    Returns
    -------
    l_D : float
        The Debye length
    """
    ep0 = 1e-9 / 36 / np.pi
    l_D = (ep0 * const.k * T / (n * const.e**2)) ** (1 / 2)
    return l_D


def gyro(p: str, B: float, m=16) -> float:
    """Calculate the gyro frequency of a particle species.

    Parameters
    ----------
    p : str
        A string specifying the particle species
        ({'e', 'i'} are recognized as electrons and ions)
    B : float
        The magnetic field strength
    m : float (default: 16)
        Particle mass in amu

    Returns
    -------
    w : float
        The gyro frequency
    """
    if p == "e":
        w = const.e * B / const.m_e
    elif p == "i":
        w = const.e * B / (m * (const.m_p + const.m_n) / 2)
    else:
        sys.exit(f"I do not know what kind of particle {p} is.")
    return w


def plot(dB=True):
    x, y = isr()
    y = 10 * np.log10(y) if dB else y

    plt.figure(figsize=(7, 4))
    plt.plot(x, y)
    plt.xlabel("Frequency $f$")
    plt.ylabel("Echo power [dB]")
    plt.grid(alpha=0.4)
    # plt.savefig('gyrolines.pdf')
    plt.show()


if __name__ == "__main__":
    print("Calculating isr spectrum...")
    plot()
