from rectification.utils import unitcheck
import numpy as np


def operating_line(x, R, xp):
    return R * x / (R + 1) + xp / (R + 1)


@unitcheck(L="kg/s", rhol="kg/m**3", Ft="m**2", res_unit="m**3/(m**2*s)")
def calc_Ucoef(L, rhol, Ft):
    return L / (rhol * Ft)


@unitcheck(Massl="g/mol", Massh="g/mol", mu_solv="Pa*s", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def calc_Diffcoef20(Massl, Massh, A , B, mu_solv, nul, nuh):
    """
    Calculates the diffusion coefficient at 20 degrees celcium.
    Parameters
    ----------
    Massl : float
    The molar mass of low-boilling component, [g/mol]
    Massh : float
    The molar mass of high-boilling component, [g/mol]
    A : float
    The correction coefficient depending on the properties solute, [dimensionless]
    B : float
    The correction coefficient depending on the properties solvent, [dimensionless]
    mu_solv : float
    The viscocity of solvent liquid, [Pa/s]
    nul : float
    The molar volume of solute, [sm**3/s]
    nuh : float
    The molar volume of solvent, [sm**3/s]
    Returns
    -------
    calc_Diffcoef20 : float
    The diffusion coefficient at 20 degrees celcium, [m**2/s]
    References
    ----------
    Романков, страница 289, формула 6.22
    """
    return 1e-6 * ((1/Massl) + (1/Massh))**0.5 / (A * B * mu_solv**0.5 * ((nul)**0.66 + (nuh)*0.66)**2)


def b(mul_mix, rhol_mix):
    """
    Calculates the temperature coefficient.
    Parameters
    ----------
    mul_mix : float
    The viscocity of low-boilling component of liquid, [Pa/s]
    rhol_mix : float
    The destinity of low-boilling component of liquid, [kg/m**3]
    Returns
    -------
    b : float
    The temperature coefficient, [dimensionless]
    References
    ----------
    Романков, страница 289, формула 6.24
    """
    return mul_mix**0.5/rhol_mix**0.66


def calc_Diffliq(calc_Diffcoef20, b, t_boil):
    """
    Calculates the diffusion coefficient of liquid phaze.
    Parameters
    ----------
    calc_Diffcoef20 : float
    The diffusion coefficient at 20 degrees celcium, [m**2/s]
    b : float
    The temperature coefficient, [dimensionless]
    t_boil : float
    The temperature of low-boilling component of liquid, [degrees celcium]
    Returns
    -------
    calc_Diffliq : float
    The diffusion coefficient of liquid phaze.
    References
    ----------
    Романков, страница 289, формула 6.23
    """
    return calc_Diffcoef20 * (1 + b * (t_boil - 20))


@unitcheck(t_boil="K", Massl="g/mol", Massh="g/mol", P_abs="Pa", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def calc_Diffvapor(t_boil, P_abs, Massl, Massh, nul, nuh):
    """
    Calculates the diffusion coefficient of vapor.
    Parameters
    ----------
    Massl : float
    The molar mass of the low-boilling component, [g/mol]
    Massh : float
    The molar mass of the high-boilling component, [g/mol]
    t_boil : float
    The boiling temperature of the low-boiling component, [K]
    P : float
    The absolute pressure of the column, [Pa]
    nul : float
    The molar volume of solute, [sm**3/s]
    nuh : float
    The molar volume of solvent, [sm**3/s]
    Returns
    -------
    calc_Diffvapor : float
    The diffusion coefficient of vapor, [m**2/s]
    References
    ----------
    Романков, страница 234, формула 6.25
    """
    return 4.22e-2 * t_boil**(3/2) * ((1/Massl) + (1/Massh))**0.5 / (P_abs * ((nul)**0.66 + (nuh)*0.66)**2)

