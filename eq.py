import math
def w_column(rliq, rvap):
    """
    Calculates working speed of column
    Parameters
    ----------
    rliq : [ kg / m**3 ]
    Mixture density of liquid
    rvap : [ kg / m**3 ]
    Mixture density of vapor
    Returns
    -------
    w : [ m / s ] 
    Working speed of column
    References
    ---------- 
    Дытнерский, формула 5.33, страница 205     
    """
    w = 0.05*((rliq/rgas)**0.5)
    return w

def mu_liq(mulow, muhigh, xlow)
    """
    Calculates  mixture viscosity of liquid
    Parameters
    ----------
    mulow : [ Pa * s ]
    Viscosity of lowboiling component
    muhigh : [ Pa * s ]
    Viscosity of highboiling component
    xlow : [ dimensionless ]
    Mass concentation of lowboiling component
    Returns
    -------
    muliq : [ Pa * s ]
    Mixture viscosity of liqud
    References:
    ----------
    Романков, формула 1.14, стр.15
    """
    muliq = math.log10(mulow) * xlow * math.log10(muhigh) * (1 - xlow)
    return muliq 