# this file contains the functions written and needed for the beam narrowing study

import numpy as np

def gauss2D(X, Y, wx, wy, x0, y0, alpha0, beta0):
    """Function that build the array of the magnetic field for a 2D Gaussian beam coming from a pinhole.
    The values are normalised to 1.

    Args:
        X (numpy array): 2D grid of the x coordinates
        Y (numpy array): 2D grid of the y coordinates
        wx (float): waist of the beam along the x-axis
        wy (float): waist of the beam along the y-axis
        x0 (float): center of the beam position on the x-axis
        y0 (float): center of the beam position on the y-axis
        alpha0 (float): center of the wave vector projection on the x-axis
        beta0 (float): center of the wave vector projection on the y-axis

    Returns:
        out (numpy array):  2D grid of the field (complex)
    """
    out = np.exp(-(((Y - y0) / wy) ** 2) - ((X - x0) / wx) ** 2) * np.exp(
        1j * (alpha0 * (X - x0) + beta0 * (Y - y0))
    )
    return out


def SO1(n1, f1, n2, SC, OA, Delta):
    """Function that simulate the first optic system of our simulation (lens + spherical diopter) and provides its magnification and geometric parameters.
    Args:
        n1 (float): refractive index of middle 1
        f1 (float): focal length of the lens
        n2 (float): refractive index of middle 2 (after de spherical diopter)
        SC (float): radius of curvature of the spherical diopter
        OA (float): distance between the lens center and the object
        Delta (float): optical interval or distance between the image focus of the lens and the object focus of the spherical diopter

    Returns:
        out (float tuple): determined quantities of the SO such as the lens to diopter and diopter to image and object to image distances; and the SO magnification.
    """
    # use of algebric distances in mm
    # lens
    OF1i = f1  # distance lens to image focus
    OF1o = -f1  # distance lens to object focus
    # spherical diopter
    SF2i = n2 * SC / (n2 - n1)  # distance top of the diopter to image focus
    SF2o = n1 * SC / (n1 - n2)  # distance top of the diopter to object focus
    # total centered optical system
    # focus
    if (Delta==0): # security to avoid a per 0 division, this configuration is automatically rejected
        return 0,0,0,0
    F1oFo = OF1i * OF1o / Delta  # distance lens object focus to total system object focus
    F2iFi = -SF2i * SF2o / Delta  # distance spherical diopter image focus to total system image focus
    # principal points
    HoFo = OF1o * SF2o / Delta  # distance object principal point to object focus (fo)
    HiFi = -SF2i * OF1i / Delta  # distance image principal point to image focus (fi)
    # intermediate quantities
    OFo = OF1o + F1oFo  # distance lens to total object focus
    FoA = OA - OFo  # distance total object focus to object
    if (FoA==0): # security to avoid a per 0 division, this configuration is automatically rejected
        return 0,0,0,0
    FiAi = HoFo * HiFi / FoA  # distance image focus to image (from Newton relation of the system)
    # results
    OS = OF1i + Delta - SF2o  # distance lens to diopter top
    SAi = SF2i + F2iFi + FiAi  # distance diopter top to image
    AAi = -OA + OF1i + Delta - SF2o + SAi  # distance object to image
    G = -HoFo / FoA  # magnification of the system
    return OS, SAi, AAi, G


def SO2(n1, f1, n2, SC, SA, Delta):
    """Function that simulate the second optic system of our simulation (lens + spherical diopter) and provides its magnification and geometric parameters.
    Args:
        n1 (float): refractive index of middle 1 (after the spherical diopter)
        f1 (float): focal length of the lens
        n2 (float): refractive index of middle 2 (before the spherical diopter)
        SC (float): radius of curvature of the spherical diopter
        SA (float): distance between the top of the spherical diopter and the object
        Delta (float): optical interval or distance between the image focus of the spherical diopter and the object focus of the lens

    Returns:
        out (float tuple): determined quantities of the SO such as the diopter to lens and lens to image and object to image distances; and the SO magnification.
    """
    # use of algebric distances in mm
    # lens
    OF2i = f1  # distance lens to image focus
    OF2o = -f1  # distance lens to object focus
    # spherical diopter
    SF1o = n2 * SC / (n2 - n1)  # distance top of the diopter to object focus
    SF1i = n1 * SC / (n1 - n2)  # distance top of the diopter to image focus
    # total centered optical system
    # focus
    if (Delta==0): # security to avoid a per 0 division, this configuration is automatically rejected
        return 0,0,0,0
    F1oFo = SF1i * SF1o / Delta  # distance spherical diopter object focus to total system object focus
    F2iFi = -OF2i * OF2o / Delta  # distance lens image focus to total system image focus
    # principal points
    HoFo = OF2o * SF1o / Delta  # distance object principal point to object focus (fo)
    HiFi = -SF1i * OF2i / Delta  # distance image principal point to image focus (fi)
    # intermediate quantities
    SFo = SF1o + F1oFo  # distance spherical diopter to total object focus
    FoA = SA - SFo  # distance total object focus to object
    if (FoA==0): # security to avoid a per 0 division, this configuration is automatically rejected
        return 0,0,0,0
    FiAi = HoFo * HiFi / FoA  # distance image focus to image (from Newton relation of the system)
    # results
    SO = SF1i + Delta - OF2o  # distance diopter top to lens
    OAi = OF2i + F2iFi + FiAi  # distance lens top to image
    AAi = -SA + SF1i + Delta - OF2o + OAi  # distance object to image
    G = -HoFo / FoA  # magnification of the system
    return SO, OAi, AAi, G