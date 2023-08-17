"""Script to convert soundings from SAM into NetCDF for use as 'constants in CAM.'"""

import os
from typing import Tuple
import time
import numpy as np
import netCDF4 as nc4


VERSION = "1.0"


def read_sounding(filename: str) -> np.ndarray:
    """
    Read SAM soundings from txt file into numpy arrays.

    Parameters
    ----------
    filename : str
        name of txt file to read

    Returns
    -------
    alt : np.ndarray
        array of altitudes in [m] as floats
    var_1 : np.ndarray
        array of floats corrresponding to altitudes, alt
        air densities [kg/m^3] or pressures [hPa]
    """
    alt, var_1, _ = np.loadtxt(
        filename,
        comments="#",
        delimiter=None,
        skiprows=1,
        unpack=True,
    )

    return alt, var_1


def read_grid(filename: str) -> np.ndarray:
    """
    Read SAM grid from grd file into numpy array.

    Parameters
    ----------
    filename : str
        name of txt file to read

    Returns
    -------
    alt : np.ndarray
        array of altitudes in [m] as floats
    """
    alt = np.loadtxt(
        filename,
    )

    return alt


def grid_calcs(alt: np.ndarray) -> Tuple[float, np.ndarray]:
    """
    Perform grid calculations as in setgrid.f90 of SAM for dz and adz.

    Parameters
    ----------
    alt : np.ndarray
        array of altitudes

    Returns
    -------
    dz : float
        dz at bottom of grid
    adz : np.ndarray
        array of adz grid values as floats
    """
    dz = 0.5 * (alt[0] + alt[1])

    # From setgrid.f90 (z = alt here):
    # z(k+1) = z(k)+(z(k)-z(k-1))
    # adz(1) = 1.
    # adz(k) = 0.5*(z(k+1)-z(k-1))/dz
    adz = np.ones(len(alt))
    adz[1:-1] = 0.5 * (alt[2:] - alt[:-2]) / dz
    adz[-1] = (alt[-1] - alt[-2]) / dz

    return dz, adz


def save_to_netcdf(
    filename: str,
    alt_r: np.ndarray,
    rho_a: np.ndarray,
    alt_p: np.ndarray,
    pres_a: np.ndarray,
    dz: float,
    adz: np.ndarray,
):
    """
    Save sounding data to NetCDF.

    Parameters
    ----------
    filename : str
        name of NetCFD file to save to
    alt_r : np.ndarray
        array of altitudes in [m] as floats correspondig to densities
    rho_a : np.ndarray
        array of air densities [kg/m^3] as floats corrresponding to altitudes, alt
    alt_p : np.ndarray
        array of altitudes in [m] as floats correspondig to pressures
    pres_a : np.ndarray
        array of air pressures [hPa] as floats corrresponding to altitudes, alt
    dz : float
        dz at bottom of grid
    adz : np.ndarray
        array of adz grid values as floats

    Returns
    -------
    """
    ncfile = nc4.Dataset(filename, mode="w", format="NETCDF4_CLASSIC")

    ncfile.title = "SAM Soundings for neural net convection parameterisation"
    ncfile.history = f"Created: {time.ctime(time.time())} by user {os.getlogin()}"
    ncfile.software = (
        f"Generated using sounding_to_netcdf.py v {VERSION} by Jack Atkinson (ICCS)"
    )
    ncfile.software = (
        f"Generated using sounding_to_netcdf.py v {VERSION} by Jack Atkinson (ICCS)"
    )

    alt_r_dim = ncfile.createDimension("z_r", len(alt_r))
    alt_p_dim = ncfile.createDimension("z_p", len(alt_p))
    alt_dz_dim = ncfile.createDimension("dz", 1)

    alt_r_var = ncfile.createVariable("z_r", np.float64, ("z_r",))
    alt_r_var.units = "metres"
    alt_r_var.long_name = "air density altitudes"
    alt_r_var.standard_name = "altitude"
    alt_r_var[:] = alt_r

    rho_a_var = ncfile.createVariable("rho", np.float64, ("z_r"))
    rho_a_var.units = "kg m-3"
    rho_a_var.long_name = "air density"
    rho_a_var.standard_name = "air_density"
    rho_a_var[:] = rho_a

    alt_p_var = ncfile.createVariable("z_p", np.float64, ("z_p",))
    alt_p_var.units = "metres"
    alt_p_var.long_name = "air pressure altitudes"
    alt_p_var.standard_name = "altitude"
    alt_p_var[:] = alt_p

    pres_var = ncfile.createVariable("pres", np.float64, ("z_p"))
    pres_var.units = "Pa"
    pres_var.long_name = "air pressure"
    pres_var.standard_name = "air_pressure"
    pres_var[:] = 100.0 * pres_a

    dz_var = ncfile.createVariable("dz", np.float64, ("dz"))
    dz_var.units = "m"
    dz_var.long_name = "delta_z for lowest grid cell"
    dz_var.standard_name = "magnitude_of_derivative_of_position_wrt_model_level_number"
    dz_var[:] = dz

    adz_var = ncfile.createVariable("adz", np.float64, ("z_p"))
    adz_var.units = "1"
    adz_var.long_name = "delta_z for all cells scaled by lowest delta_z (dz)"
    adz_var[:] = adz

    ncfile.close()


if __name__ == "__main__":
    # File to read
    RHO_SOUNDING_FILE = "sounding_z_rho_rhow.txt"
    PRES_SOUNDING_FILE = "sounding_z_pres0_tabs0.txt"
    GRD_FILE = "grd"
    NC_FILE = "SAM_sounding.nc"

    # Extract data from files as columns of z, rho, pres
    altitude_r, rho = read_sounding(RHO_SOUNDING_FILE)
    altitude_p, pres = read_sounding(PRES_SOUNDING_FILE)
    altitude_g = read_grid(GRD_FILE)

    dz, adz = grid_calcs(altitude_g)

    # Save data to NetCDF file
    save_to_netcdf(NC_FILE, altitude_r, rho, altitude_p, pres, dz, adz)
