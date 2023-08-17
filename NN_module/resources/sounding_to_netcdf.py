"""Script to convert soundings from SAM into NetCDF for use as 'constants in CAM.'"""

import os
from typing import Tuple
import time
import numpy as np
import netCDF4 as nc4


VERSION = "1.1"


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
        may be air densities [kg/m^3] or pressures [hPa] depending on sounding
    """
    alt, var_1, _ = np.loadtxt(
        filename,
        comments="#",
        delimiter=None,
        skiprows=1,
        unpack=True,
    )

    # Flip the arrays before returning as soundings are given in descending altitude:
    return np.flip(alt), np.flip(var_1)


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
    alt: np.ndarray,
    rho_a: np.ndarray,
    dz: float,
    adz: np.ndarray,
):
    """
    Save sounding data to NetCDF.

    Parameters
    ----------
    filename : str
        name of NetCFD file to save to
    alt : np.ndarray
        array of altitudes in [m] as floats corresponding to variables
    rho_a : np.ndarray
        array of air densities [kg/m^3] as floats corrresponding to altitudes, alt
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
        f"Generated using sounding_to_netcdf.py v{VERSION} by Jack Atkinson (ICCS)"
    )

    alt_dim = ncfile.createDimension("z", len(alt))
    dz_dim = ncfile.createDimension("dz", 1)

    alt_var = ncfile.createVariable("z", np.float64, ("z",))
    alt_var.units = "m"
    alt_var.long_name = "altitudes"
    alt_var.standard_name = "altitude"
    alt_var[:] = alt

    rho_a_var = ncfile.createVariable("rho", np.float64, ("z"))
    rho_a_var.units = "kg m-3"
    rho_a_var.long_name = "air density"
    rho_a_var.standard_name = "air_density"
    rho_a_var[:] = rho_a

    dz_var = ncfile.createVariable("dz", np.float64, ("dz"))
    dz_var.units = "m"
    dz_var.long_name = "delta_z for lowest grid cell"
    dz_var.standard_name = "magnitude_of_derivative_of_position_wrt_model_level_number"
    dz_var[:] = dz

    adz_var = ncfile.createVariable("adz", np.float64, ("z"))
    adz_var.units = "1"
    adz_var.long_name = "delta_z for all cells scaled by lowest delta_z (dz)"
    adz_var[:] = adz

    ncfile.close()


if __name__ == "__main__":
    # Files to read/save
    RHO_SOUNDING_FILE = "sounding_z_rho_rhow.txt"
    GRD_FILE = "grd"
    NC_FILE = "SAM_sounding.nc"

    # Extract data from files as columns of z, rho, pres
    altitude_r, rho = read_sounding(RHO_SOUNDING_FILE)
    altitude_g = read_grid(GRD_FILE)

    # Perform grid calculations
    dz_g, adz_g = grid_calcs(altitude_g)

    # Save data to NetCDF file
    save_to_netcdf(NC_FILE, altitude_g, rho, dz_g, adz_g)
