"""Script to convert soundings from SAM into NetCDF for use as 'constants in CAM.'"""

import numpy as np
import netCDF4 as nc4


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
    rho_a : np.ndarray
        array of air densities [kg/m^3] as floats corrresponding to altitudes, alt
    """
    alt, rho_a, rho_w = np.loadtxt(
        filename,
        comments="#",
        delimiter=None,
        skiprows=1,
        unpack=True,
    )
    return alt, rho_a


def save_to_netcdf(filename: str, alt: np.ndarray, rho_a: np.ndarray):
    """
    Save sounding data to NetCDF.

    Parameters
    ----------
    filename : str
        name of NetCFD file to save to
    alt : np.ndarray
        array of altitudes in [m] as floats
    rho_a : np.ndarray
        array of air densities [kg/m^3] as floats corrresponding to altitudes, alt

    Returns
    -------
    """
    ncfile = nc4.Dataset(filename, mode="w", format="NETCDF4_CLASSIC")

    ncfile.title = "SAM Soundings for neural net convection parameterisation"

    alt_dim = ncfile.createDimension("z", len(alt))

    alt_var = ncfile.createVariable("z", np.float64, ("z",))
    alt_var.units = "metres"
    alt_var.long_name = "altitude"
    alt_var.standard_name = "altitude"
    alt_var[:] = alt

    rho_a_var = ncfile.createVariable("rho", np.float64, ("z"))
    rho_a_var.units = "kg m-3"
    alt_var.long_name = "air density"
    rho_a_var.standard_name = "air_density"
    rho_a_var[:] = rho_a

    ncfile.close()


if __name__ == "__main__":
    # File to read
    SOUNDING_FILE = "sounding_z_rho_rhow.txt"
    NC_FILE = "SAM_sounding.nc"

    # Extract data from file as columns of z, rho, rhow
    altitude, rho = read_sounding(SOUNDING_FILE)

    # Save data to NetCDF file
    save_to_netcdf(NC_FILE, altitude, rho)
