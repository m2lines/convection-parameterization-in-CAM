"""Script to plot the results from running the YOG parameterisation on CAM profiles."""

import netCDF4 as nc4
import matplotlib.pyplot as plt
import numpy as np


def get_ncvals(data, varname):
    """Extract values from netCDF dataset for a particular variable."""
    return data[varname][:]


def get_ncvar(data, varname, coordname, varlabel="", coordlabel=r"$\hat p$ [-]"):
    """Extract variable from netCDF dataset as a dict of info."""
    return {
        "name": varname,
        "varlabel": varlabel,
        "values": get_ncvals(data, varname),
        "coordlabel": coordlabel,
        "coords": get_ncvals(data, coordname),
    }


def get_coord(data, coordname):
    """Extract the normalised pressure coordinates"""
    # return {"name": coordname, "values": get_ncvals(data, coordname), "coords": np.linspace(0, 1, len(get_ncvals(data, coordname)))}
    return {
        "name": coordname,
        "values": get_ncvals(data, coordname),
        "coords": np.arange(len(get_ncvals(data, coordname))),
    }


def comparison_plot(vars, xlab="", ylab="", title="", n=0):
    """
    Plot profiles for comparison based on inputs

    This assumes we are plotting soundings against pressure, so pressure is on y-axis
    and variable is on x-axis. We also flip the y-axis so highest pressure is at the
    bottom.
    """
    for var in vars:
        plt.scatter(var["values"][n,:], var["coords"][n,:], label=var["name"])

    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.gca().invert_yaxis()

    plt.title(title)

    plt.show()


def norm_comparison_plot(vars, xlab="", ylab="", title="", n=0):
    """
    Plot normalised profiles for shape comparison based on inputs

    This assumes we are plotting soundings against pressure, so pressure is on y-axis
    and variable is on x-axis. We also flip the y-axis so highest pressure is at the
    bottom.
    """
    for var in vars:
        plt.scatter(
            var["values"][n,:] / np.max(var["values"][n,:]), var["coords"][n,:], label=var["name"]
        )

    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.gca().invert_yaxis()

    plt.title(title)

    plt.show()


def conversion_plot(vars, xlab="", ylab="", title="", n=0):
    """
    Plot profiles for checking variable conversion based on inputs.

    This assumes we are plotting soundings against pressure, so pressure is on y-axis
    and variable is on x-axis. We also flip the y-axis so highest pressure is at the
    bottom.
    """

    # Set up a plots with shared yaxes
    # Plot each of the variables in the list.
    fig, ax = plt.subplots(1, len(vars), sharey=True)

    for ival, var in enumerate(vars):
        ax[ival].scatter(var["values"][n,:], var["coords"][n,:], label=var["name"])
        ax[ival].set_xlabel(var["varlabel"])
        ax[ival].legend()
    plt.gca().invert_yaxis()

    ax[0].set_ylabel(ylab)

    fig.suptitle(title)

    plt.show()


def coord_plot(vars, xlab="", ylab="", title=""):
    """
    Plot coordinates for comparison.

    This assumes coordinates are pressures, so we flip the x-axis so highest pressure
    is at the left.
    """
    for var in vars:
        plt.scatter(var["coords"][n,:], var["values"][n,:], label=var["name"])

    plt.legend()
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.gca().invert_yaxis()

    plt.title(title)

    plt.show()


if __name__ == "__main__":
    data = nc4.Dataset("NN_test_output.nc")

    pnorm_cam = get_coord(data, "PNORM_CAM")
    pnorm_sam = get_coord(data, "PNORM_SAM")

    # Check Coordinates
    #    coord_plot([pnorm_cam, pnorm_sam], ylab=r"$\hat p$ [-]", xlab=r"$i$ [-]")

    # Check interpolation of inputs from CAM grid to SAM grid
    tabs_cam_in = get_ncvar(data, "TABS_CAM_IN", "PNORM_CAM", varlabel=r"$T$ [-]")
    tabs_sam_in = get_ncvar(data, "TABS_SAM_IN", "PNORM_SAM", varlabel=r"$T$ [-]")

    qv_cam_in = get_ncvar(data, "QV_CAM_IN", "PNORM_CAM", varlabel=r"$q_v$ [-]")
    qv_sam_in = get_ncvar(data, "QV_SAM_IN", "PNORM_SAM", varlabel=r"$q_v$ [-]")

    qc_cam_in = get_ncvar(data, "QC_CAM_IN", "PNORM_CAM", varlabel=r"$q_c$ [-]")
    qc_sam_in = get_ncvar(data, "QC_SAM_IN", "PNORM_SAM", varlabel=r"$q_c$ [-]")

    qi_cam_in = get_ncvar(data, "QI_CAM_IN", "PNORM_CAM", varlabel=r"$q_i$ [-]")
    qi_sam_in = get_ncvar(data, "QI_SAM_IN", "PNORM_SAM", varlabel=r"$q_i$ [-]")

    comparison_plot(
        [tabs_sam_in, tabs_cam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$T$ [K]",
        title="Comparison of T [K] for CAM and SAM grids to check interpolation.",
    )
    comparison_plot(
        [qv_sam_in, qv_cam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q_{v}$ [-]",
        title=r"Comparison of $q_v$ [-] for CAM and SAM grids to check interpolation.",
    )
    comparison_plot(
        [qc_sam_in, qc_cam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q_{v}$ [-]",
        title=r"Comparison of $q_c$ [-] for CAM and SAM grids to check interpolation.",
    )
    comparison_plot(
        [qi_sam_in, qi_cam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q_{v}$ [-]",
        title=r"Comparison of $q_i$ [-] for CAM and SAM grids to check interpolation.",
    )

    # Check variable conversion from CAM variables to SAM variables
    t_sam_in = get_ncvar(data, "T_SAM_IN", "PNORM_SAM", varlabel=r"$t$ [K]")
    q_sam_in = get_ncvar(data, "Q_SAM_IN", "PNORM_SAM", varlabel=r"$q_{tot}$ [-]")
    gamaz_sam = get_ncvar(
        data, "GAMAZ_SAM", "PNORM_SAM", varlabel=r"$g \cdot z / c_p$ [K]"
    )

    comparison_plot(
        [qv_sam_in, q_sam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q$ [-]",
        title=r"Comparison of $q_v$ [-] from CAM and $q$ [-] from SAM grids to check variable conversion. Should match as $q_i$ and $q_c$ are initially 0.",
    )
    comparison_plot(
        [tabs_sam_in, gamaz_sam, t_sam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$T$ [K]",
        title=r"Comparison of $T$ [K] from CAM, $gz/c_p$ component of SAM, and $t$ in SAM to check variable conversion.",
    )
    conversion_plot(
        [tabs_sam_in, gamaz_sam, t_sam_in],
        ylab=r"$\hat p$ [-]",
        title=r"Comparison of $T$ [K] from CAM, $gz/c_p$ component of SAM, and $t$ in SAM to check variable conversion.",
    )

    comparison_plot(
        [q_sam_in, qv_sam_in, qc_sam_in, qi_sam_in],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q$ [-]",
        title=r"Comparison of $q$ components [-] from CAM, and $q$ value in SAM (sent to NN) to check variable conversion.",
    )

    # From the NN
    # Fluxes
    t_flx_adv_sam = get_ncvar(data, "T_FLX_ADV", "PNORM_SAM", varlabel=r"$dt$ [K]")
    t_tend_auto_sam = get_ncvar(data, "T_TEND_AUTO", "PNORM_SAM", varlabel=r"$dt$ [K]")
    t_flx_sed_sam = get_ncvar(data, "T_FLX_SED", "PNORM_SAM", varlabel=r"$dt$ [K]")
    q_flx_adv_sam = get_ncvar(data, "Q_FLX_ADV", "PNORM_SAM", varlabel=r"$dq$ [-]")
    q_tend_auto_sam = get_ncvar(data, "Q_TEND_AUTO", "PNORM_SAM", varlabel=r"$dq$ [-]")
    q_flx_sed_sam = get_ncvar(data, "Q_FLX_SED", "PNORM_SAM", varlabel=r"$dq$ [-]")
    # Deltas
    dt_rad_sam = get_ncvar(data, "DT_RAD", "PNORM_SAM", varlabel=r"$dt$ [K]")
    dt_adv_sam = get_ncvar(data, "DT_ADV", "PNORM_SAM", varlabel=r"$dt$ [K]")
    dt_auto_sam = get_ncvar(data, "DT_AUTO", "PNORM_SAM", varlabel=r"$dt$ [K]")
    dt_sed_sam = get_ncvar(data, "DT_SED", "PNORM_SAM", varlabel=r"$dt$ [K]")
    dq_adv_sam = get_ncvar(data, "DQ_ADV", "PNORM_SAM", varlabel=r"$dq$ [-]")
    dq_auto_sam = get_ncvar(data, "DQ_AUTO", "PNORM_SAM", varlabel=r"$dq$ [-]")
    dq_sed_sam = get_ncvar(data, "DQ_SED", "PNORM_SAM", varlabel=r"$dq$ [-]")

    # Plot fluxes and tendencies for t and q as they come out of the NN
    comparison_plot(
        [q_flx_adv_sam, q_tend_auto_sam, q_flx_sed_sam],
        ylab=r"$\hat p$ [-]",
        xlab=r"Flux $q$ [-]",
        title=r"Comparison of $q$-based outputs from the NN.",
    )
    comparison_plot(
        [t_flx_adv_sam, t_tend_auto_sam, t_flx_sed_sam],
        ylab=r"$\hat p$ [-]",
        xlab=r"Flux $T$ [K]",
        title=r"Comparison of $t$-based outputs from the NN.",
    )

    # Plot deltas for T and q as calculated by the code
    comparison_plot(
        [dq_adv_sam, dq_auto_sam, dq_sed_sam],
        ylab=r"$\hat p$ [-]",
        xlab=r"$dq$ [-]",
        title=r"Comparison of $q$ delta components [-] from the parameterisation.",
    )
    comparison_plot(
        [dt_rad_sam, dt_adv_sam, dt_auto_sam, dt_sed_sam],
        ylab=r"$\hat p$ [-]",
        xlab=r"$dT$ [K]",
        title=r"Comparison of $t$ delta components [K] outputs from the parameterisation.",
    )

    # After NN parameterisation
    t_sam_out = get_ncvar(data, "T_SAM_OUT", "PNORM_SAM", varlabel=r"$t$ [K]")
    q_sam_out = get_ncvar(data, "Q_SAM_OUT", "PNORM_SAM", varlabel=r"$q_{tot}$ [-]")
    # Plot comparison to T/q before and after
    comparison_plot(
        [q_sam_in, q_sam_out],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q$ [-]",
        title=r"Comparison of $q$ [-] on SAM grid before and after parameterisation.",
    )
    comparison_plot(
        [t_sam_in, t_sam_out],
        ylab=r"$\hat p$ [-]",
        xlab=r"$T$ [K]",
        title=r"Comparison of $T$ [K] on SAM grid before and after parameterisation.",
    )

    # Compare to ZM tendencies output by CAM
    dt_yog = get_ncvar(data, "YOGDT", "PNORM_CAM", varlabel=r"$dt$ [K]")
    dq_yog = get_ncvar(data, "YOGDQ", "PNORM_CAM", varlabel=r"$dq$ [-]")
    dt_zm = get_ncvar(data, "ZMDT", "PNORM_CAM", varlabel=r"$dt$ [K]")
    dq_zm = get_ncvar(data, "ZMDQ", "PNORM_CAM", varlabel=r"$dq$ [-]")
    comparison_plot(
        [dq_zm, dq_yog],
        ylab=r"$\hat p$ [-]",
        xlab=r"$q$ [-]",
        title=r"Comparison of $q$ tendencies [-] from ZM and this routine on CAM grid.",
    )
    comparison_plot(
        [dt_zm, dt_yog],
        ylab=r"$\hat p$ [-]",
        xlab=r"$T$ [K]",
        title=r"Comparison of $T$ tendencies [K] from ZM and this routine on CAM grid.",
    )
    norm_comparison_plot(
        [dq_zm, dq_yog],
        ylab=r"$\hat p$ [-]",
        xlab=r"$\hat q$ [-]",
        title=r"Comparison of normalised $q$ tendencies from ZM and this routine on CAM grid.",
    )
    norm_comparison_plot(
        [dt_zm, dt_yog],
        ylab=r"$\hat p$ [-]",
        xlab=r"$\hat T$ [-]",
        title=r"Comparison of normalised $T$ tendencies from ZM and this routine on CAM grid.",
    )