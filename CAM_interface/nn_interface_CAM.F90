module nn_interface_CAM
    !! Interface to convection parameterisation for the CAM model
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use netcdf
use precision, only: dp
use nn_convection_flux_mod, only: nn_convection_flux, &
                                  nn_convection_flux_init, nn_convection_flux_finalize, &
                                  esati, rsati, esatw, rsatw, dtrsatw, dtrsati
use SAM_consts_mod, only: nrf, ggr, cp, tbgmax, tbgmin, tprmax, tprmin, &
                          fac_cond, fac_sub, fac_fus, &
                          a_bg, a_pr, an, bn, ap, bp, &
                          omegan, check
#ifdef CAM_PROFILE
use nf_out, only: nf_write_sam, nf_write_cam, nf_write_scalar
#endif

implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux_CAM, &
        nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize

! Make these routines public for purposes of testing,
! otherwise not required outside this module.
public interp_to_sam, interp_to_cam
public SAM_var_conversion, CAM_var_conversion
public fetch_sam_data
public calculate_dry_mixing_ratio, calculate_moist_mixing_ratio

!---------------------------------------------------------------------
! local/private data

!= unit 1 :: nz_sam
integer :: nz_sam
    !! number of vertical values in the SAM sounding profiles

!= unit m :: z
!= unit hPa :: pres, presi
!= unit kg m-3 :: rho
!= unit 1 :: adz
!= unit K :: gamaz
real(dp), allocatable, dimension(:) :: z, pres, presi, rho, adz, gamaz
    !! SAM sounding variables

!= unit m :: dz
real(dp) :: dz
    !! grid spacing in z direction for the lowest grid layer

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine fetch_sam_data(pressam, presisam, gamazsam, rhosam, zsam)
        !! Temporary subroutine to extract SAM pressure profile
        real(dp), dimension(48), intent(out) :: pressam, presisam, gamazsam, rhosam, zsam
        pressam(:) = pres(:)
        presisam(:) = presi(:)
        gamazsam(:) = gamaz(:)
        rhosam(:) = rho(:)
        zsam(:) = z(:)

    end subroutine fetch_sam_data


    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_CAM_init(nn_filename, sounding_filename)
        !! Initialise the NN module

        character(len=136), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights
        character(len=136), intent(in) :: sounding_filename
            !! NetCDF filename from which to read SAM sounding and grid data

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename)

        ! Read in the SAM sounding and grid data required
        call sam_sounding_init(sounding_filename)

    end subroutine nn_convection_flux_CAM_init


    subroutine nn_convection_flux_CAM(pres_cam, pres_int_cam, pres_sfc_cam, &
                                      tabs_cam, qv_cam, qc_cam, qi_cam, &
                                      cp_cam, &
                                      dtn, &
                                      nx, nz, &
                                      precsfc, &
                                      dqi, dqv, dqc, ds)
        !! Interface to the nn_convection parameterisation for the CAM model

        real(dp), intent(in) :: dtn    ! Seconds
        integer, intent(in) :: nx, nz

        real(dp), intent(in) :: cp_cam
            !! specific heat capacity of dry air from CAM [J/kg/K]
        real(dp), dimension(:,:), intent(in) :: pres_cam
            !! pressure [Pa] from the CAM model
        real(dp), dimension(:,:), intent(in) :: pres_int_cam
            !! interface pressure [Pa] from the CAM model
        real(dp), dimension(:), intent(in) :: pres_sfc_cam
            !! surface pressure [Pa] from the CAM model
        real(dp), dimension(:,:), intent(in) :: tabs_cam
            !! absolute temperature [K] from the CAM model
        real(dp), dimension(:,:), intent(in) :: qv_cam, qc_cam, qi_cam
            !! moisture content [kg/kg] from the CAM model

        real(dp), dimension(nx, nz), intent(out) :: dqi, dqv, dqc, ds
            !! Output tendencies - CAM variables on the CAM grid
        real(dp), dimension(:), intent(out) :: precsfc
            !! Output surface precipitation in CAM

! Local input variables on the SAM grid
        real(dp), dimension(nx, nrf) :: tabs_sam
            !! absolute temperature [K] on the SAM grid
        real(dp), dimension(nx, nrf) :: r_sam
            !! non-precipitating water mixing ratio [kg/kg] on the SAM grid
        real(dp), dimension(nx, nrf) :: t_sam
            !! static energy [J/kg] on the SAM grid
        real(dp), dimension(nx, nrf) :: qi_sam, qv_sam, qc_sam
            !! mixing ratios [kg/kg] on the SAM grid
        real(dp), dimension(nx, nrf) :: qi0_sam, qv0_sam, qc0_sam, tabs0_sam, r0_sam
            !! Variables at the start of the parameterisation on the SAM grid - used to calculate tendencies
        real(dp), dimension(nx, nrf) :: dqi_sam, dqv_sam, dqc_sam, ds_sam
            !! CAM variable tendencies calculated on the SAM grid
        real(dp), dimension(nx)      :: qi_surf, qv_surf, qc_surf, tabs_surf
            !! Surface values
        real(dp) :: y_in(nx)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        real(dp), dimension(nx)      :: precsfc_i
            !! precipitation at surface from one call to parameterisation

#ifdef CAM_PROFILE
        real(dp), dimension(nx, nz) :: rh_cam
        real(dp), dimension(nx, nrf) :: rh_sam
        real(dp) :: vap_pres
#endif
        integer :: k

        ! Initialise precipitation to 0 if required and at start of cycle if subcycling
        precsfc(:)=0.

        ! distance to the equator
        ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
        ! TODO: Set y_in as appropriate for CAM
        y_in(:) = 0.0

        !-----------------------------------------------------

#ifdef CAM_PROFILE
        call nf_write_cam(tabs_cam(1,:), "TABS_CAM_IN")
        call nf_write_cam(qi_cam(1,:), "QI_CAM_IN")
        call nf_write_cam(qc_cam(1,:), "QC_CAM_IN")
        call nf_write_cam(qv_cam(1,:), "QV_CAM_IN")
#endif

        ! Interpolate CAM variables to the SAM pressure levels
        ! TODO Interpolate all variables in one call
        ! Set surface values
        call extrapolate_to_surface(pres_cam, qi_cam, pres_sfc_cam, qi_surf)
        where (qi_surf>=0)
            qi_surf = qi_surf
        elsewhere
            qi_surf = 0
        end where
        call extrapolate_to_surface(pres_cam, qc_cam, pres_sfc_cam, qc_surf)
        where (qc_surf>=0)
            qc_surf = qc_surf
        elsewhere
            qc_surf = 0
        end where
        call extrapolate_to_surface(pres_cam, qv_cam, pres_sfc_cam, qv_surf)
        where (qv_surf>=0)
            qv_surf = qv_surf
        elsewhere
            qv_surf = 0
        end where
        call extrapolate_to_surface(pres_cam, tabs_cam, pres_sfc_cam, tabs_surf)

        call interp_to_sam(pres_cam, pres_sfc_cam, &
                           tabs_cam, tabs_sam, tabs_surf)
        call interp_to_sam(pres_cam, pres_sfc_cam, &
                           qi_cam, qi_sam, qi_surf)
        call interp_to_sam(pres_cam, pres_sfc_cam, &
                           qc_cam, qc_sam, qc_surf)
        call interp_to_sam(pres_cam, pres_sfc_cam, &
                           qv_cam, qv_sam, qv_surf)

#ifdef CAM_PROFILE
        call nf_write_sam(tabs_sam(1,:), "TABS_SAM_IN")
        call nf_write_sam(qi_sam(1,:), "QI_SAM_IN")
        call nf_write_sam(qc_sam(1,:), "QC_SAM_IN")
        call nf_write_sam(qv_sam(1,:), "QV_SAM_IN")

        call nf_write_cam(pres_cam(1,:), "P_CAM")
        call nf_write_sam(pres(1:nrf) * 100D0, "P_SAM")
        call nf_write_cam(pres_cam(1,:) / pres_sfc_cam(1), "PNORM_CAM")
        call nf_write_sam(pres(1:nrf) / presi(1), "PNORM_SAM")
        call nf_write_sam(gamaz(1:nrf), "GAMAZ_SAM")

        !-----------------------------------------------------

        ! Check relative humidities at inputs in SAM
        do k = 1, nrf
            vap_pres = qv_sam(1,k) * pres(k) / (qv_sam(1,k) + 0.622 * (1.0 - qv_sam(1,k)))
            rh_sam(1,k) = vap_pres / esatw(tabs_sam(1,k))
        end do
        call nf_write_sam(rh_sam(1,:), "RH_SAM_IN")

#endif
        !-----------------------------------------------------

        ! Convert CAM Moistures and tabs to SAM q and t
        call  CAM_var_conversion(qv_sam, qc_sam, qi_sam, r_sam, tabs_sam, t_sam)

        ! Store variables on SAM grid at the start of the timestep
        ! for calculating tendencies later
        qv0_sam = qv_sam
        qc0_sam = qc_sam
        qi0_sam = qi_sam
        tabs0_sam = tabs_sam

#ifdef CAM_PROFILE
        call nf_write_sam(r_sam(1,:), "R_SAM_IN")
        call nf_write_sam(t_sam(1,:), "T_SAM_IN")
#endif
        !-----------------------------------------------------

        tabs0_sam = tabs_sam
        r0_sam = r_sam
        ! Run the neural net parameterisation
        ! Updates q and t with delta values
        ! advective, autoconversion (dt = -dq*(latent_heat/cp)),
        ! sedimentation (dt = -dq*(latent_heat/cp)),
        ! radiation rest tendency (multiply by dtn to get dt)
        call nn_convection_flux(tabs0_sam(:,1:nrf), r0_sam(:,1:nrf), y_in, &
                                tabs_sam(:,1:nrf), &
                                t_sam(:,1:nrf), r_sam(:,1:nrf), &
                                rho, adz, dz, dtn, &
                                precsfc_i)
        ! Update precsfc with prec from this timestep
        precsfc = precsfc + precsfc_i

#ifdef CAM_PROFILE
        call nf_write_sam(t_sam(1,:), "T_SAM_OUT")
        call nf_write_sam(r_sam(1,:), "R_SAM_OUT")
        call nf_write_scalar(precsfc_i(1), "PREC_SAM_OUT")
#endif
        !-----------------------------------------------------

        ! Formulate the output variables to CAM as required.
        call SAM_var_conversion(t_sam, r_sam, tabs_sam, qv_sam, qc_sam, qi_sam)
        ! Convert precipitation from kg/m^2 to m by dividing by density (1000)
        precsfc = precsfc * 1.0D-3

#ifdef CAM_PROFILE
        call nf_write_sam(tabs_sam(1,:), "TABS_SAM_OUT")
        call nf_write_sam(qi_sam(1,:), "QI_SAM_OUT")
        call nf_write_sam(qc_sam(1,:), "QC_SAM_OUT")
        call nf_write_sam(qv_sam(1,:), "QV_SAM_OUT")
#endif
        !-----------------------------------------------------

        ! Convert back into CAM variable tendencies (diff div by dtn) on SAM grid
        ! q into components and tabs to s (mult by cp)
        dqv_sam = (qv_sam - qv0_sam) / dtn
        dqc_sam = (qc_sam - qc0_sam) / dtn
        dqi_sam = (qi_sam - qi0_sam) / dtn
        ds_sam = cp_cam * (tabs_sam - tabs0_sam) / dtn
        ! Convert precipitation from total in m to date in m / s by div by dtn
        precsfc = precsfc / dtn

#ifdef CAM_PROFILE
        call nf_write_sam(ds_sam(1,:), "DS_SAM_OUT")
        call nf_write_sam(dqi_sam(1,:), "DQI_SAM_OUT")
        call nf_write_sam(dqc_sam(1,:), "DQC_SAM_OUT")
        call nf_write_sam(dqv_sam(1,:), "DQV_SAM_OUT")
#endif
        !-----------------------------------------------------

        ! Interpolate SAM variables to the CAM pressure levels
        ! setting tendencies above the SAM grid to 0.0
        ! TODO Interpolate all variables in one call
        call interp_to_cam(pres_cam, pres_int_cam, pres_sfc_cam, dqv_sam, dqv)
        call interp_to_cam(pres_cam, pres_int_cam, pres_sfc_cam, dqc_sam, dqc)
        call interp_to_cam(pres_cam, pres_int_cam, pres_sfc_cam, dqi_sam, dqi)
        call interp_to_cam(pres_cam, pres_int_cam, pres_sfc_cam, ds_sam, ds)

#ifdef CAM_PROFILE
        call nf_write_cam(ds(1,:), "DS_CAM_OUT")
        call nf_write_cam(dqi(1,:), "DQI_CAM_OUT")
        call nf_write_cam(dqc(1,:), "DQC_CAM_OUT")
        call nf_write_cam(dqv(1,:), "DQV_CAM_OUT")
#endif
    end subroutine nn_convection_flux_CAM


    subroutine nn_convection_flux_CAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()
        call sam_sounding_finalize()

    end subroutine nn_convection_flux_CAM_finalize


    !-----------------------------------------------------------------
    subroutine calculate_dry_mixing_ratio(qv, qc, qi, r)
        !! Calculate dry mixing ratio as required by SAM from moist variables as
        !! provided by CAM
        real(dp), intent(in) :: qv, qc, qi
        real(dp), intent(out) :: r
        real(dp) :: rv, rc, ri
            rv = qv / (1. - qv)
            rc = qc * (1.+rv)
            ri = qi * (1.+rv)
            r = rv + rc + ri
    end subroutine

    !-----------------------------------------------------------------
    subroutine calculate_moist_mixing_ratio(rv, rn, omegan, qc, qi, qv)
        !! Code for calculating qc and qi from rn.
        !! Adapted from statistics.f90 in SAM assuming dokruegermicro=.false.
        !! Also implements conversion from SAM dry mixing ratios to CAM moist mixing
        !! ratios
        real(dp), intent(in) :: rv, rn, omegan
        real(dp), intent(out) :: qc, qi, qv
        qc = rn*omegan/(1.+rv)
        qi = rn*(1.-omegan)/(1.+rv)
        qv = rv/(1.+rv)
    end subroutine

    !-----------------------------------------------------------------
    ! Private Subroutines

    subroutine interp_to_sam(p_cam, p_surf_cam, var_cam, var_sam, var_cam_surface)
        !! Interpolate from the CAM pressure grid to the SAM pressure grid.
        !! Uses linear interpolation between nearest grid points on domain (CAM) grid.

        !! The CAM and SAM grids have pressures measured at cell centres which is where
        !! the variables that we are interpolating are stored.

        !! It can be expected, but not guaranteed, that the domain (CAM) grid is
        !! coarser than the target (SAM) grid.

        !! NOTE: When reading this code remember that pressure DECREASES monotonically
        !!       with altitude, so comparisons are not intuitive: Pa > Pb => a is lower altitude than b

        != unit hPa :: p_cam, p_surf_cam
        real(dp), dimension(:,:), intent(in) :: p_cam
        real(dp), dimension(:), intent(in) :: p_surf_cam
            !! pressure [Pa] from the CAM model
        != unit 1 :: var_cam, var_sam
        real(dp), dimension(:,:), intent(in) :: var_cam
            !! variable from CAM grid to be interpolated from
        real(dp), dimension(:), intent(in) :: var_cam_surface
            !! Surface value of CAM variable for interpolation/boundary condition
        real(dp), dimension(:,:), intent(out) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real(dp), dimension(nrf) :: p_mid_sam, p_norm_sam
        real(dp), dimension(:,:), allocatable :: p_norm_cam
            !! Normalised pressures (using surface pressure) as a sigma-coordinate

        integer :: i, k, nz_cam, ncol_cam

        integer :: kc_u, kc_l
        real(dp) :: pc_u, pc_l
            !! Pressures of the CAM cell above and below the target SAM cell

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(ncol_cam, nz_cam))

        ! Normalise pressures by surface value
        ! Note - We work in pressure space but use a pressure 'midpoint' based on altitude from
        !        SAM as this is concurrent to the variable value at this location.
        p_norm_sam(:) = pres(:) / presi(1)
        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_surf_cam(:)
        end do

        ! Loop over columns and SAM levels and interpolate each column
        do k = 1, nrf
        do i = 1, ncol_cam

            ! If SAM cell centre is below lowest CAM centre - interpolate to surface value
            if (p_norm_sam(k) > p_norm_cam(i, 1)) then
                var_sam(i, k) = var_cam_surface(i) &
                                + (p_norm_sam(k)-1.0) &
                                * (var_cam(i, 1)-var_cam_surface(i))/(p_norm_cam(i, 1)-1.0)

            ! Check CAM grid top cell is above SAM grid top
            elseif (p_norm_cam(i, nz_cam) > p_norm_sam(nrf)) then
                ! This should not happen as CAM grid extends to higher altitudes than SAM
                ! TODO This check is run on every iteration - move outside
                 write(*,*) "CAM upper pressure level is lower than that of SAM: Stopping."
                stop

            ! Locate the neighbouring CAM indices to interpolate between
            else
                ! TODO - this will be slow - speed up later
                kc_u = 1
                pc_u = p_norm_cam(i, kc_u)
                do while (pc_u > p_norm_sam(k))
                    kc_u = kc_u + 1
                    pc_u = p_norm_cam(i, kc_u)
                end do
                kc_l = kc_u - 1
                pc_l = p_norm_cam(i, kc_l)

                ! interpolate variables - Repeat for as many variables as need interpolating
                var_sam(i, k) = var_cam(i, kc_l) &
                                + (p_norm_sam(k)-pc_l) &
                                * (var_cam(i, kc_u)-var_cam(i, kc_l))/(pc_u-pc_l)
            endif
        end do
        end do

        ! Clean up
        deallocate(p_norm_cam)

    end subroutine interp_to_sam

    subroutine interp_to_cam(p_cam, p_int_cam, p_surf_cam, var_sam, var_cam)
      ! TODO This would seriously benefit from a refactor to add the top interface to p_sam_int and p_cam_int
        !! Interpolate from the SAM NN pressure grid to the CAM pressure grid.
        !! Uses conservative regridding between grids.

        !! Requires dealing in interfaces rather than centres.
        !! These  start at the surface of both grids
        !! and extend to grid-top for CAM, but stop at the base of the top cell
        !! (i.e. no grid-top value) for SAM.

        != unit Pa :: p_cam, p_int_cam
        real(dp), dimension(:,:), intent(in) :: p_cam
            !! pressure [Pa] from the CAM model
        real(dp), dimension(:,:), intent(in) :: p_int_cam
            !! interface pressures [Pa] for each column in CAM
        real(dp), dimension(:), intent(in) :: p_surf_cam
            !! surface pressures [Pa] for each column in CAM
        != unit 1 :: var_cam, var_sam
        real(dp), dimension(:,:), intent(out) :: var_cam
            !! variable from CAM grid to be interpolated from
        real(dp), dimension(:,:), intent(in) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real(dp), dimension(nrf) :: p_norm_sam
        real(dp), dimension(nrf+1) :: p_int_norm_sam
        real(dp), dimension(:,:), allocatable :: p_norm_cam
        real(dp), dimension(:,:), allocatable :: p_int_norm_cam

        integer :: i, k, c, nz_cam, ncol_cam, ks, ks_u, ks_l
        real(dp) :: ps_u, ps_l, pc_u, pc_l
        real(dp) :: p_norm_sam_min

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(ncol_cam, nz_cam))
        allocate(p_int_norm_cam(ncol_cam, nz_cam+1))

        ! Normalise pressures by surface value (extending SAM to add top interface)
        p_norm_sam(:) = pres(1:nrf) / presi(1)
        p_int_norm_sam(1:nrf+1) = presi(1:nrf+1) / presi(1)
        p_int_norm_sam(1) = 1.0

        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_surf_cam(:)
            p_int_norm_cam(:,k) = p_int_cam(:,k) / p_surf_cam(:)
        end do
        p_int_norm_cam(:,nz_cam+1) = p_int_cam(:,nz_cam+1) / p_surf_cam(:)
        p_int_norm_cam(:,1) = 1.0

        ! Loop over columns and SAM levels and interpolate each column
        ! TODO Revert indices to i,k for speed
        do k = 1, nz_cam
        do i = 1, ncol_cam
            ! Check we are within array overlap:
            ! No need to compare bottom end as pressure is normalised to 1.0 for both.

            if (p_int_norm_cam(i, k) < p_int_norm_sam(nrf+1)) then
                ! Anything above the SAM NN grid should be set to 0.0
!                 write(*,*) "CAM pressure is lower than SAM top: set tend to 0.0."
                var_cam(i, k) = 0.0

            else
                ! Get the pressures at the top and bottom of the CAM cell
                pc_u = p_int_norm_cam(i, k+1)
                pc_l = p_int_norm_cam(i,k)

                ! Locate the SAM indices to interpolate between
                ! This will be slow - speed up later
                ! ks_u is the SAM block with an upper interface above the CAM block
                ks_u = 1
                ps_u = p_int_norm_sam(ks_u+1)
                do while ((ps_u > pc_u) .and. (ks_u < nrf))
                    ks_u = ks_u + 1
                    ps_u = p_int_norm_sam(ks_u+1)
                end do
                ! ks_l is the SAM block with a lower interface below the CAM block
                ks_l = ks_u
                ps_l = p_int_norm_sam(ks_l)
                do while ((ps_l < pc_l) .and. (ks_l > 1))
                    ks_l = ks_l - 1
                    ps_l = p_int_norm_sam(ks_l)
                end do

                ! Combine all SAM cells that this CAM_CELL(i, k) overlaps
                var_cam(i,k) = 0.0
                if (ks_u == ks_l) then
                    ! CAM cell fully enclosed by SAM cell => match 'density' of SAM cell
                    var_cam(i,k) = var_cam(i,k) + var_sam(i,ks_u)
                else
                  do c = ks_l, ks_u
                    if (c == ks_l) then
                      ! Bottom cell
                      var_cam(i,k) = var_cam(i,k) + (pc_l-p_int_norm_sam(ks_l+1))*var_sam(i,ks_l) / (pc_l - pc_u)
                    elseif (c == ks_u) then
                      ! Top cell
                      var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(ks_u)-pc_u)*var_sam(i, ks_u) / (pc_l - pc_u)
                    else
                      ! Intermediate cell (SAM Cell fully enclosed by CAM Cell => Absorb all)
                      var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(c)-p_int_norm_sam(c+1))*var_sam(i, c) / (pc_l - pc_u)

                    endif

                  enddo
                endif
            endif
        end do
        end do

        ! Clean up
        deallocate(p_norm_cam)
        deallocate(p_int_norm_cam)

    end subroutine interp_to_cam

    subroutine extrapolate_to_surface(p_cam, var_cam, p_surf_cam, var_surf_cam)
        !! Extrapolate variable profile to the surface, required for interpolation.

        != unit Pa :: p_cam
        real(dp), dimension(:,:), intent(in) :: p_cam
            !! pressure [Pa] from the CAM model
        real(dp), dimension(:), intent(in) :: p_surf_cam
            !! surface pressures [Pa] for each column in CAM
        != unit 1 :: var_cam
        real(dp), dimension(:,:), intent(in) :: var_cam
            !! variable from CAM grid to be interpolated from
        real(dp), dimension(:), intent(out) :: var_surf_cam
            !! variable from CAM grid to be interpolated from

        integer :: i, k
        real(dp) :: gdt(size(var_surf_cam))

        gdt(:) = (var_cam(:,2) - var_cam(:,1)) / (p_cam(:,2) - p_cam(:,1))
        var_surf_cam(:) = var_cam(:,1) + gdt(:) * (p_surf_cam(:) - p_cam(:,1))

    end subroutine extrapolate_to_surface

    subroutine sam_sounding_init(filename)
        !! Read various profiles in from SAM sounding file

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, k
        integer :: z_dimid, dz_dimid
        integer :: z_varid, dz_varid, pres_varid, presi_varid, rho_varid, adz_varid

        character(len=136), intent(in) :: filename
            !! NetCDF filename from which to read data


        !-------------allocate arrays and read data-------------------

        ! Open the file. NF90_NOWRITE tells netCDF we want read-only
        ! access
        ! Get the varid or dimid for each variable or dimension based
        ! on its name.

        call check( nf90_open(trim(filename),NF90_NOWRITE,ncid ))

        call check( nf90_inq_dimid(ncid, 'z', z_dimid))
        call check( nf90_inquire_dimension(ncid, z_dimid, len=nz_sam))
        ! call check( nf90_inq_dimid(ncid, 'dz', dz_dimid))
        ! call check( nf90_inquire_dimension(ncid, dz_dimid, len=ndz))

        ! Note that nz of sounding may be longer than nrf
        allocate(z(nz_sam))
        allocate(pres(nz_sam))
        allocate(presi(nz_sam))
        allocate(rho(nz_sam))
        allocate(adz(nz_sam))
        allocate(gamaz(nz_sam))

        ! Read data in from nc file - convert pressures to hPa
        call check( nf90_inq_varid(ncid, "z", z_varid))
        call check( nf90_get_var(ncid, z_varid, z))
        call check( nf90_inq_varid(ncid, "pressure", pres_varid))
        call check( nf90_get_var(ncid, pres_varid, pres))
        pres(:) = pres / 100.0
        call check( nf90_inq_varid(ncid, "interface_pressure", presi_varid))
        call check( nf90_get_var(ncid, presi_varid, presi))
        presi(:) = presi(:) / 100.0
        call check( nf90_inq_varid(ncid, "rho", rho_varid))
        call check( nf90_get_var(ncid, rho_varid, rho))
        call check( nf90_inq_varid(ncid, "adz", adz_varid))
        call check( nf90_get_var(ncid, adz_varid, adz))
        call check( nf90_inq_varid(ncid, "dz", dz_varid))
        call check( nf90_get_var(ncid, dz_varid, dz))

        ! Close the nc file
        call check( nf90_close(ncid))

        ! Calculate gamaz required elsewhere
        do k = 1, nz_sam
            gamaz(k) = ggr/cp*z(k)
        end do

        write(*,*) 'Finished reading SAM sounding file.'

    end subroutine sam_sounding_init


    subroutine sam_sounding_finalize()
        !! Deallocate module variables read from sounding

        deallocate(z, pres, presi, rho, adz, gamaz)

    end subroutine sam_sounding_finalize


    subroutine CAM_var_conversion(qv, qc, qi, r, tabs, t)
        !! Convert CAM moist mixing ratios for species qv, qc, qi to
        !! dry mixing ratio r to used by SAM parameterisation
        !! q is total water qv/c/i is cloud vapor/liquid/ice
        !! Convert CAM absolute temperature to moist static energy t used by SAM
        !! by inverting lines 49-53 of diagnose.f90 from SAM where

        !! WARNING: This routine uses gamaz(k) which is defined on the SAM grid.
        !!          Using a grid that does not match the SAM grid for input/output
        !!          data will produce unphysical values.

        integer :: nx, nz
            !! array sizes
        integer :: i, k
            !! Counters
        real(dp) :: omn
            !! intermediate omn factor used in variable conversion
        real(dp) :: rv, rc, ri
            !! intermediate r values to hold dry mixing ratios of vapour, cld, ice

        ! ---------------------
        ! Fields from CAM/SAM
        ! ---------------------
        != unit 1 :: r
        real(dp), intent(out) :: r(:, :)
            !! Total non-precipitating water mixing ratio as required by SAM NN
        real(dp), intent(in) :: qv(:, :)
            !! Cloud water vapour from CAM
        real(dp), intent(in) :: qc(:, :)
            !! Cloud water from CAM
        real(dp), intent(in) :: qi(:, :)
            !! Cloud ice from CAM

        real(dp), intent(out) :: t(:, :)
            !! Static energy as required by SAM NN
        real(dp), intent(in) :: tabs(:, :)
            !! Absolute temperature from CAM


        nx = size(r, 1)
        nz = size(r, 2)

        do k = 1, nz
          do i = 1, nx
            ! Calculate dry mixing ratio from moist variables
            call calculate_dry_mixing_ratio(qv(i,k), qc(i,k), qi(i,k), r(i,k))

            ! omp  = max(0.,min(1.,(tabs(i,k)-tprmin)*a_pr))  ! There is no qp in CAM
            omn  = max(0.,min(1.,(tabs(i,k)-tbgmin)*a_bg))
            t(i,k) = tabs(i,k) &
                   ! - (fac_cond+(1.-omp)*fac_fus)*qp(i,k) &  ! There is no qp in CAM
                     - (fac_cond+(1.-omn)*fac_fus)*(rc + ri) &
                     + gamaz(k)
          end do
        end do

    end subroutine CAM_var_conversion


    subroutine SAM_var_conversion(t, r, tabs, qv, qc, qi)
        !! Convert SAM t and r to tabs, qv, qc, qi used by CAM
        !! t is normalised liquid ice static energy, r is a dry mixing ratio
        !! tabs is absolute temperature, q is cloud vapor/liquid/ice,

        integer :: nx, nz
            !! array sizes
        integer :: i, k, niter
            !! Counters

        ! ---------------------
        ! Fields from SAM/CAM
        ! ---------------------
        != unit K :: tabs, tabs1
        real(dp), intent(inout) :: tabs(:, :)
            !! absolute temperature
        real(dp) :: tabs1
            !! Temporary variable for tabs

        != unit kg/kg :: q, qv, qc, qi
        real(dp), intent(in) :: r(:, :)
            !! Total non-precipitating water mixing ratio from SAM
        real(dp), intent(out) :: qv(:, :)
            !! Cloud water vapour in CAM
        real(dp), intent(out) :: qc(:, :)
            !! Cloud water (liquid) in CAM
        real(dp), intent(out) :: qi(:, :)
            !! Cloud ice in CAM

        != K :: t
        real(dp), intent(in) :: t(:, :)
            !! normalised liquid ice static energy

        ! Intermediate variables
        real(dp) :: rsat, om, omn, dtabs, drsat, lstarn, dlstarn, fff, dfff, rv, rn, r_temp

        nx = size(tabs, 1)
        nz = size(tabs, 2)

        ! Loop over all cells and iterate until equilibrium of condensed species is achieved.
        ! This code is adapted from cloud.f90 in SAM
        do k = 1, nz
        do i = 1, nx

            ! Enforce r >= 0.0
            r_temp=max(0.,r(i,k))

            ! Initial guess for temperature assuming no cloud water/ice:
            tabs(i,k) = t(i,k)-gamaz(k)
            tabs1=tabs(i,k)

            ! Set saturation vapour based on cloud type
            ! Warm cloud:
            if(tabs1.ge.tbgmax) then
                rsat = rsatw(tabs1,pres(k))
            ! Ice cloud:
            elseif(tabs1.le.tbgmin) then
                rsat = rsati(tabs1,pres(k))
            ! Mixed-phase cloud:
            else
                om = an*tabs1-bn
                rsat = om*rsatw(tabs1,pres(k))+(1.-om)*rsati(tabs1,pres(k))
            endif

            ! Test if condensation is possible (humidity is above saturation) and iterate:
            if(r_temp .gt. rsat) then
                niter=0
                dtabs = 100.
                do while(abs(dtabs).gt.0.01.and.niter.lt.10)
                ! Warm cloud regime
                    if(tabs1.ge.tbgmax) then
                        om=1.
                        lstarn=fac_cond
                        dlstarn=0.
                        rsat=rsatw(tabs1,pres(k))
                        drsat=dtrsatw(tabs1,pres(k))
                    ! Ice cloud regime
                    else if(tabs1.le.tbgmin) then
                        om=0.
                        lstarn=fac_sub
                        dlstarn=0.
                        rsat=rsati(tabs1,pres(k))
                        drsat=dtrsati(tabs1,pres(k))
                    ! Mixed cloud regime
                    else
                        om=an*tabs1-bn
                        lstarn=fac_cond+(1.-om)*fac_fus
                        dlstarn=an
                        rsat=om*rsatw(tabs1,pres(k))+(1.-om)*rsati(tabs1,pres(k))
                        drsat=om*dtrsatw(tabs1,pres(k))+(1.-om)*dtrsati(tabs1,pres(k))
                    endif

                    ! Update thermodynamics and check for convergence
                    fff = tabs(i,k)-tabs1+lstarn*(r_temp-rsat)
                    dfff=dlstarn*(r_temp-rsat)-lstarn*drsat-1.
                    dtabs=-fff/dfff
                    niter=niter+1
                    tabs1=tabs1+dtabs
               end do

               ! Update saturation point and then calculate water and residual vapour content
               rsat = rsat + drsat * dtabs
               rn = max(0., r_temp-rsat)
               rv = max(0., r_temp-rn)

            ! If condensation not possible rn is 0.0
            else
              rn = 0.
              rv = r_temp

            endif

            ! Set tabs to iterated tabs after convection
            tabs(i,k) = tabs1

            ! Calculating qc and qi from rn, converting from SAM dry
            ! mixing ratios to CAM moist mixing ratios
            omn = omegan(tabs(i,k))
            call calculate_moist_mixing_ratio(rv, rn, omn, qc(i,k), qi(i, k), qv(i, k))

        end do
        end do

    end subroutine SAM_var_conversion

end module nn_interface_CAM
