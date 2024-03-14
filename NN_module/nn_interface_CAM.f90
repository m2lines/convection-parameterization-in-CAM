module nn_interface_CAM
    !! Interface to convection parameterisation for the CAM model
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

    ! TODO: Check if CAM needs to allocate MPI Subdomains in some way
    
    ! TODO: Check for redundant variables once refactored (icycle, nstatis, precip etc.)

    ! TODO: Remove commented lines for variables not required from SAM Sounding
    ! TODO: Check precisions when defined in python, saved as netcdf, and read to Fortran

!---------------------------------------------------------------------
! Libraries to use
use netcdf
use nn_convection_flux_mod, only: nn_convection_flux, &
                                  nn_convection_flux_init, nn_convection_flux_finalize, &
                                  esati, qsati, qsatw, dtqsatw, dtqsati
use SAM_consts_mod, only: nrf, ggr, cp, tbgmax, tbgmin, tprmax, tprmin, &
                          fac_cond, fac_sub, fac_fus, &
                          a_bg, a_pr, an, bn, ap, bp, &
                          omegan, check
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux_CAM, &
        nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize

! TODO Make routines public for purposes of testing
public interp_to_sam, interp_to_cam, fetch_sam_data
public SAM_var_conversion, CAM_var_conversion

!---------------------------------------------------------------------
! local/private data

!= unit 1 :: nz_sam
integer :: nz_sam
    !! number of vertical values in the SAM sounding profiles
!= unit m :: z
!= unit hPa :: pres, presi
!= unit kg m-3 :: rho
!= unit 1 :: adz
real, allocatable, dimension(:) :: z, pres, presi, pres_mid, rho, adz, gamaz
    !! SAM sounding variables
!= unit m :: dz
real :: dz
    !! grid spacing in z direction for the lowest grid layer

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine fetch_sam_data(pressam, presisam, gamazsam, rhosam, zsam)
        !! Temporary subroutine to extract SAM pressure profile
        real, dimension(48), intent(out) :: pressam, presisam, gamazsam, rhosam, zsam
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

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights
        character(len=1024), intent(in) :: sounding_filename
            !! NetCDF filename from which to read SAM sounding and grid data

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename)

        ! Read in the SAM sounding and grid data required
        call sam_sounding_init(sounding_filename)

    end subroutine nn_convection_flux_CAM_init


    subroutine nn_convection_flux_CAM(pres_cam, pres_int_cam, &
                                      tabs_cam, qv_cam, qc_cam, qi_cam, &
                                      cp_cam, &
                                      dtn, &
                                      nx, nz, &
                                      nstep, nstatis, icycle, &
                                      precsfc, &
                                      dqi, dqv, dqc, ds)
        !! Interface to the nn_convection parameterisation for the CAM model

        integer :: k

        real, intent(in) :: cp_cam
            !! specific heat capacity of dry air from CAM
        real, dimension(:,:) :: pres_cam
            !! pressure [hPa] from the CAM model
        real, dimension(:,:) :: pres_int_cam
            !! interface pressure [hPa] from the CAM model (element 1 is surface pressure)
        real, dimension(:,:) :: tabs_cam
            !! absolute temperature [K] from the CAM model
        real, dimension(:,:) :: qv_cam, qc_cam, qi_cam
            !! moisture content [-] from the CAM model
        real, dimension(:) :: precsfc
        real, intent(in) :: dtn
        integer, intent(in) :: nx, nz, nstep, nstatis, icycle

        real :: y_in(nx)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        real, dimension(nx)      :: prec_sed
            !! Sedimenting precipitation at surface

        ! Variables on the SAM grid
        real, dimension(nx, nrf) :: q0_sam, tabs0_sam
        real, dimension(nx, nrf) :: q_sam, t_sam, tabs_sam
        real, dimension(nx, nrf) :: qi_sam, qv_sam, qc_sam
        real, dimension(nx, nrf) :: qi0_sam, qv0_sam, qc0_sam
        real, dimension(nx, nrf) :: dqi_sam, dqv_sam, dqc_sam, ds_sam

        real, dimension(nx, nz), intent(out) :: dqi, dqv, dqc, ds

        real, dimension(nx) :: qi_surf, qv_surf, qc_surf, tabs_surf

        ! TODO: CAM requires surface precipitation
        ! Initialise precipitation to 0 if required and at start of cycle if subcycling
        if(mod(nstep-1,nstatis).eq.0 .and. icycle.eq.1) then
            precsfc(:)=0.
        end if

        ! distance to the equator
        ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
        ! TODO: Set y_in as appropriate for CAM
        y_in(:) = 0.0

        !-----------------------------------------------------
        
        ! Interpolate CAM variables to the SAM pressure levels
        ! TODO Interpolate all variables in one call
        ! TODO Check Boundary Condition
        qi_surf = 0.0
        qc_surf = 0.0
        qv_surf = 0.0
        call interp_to_sam(pres_cam, pres_int_cam(:,1), &
                           tabs_cam, tabs_sam, tabs_surf)
        call interp_to_sam(pres_cam, pres_int_cam(:,1), &
                           qi_cam, qi_sam, qi_surf)
        call interp_to_sam(pres_cam, pres_int_cam(:,1), &
                           qc_cam, qc_sam, qc_surf)
        call interp_to_sam(pres_cam, pres_int_cam(:,1), &
                           qv_cam, qv_sam, qv_surf)

        !-----------------------------------------------------
        
        ! Convert CAM Moistures and tabs to SAM q and t
        call  CAM_var_conversion(qv_sam, qc_sam, qi_sam, q_sam, t_sam, tabs_sam)

        ! Store variables on SAM grid at the start of the timestep
        ! for calculating tendencies later
        qv0_sam = qv_sam
        qc0_sam = qc_sam
        qi0_sam = qi_sam
        tabs0_sam = tabs_sam

        !-----------------------------------------------------

        tabs0_sam = tabs_sam
        q0_sam = q_sam
        ! Run the neural net parameterisation
        ! Updates q and t with delta values
        ! advective, autoconversion (dt = -dq*(latent_heat/cp)),
        ! sedimentation (dt = -dq*(latent_heat/cp)),
        ! radiation rest tendency (multiply by dtn to get dt)
        call nn_convection_flux(tabs0_sam(:,1:nrf), q0_sam(:,1:nrf), y_in, &
                                tabs_sam(:,1:nrf), &
                                t_sam(:,1:nrf), q_sam(:,1:nrf), &
                                rho, adz, dz, dtn, &
                                prec_sed)

        ! TODO: Update precipitation if required

        !-----------------------------------------------------
        
        ! Formulate the output variables to CAM as required.
        call SAM_var_conversion(t_sam, q_sam, tabs_sam, qv_sam, qc_sam, qi_sam)

        !-----------------------------------------------------

        ! Convert back into CAM variable tendencies (diff div by dtn) on SAM grid
        ! q into components and tabs to s (mult by cp)
        dqv_sam = (qv_sam - qv0_sam) / dtn
        dqc_sam = (qc_sam - qc0_sam) / dtn
        dqi_sam = (qi_sam - qi0_sam) / dtn
        ds_sam = cp_cam * (tabs_sam - tabs0_sam) / dtn

        !-----------------------------------------------------

        ! Interpolate SAM variables to the CAM pressure levels
        ! setting tendencies above the SAM grid to 0.0
        ! TODO Interpolate all variables in one call
        call interp_to_cam(pres_cam, pres_int_cam, dqv_sam, dqv)
        call interp_to_cam(pres_cam, pres_int_cam, dqc_sam, dqc)
        call interp_to_cam(pres_cam, pres_int_cam, dqi_sam, dqi)
        call interp_to_cam(pres_cam, pres_int_cam, ds_sam, ds)

    end subroutine nn_convection_flux_CAM


    subroutine nn_convection_flux_CAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()
        call sam_sounding_finalize()

    end subroutine nn_convection_flux_CAM_finalize


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
        !!       with altitude, so comparisons are not intuitive: Pa > Pb => a is below b

        != unit hPa :: p_cam, p_surf_cam
        real, dimension(:,:), intent(in) :: p_cam
        real, dimension(:), intent(in) :: p_surf_cam
            !! pressure [hPa] from the CAM model
            !! p_surf_cam is from pint(1)
        != unit 1 :: var_cam, var_sam
        real, dimension(:,:), intent(in) :: var_cam
            !! variable from CAM grid to be interpolated from
        real, dimension(:), intent(in) :: var_cam_surface
            !! Surface value of CAM variable for interpolation/boundary condition
        real, dimension(:,:), intent(out) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real, dimension(nrf) :: p_mid_sam, p_norm_sam
        real, dimension(:,:), allocatable :: p_norm_cam
            !! Normalised pressures (using surface pressure) as a sigma-coordinate

        integer :: i, k, nz_cam, ncol_cam

        integer :: kc_u, kc_l
        real :: pc_u, pc_l
            !! Pressures of the CAM cell above and below the target SAM cell

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(ncol_cam, nz_cam))

        ! Normalise pressures by surface value
        ! Note - We work in pressure space but use a 'midpoint' based on altitude from
        !        SAM as this is concurrent to the variable value at this location.
        p_norm_sam(:) = pres(:) / presi(1)
        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_surf_cam(:)
        end do

        ! Loop over columns and SAM levels and interpolate each column
        ! TODO - this can be made more efficient if the CAM grid is the same for all columns
        !        as currently interpolation is run for every cell.
        do k = 1, nrf
        do i = 1, ncol_cam
            ! Check we are within array

            ! If SAM cell centre is below lowest CAM centre - interpolate to surface value
            if (p_norm_sam(k) > p_norm_cam(i, 1)) then
                ! TODO Remove warning and interpolate to boundary
                ! TODO Add boundary conditions for any input variables as required.
!                write(*,*) "Interpolating to surface."
                var_sam(i, k) = var_cam_surface(i) &
                                + (p_norm_sam(k)-1.0) &
                                * (var_cam(i, 1)-var_cam_surface(i))/(p_norm_cam(i, 1)-1.0)
            ! Check CAM grid top cell is above SAM grid top
            elseif (p_norm_cam(i, nz_cam) > p_norm_sam(nrf)) then
                ! This should not happen as CAM grid extends to higher altitudes than SAM
                ! TODO This check is run on every iteration - move outside
                write(*,*) "CAM upper pressure level is lower than that of SAM: Stopping."
                stop
            else
                ! Locate the neighbouring CAM indices to interpolate between
                ! TODO - this will be slow - speed up later
                kc_u = 1
                pc_u = p_norm_cam(i, kc_u)
                do while (pc_u > p_norm_sam(k))
                    kc_u = kc_u + 1
                    pc_u = p_norm_cam(i, kc_u)
                end do
                kc_l = kc_u - 1
                pc_l = p_norm_cam(i, kc_l)
                ! Redundant following the lines above
                ! do while (pc_l < p_norm_sam(k))
                !     kc_l = kc_l - 1
                !     pc_l = p_norm_cam(i, kc_l)
                ! end do

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

    subroutine interp_to_cam(p_cam, p_int_cam, var_sam, var_cam)
      ! TODO This would seriously benefit from a refactor to add the top interface to p_sam_int and p_cam_int
        !! Interpolate from the SAM NN pressure grid to the CAM pressure grid.
        !! Uses conservative regridding between grids.

        !! Requires dealing in interfaces rather than centres.
        !! These  start at the surface of both grids
        !! and extend to grid-top for CAM, but stop at the base of the top cell
        !! (i.e. no grid-top value) for SAM.

        != unit hPa :: p_cam, p_int_cam
        real, dimension(:,:), intent(in) :: p_cam
            !! pressure [hPa] from the CAM model
        real, dimension(:,:), intent(in) :: p_int_cam
            !! interface pressures [hPa] for each column in CAM
        != unit 1 :: var_cam, var_sam
        real, dimension(:,:), intent(out) :: var_cam
            !! variable from CAM grid to be interpolated from
        real, dimension(:,:), intent(in) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real, dimension(nrf) :: p_norm_sam
        real, dimension(nrf+1) :: p_int_norm_sam
        real, dimension(:,:), allocatable :: p_norm_cam
        real, dimension(:,:), allocatable :: p_int_norm_cam

        integer :: i, k, c, nz_cam, ncol_cam, ks, ks_u, ks_l
        real :: ps_u, ps_l, pc_u, pc_l
        real :: p_norm_sam_min

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(ncol_cam, nz_cam))
        allocate(p_int_norm_cam(ncol_cam, nz_cam+1))

        ! Normalise pressures by surface value (extending SAM to add top interface)
        p_norm_sam(:) = pres_mid(1:nrf) / presi(1)
        p_int_norm_sam(1:nrf+1) = presi(1:nrf+1) / presi(1)
        p_int_norm_sam(1) = 1.0

        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_int_cam(:,1)
            p_int_norm_cam(:,k) = p_int_cam(:,k) / p_int_cam(:,1)
        end do
        p_int_norm_cam(:,nz_cam+1) = p_int_cam(:,nz_cam+1) / p_int_cam(:,1)
        p_int_norm_cam(:,1) = 1.0

        ! Loop over columns and SAM levels and interpolate each column
        ! TODO Revert indices to i,k for speed
        do k = 1, nz_cam
        do i = 1, ncol_cam
            ! Check we are within array overlap:
            ! No need to compare bottom end as pressure is normalised to 1.0 for both.

            if (p_int_norm_cam(i, k) < p_int_norm_sam(nrf+1)) then
                ! Anything above the SAM NN grid should be set to 0.0
                write(*,*) "CAM pressure is lower than SAM top: set tend to 0.0."
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


                ! write(*,*) ks_l, ks_u
                if (ks_u == ks_l) then
                    ! CAM cell fully enclosed by SAM cell => match 'density' of SAM cell
                    var_cam(i,k) = var_cam(i,k) + (pc_u - pc_l)*var_sam(i, ks_u)/(ps_u-ps_l)
                else
                  do c = ks_l, ks_u
                    if (c == ks_l) then
                      ! Bottom cell
                      var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(ks_l+1)-pc_l)*var_sam(i,ks_l)/(p_int_norm_sam(ks_l+1)-ps_l)
                    elseif (c == ks_u) then
                      ! Top cell
                      var_cam(i,k) = var_cam(i,k) + (pc_u - p_int_norm_sam(ks_u))*var_sam(i, ks_u)/(ps_u-p_int_norm_sam(ks_u))
                    else
                      ! Intermediate cell (SAM Cell fully enclosed by CAM Cell => Absorb all)
                      var_cam(i,k) = var_cam(i,k) + var_sam(i,c)

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

    subroutine sam_sounding_init(filename)
        !! Read various profiles in from SAM sounding file

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, k
        integer :: z_dimid, dz_dimid
        integer :: z_varid, dz_varid, pres_varid, presi_varid, rho_varid, adz_varid

        character(len=1024), intent(in) :: filename
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
        allocate(pres_mid(nz_sam))
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

        ! Calculate the pressure centre of the cell (as opposed to the altitude centre)
        do k = 1, nz_sam-1
            pres_mid(k) = (presi(k) + presi(k+1)) / 2.0
        end do
        pres_mid(nz_sam) = -9999.0

        write(*,*) 'Finished reading SAM sounding file.'

    end subroutine sam_sounding_init


    subroutine sam_sounding_finalize()
        !! Deallocate module variables read from sounding

        deallocate(z, pres, presi, pres_mid, rho, adz, gamaz)

    end subroutine sam_sounding_finalize


    subroutine CAM_var_conversion(qv, qc, qi, q, tabs, t)
        !! Convert CAM qv, qc, qi to q to used by SAM parameterisation
        !! q is total water qv/c/i is cloud vapor/liquid/ice
        !! Convert CAM absolute temperature to moist static energy t used by SAM
        !! by inverting lines 49-53 of diagnose.f90 from SAM where
        !! qn is cloud water + ice, qp is precipitable water (0 here)

        !! WARNING: This routine uses gamaz(k) which is defined on the SAM grid.
        !!          Using a grid that does not match the SAM grid for input/output
        !!          data will produce unphysical values.

        integer :: nx, nz
            !! array sizes
        integer :: i, k
            !! Counters
        real :: omn
            !! intermediate omn factor used in variable conversion

        ! ---------------------
        ! Fields from CAM/SAM
        ! ---------------------
        != unit 1 :: q, qv, qc, qi
        real, intent(out) :: q(:, :)
            !! Total non-precipitating water mixing ratio as required by SAM NN
        real, intent(in) :: qv(:, :)
            !! Cloud water vapour from CAM
        real, intent(in) :: qc(:, :)
            !! Cloud water from CAM
        real, intent(in) :: qi(:, :)
            !! Cloud ice from CAM

        real, intent(out) :: t(:, :)
            !! Static energy as required by SAM NN
        real, intent(in) :: tabs(:, :)
            !! Absolute temperature from CAM


        nx = size(q, 1)
        nz = size(q, 2)

        do k = 1, nz
          do i = 1, nx
            q(i,k) = qv(i,k) + qc(i,k) + qi(i,k)

            ! omp  = max(0.,min(1.,(tabs(i,k)-tprmin)*a_pr))
            omn  = max(0.,min(1.,(tabs(i,k)-tbgmin)*a_bg))
            t(i,k) = tabs(i,k) &
                     ! - (fac_cond+(1.-omp)*fac_fus)*qp(i,k) &
                     - (fac_cond+(1.-omn)*fac_fus)*(qc(i,k) + qi(i,k)) &
                     + gamaz(k)
          end do
        end do

    end subroutine CAM_var_conversion
    

    subroutine SAM_var_conversion(t, q, tabs, qv, qc, qi)
        !! Convert SAM t and q to tabs, qv, qc, qi used by CAM
        !! t is normalised liquid ice static energy, q is total water
        !! tabs is absolute temperature, q is cloud vapor/liquid/ice,

        integer :: nx, nz
            !! array sizes
        integer :: i, k, niter
            !! Counters

        ! ---------------------
        ! Fields from SAM/CAM
        ! ---------------------
        != unit K :: tabs, tabs1
        real, intent(inout) :: tabs(:, :)
            !! absolute temperature
        real, allocatable :: tabs1
            !! Temporary variable for tabs

        != unit 1 :: q, qn, qv, qc, qi
        real :: q(:, :)
            !! Total non-precipitating water mixing ratio from SAM
        real, intent(out) :: qv(:, :)
            !! Cloud water vapour
        real, intent(out) :: qc(:, :)
            !! Cloud water
        real, intent(out) :: qi(:, :)
            !! Cloud ice
        real :: qn
            !! Cloud liquid + cloud ice

        != unit  :: t
        real, intent(in) :: t(:, :)
            !! normalised liquid ice static energy

        real :: qsat, om, omn, dtabs, dqsat, lstarn, dlstarn, fff, dfff

        nx = size(tabs, 1)
        nz = size(tabs, 2)

        do k = 1, nz
        do i = 1, nx
        
            ! Enforce q >= 0.0
            q(i,k)=max(0.,q(i,k))
        
        
            ! Initial guess for temperature assuming no cloud water/ice:
            tabs(i,k) = t(i,k)-gamaz(k)
            tabs1=tabs(i,k)
        
            ! Warm cloud:
            if(tabs1.ge.tbgmax) then
                qsat = qsatw(tabs1,pres(k))
        
            ! Ice cloud:
            elseif(tabs1.le.tbgmin) then
                qsat = qsati(tabs1,pres(k))
        
            ! Mixed-phase cloud:
            else
                om = an*tabs1-bn
                qsat = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))

            endif

            !  Test if condensation is possible and iterate:
            if(q(i,k) .gt. qsat) then
                niter=0
                dtabs = 100.
                do while(abs(dtabs).gt.0.01.and.niter.lt.10)
                    if(tabs1.ge.tbgmax) then
                        om=1.
                        lstarn=fac_cond
                        dlstarn=0.
                        qsat=qsatw(tabs1,pres(k))
                        dqsat=dtqsatw(tabs1,pres(k))
                           else if(tabs1.le.tbgmin) then
                        om=0.
                        lstarn=fac_sub
                        dlstarn=0.
                        qsat=qsati(tabs1,pres(k))
                        dqsat=dtqsati(tabs1,pres(k))
                    else
                        om=an*tabs1-bn
                        lstarn=fac_cond+(1.-om)*fac_fus
                        dlstarn=an
                        qsat=om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))
                        dqsat=om*dtqsatw(tabs1,pres(k))+(1.-om)*dtqsati(tabs1,pres(k))
                    endif

                    fff = tabs(i,k)-tabs1+lstarn*(q(i,k)-qsat)
                    dfff=dlstarn*(q(i,k)-qsat)-lstarn*dqsat-1.
                    dtabs=-fff/dfff
                    niter=niter+1
                    tabs1=tabs1+dtabs
               end do

               qsat = qsat + dqsat * dtabs
               qn = max(0., q(i,k)-qsat)
               qv(i,k) = max(0., q(i,k)-qn)

            ! If condensation not possible qn is 0.0
            else

              qn = 0.
              qv(i,k) = q(i,k)

            endif

            ! Set tabs to iterated tabs after convection
            tabs(i,k) = tabs1

            !! Code for calculating qcc and qii from qn.
            !! Assumes dokruegermicro=.false. in SAM.
            !! Taken from statistics.f90
            omn = omegan(tabs(i,k))
            qc(i,k) = qn*omn
            qi(i,k) = qn*(1.-omn)

        end do
        end do

    end subroutine SAM_var_conversion

end module nn_interface_CAM
