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
                                  esati, qsati, qsatw
use SAM_consts_mod, only: nrf, ggr, cp, tbgmax, tbgmin, &
                          fac1, fac2, fac_cond, fac_sub, &
                          an, bn
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux_CAM, &
        nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize

! TODO Make routines public for purposes of testing
public interp_to_sam, interp_to_cam, fetch_sam_data

!---------------------------------------------------------------------
! local/private data

!= unit 1 :: nz_sam
integer :: nz_sam
    !! number of vertical values in the SAM sounding profiles
!= unit m :: z
!= unit hPa :: pres, presi
!= unit kg m-3 :: rho
!= unit 1 :: adz
real, allocatable, dimension(:) :: z, pres, presi, rho, adz, gamaz
    !! SAM sounding variables
!= unit m :: dz
real :: dz
    !! grid spacing in z direction for the lowest grid layer

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine fetch_sam_data(pressam, presisam)
        !! Temporary subroutine to extract SAM pressure profile
        real, dimension(48), intent(out) :: pressam, presisam
        pressam(:) = pres(:)
        presisam(:) = presi(:)

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
                                      tabs_i, q_i, &
                                      tabs, &
                                      dtn, dy, &
                                      nx, ny, &
                                      nstep, nstatis, icycle, &
                                      precsfc, prec_xy)
        !! Interface to the nn_convection parameterisation for the SAM model

        integer :: j, k

        real, dimension(:,:) :: pres_cam
            !! pressure [hPa] from the CAM model
        real, dimension(:) :: pres_int_cam
            !! interface pressure [hPa] from the CAM model (element 1 is surface pressure)
        real, dimension(:,:) :: tabs_i, tabs!, t
        real, dimension(:,:,:) :: q_i!, q
        real, dimension(:, :) :: precsfc, prec_xy
        real, intent(in) :: dtn
        real, intent(in) :: dy
        integer, intent(in) :: nx, ny, nstep, nstatis, icycle

        real :: y_in(nx, ny)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        real, dimension(nx, ny, nrf) :: t_rad_rest_tend, t_delta_adv, q_delta_adv, &
                                        t_delta_auto, t_delta_sed, &
                                        q_delta_auto, q_delta_sed
            !! deltas/tendencies returned by the parameterisation:
            !! radiation rest tendency, advective, autoconversion, sedimentation

        real, dimension(nx, ny)      :: prec_sed
            !! Sedimenting precipitation at surface
       
        ! TODO: Does CAM require precipitation?
        ! Initialise precipitation to 0 if required and at start of cycle
        if(mod(nstep-1,nstatis).eq.0 .and. icycle.eq.1) then
            precsfc(:,:)=0.
        end if

        ! distance to the equator
        ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
        do j=1, ny
            ! TODO: Set y_in as appropriate for CAM
            y_in(:,j) = 0.0
        enddo

        !-----------------------------------------------------
        
        ! TODO: Formulate the input variables to the parameterisation as required.

        !-----------------------------------------------------
        
        ! TODO
        ! Interpolate CAM variables to the SAM pressure levels
        ! call interp_to_sam(pres_cam, pres_int_cam(1), var_cam, var_sam)

        !-----------------------------------------------------
        
        ! Run the neural net parameterisation
!        call nn_convection_flux(tabs_i(:,:,1:nrf), q_i(:,:,1:nrf), y_in, &
!                                tabs(:,:,1:nrf), &
!                                t(:,:,1:nrf), q(:,:,1:nrf), &
!                                rho, adz, dz, dtn, &
!                                t_rad_rest_tend, &
!                                t_delta_adv, q_delta_adv, &
!                                t_delta_auto, q_delta_auto, &
!                                t_delta_sed, q_delta_sed, prec_sed)

        !-----------------------------------------------------
        
        ! TODO: Update CAM variables with NN outputs
        ! q_delta_adv(:,:,:)
        ! q_delta_auto(:,:,:)
        ! q_delta_sed(:,:,:)
        ! t_delta_adv(:,:,:)
        ! t_delta_auto(:,:,:)
        ! t_delta_sed(:,:,:)
        ! t_rad_rest_tend(:,:,:)*dtn

        ! TODO: Update precipitation if required

        
        ! TODO: Formulate the output variables to CAM as required.

        !-----------------------------------------------------
        
        ! TODO
        ! Interpolate SAM variables to the CAM pressure levels
        ! call interp_to_cam(pres_cam, pres_int_cam, var_sam, var_cam)


    end subroutine nn_convection_flux_CAM


    subroutine nn_convection_flux_CAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()
        call sam_sounding_finalize()

    end subroutine nn_convection_flux_CAM_finalize
    

    !-----------------------------------------------------------------
    ! Private Subroutines
    
    subroutine interp_to_sam(p_cam, ps_cam, var_cam, var_sam, var_cam_surface)
        !! Interpolate from the CAM pressure grid to the SAM pressure grid.
        !! Uses linear interpolation between grid points.

        != unit hPa :: p_cam, ps_cam
        real, dimension(:,:), intent(in) :: p_cam
            !! pressure [hPa] from the CAM model
        real, dimension(:), intent(in) :: ps_cam
            !! surface pressure [hPa] for each column in CAM
        != unit 1 :: var_cam, var_sam
        real, dimension(:,:), intent(in) :: var_cam
            !! variable from CAM grid to be interpolated from
        real, dimension(:), intent(in) :: var_cam_surface
            !! Surface value of CAM variable for interpolation/boundary condition
        real, dimension(:,:), intent(out) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real, dimension(nrf) :: p_norm_sam
        real, dimension(:,:), allocatable :: p_norm_cam

        integer :: i, k, nz_cam, ncol_cam, kc_u, kc_l
        real :: pc_u, pc_l

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(size(p_cam, 1), nz_cam))

        ! Normalise pressures by surface value
        p_norm_sam(:) = pres(:) / presi(1)
        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / ps_cam(:)
        end do

        ! Loop over columns and SAM levels and interpolate each column
        do k = 1, nrf
        do i = 1, ncol_cam
            ! Check we are within array

            ! If lowest SAM is below lowest CAM - interpolate to surface value
            if (p_norm_sam(k) > p_norm_cam(i, 1)) then
                ! TODO Remove warning and interpolate to boundary
                ! TODO Add boundary conditions for any input variables as required.
                write(*,*) "Interpolating to surface."
                var_sam(i, k) = var_cam_surface(i) &
                                + (p_norm_sam(k)-1.0) &
                                *(var_cam(i, 1)-var_cam_surface(i))/(p_norm_cam(i, 1)-1.0)
            elseif (p_norm_cam(i, nz_cam) > p_norm_sam(nrf)) then
                ! Also check the top
                ! This should not happen as CAM grid extends to higher altitudes than SAM
                write(*,*) "CAM upper pressure level is lower than that of SAM: Stopping."
                stop
            else
                ! Locate the neighbouring CAM indices to interpolate between
                ! This will be slow - speed up later
                kc_u = 1
                pc_u = p_norm_cam(i, kc_u)
                do while (pc_u > p_norm_sam(k))
                    kc_u = kc_u + 1
                    pc_u = p_norm_cam(i, kc_u)
                end do
                kc_l = kc_u - 1
                pc_l = p_norm_cam(i, kc_l)
                do while (pc_l < p_norm_sam(k))
                    kc_l = kc_l - 1
                    pc_l = p_norm_cam(i, kc_l)
                end do

                ! interpolate variables - Repeat for as many variables as need interpolating
                var_sam(i, k) = var_cam(i, kc_l) &
                                + (p_norm_sam(k)-pc_l) &
                                *(var_cam(i, kc_u)-var_cam(i, kc_l))/(pc_u-pc_l)
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
        real, dimension(nrf) :: p_int_norm_sam
        real, dimension(:,:), allocatable :: p_norm_cam
        real, dimension(:,:), allocatable :: p_int_norm_cam

        integer :: i, k, nz_cam, ncol_cam, ks, ks_u, ks_l
        real :: ps_u, ps_l, pc_u, pc_l
        real :: p_sam_min

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(size(p_cam, 1), nz_cam))
        allocate(p_int_norm_cam(size(p_int_cam, 1), nz_cam))

        ! Normalise pressures by surface value
        p_norm_sam(:) = pres(:) / presi(1)
        p_int_norm_sam(:) = presi(:) / presi(1)
        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_int_cam(:, 1)
            !! TODO Check size of CAM column interfaces.
            p_int_norm_cam(:,k) = p_int_cam(:,k) / p_int_cam(:, 1)
        end do

        ! Loop over columns and SAM levels and interpolate each column
        ! TODO Revert indices
        do k = 1, nz_cam
        do i = 1, ncol_cam
            ! Check we are within array overlap:
            ! No need to compare bottom end as pressure is normalised.

            ! Set top of column values
            p_sam_min = 2.0*p_norm_sam(nrf) - p_int_norm_sam(nrf)
!            write(*,*) p_norm_sam(nrf), p_int_norm_sam(nrf), p_sam_min

            if (p_int_norm_cam(i, k) < p_sam_min) then
                ! Anything above the SAM NN grid should be set to 0.0
                write(*,*) "CAM pressure is lower than SAM top: set tend to 0.0."
                var_cam(i, k) = 0.0

            else
                ! Get the pressures at the top and bottom of the CAM cell
                if (k == nz_cam) then
                    pc_u = 2.0*p_norm_cam(i, k) - p_int_norm_cam(i, k)
                else
                    pc_u = p_int_norm_cam(i,k+1)
                endif
                pc_l = p_int_norm_cam(i,k)

                ! Locate the SAM indices to interpolate between
                ! This will be slow - speed up later
                ! ks_u is the SAM block with an upper interface above the CAM block
                ks_u = 1
                ps_u = p_int_norm_sam(ks_u+1)
                do while ((ps_u > pc_u) .and. (ks_u < nrf))
                    ks_u = ks_u + 1
                    if (ks_u == nrf) then
                        ps_u = p_sam_min
                    else
                        ps_u = p_int_norm_sam(ks_u+1)
                    endif
                end do
                ! ks_l is the SAM block with a lower interface below the CAM block
                ks_l = ks_u
                ps_l = p_int_norm_sam(ks_l)
                do while ((ps_l < pc_l) .and. (ks_l > 1))
                    ks_l = ks_l - 1
                    ps_l = p_int_norm_sam(ks_l)
                end do
!                write(*,*) "ks_u = ", ks_u, ps_u, "ks_l = ", ks_l, ps_l, "CAM = ", pc_u, pc_l

                ! Combine all SAM cells that this CAM_CELL(i, k) overlaps
                var_cam(i,k) = 0.0
!                write(*,*) "pre: Var_Cam = ", var_cam(i,k)/(pc_u-pc_l)
                
                if (pc_l < p_int_norm_sam(ks_u)) then
                    ! CAM cell fully enclosed by SAM cell
                    var_cam(i,k) = var_cam(i,k) + (max(pc_u, ps_u) - pc_l)*var_sam(i, ks_u)/(ps_u-p_int_norm_sam(ks_u))
                else
                    ! Top cell if CAM and SAM cells intermesh
                    var_cam(i,k) = var_cam(i,k) + (pc_u - p_int_norm_sam(ks_u))*var_sam(i, ks_u)/(ps_u-p_int_norm_sam(ks_u))
                endif
!                write(*,*) "upp: Var_Cam = ", var_cam(i,k)/(pc_u-pc_l)
                ! Lower cell
                if (ks_l /= ks_u) then
                    var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(ks_l+1)-pc_l)*var_sam(i,ks_l)/(p_int_norm_sam(ks_l+1)-ps_l)
                endif
!                write(*,*) "low: Var_Cam = ", var_cam(i,k)/(pc_u-pc_l)
                ! Intermediate cells
                if (ks_u > ks_l+1) then
                    ! TODO Look at index labelling as currently moving 2, 0 here
                    do ks = ks_l+1, ks_u-1
                        var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(ks+1)-p_int_norm_sam(ks))*var_sam(i,ks)/(p_int_norm_sam(ks+1)-p_int_norm_sam(ks))
                    enddo
                endif
!                write(*,*) "int: Var_Cam = ", var_cam(i,k)/(pc_u-pc_l)
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


    ! TODO Duplicated from nn_cf_net.f90 - consider placing 1 copy in shared location?
    subroutine check(err_status)
        !! Check error status after netcdf call and print message for
        !! error codes.

        integer, intent(in) :: err_status
            !! error status from nf90 function

        if(err_status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(err_status))
        end if

    end subroutine check

    subroutine t_q_conversion(tabs, qn, qp, t)

        integer :: nx, ny
            !! array sizes
        integer :: nzm
            !! Number of z points in a subdomain - 1
        integer :: i, j, k
            !! Counters

        ! ---------------------
        ! Fields from SAM
        ! ---------------------
        != unit K :: tabs, tabs1
        real, intent(inout) :: tabs(:, :, :)
            !! absolute temperature
        real, allocatable :: tabs1
            !! Temporary variable for tabs

        != unit 1 :: q, qp, qn0, q
        real, intent(in) :: qn(:, :, :)
            !! Total water
        real, intent(in) :: qp(:, :, :)
            !! Precipitable water (rain+snow)
        real :: qn0
            !! Temporary variable for q
        real, allocatable :: q(:, :, :)
            !! Copy

        != unit  :: t
        real, intent(in) :: t(:, :, :)
            !! t

        real :: qsat, om

        nx = size(tabs, 1)
        ny = size(tabs, 2)
        nzm = size(tabs, 3)

        allocate(q(nx,ny,nrf))

        do k = 1, nzm
        do j = 1, ny
        do i = 1, nx
        
!              if(domicroscaling) dtn = dtn_scaled(i, j, k)
        
            qn0 = qn(i,j,k)
        
            q(i,j,k)=max(0.,q(i,j,k))
        
        
            ! Initial guess for temperature assuming no cloud water/ice:
            tabs(i,j,k) = t(i,j,k)-gamaz(k)
            tabs1=(tabs(i,j,k)+fac1*qp(i,j,k))/(1.+fac2*qp(i,j,k))
        
            ! Warm cloud:
            if(tabs1.ge.tbgmax) then
                tabs1=tabs(i,j,k)+fac_cond*qp(i,j,k)
                qsat = qsatw(tabs1,pres(k))
        
            ! Ice cloud:
            elseif(tabs1.le.tbgmin) then
                tabs1=tabs(i,j,k)+fac_sub*qp(i,j,k)
                qsat = qsati(tabs1,pres(k))
        
            ! Mixed-phase cloud:
            else
                om = an*tabs1-bn
                qsat = om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))

            endif
!         
!         
!         !  Test if condensation is possible:
!         
!         
!             if(q(i,j,k) .gt. qsat) then
!         
!               niter=0
!               dtabs = 100.
!               do while(abs(dtabs).gt.0.01.and.niter.lt.10)
!         	if(tabs1.ge.tbgmax) then
!         	   om=1.
!         	   lstarn=fac_cond
!         	   dlstarn=0.
!         	   qsat=qsatw(tabs1,pres(k))
!         	   dqsat=dtqsatw(tabs1,pres(k))
!                 else if(tabs1.le.tbgmin) then
!         	   om=0.
!         	   lstarn=fac_sub
!         	   dlstarn=0.
!         	   qsat=qsati(tabs1,pres(k))
!         	   dqsat=dtqsati(tabs1,pres(k))
!         	else
!         	   om=an*tabs1-bn
!         	   lstarn=fac_cond+(1.-om)*fac_fus
!         	   dlstarn=an
!         	   qsat=om*qsatw(tabs1,pres(k))+(1.-om)*qsati(tabs1,pres(k))
!         	   dqsat=om*dtqsatw(tabs1,pres(k))+(1.-om)*dtqsati(tabs1,pres(k))
!         	endif
!         	if(tabs1.ge.tprmax) then
!         	   omp=1.
!         	   lstarp=fac_cond
!         	   dlstarp=0.
!                 else if(tabs1.le.tprmin) then
!         	   omp=0.
!         	   lstarp=fac_sub
!         	   dlstarp=0.
!         	else
!         	   omp=ap*tabs1-bp
!         	   lstarp=fac_cond+(1.-omp)*fac_fus
!         	   dlstarp=ap
!         	endif
!         	fff = tabs(i,j,k)-tabs1+lstarn*(q(i,j,k)-qsat)+lstarp*qp(i,j,k)
!         	dfff=dlstarn*(q(i,j,k)-qsat)+dlstarp*qp(i,j,k)-lstarn*dqsat-1.
!         	dtabs=-fff/dfff
!         	niter=niter+1
!         	tabs1=tabs1+dtabs
!               end do   
!         
!               qsat = qsat + dqsat * dtabs
!               qn(i,j,k) = max(0.,q(i,j,k)-qsat)
!         
!             else
!         
!               qn(i,j,k) = 0.
!         
!             endif
!         
!             tabs(i,j,k) = tabs1
!             qp(i,j,k) = max(0.,qp(i,j,k)) ! just in case
!         
        end do
        end do
        end do

        deallocate(q)

    end subroutine t_q_conversion

end module nn_interface_CAM
