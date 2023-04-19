! neural_net convection emulator

module nn_convection_flux_mod

!----------------------------------------------------------------------
! Libraries to use

use netcdf
use vars
use grid
use params , only: fac_cond, fac_fus, tprmin, a_pr, cp, a_bg, tbgmin
implicit none
private

!---------------------------------------------------------------------
!  ---- public interfaces ----
public  relu, nn_convection_flux_init, nn_convection_flux

!-----------------------------------------------------------------------
!   ---- local/private data ----
! Parameters from that are used here
!! From domain.f90
!integer, parameter :: YES3D = 1   ! Domain dimensionality: 1 - 3D, 0 - 2D
!integer, parameter :: nx_gl = 72  ! Number of grid points in X - Yani changed to 36 from 32
!integer, parameter :: ny_gl = 180 ! Number of grid points in Y
!integer, parameter :: nz_gl = 48  ! Number of pressure (scalar) levels
!integer, parameter :: nsubdomains_x  = 3  ! No of subdomains in x
!integer, parameter :: nsubdomains_y  = 12 ! No of subdomains in y
!
!! From grid.f90
!integer, parameter :: nx = nx_gl/nsubdomains_x   ! Number of x points in a subdomain
!integer, parameter :: ny = ny_gl/nsubdomains_y   !           y
!integer, parameter :: nz = nz_gl+1               !           z
!integer, parameter :: nzm = nz-1                 ! ???


logical :: do_init=.true.

  ! Neural Net Parameters

  integer :: n_in ! Input dim features
  integer :: n_h1 ! hidden dim
  integer :: n_h2 ! hidden dim
  integer :: n_h3 ! hidden dim
  integer :: n_h4 ! hidden dim
  integer :: n_out ! outputs dim
  integer :: nrf  ! number of vertical levels the NN uses
  integer :: nrfq ! number of vertical levels the NN uses
  integer :: it,jt
  integer :: o_var_dim ! number of output variables  (different types)

  real(4), allocatable, dimension(:,:)     :: r_w1
  real(4), allocatable, dimension(:,:)     :: r_w2
  real(4), allocatable, dimension(:,:)     :: r_w3
  real(4), allocatable, dimension(:,:)     :: r_w4
  real(4), allocatable, dimension(:,:)     :: r_w5
  real(4), allocatable, dimension(:)       :: r_b1
  real(4), allocatable, dimension(:)       :: r_b2
  real(4), allocatable, dimension(:)       :: r_b3
  real(4), allocatable, dimension(:)       :: r_b4
  real(4), allocatable, dimension(:)       :: r_b5
  real(4), allocatable, dimension(:)       :: xscale_mean
  real(4), allocatable, dimension(:)       :: xscale_stnd

  real(4), allocatable, dimension(:)       :: z1
  real(4), allocatable, dimension(:)       :: z2
  real(4), allocatable, dimension(:)       :: z3
  real(4), allocatable, dimension(:)       :: z4
  real(4), allocatable, dimension(:)       :: z5

  real(4), allocatable, dimension(:)       :: yscale_mean
  real(4), allocatable, dimension(:)       :: yscale_stnd


!-----------------------------------------------------------------------



contains

  subroutine relu(logits)

     real(4), dimension(:), intent(inout)   :: logits
     where (logits .lt. 0.0)  logits = 0.0

  end subroutine relu

  subroutine net_forward(features, logits)

    real(4), dimension(:), intent(in)   :: features
    real(4), dimension(:), intent(out)  :: logits

    z1 = matmul( features, r_w1) + r_b1
    call relu(z1)

    z2 = matmul( z1,r_w2) + r_b2
    call relu(z2)

    z3 = matmul( z2,r_w3) + r_b3
    call relu(z3)

    z4 = matmul( z3,r_w4) + r_b4
    call relu(z4)

    logits = matmul( z4,r_w5) + r_b5

  end subroutine


!#######################################################################

  subroutine nn_convection_flux_init(nn_filename)

    !-----------------------------------------------------------------------
    !
    !        initialization for nn convection
    !
    !-----------------------------------------------------------------------
    integer  unit,io,ierr
    
    ! This will be the netCDF ID for the file and data variable.
    integer :: ncid
    integer :: in_dimid, h1_dimid, out_dimid, single_dimid
    integer :: h2_dimid, h3_dimid, h4_dimid
    integer :: r_w1_varid, r_w2_varid, r_b1_varid, r_b2_varid
    integer :: r_w3_varid, r_w4_varid, r_b3_varid, r_b4_varid
    integer :: r_w5_varid, r_b5_varid
    
    
    integer :: xscale_mean_varid, xscale_stnd_varid
    integer :: yscale_mean_varid, yscale_stnd_varid
    
    character(len=1024), intent(in) :: nn_filename
    
    call task_rank_to_index(rank,it,jt)
    
    
    !-------------allocate arrays and read data-------------------------
    
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
    ! Get the varid or dimid for each variable or dimension based on its name.

    call check( nf90_open(     trim(nn_filename),NF90_NOWRITE,ncid ))

    call check( nf90_inq_dimid(ncid, 'N_in', in_dimid))
    call check( nf90_inquire_dimension(ncid, in_dimid, len=n_in))

    call check( nf90_inq_dimid(ncid, 'N_h1', h1_dimid))
    call check( nf90_inquire_dimension(ncid, h1_dimid, len=n_h1))
    call check( nf90_inq_dimid(ncid, 'N_h2', h2_dimid))
    call check( nf90_inquire_dimension(ncid, h2_dimid, len=n_h2))

    call check( nf90_inq_dimid(ncid, 'N_h3', h3_dimid))
    call check( nf90_inquire_dimension(ncid, h3_dimid, len=n_h3))

    call check( nf90_inq_dimid(ncid, 'N_h4', h4_dimid))
    call check( nf90_inquire_dimension(ncid, h4_dimid, len=n_h4))
    call check( nf90_inq_dimid(ncid, 'N_out', out_dimid))
    call check( nf90_inquire_dimension(ncid, out_dimid, len=n_out))

    call check( nf90_inq_dimid(ncid, 'N_out_dim', out_dimid))
    call check( nf90_inquire_dimension(ncid, out_dimid, len=o_var_dim))

    print *, 'size of features', n_in
    print *, 'size of outputs', n_out

    nrf = 30 ! Size in the vertical
    nrfq = 29 !Size in the vertical  for advection

    call check( nf90_open(     trim(nn_filename),NF90_NOWRITE,ncid ))

    ! Allocate arrays for NN based on sizes read from nc file
    allocate(r_w1(n_in,n_h1))
    allocate(r_w2(n_h1,n_h2))
    allocate(r_w3(n_h2,n_h3))
    allocate(r_w4(n_h3,n_h4))
    allocate(r_w5(n_h4,n_out))

    allocate(r_b1(n_h1))
    allocate(r_b2(n_h2))
    allocate(r_b3(n_h3))
    allocate(r_b4(n_h4))
    allocate(r_b5(n_out))
    allocate(z1(n_h1))
    allocate(z2(n_h2))
    allocate(z3(n_h3))
    allocate(z4(n_h4))
    allocate(z5(n_out))

    allocate(xscale_mean(n_in))
    allocate(xscale_stnd(n_in))

    allocate(yscale_mean(o_var_dim))
    allocate(yscale_stnd(o_var_dim))

    ! Read NN data in from nc file
    call check( nf90_inq_varid(ncid, "w1", r_w1_varid))
    call check( nf90_get_var(ncid, r_w1_varid, r_w1))
    call check( nf90_inq_varid(ncid, "w2", r_w2_varid))
    call check( nf90_get_var(ncid, r_w2_varid, r_w2))

    call check( nf90_inq_varid(ncid, "w3", r_w3_varid))
    call check( nf90_get_var(ncid, r_w3_varid, r_w3))
    call check( nf90_inq_varid(ncid, "w4", r_w4_varid))
    call check( nf90_get_var(ncid, r_w4_varid, r_w4))

    call check( nf90_inq_varid(ncid, "w5", r_w5_varid))
    call check( nf90_get_var(ncid, r_w5_varid, r_w5))

    call check( nf90_inq_varid(ncid, "b1", r_b1_varid))
    call check( nf90_get_var(ncid, r_b1_varid, r_b1))
    call check( nf90_inq_varid(ncid, "b2", r_b2_varid))
    call check( nf90_get_var(ncid, r_b2_varid, r_b2))

    call check( nf90_inq_varid(ncid, "b3", r_b3_varid))
    call check( nf90_get_var(ncid, r_b3_varid, r_b3))
    call check( nf90_inq_varid(ncid, "b4", r_b4_varid))
    call check( nf90_get_var(ncid, r_b4_varid, r_b4))
    call check( nf90_inq_varid(ncid, "b5", r_b5_varid))
    call check( nf90_get_var(ncid, r_b5_varid, r_b5))

    call check( nf90_inq_varid(ncid,"fscale_mean",     xscale_mean_varid))
    call check( nf90_get_var(  ncid, xscale_mean_varid,xscale_mean      ))
    call check( nf90_inq_varid(ncid,"fscale_stnd",     xscale_stnd_varid))
    call check( nf90_get_var(  ncid, xscale_stnd_varid,xscale_stnd      ))

    call check( nf90_inq_varid(ncid,"oscale_mean",     yscale_mean_varid))
    call check( nf90_get_var(  ncid, yscale_mean_varid,yscale_mean      ))
    call check( nf90_inq_varid(ncid,"oscale_stnd",     yscale_stnd_varid))
    call check( nf90_get_var(  ncid, yscale_stnd_varid,yscale_stnd      ))

    ! Close the nc file
    call check( nf90_close(ncid))

    write(*, *) 'Finished reading NN regression file.'

    do_init=.false.

  end subroutine nn_convection_flux_init
  !#######################################################################

  subroutine nn_convection_flux(tabs)

    !-----------------------------------------------------------------------
    !
    !  NN subgrid parameterization
    !
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    !
    !   input:  tabs               absolute temperature
    !           q                  total non-precipitating water
    !           distance from equator
    !   changes: t                 liquid static energy as temperature
    !            q                 total non-precipitating water
    !
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
  
  
    !-----------------------------------------------------------------------
    !---------------------- inputs -------------------------------------
    real, intent(in), dimension(:, :, :)  :: tabs                ! temperature
  
    !-----------------------------------------------------------------------
    !---------------------- local data -------------------------------------
    real,   dimension(nrf) :: t_tendency_adv, q_tendency_adv, q_tendency_auto, q_tendency_sed,t_tendency_auto
    real,   dimension(nrf) :: q_flux_sed, qp_flux_fall,t_tendency_sed, q_tend_tot
    real,   dimension(nrf) :: t_flux_adv, q_flux_adv, t_sed_flux, t_rad_rest_tend, omp,fac ! Do not predict surface adv flux
    real,   dimension(nzm) :: qsat, irhoadz, irhoadzdz, irhoadzdz2
    real(4), dimension(n_in) :: features
    real(4), dimension(n_out) :: outputs
    real    omn, rev_dz, rev_dz2
    integer  i, j, k,dd, dim_counter, out_dim_counter, out_var_control
  
    ! Check that we have initialised all of the variables.
    if (do_init) call error_mesg('nn_convection_flux_init has not been called.')
  
    ! Initialise precipitation if required
    if (.not. rf_uses_qp) then
      if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) precsfc(:,:)=0.
    end if
  
    ! Define TODO Magical Mystery variables
    rev_dz = 1/dz
    rev_dz2 = 1/(dz*dz)
    do k=1,nzm
      irhoadz(k) = dtn/(rho(k)*adz(k)) ! Useful factor
      irhoadzdz(k) = irhoadz(k)/dz ! Note the time step
      irhoadzdz2(k) = irhoadz(k)/(dz*dz) ! Note the time step
    end do
  
    ! The NN operates on columns, so loop over x and y coordinates in turn
    do j=1,ny
      do i=1,nx
  
        ! Initialize variables
        features = 0.
        outputs = 0.
        t_tendency_adv = 0.
        q_tendency_adv = 0.
        q_tendency_auto = 0.
        t_tendency_auto = 0.
        q_tendency_sed = 0.
        t_tendency_sed = 0.
        t_rad_rest_tend = 0.
        q_tend_tot = 0.
        t_flux_adv = 0.
        q_flux_adv = 0.
        q_flux_sed = 0.
        z1 = 0.
        z2 = 0.
        z3 = 0.
        z4 = 0.
        z5 = 0.
        dim_counter = 0
        omp = 0.
        fac = 0.
  
        !-----------------------------------------------------------------------
        ! Combine all features into one vector
  
        ! If using temperature then add as input feature
        if (Tin_feature_rf) then
          features(dim_counter+1:dim_counter + input_ver_dim) = real(t_i(i,j,1:input_ver_dim),4)
          dim_counter = dim_counter + input_ver_dim
        endif
  
        ! If using water content then add as input feature
        if (qin_feature_rf) then
          ! If using generalised relative humidity to estimate water content
          if (rf_uses_rh) then
            do k=1,nzm
              omn = omegan(tabs(i,j,k))
              qsat(k) = omn*qsatw(tabs(i,j,k),pres(k))+(1.-omn)*qsati(tabs(i,j,k),pres(k))
            end do
            features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim)/qsat(1:input_ver_dim),4)
            dim_counter = dim_counter + input_ver_dim
          ! if using non-precipitating water as water content
          else
            features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim),4)
            dim_counter =  dim_counter + input_ver_dim
          endif
        endif
        
        ! If using TODO?? then add as input feature
        if (rf_uses_qp) then
          features(dim_counter+1:dim_counter+input_ver_dim) = real(qp_i(i,j,1:input_ver_dim),4)
          dim_counter = dim_counter + input_ver_dim
        endif
   
        ! If using TODO - Some magical mystery input?? then add as input feature
        if(do_yin_input) then
          features(dim_counter+1) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
          dim_counter = dim_counter+1
        endif

        !-----------------------------------------------------------------------
        !Normalize features
        features = (features - xscale_mean) / xscale_stnd
        
        call net_forward(features, outputs)

        !-----------------------------------------------------------------------
        ! Separate out outputs into heating and moistening tendencies
        out_var_control =1
        t_rad_rest_tend(1:nrf) = (outputs(1:nrf) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)
        out_dim_counter = nrf

        out_var_control = out_var_control + 1
        t_flux_adv(2:nrf) = (outputs(out_dim_counter+1:out_dim_counter+nrfq)* yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)
        out_dim_counter = out_dim_counter + nrfq

        out_var_control = out_var_control + 1
        q_flux_adv(2:nrf) = (outputs(out_dim_counter+1:out_dim_counter+nrfq) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)
        out_dim_counter = out_dim_counter + nrfq

        out_var_control = out_var_control + 1
        q_tendency_auto(1:nrf) = (outputs(out_dim_counter+1:out_dim_counter+nrf) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)
        out_dim_counter = out_dim_counter + nrf

        out_var_control = out_var_control + 1
        q_flux_sed(1:nrf) = (outputs(out_dim_counter+1:out_dim_counter+nrf) * yscale_stnd(out_var_control))  +  yscale_mean(out_var_control)
        out_dim_counter = out_dim_counter + nrf

        out_var_control = out_var_control + 1

        ! advection surface flux is zero
        t_flux_adv(1) = 0.0
        q_flux_adv(1) = 0.0
  
        do k=2,nrf
          if (q_flux_adv(k).lt.0) then
            if ( q(i,j,k).lt.-q_flux_adv(k)* irhoadzdz(k)) then
              q_flux_adv(k) = -q(i,j,k)/irhoadzdz(k)
            end if
          else
            if (q(i,j,k-1).lt.q_flux_adv(k)* irhoadzdz(k)) then
              q_flux_adv(k) = q(i,j,k-1)/irhoadzdz(k)
            end if
          end if
        end do
  
        do k=1,nrf-1
           t_tendency_adv(k) = - (t_flux_adv(k+1) - t_flux_adv(k)) * irhoadzdz(k)
           q_tendency_adv(k) = - (q_flux_adv(k+1) - q_flux_adv(k)) * irhoadzdz(k)
        end do

        k = nrf
        t_tendency_adv(k) = - (0.0 - t_flux_adv(k)) * irhoadzdz(k)
        q_tendency_adv(k) = - (0.0 - q_flux_adv(k)) * irhoadzdz(k)

        do k=1,nrf
          if (q(i,j,k).lt.-q_tendency_adv(k)) then
            q_tendency_adv(k) = -q(i,j,k)
          end if
        end do
  
        t(i,j,1:nrf) = t(i,j,1:nrf) + t_tendency_adv(1:nrf)
        q(i,j,1:nrf) = q(i,j,1:nrf) + q_tendency_adv(1:nrf)
  
        do k=1,nrf
          omp(k) = max(0.,min(1.,(tabs(i,j,k)-tprmin)*a_pr))
          fac(k) = (fac_cond + fac_fus * (1.0 - omp(k)))
          if (q_tendency_auto(k).lt.0) then
            q_tend_tot(k) = min(-q_tendency_auto(k) * dtn, q(i,j,k))!q_tendency_auto(1:nrf) * dtn + q_tendency_adv(1:nrf) + q_tendency_sed(1:nrf)
            q_tend_tot(k) = -q_tend_tot(k)
          else
            q_tend_tot(k) = q_tendency_auto(k) * dtn
          endif
        end do

        q(i,j,1:nrf) = q(i,j,1:nrf) +q_tend_tot(1:nrf)
        t(i,j,1:nrf) = t(i,j,1:nrf)  -q_tend_tot(1:nrf)*fac(1:nrf)

        do k=2,nrf
          if (q_flux_sed(k).lt.0) then
            if ( q(i,j,k).lt.-q_flux_sed(k)* irhoadzdz(k)) then
              q_flux_sed(k) = -q(i,j,k)/irhoadzdz(k)
            end if
          else
            if (q(i,j,k-1).lt.q_flux_sed(k)* irhoadzdz(k)) then
              q_flux_sed(k) = q(i,j,k-1)/irhoadzdz(k)
            end if
          end if
        end do

        do k=1,nrf-1 ! One level less than I actually use
          q_tendency_sed(k) = - (q_flux_sed(k+1) - q_flux_sed(k)) * irhoadzdz(k)
        end do

        k = nrf
        q_tendency_sed(k) = - (0.0 - q_flux_sed(k)) * irhoadzdz(k)

        do k=1,nrf
          if (q_tendency_sed(k).lt.0) then
            q_tendency_sed(k) = min(-q_tendency_sed(k), q(i,j,k))
            q_tendency_sed(k) = -q_tendency_sed(k)
          end if
        end do

        t(i,j,1:nrf) = t(i,j,1:nrf) - q_tendency_sed(1:nrf)*(fac_fus+fac_cond)
        q(i,j,1:nrf) = q(i,j,1:nrf) +q_tendency_sed(1:nrf)

        t(i,j,1:nrf) = t(i,j,1:nrf) + t_rad_rest_tend(1:nrf)*dtn

        precsfc(i,j) = precsfc(i,j)  - q_flux_sed(1)*dtn*rev_dz ! For statistics
        prec_xy(i,j) = prec_xy(i,j)  - q_flux_sed(1)*dtn*rev_dz ! For 2D output

        do k=1, nrf
          precsfc(i,j) = precsfc(i,j)-q_tend_tot(k)*adz(k)*dz*rho(k)*(1/dz)! removed the time step mult because q_tend_tot is already mult
          prec_xy(i,j) = prec_xy(i,j)-q_tend_tot(k)*adz(k)*dz*rho(k)*(1/dz)
        end do

        do k = 1,nrf
          q(i,j,k)=max(0.,q(i,j,k))
        end do

        where (qn(i,j,1:nrf).gt.0.0)
          qn(i,j,1:nrf) = qn(i,j,1:nrf)+ q_tend_tot(1:nrf) + q_tendency_adv(1:nrf) + q_tendency_sed(1:nrf)
        end where

        where (qn(i,j,:).lt.0.0)
          qn(i,j,:) = 0.0
        end where

      end do
    end do
    ! End of loop over x, y, columns

  end subroutine nn_convection_flux


  !##############################################################################
  subroutine check(status)

    ! checks error status after each netcdf, prints out text message each time
    !   an error code is returned.

    integer, intent(in) :: status

    if(status /= nf90_noerr) then
    !   write(*, *) trim(nf90_strerror(status))
    end if
  end subroutine check


  !#######################################################################
  subroutine error_mesg (message)
    character(len=*), intent(in) :: message

    !  input:
    !      message   message written to output   (character string)

    if(masterproc) print*, 'Neural network  module: ', message
    stop

  end subroutine error_mesg


  !#######################################################################
  subroutine task_rank_to_index (rank,i,j)

    ! returns the pair of  beginning indices for the subdomain on the  
    ! global grid given the subdomain's rank.
    integer ::  rank, i, j

    j = rank/nsubdomains_x 
    i = rank - j*nsubdomains_x

    i = i * (nx_gl/nsubdomains_x) 
    j = j * (ny_gl/nsubdomains_y) 

  end subroutine task_rank_to_index

  !#######################################################################
  ! Need omegan function to run, ripped from https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/omega.f90
  ! 
  real function omegan(tabs)
    real :: tabs
    omegan = max(0.,min(1.,(tabs-tbgmin)*a_bg))
    return
  end function omegan

  !#######################################################################
  ! Need qsatw functions to run, ripped from https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/sat.f90
  ! 
  ! Saturation vapor pressure and mixing ratio. 
  ! Based on Flatau et.al, (JAM, 1992:1507)
  !

  real function esatw(t)
    implicit none
    real :: t  ! temperature (K)
    real :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
            6.105851, 0.4440316, 0.1430341e-1, &
            0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
            0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
    !       6.11239921, 0.443987641, 0.142986287e-1, &
    !       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
    !       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real :: dt
    dt = max(-80.,t-273.16)
    esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  end function esatw

  real function qsatw(t,p)
    implicit none
    real :: t  ! temperature (K)
    real :: p  ! pressure    (mb)
    real :: esat
    esat = esatw(t)
    qsatw = 0.622 * esat/max(esat, p-esat)
  end function qsatw

  ! real function dtesatw(t)
  !   implicit none
  !   real :: t  ! temperature (K)
  !   real :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
  !   data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
  !             0.443956472, 0.285976452e-1, 0.794747212e-3, &
  !             0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
  !            -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
  !   real :: dt
  !   dt = max(-80.,t-273.16)
  !   dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  ! end function dtesatw
  ! 
  ! real function dtqsatw(t,p)
  !   implicit none
  !   real :: t  ! temperature (K)
  !   real :: p  ! pressure    (mb)
  !   real :: dtesatw
  !   dtqsatw=0.622*dtesatw(t)/p
  ! end function dtqsatw

  real function esati(t)
    implicit none
    real :: t  ! temperature (K)
    real :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
            6.11147274, 0.503160820, 0.188439774e-1, &
            0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
            0.385852041e-9, 0.146898966e-11, 0.252751365e-14/       
    real :: dt
    dt = max(-80.0, t-273.16)
    esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  end function esati

  real function qsati(t,p)
    implicit none
    real :: t  ! temperature (K)
    real :: p  ! pressure    (mb)
    real :: esat
    esat = esati(t)
    qsati = 0.622 * esat/max(esat,p-esat)
  end function qsati
  !         
  !         
  ! function dtesati(t)
  !   implicit none
  !   real :: t  ! temperature (K)
  !   real :: a0,a1,a2,a3,a4,a5,a6,a7,a8 
  !   data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
  !           0.503223089, 0.377174432e-1,0.126710138e-2, &
  !           0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
  !           0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
  !   real :: dtesati, dt
  !   dt = max(-800. ,t-273.16)
  !   dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt))))))) 
  !   return
  ! end function dtesati
  !         
  !         
  ! real function dtqsati(t,p)
  !   implicit none
  !   real :: t  ! temperature (K)
  !   real :: p  ! pressure    (mb)
  !   real :: dtesati
  !   dtqsati = 0.622 * dtesati(t) / p
  ! end function dtqsati

end module nn_convection_flux_mod
