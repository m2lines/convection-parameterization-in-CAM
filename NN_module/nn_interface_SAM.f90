module nn_interface_SAM_mod
    !! Interface to convection NN for the SAM model
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use nn_convection_flux_mod, only: nn_convection_flux_init, net_forward
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux_SAM, nn_convection_flux_SAM_init


!---------------------------------------------------------------------
! local/private data

integer :: it,jt
    !! indices corresponding to the start of the grid domain for
    !! current MPI rank

! Neural Net parameters
integer :: n_inputs, n_outputs
    !! Length of input/output vector to the NN
integer :: nrf
    !! number of vertical levels the NN uses
integer :: nrfq
    !! number of vertical levels the NN uses
logical :: do_init=.true.
    !! model initialisation is yet to be performed


! Parameters from SAM that are used here
! From domain.f90
integer, parameter :: YES3D = 1
    !! Domain dimensionality: 1 - 3D, 0 - 2D
integer, parameter :: nx_gl = 72
    !! Number of grid points in X - Yani changed to 36 from 32
integer, parameter :: ny_gl = 180
    !! Number of grid points in Y
integer, parameter :: nz_gl = 48
    !! Number of pressure (scalar) levels
integer, parameter :: nsubdomains_x  = 3
    !! No of subdomains in x
integer, parameter :: nsubdomains_y  = 12
    !! No of subdomains in y

! From grid.f90
integer input_ver_dim
integer, parameter :: nx = nx_gl/nsubdomains_x
    !! Number of x points in a subdomain
integer, parameter :: ny = ny_gl/nsubdomains_y
    !! Number of y points in a subdomain
integer, parameter :: nz = nz_gl+1
    !! Number of z points in a subdomain
! Store useful variations on these values
integer, parameter :: nzm = nz-1
integer, parameter :: nxp3 = nx + 3
integer, parameter :: nyp3 = ny + 3 * YES3D
integer, parameter :: dimx1_s = -2
integer, parameter :: dimx2_s = nxp3
integer, parameter :: dimy1_s = 1-3*YES3D
integer, parameter :: dimy2_s = nyp3
real dy
    !! grid spacing in y direction
real dz
    !! grid spacing in z direction for the lowest grid layer
real dtn
    !! current dynamical timestep (can be smaller than dt)
real pres(nzm)
    !! pressure,mb at scalar levels
real adz(nzm)
    !! ratio of the grid spacing to dz for pressure levels
integer nstep
    !! current number of performed time steps 
integer icycle
    !! current subcycle 
integer nstatis
    !! the interval between substeps to compute statistics
logical tin_feature_rf, rf_uses_rh, rf_uses_qp, qin_feature_rf, do_yin_input
!! Logical switches and flags:
! Multitasking stuff
integer rank
    !! rank of the current subdomain task (default 0) 
logical masterproc
    !! logical variable .true. if rank .eq. 0 

! From params.f90
! Constants:
real, parameter :: lfus = 0.3336e+06
    !! Latent heat of fusion, J/kg
real, parameter :: lcond = 2.5104e+06
    !! Latent heat of condensation, J/kg
real, parameter :: cp = 1004.
    !! Specific heat of air, J/kg/K
real, parameter :: fac_cond = lcond/cp
    !! 
real, parameter :: fac_fus = lfus/cp
    !! 
! Temperatures limits for various hydrometeors
real, parameter :: tprmin = 268.16
    !! Minimum temperature for rain, K
real, parameter :: tbgmin = 253.16
    !! Minimum temperature for cloud water., K

real :: a_pr, a_bg
    !! Misc. microphysics variables

! From vars.f90
! prognostic variables:
real t(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    !! moist static energy
real q(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
    !! total water
! diagnostic variables:
real qn(nx, ny, nzm)
    !! cloud water+cloud ice
! fluxes at the top and bottom of the domain:
real precsfc(nx,ny)
    !! surface precip. rate
!  Horizontally varying stuff (as a function of xy)
real prec_xy(nx,ny)
    !! surface precipitation rate
! reference vertical profiles:
real rho(nzm)
    !! air density at pressure levels,kg/m3 

! Fields from beginning of time step used as NN inputs
real t_i(nx,ny,nzm)
    !! Temperature
real q_i(nx,ny,nzm)
    !! Non-precipitating water mixing ratio
real qp_i (nx,ny,nzm)
    !! Precipitating water mixing ratio

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine nn_convection_flux_SAM_init(nn_filename)
        !! Initialise the NN module for use in SAM

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights
        
        ! Get indices for grid section on this MPI rank
        call task_rank_to_index(rank,it,jt)

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename, n_inputs, n_outputs)

        ! Set init flag as complete
        do_init=.false.

        ! TODO Do we really need these variables, or are they defined through something else somewhere?
        ! Hardcoded magic numbers
        ! nrf is the fact results are supplied at lowest 30 half-model levels for sedimentation fluxes, and at 29 levels for fluxes
        ! (as flux at bottom boundary is zero).
        nrf = 30
        nrfq = 29
  
    end subroutine nn_convection_flux_SAM_init

    subroutine nn_convection_flux_SAM()
        ! Inputs?: t_i, q_i, qp_i
        !! Interface to the nn_convection for the SAM model
  
        ! real, intent(in), dimension(:, :, :)  :: tabs
            !! absolute temperature
        real, dimension(nx, ny, nzm)  :: tabs
            !! absolute temperature

        ! Local Variables
        real,   dimension(nrf) :: t_tendency_adv, q_tendency_adv, q_tendency_auto, q_tendency_sed, t_tendency_auto
        real,   dimension(nrf) :: q_flux_sed, qp_flux_fall, t_tendency_sed, q_tend_tot
        real,   dimension(nrf) :: t_flux_adv, q_flux_adv, t_sed_flux, t_rad_rest_tend, omp, fac ! Do not predict surface adv flux
        real,   dimension(nzm) :: qsat, irhoadz, irhoadzdz, irhoadzdz2

        ! NN parameters
        real(4), dimension(n_inputs) :: features
            !! Vector of input features for the NN
        real(4), dimension(n_outputs) :: outputs
            !! vector of output features from the NN

        real :: omn
            !! variable to store result of omegan function
        real :: rev_dz
            !! variable to store 1/dz factor
        integer  i, j, k,dd, dim_counter, out_dim_counter, out_var_control

        ! Check that we have initialised all of the variables.
        if (do_init) call error_mesg('NN has not yet been initialised using nn_convection_flux_init.')

        ! Initialise precipitation to 0 if required and at start of cycle
        if (.not. rf_uses_qp) then
            if(mod(nstep-1,nstatis).eq.0 .and. icycle.eq.1) then
                precsfc(:,:)=0.
            end if
        end if

        ! Define TODO Magical Mystery variables
        rev_dz = 1/dz
        do k=1,nzm
            irhoadz(k) = dtn/(rho(k)*adz(k)) ! Useful factor
            irhoadzdz(k) = irhoadz(k)/dz ! Note the time step
            irhoadzdz2(k) = irhoadz(k)/(dz*dz) ! Note the time step
        end do

        ! The NN operates on atmospheric columns, so loop over x and y coordinates in turn
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
                dim_counter = 0
                omp = 0.
                fac = 0.

                !-----------------------------------------------------
                ! Combine all features into one vector

                ! If using temperature then add as input feature
                if (Tin_feature_rf) then
                    features(dim_counter+1:dim_counter + input_ver_dim) = real(t_i(i,j,1:input_ver_dim),4)
                    dim_counter = dim_counter + input_ver_dim
                endif

                ! If using non-precipitating water mixing ratio then add as input feature using random forest (rf) approach from earlier paper
                if (qin_feature_rf) then
                    if (rf_uses_rh) then
                    ! If using generalised relative humidity convert non-precip. water content to rel. hum
                        do k=1,nzm
                            omn = omegan(tabs(i,j,k))
                            qsat(k) = omn * qsatw(tabs(i,j,k),pres(k)) + (1.-omn) * qsati(tabs(i,j,k),pres(k))
                        end do
                        features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim)/qsat(1:input_ver_dim),4)
                        dim_counter = dim_counter + input_ver_dim
                    else
                    ! if using non-precipitating water as water content
                        features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim),4)
                        dim_counter =  dim_counter + input_ver_dim
                    endif
                endif

                ! If using precipitating water (rain+snow) mixing ration content then add as input feature
                if (rf_uses_qp) then
                    features(dim_counter+1:dim_counter+input_ver_dim) = real(qp_i(i,j,1:input_ver_dim),4)
                    dim_counter = dim_counter + input_ver_dim
                endif

                ! If using distance to the equator then add as input feature
                ! y is a proxy for insolation and surface albedo as both are only a function of |y| in SAM
                ! TODO Will this be true in CAM/SCAM?
                if(do_yin_input) then
                    features(dim_counter+1) = real(abs(dy*(j+jt-(ny_gl+YES3D-1)/2-0.5)))
                    dim_counter = dim_counter+1
                endif

                !-----------------------------------------------------
                ! Call the forward method of the NN on the input features
                ! Scaling and normalisation done as layers in NN

                call net_forward(features, outputs, nrf)

                !-----------------------------------------------------
                ! Separate physical outputs from NN output vector and apply physical constraints

                ! TODO: The following code applies various physical constraints to the output of the NN
                ! These are likely to be required for any application, so should they be in the NN rather than the interface?
                ! If in the NN can they be generalised?

                ! Temperature rest tendency
                out_var_control =1
                t_rad_rest_tend(1:nrf) = outputs(1:nrf)
                out_dim_counter = nrf

                ! Temperature flux
                ! BC: advection surface flux is zero
                t_flux_adv(1) = 0.0
                out_var_control = out_var_control + 1
                t_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
                out_dim_counter = out_dim_counter + nrfq

                ! Non-precip. water flux
                ! BC: advection surface flux is zero
                q_flux_adv(1) = 0.0
                out_var_control = out_var_control + 1
                q_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
                out_dim_counter = out_dim_counter + nrfq

                ! Non-precip. water must be >= 0, so ensure advective flux will not reduce it below 0
                do k=2,nrf
                    if (q_flux_adv(k).lt.0) then
                        if ( q(i,j,k).lt.-q_flux_adv(k)* irhoadzdz(k)) then
                            q_flux_adv(k) = -q(i,j,k)/irhoadzdz(k)
                        end if
                    else
                        ! If gaining water content ensure that we are not gaining more than is in cell below
                        ! TODO: Not sure I fully understant this constraint
                        if (q(i,j,k-1).lt.q_flux_adv(k)* irhoadzdz(k)) then
                            q_flux_adv(k) = q(i,j,k-1)/irhoadzdz(k)
                        end if
                    end if
                end do

                ! Calculate tendencies
                out_var_control = out_var_control + 1
                q_tendency_auto(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
                out_dim_counter = out_dim_counter + nrf

                ! Convert flux to advective tendency via finite difference
                do k=1,nrf-1
                    t_tendency_adv(k) = - (t_flux_adv(k+1) - t_flux_adv(k)) * irhoadzdz(k)
                    q_tendency_adv(k) = - (q_flux_adv(k+1) - q_flux_adv(k)) * irhoadzdz(k)
                end do
                ! Apply finite difference boundary condition to advective tendencies
                t_tendency_adv(nrf) = - (0.0 - t_flux_adv(nrf)) * irhoadzdz(nrf)
                q_tendency_adv(nrf) = - (0.0 - q_flux_adv(nrf)) * irhoadzdz(nrf)
                ! q must be >= 0 so ensure tendency won't reduce below zero
                do k=1,nrf
                    if (q(i,j,k) .lt. -q_tendency_adv(k)) then
                        q_tendency_adv(k) = -q(i,j,k)
                    end if
                end do
                ! Apply advective tendencies to variables
                t(i,j,1:nrf) = t(i,j,1:nrf) + t_tendency_adv(1:nrf)
                q(i,j,1:nrf) = q(i,j,1:nrf) + q_tendency_adv(1:nrf)

                ! TODO ???
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

                ! Update q and t
                q(i,j,1:nrf) = q(i,j,1:nrf) + q_tend_tot(1:nrf)
                t(i,j,1:nrf) = t(i,j,1:nrf) - q_tend_tot(1:nrf)*fac(1:nrf)

                ! TODO ???
                out_var_control = out_var_control + 1
                q_flux_sed(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
                ! q_flux must be >= 0
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

                ! Calculate sed tendency via finite difference
                do k=1,nrf-1 ! One level less than I actually use
                    q_tendency_sed(k) = - (q_flux_sed(k+1) - q_flux_sed(k)) * irhoadzdz(k)
                end do
                ! Set value at top of nrf layr
                q_tendency_sed(nrf) = - (0.0 - q_flux_sed(nrf)) * irhoadzdz(nrf)

                ! If q sed tendency < 0 ensure it will not reduce q below 0
                do k=1,nrf
                    if (q_tendency_sed(k).lt.0) then
                        q_tendency_sed(k) = min(-q_tendency_sed(k), q(i,j,k))
                        q_tendency_sed(k) = -q_tendency_sed(k)
                    end if
                end do

                ! Apply sed tendency to variables
                t(i,j,1:nrf) = t(i,j,1:nrf) - q_tendency_sed(1:nrf)*(fac_fus+fac_cond)
                q(i,j,1:nrf) = q(i,j,1:nrf) +q_tendency_sed(1:nrf)

                ! Apply rest tendency to variables
                t(i,j,1:nrf) = t(i,j,1:nrf) + t_rad_rest_tend(1:nrf)*dtn

                ! Calculate surface precipitation
                ! Apply sed flux at surface
                precsfc(i,j) = precsfc(i,j)  - q_flux_sed(1)*dtn*rev_dz ! For statistics
                prec_xy(i,j) = prec_xy(i,j)  - q_flux_sed(1)*dtn*rev_dz ! For 2D output
                ! 
                do k=1, nrf
                    ! TODO: Both these lines have dz*(1/dz) can we just remove this?
                    precsfc(i,j) = precsfc(i,j) - q_tend_tot(k)*adz(k)*dz*rho(k)*rev_dz! removed the time step mult because q_tend_tot is already mult
                    prec_xy(i,j) = prec_xy(i,j) - q_tend_tot(k)*adz(k)*dz*rho(k)*rev_dz
                end do

                ! As a final ckeck q must be >= 0.0, if not then set to 0.0, otherwise add tendencies
                do k = 1,nrf
                    q(i,j,k) = max(0.,q(i,j,k))
                end do

                ! qn (precip.) must be >= 0.0, if not then set to 0.0, otherwise add tendencies
                where (qn(i,j,1:nrf) .gt. 0.0)
                    qn(i,j,1:nrf) = qn(i,j,1:nrf) + q_tend_tot(1:nrf) + q_tendency_adv(1:nrf) + q_tendency_sed(1:nrf)
                end where
                where (qn(i,j,:) .lt. 0.0)
                    qn(i,j,:) = 0.0
                end where

            end do
        end do
        ! End of loop over x, y, columns
  
    end subroutine nn_convection_flux_SAM
  
    !###################################################################
  
    subroutine error_mesg (message)
      character(len=*), intent(in) :: message
          !! message to be written to output   (character string)
  
      if(masterproc) print*, 'Neural network  module: ', message
      stop
  
    end subroutine error_mesg
 
    !###################################################################
  
    subroutine task_rank_to_index (rank,i,j)
  
      ! returns the pair of  beginning indices for the subdomain on the  
      ! global grid given the subdomain's rank.
      integer, intent(in) :: rank
          !! rank of MPI process
      integer :: i, j
          !! indices at which subdomain starts
  
      j = rank/nsubdomains_x 
      i = rank - j*nsubdomains_x
  
      i = i * (nx_gl/nsubdomains_x) 
      j = j * (ny_gl/nsubdomains_y) 
  
    end subroutine task_rank_to_index
  
    !###################################################################
    ! Need omegan function to run, ripped from https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/omega.f90
  
    real function omegan(tabs)

        real, intent(in) :: tabs
            !! Absolute temperature

        omegan = max(0., min(1., (tabs-tbgmin)*a_bg))

        return

    end function omegan
  
    !###################################################################
    ! Need qsatw functions to run, ripped from https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/sat.f90
  
    ! Saturation vapor pressure and mixing ratio. 
    ! Based on Flatau et.al, (JAM, 1992:1507)
  
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

end module nn_interface_SAM_mod
