! neural_net convection emulator

module nn_convection_flux_mod

!---------------------------------------------------------------------
! Libraries to use
use netcdf
implicit none
private

!---------------------------------------------------------------------
! public interfaces
public  relu, net_forward, nn_convection_flux_init

!---------------------------------------------------------------------
! local/private data

! Neural Net Parameters
integer :: n_in
    !! Input dim features
! Dimension of each hidden layer
integer :: n_h1
integer :: n_h2
integer :: n_h3
integer :: n_h4
integer :: n_out
integer :: out_var_dim
    !! number of output variables  (different types)

! Weights at each hidden layer of the NN
real(4), allocatable, dimension(:,:)     :: r_w1
real(4), allocatable, dimension(:,:)     :: r_w2
real(4), allocatable, dimension(:,:)     :: r_w3
real(4), allocatable, dimension(:,:)     :: r_w4
real(4), allocatable, dimension(:,:)     :: r_w5

! Biases at each hidden layer of the NN
real(4), allocatable, dimension(:)       :: r_b1
real(4), allocatable, dimension(:)       :: r_b2
real(4), allocatable, dimension(:)       :: r_b3
real(4), allocatable, dimension(:)       :: r_b4
real(4), allocatable, dimension(:)       :: r_b5

! Scale factors for inputs
real(4), allocatable, dimension(:)       :: xscale_mean
real(4), allocatable, dimension(:)       :: xscale_stnd

! vectors at each hidden layer of the NN
real(4), allocatable, dimension(:)       :: z1
real(4), allocatable, dimension(:)       :: z2
real(4), allocatable, dimension(:)       :: z3
real(4), allocatable, dimension(:)       :: z4

! Scale factors for outputs
real(4), allocatable, dimension(:)       :: yscale_mean
real(4), allocatable, dimension(:)       :: yscale_stnd


!---------------------------------------------------------------------
! Functions and Subroutines

contains

    subroutine relu(logits)
        !! Applies ReLU to a vector.
  
         real(4), dimension(:), intent(inout)   :: logits
             !! vector to which ReLU will be applied
    
         where (logits .lt. 0.0)  logits = 0.0
  
    end subroutine relu

    !#################################################################

    subroutine net_forward(features, logits, nrf)
        !! Run forward method of the Neural Net.
  
        real(4), dimension(:) :: features
            !! Vector of input features
        real(4), dimension(:), intent(out)  :: logits
            !! Output vector
        integer :: i
            !! Loop counter
        integer, intent(in) :: nrf
            !! number of layers to apply 
        integer :: out_dim_counter
    
        ! Initialise hidden layer values
        ! TODO Is this really neccessary...???
        z1 = 0.
        z2 = 0.
        z3 = 0.
        z4 = 0.
        
        !Normalize features
        features = (features - xscale_mean) / xscale_stnd

        z1 = matmul(features, r_w1) + r_b1
        call relu(z1)
    
        z2 = matmul(z1, r_w2) + r_b2
        call relu(z2)
    
        z3 = matmul(z2, r_w3) + r_b3
        call relu(z3)
    
        z4 = matmul(z3, r_w4) + r_b4
        call relu(z4)
    
        logits = matmul(z4, r_w5) + r_b5
  
        
        ! Apply scaling and normalisation of output features
        ! TODO: This is ugly, but not sure much can be done now.
        i =1
        logits(1:nrf) = (logits(1:nrf) * yscale_stnd(i))  +  yscale_mean(i)
        out_dim_counter = nrf

        i = i + 1
        logits(out_dim_counter+1:out_dim_counter+nrf-1) = (logits(out_dim_counter+1:out_dim_counter+nrf-1)* yscale_stnd(i))  +  yscale_mean(i)
        out_dim_counter = out_dim_counter + nrf-1

        i = i + 1
        logits(out_dim_counter+1:out_dim_counter+nrf-1) = (logits(out_dim_counter+1:out_dim_counter+nrf-1) * yscale_stnd(i))  +  yscale_mean(i)
        out_dim_counter = out_dim_counter + nrf-1

        i = i + 1
        logits(out_dim_counter+1:out_dim_counter+nrf) = (logits(out_dim_counter+1:out_dim_counter+nrf) * yscale_stnd(i))  +  yscale_mean(i)
        out_dim_counter = out_dim_counter + nrf

        i = i + 1
        logits(out_dim_counter+1:out_dim_counter+nrf) = (logits(out_dim_counter+1:out_dim_counter+nrf) * yscale_stnd(i))  +  yscale_mean(i)
        out_dim_counter = out_dim_counter + nrf

    end subroutine

    !#################################################################

    subroutine nn_convection_flux_init(nn_filename, n_inputs, n_outputs)
        !! Initialise the neural net

        integer  unit,io,ierr
        integer, intent(out) :: n_inputs, n_outputs
        
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
            !! NetCDF filename from which to read model weights
        
 
        !-------------allocate arrays and read data-------------------
        
        ! Open the file. NF90_NOWRITE tells netCDF we want read-only
        !   access
        ! Get the varid or dimid for each variable or dimension based
        !   on its name.
  
        call check( nf90_open(trim(nn_filename),NF90_NOWRITE,ncid ))
  
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
        call check( nf90_inquire_dimension(ncid, out_dimid, len=out_var_dim))
  
        print *, 'size of features to NN', n_in
        print *, 'size of outputs from NN', n_out
  
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
  
        allocate(xscale_mean(n_in))
        allocate(xscale_stnd(n_in))
  
        allocate(yscale_mean(out_var_dim))
        allocate(yscale_stnd(out_var_dim))
  
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

        n_inputs = n_in
        n_outputs = n_out
  
    end subroutine nn_convection_flux_init

    !#################################################################
  
    subroutine check(status) 
        !! check error status after netcdf call and print message for
        !! error codes.
  
        integer, intent(in) :: status
            !! error status from nf90 function
  
        if(status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(status))
        end if
    end subroutine check
  
end module nn_convection_flux_mod
