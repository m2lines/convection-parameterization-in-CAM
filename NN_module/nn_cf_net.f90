module nn_cf_net_mod
    !! neural_net convection emulator
    !! Module containing code pertaining only to the Neural Net functionalities of
    !! the M2LiNES Convection Flux parameterisation.

!---------------------------------------------------------------------
! Libraries to use
use netcdf
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  relu, nn_cf_net_forward, nn_cf_net_init, nn_cf_net_finalize


!---------------------------------------------------------------------
! local/private data

! Neural Net Parameters
integer :: n_in
    !! Combined length of input features
integer :: n_features_out
    !! Number of output features (variables)
integer, allocatable, dimension(:) :: feature_out_sizes
    !! Vector storing the length of each of the input features
integer :: n_lev
    !! number of atmospheric layers (model levels) output is supplied for
! Dimension of each hidden layer
integer :: n_h1
integer :: n_h2
integer :: n_h3
integer :: n_h4
integer :: n_out

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

    !-----------------------------------------------------------------
    ! Public Subroutines
    
    subroutine relu(logits)
        !! Applies ReLU to a vector.

         real(4), dimension(:), intent(inout)   :: logits
             !! vector to which ReLU will be applied

         where (logits .lt. 0.0)  logits = 0.0

    end subroutine relu


    subroutine nn_cf_net_forward(features, logits)
        !! Run forward method of the Neural Net.

        real(4), dimension(:), target :: features
            !! Vector of input features
        real(4), dimension(:), target, intent(out)  :: logits
            !! Output vector
        integer :: i
            !! Loop counter
        integer :: out_dim_counter
            !!
        integer :: out_pos, feature_size, f
            !! 

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

        ! Apply scaling and normalisation of each output feature in logits
        out_pos = 0
        do f = 1, n_features_out
           feature_size = feature_out_sizes(f)
           logits(out_pos+1:out_pos+feature_size) = (logits(out_pos+1:out_pos+feature_size)*yscale_stnd(f)) + yscale_mean(f)
           out_pos = out_pos + feature_size
        end do

    end subroutine


    subroutine nn_cf_net_init(nn_filename, n_inputs, n_outputs, nrf)
        !! Initialise the neural net

        integer, intent(out) :: n_inputs, n_outputs
        integer, intent(in) :: nrf
            !! number of atmospheric layers in each input

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid
        integer :: in_dimid, out_dimid, single_dimid
        integer :: h1_dimid, h2_dimid, h3_dimid, h4_dimid
        integer :: r_w1_varid, r_w2_varid, r_b1_varid, r_b2_varid
        integer :: r_w3_varid, r_w4_varid, r_b3_varid, r_b4_varid
        integer :: r_w5_varid, r_b5_varid

        integer :: xscale_mean_varid, xscale_stnd_varid
        integer :: yscale_mean_varid, yscale_stnd_varid

        character(len=1024), intent(in) :: nn_filename
            !! NetCDF filename from which to read model weights


        !-------------allocate arrays and read data-------------------

        ! Open the file. NF90_NOWRITE tells netCDF we want read-only
        ! access
        ! Get the varid or dimid for each variable or dimension based
        ! on its name.

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
        call check( nf90_inquire_dimension(ncid, out_dimid, len=n_features_out))

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

        allocate(yscale_mean(n_features_out))
        allocate(yscale_stnd(n_features_out))
        allocate(feature_out_sizes(n_features_out))

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

        ! Set sizes of output features based on nrf
        n_lev = nrf
        
        ! Set sizes of the feature groups
        feature_out_sizes(:)   = n_lev
        feature_out_sizes(2:3) = n_lev-1

    end subroutine nn_cf_net_init


    subroutine nn_cf_net_finalize()
        !! Clean up NN space by deallocating any arrays.

        ! Deallocate arrays
        deallocate(r_w1, r_w2, r_w3, r_w4, r_w5)
        deallocate(r_b1, r_b2, r_b3, r_b4, r_b5)
        deallocate(z1, z2, z3, z4)
        deallocate(xscale_mean, xscale_stnd, yscale_mean, yscale_stnd)
        deallocate(feature_out_sizes)

    end subroutine nn_cf_net_finalize


    subroutine check(err_status)
        !! Check error status after netcdf call and print message for
        !! error codes.

        integer, intent(in) :: err_status
            !! error status from nf90 function

        if(err_status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(err_status))
        end if

    end subroutine check


end module nn_cf_net_mod
