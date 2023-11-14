module nn_cf_net_torch
    !! neural_net convection emulator
    !! Module containing code pertaining only to the Neural Net functionalities of
    !! the M2LiNES Convection Flux parameterisation.

!---------------------------------------------------------------------
! Libraries to use
use netcdf
! Import our library for interfacing with PyTorch
use :: ftorch
use, intrinsic :: iso_c_binding, only: c_int, c_int64_t, c_loc, c_wp=>c_float

implicit none
private


!---------------------------------------------------------------------
! public interfaces
public nn_cf_net_torch_init, nn_cf_net_torch_forward, nn_cf_net_torch_finalize


!---------------------------------------------------------------------
! local/private data

! Neural Net Parameters
integer :: n_in
    !! Combined length of input features
integer :: n_out
    !! Number of output features
integer :: n_features_out
    !! Number of output features (variables)
integer, allocatable, dimension(:) :: feature_out_sizes
    !! Vector storing the length of each of the input features
integer :: n_lev
    !! number of atmospheric layers (model levels) output is supplied for

! Scale factors for inputs
real(4), allocatable, dimension(:)       :: xscale_mean
real(4), allocatable, dimension(:)       :: xscale_stnd

! Scale factors for outputs
real(4), allocatable, dimension(:)       :: yscale_mean
real(4), allocatable, dimension(:)       :: yscale_stnd

! FTorch variables
type(torch_module) :: model
type(torch_tensor), dimension(1) :: in_tensor
type(torch_tensor) :: out_tensor

integer(c_int), parameter :: n_tensor_inputs = 1

integer(c_int), parameter :: in_dims = 1
integer(c_int) :: in_layout(in_dims) = [1]
integer(c_int), parameter :: out_dims = 1
integer(c_int) :: out_layout(out_dims) = [1]

integer(c_int64_t) :: in_shape(in_dims)
integer(c_int64_t) :: out_shape(out_dims)

integer, parameter :: torch_wp = torch_kFloat32

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_cf_net_torch_forward(features, logits)
        !! Run forward method of the Neural Net.

        real(4), dimension(:), target :: features
            !! Vector of input features
        real(4), dimension(:), target, intent(out)  :: logits
            !! Output vector
        integer :: out_pos, feature_size, f, i
            !!

        real(c_wp), dimension(:), allocatable, target :: in_data
        real(c_wp), dimension(:), allocatable, target :: out_data

        ! Normalize features
        features = (features - xscale_mean) / xscale_stnd

        allocate(in_data(in_shape(1)))
        allocate(out_data(out_shape(1)))

        ! Cast from input type to C type required
        in_data = real(features, kind=c_wp)

        ! Create tensors
        in_tensor(1) = torch_tensor_from_blob(c_loc(in_data), in_dims, in_shape, torch_wp, torch_kCPU, in_layout)
        out_tensor = torch_tensor_from_blob(c_loc(out_data), out_dims, out_shape, torch_wp, torch_kCPU, out_layout)

        ! Run forward pass
        call torch_module_forward(model, in_tensor, n_tensor_inputs, out_tensor)

        ! Cast back from C type to Fortran type
        logits = real(out_data, kind=4)

        deallocate(in_data)
        deallocate(out_data)

        ! Apply scaling and normalisation of each output feature in logits
        out_pos = 0
        do f = 1, n_features_out
           feature_size = feature_out_sizes(f)
           logits(out_pos+1:out_pos+feature_size) = (logits(out_pos+1:out_pos+feature_size)*yscale_stnd(f)) + yscale_mean(f)
           out_pos = out_pos + feature_size
        end do

        ! Clean up tensors
        call torch_tensor_delete(out_tensor)
        do i = 1, n_tensor_inputs
            call torch_tensor_delete(in_tensor(i))
        end do

    end subroutine nn_cf_net_torch_forward


    subroutine nn_cf_net_torch_init(nn_filename, n_inputs, n_outputs, nrf)
        !! Initialise the neural net

        integer(c_int), intent(out) :: n_inputs, n_outputs
        integer, intent(in) :: nrf
            !! number of atmospheric layers in each input

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid
        integer :: in_dimid, out_dimid

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
        call check( nf90_inq_dimid(ncid, 'N_out', out_dimid))
        call check( nf90_inquire_dimension(ncid, out_dimid, len=n_out))
        call check( nf90_inq_dimid(ncid, 'N_out_dim', out_dimid))
        call check( nf90_inquire_dimension(ncid, out_dimid, len=n_features_out))

        print *, 'size of features to NN', n_in
        print *, 'size of outputs from NN', n_out

        ! Allocate arrays for NN based on sizes read from nc file
        allocate(xscale_mean(n_in))
        allocate(xscale_stnd(n_in))

        allocate(yscale_mean(n_features_out))
        allocate(yscale_stnd(n_features_out))
        allocate(feature_out_sizes(n_features_out))

        ! Read NN data in from nc file
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

        ! Load model
        model = torch_module_load("../../torch_nets/saved_model.pt")

        in_shape(1) = n_inputs
        out_shape(1) = n_outputs

    end subroutine nn_cf_net_torch_init


    subroutine nn_cf_net_torch_finalize()
        !! Clean up NN space by deallocating any arrays.

        !  Delete model
        call torch_module_delete(model)

        ! Deallocate arrays
        deallocate(xscale_mean, xscale_stnd, yscale_mean, yscale_stnd)
        deallocate(feature_out_sizes)

    end subroutine nn_cf_net_torch_finalize


    subroutine check(err_status)
        !! Check error status after netcdf call and print message for
        !! error codes.

        integer, intent(in) :: err_status
            !! error status from nf90 function

        if(err_status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(err_status))
        end if

    end subroutine check


end module nn_cf_net_torch
