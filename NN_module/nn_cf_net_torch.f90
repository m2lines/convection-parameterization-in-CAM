module nn_cf_net_torch
    !! neural_net convection emulator
    !! Module containing code pertaining only to the Neural Net functionalities of
    !! the M2LiNES Convection Flux parameterisation.

!---------------------------------------------------------------------
! Libraries to use
use netcdf
! Import our library for interfacing with PyTorch
use :: ftorch

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

! FTorch variables
type(torch_module) :: model
type(torch_tensor), dimension(1) :: in_tensor
type(torch_tensor) :: out_tensor

integer, parameter :: n_tensor_inputs = 1

integer, parameter :: in_dims = 1
integer :: in_layout(in_dims) = [1]
integer, parameter :: out_dims = 1
integer :: out_layout(out_dims) = [1]

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

        ! Create tensors
        in_tensor(1) = torch_tensor_from_array(features, in_layout, torch_kCPU)
        out_tensor = torch_tensor_from_array(logits, out_layout, torch_kCPU)

        ! Run forward pass
        call torch_module_forward(model, in_tensor, n_tensor_inputs, out_tensor)

        ! Clean up tensors
        call torch_tensor_delete(out_tensor)
        do i = 1, n_tensor_inputs
            call torch_tensor_delete(in_tensor(i))
        end do

    end subroutine nn_cf_net_torch_forward


    subroutine nn_cf_net_torch_init(nn_filename, n_inputs, n_outputs, nrf)
        !! Initialise the neural net

        integer, intent(out) :: n_inputs, n_outputs
        integer, intent(in) :: nrf
            !! number of atmospheric layers in each input

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid
        integer :: in_dimid, out_dimid

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

        ! Close the nc file
        call check( nf90_close(ncid))

        write(*, *) 'Finished reading NN regression file.'

        n_inputs = n_in
        n_outputs = n_out

        ! Load model
        model = torch_module_load("../../torch_nets/saved_model.pt")

    end subroutine nn_cf_net_torch_init


    subroutine nn_cf_net_torch_finalize()
        !! Clean up NN space by deallocating any arrays.

        !  Delete model
        call torch_module_delete(model)

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