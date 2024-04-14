module nf_out
  !! module to write out data to NetCDF for analysis

!---------------------------------------------------------------------
! Libraries to use
use netcdf
use SAM_consts_mod, only: check, input_ver_dim, nrf, nrfq, dt_sam


implicit none
private

!---------------------------------------------------------------------
! public interfaces
public  nf_setup, nf_close, nf_write_sam, nf_write_cam, nf_write_scalar

!---------------------------------------------------------------------
! local/private data

integer, parameter :: cmode = 0  !! Overwrite existinf NC file
integer :: nfid  !! Filename index
integer :: n_sam_lev_id, n_sam_int_id, n_cam_lev_id, n_cam_int_id, n_surf_id  !! Dimension indices

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nf_setup(lev_sam, int_sam, lev_cam, int_cam)
        !! Create a new NetCDF file

        real(8), dimension(:), intent(in) :: lev_sam, int_sam, lev_cam, int_cam
        integer :: nlev_cam

        nlev_cam = size(lev_cam)

        if (size(lev_sam) /= nrf) then
            write(*,*) "length of lev_sam (",size(lev_sam),") != nrf (",nrf,"). Stopping"
            stop
        else if (size(int_sam) /= nrf+1) then
            write(*,*) "length of int_sam (",size(int_sam),") != nrf+1 (",nrf+1,"). Stopping"
            stop
        else if (size(lev_cam) /= nlev_cam) then
            write(*,*) "length of lev_cam (",size(lev_cam),") != nlev_cam (",nlev_cam,"). Stopping"
            stop
        else if (size(int_cam) /= nlev_cam+1) then
            write(*,*) "length of int_cam (",size(int_cam),") != nlev_cam+1 (",nlev_cam+1,"). Stopping"
            stop
        end if

        call check( NF90_CREATE('NN_test_output.nc', cmode, nfid) )

        call check( NF90_DEF_DIM(nfid, "n_sam_lev", nrf, n_sam_lev_id) )
        call check( NF90_DEF_DIM(nfid, "n_sam_int", nrf+1, n_sam_int_id) )
        call check( NF90_DEF_DIM(nfid, "n_cam_lev", nlev_cam, n_cam_lev_id) )
        call check( NF90_DEF_DIM(nfid, "n_cam_int", nlev_cam+1, n_cam_int_id) )
        call check( NF90_DEF_DIM(nfid, "n_surface", 1, n_surf_id) )

        call check( NF90_ENDDEF(nfid) )

    end subroutine nf_setup

    subroutine nf_close()
        !! Close the NetCDF file

        call check( NF90_CLOSE(nfid) )
    end subroutine nf_close

    subroutine nf_write_cam(var, var_name)
        !! Write a column from the CAM grid to file

        integer :: var_id
        character(*) :: var_name
        real(8), intent(in) :: var(:)

        ! Define variable
        call check( NF90_REDEF(nfid) )
        call check( NF90_DEF_VAR(nfid, var_name, NF90_DOUBLE, n_cam_lev_id, var_id) )
        call check( NF90_ENDDEF(nfid) )

        ! Assign values
        call check( NF90_PUT_VAR(nfid, var_id, var(:)) )

    end subroutine nf_write_cam

    subroutine nf_write_sam(var, var_name)
        !! Write a column from the CAM grid to file

        integer :: var_id
        character(*) :: var_name
        real(8), intent(in) :: var(:)

        ! Define variable
        call check( NF90_REDEF(nfid) )
        call check( NF90_DEF_VAR(nfid, var_name, NF90_DOUBLE, n_sam_lev_id, var_id) )
        call check( NF90_ENDDEF(nfid) )

        ! Assign values
        call check( NF90_PUT_VAR(nfid, var_id, var(:)) )

    end subroutine nf_write_sam

    subroutine nf_write_scalar(var, var_name)
        !! Write a column from the CAM grid to file

        integer :: var_id
        character(*) :: var_name
        real(8), intent(in) :: var

        ! Define variable
        call check( NF90_REDEF(nfid) )
        call check( NF90_DEF_VAR(nfid, var_name, NF90_DOUBLE, n_surf_id, var_id) )
        call check( NF90_ENDDEF(nfid) )

        ! Assign values
        call check( NF90_PUT_VAR(nfid, var_id, var) )

    end subroutine nf_write_scalar

end module nf_out
