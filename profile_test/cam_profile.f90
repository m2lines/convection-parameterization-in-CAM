module cam_profile
  !! Module containing

  !--------------------------------------------------------------------------
  ! Libraries to use
  use netcdf

  implicit none
  private

  !---------------------------------------------------------------------
  ! public interfaces
  public read_cam_profile, read_cam_outputs


  contains

    subroutine read_cam_profile(cam_filename, t, plev, qv, qc, qi, pint, ps)

        character(len=136), intent(in) :: cam_filename
        
        real(8), dimension(1,1,32,101) :: t, qv, qc, qi
        real(8), dimension(32) :: plev
        real(8), dimension(33) :: pint
        real(8), dimension(1,1,101) :: ps


        integer :: nx = 1
        integer :: nlev = 32

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid
        integer :: t_id, plev_id, pint_id, ps_id
        integer :: qv_id, qc_id, qi_id

        !-------------allocate arrays and read data-------------------

        call check( nf90_open(trim(cam_filename),NF90_NOWRITE,ncid ))

        ! Read data in from nc file
        call check( nf90_inq_varid(ncid, "T", t_id))
!        write(*,*) "T id = ", t_id
        call check( nf90_get_var(ncid, t_id, t))
!        write(*,*) "T read in as ", t
        call check( nf90_inq_varid(ncid, "Q", qv_id))
!        write(*,*) "Q id = ", qv_id
        call check( nf90_get_var(ncid, qv_id, qv))
        call check( nf90_inq_varid(ncid, "CLDLIQ", qc_id))
!        write(*,*) "QC id = ", qc_id
        call check( nf90_get_var(ncid, qc_id, qc))
        call check( nf90_inq_varid(ncid, "CLDICE", qi_id))
        call check( nf90_get_var(ncid, qi_id, qi))

        call check( nf90_inq_varid(ncid, "lev", plev_id))
        call check( nf90_get_var(ncid, plev_id, plev))
        call check( nf90_inq_varid(ncid, "ilev", pint_id))
        call check( nf90_get_var(ncid, pint_id, pint))
        call check( nf90_inq_varid(ncid, "PS", ps_id))
        call check( nf90_get_var(ncid, ps_id, ps))

        ! Close the nc file
        call check( nf90_close(ncid))

        write(*, *) 'Finished reading profile file.'

    end subroutine read_cam_profile

    subroutine read_cam_outputs(cam_out_file, yogdt_nc, yogdq_nc, yogdqi_nc, yogdqc_nc, &
                                zmdt_nc, zmdq_nc, zmdqi_nc, zmdqc_nc,zmdtevap_nc, zmdqevap_nc, &
                                clubbdq, clubbds, clubbdqc, clubbdqi, clubbdtadj, clubbdqadj, clubbdqiadj, &
                                mpdt, mpdq, mpdqc, mpdqi, &
                                prec)

        character(len=136), intent(in) :: cam_out_file
        
        real(8), dimension(1,1,32,101) :: yogdt_nc, yogdq_nc, yogdqi_nc, yogdqc_nc
        real(8), dimension(1,1,32,101) :: zmdt_nc, zmdq_nc, zmdqi_nc, zmdqc_nc, zmdtevap_nc, zmdqevap_nc
        real(8), dimension(1,1,101) :: prec
        real(8), dimension(1,1,32,101) :: clubbdq, clubbds, clubbdqc, clubbdqi, clubbdtadj, clubbdqadj, clubbdqiadj
        real(8), dimension(1,1,32,101) :: mpdt, mpdq, mpdqc, mpdqi

        integer :: nx = 1
        integer :: nlev = 32

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid
        integer :: yogdt_id, yogdq_id, yogdqi_id, yogdqc_id, prec_id
        integer :: zmdt_id, zmdq_id, zmdqi_id, zmdqc_id, zmdtevap_id, zmdqevap_id
        integer :: clubbdq_id, clubbds_id, clubbdqc_id, clubbdqi_id, clubbdtadj_id, clubbdqadj_id, clubbdqiadj_id
        integer :: mpdt_id, mpdq_id, mpdqc_id, mpdqi_id

        !-------------allocate arrays and read data-------------------

        call check( nf90_open(trim(cam_out_file),NF90_NOWRITE,ncid ))

        ! Read data in from nc file
        call check( nf90_inq_varid(ncid, "YOGDT", yogdt_id))
        call check( nf90_get_var(ncid, yogdt_id, yogdt_nc))
        call check( nf90_inq_varid(ncid, "YOGDQ", yogdq_id))
        call check( nf90_get_var(ncid, yogdq_id, yogdq_nc))
        call check( nf90_inq_varid(ncid, "YOGDICE", yogdqi_id))
        call check( nf90_get_var(ncid, yogdqi_id, yogdqi_nc))
        call check( nf90_inq_varid(ncid, "YOGDLIQ", yogdqc_id))
        call check( nf90_get_var(ncid, yogdqc_id, yogdqc_nc))
        call check( nf90_inq_varid(ncid, "ZMDT", zmdt_id))
        call check( nf90_get_var(ncid, zmdt_id, zmdt_nc))
        call check( nf90_inq_varid(ncid, "ZMDQ", zmdq_id))
        call check( nf90_get_var(ncid, zmdq_id, zmdq_nc))
        call check( nf90_inq_varid(ncid, "ZMDICE", zmdqi_id))
        call check( nf90_get_var(ncid, zmdqi_id, zmdqi_nc))
        call check( nf90_inq_varid(ncid, "ZMDLIQ", zmdqc_id))
        call check( nf90_get_var(ncid, zmdqc_id, zmdqc_nc))
        call check( nf90_inq_varid(ncid, "PRECC", prec_id))
        call check( nf90_get_var(ncid, prec_id, prec))

        call check( nf90_inq_varid(ncid, "EVAPTZM", zmdtevap_id))
        call check( nf90_get_var(ncid, zmdtevap_id, zmdtevap_nc))
        call check( nf90_inq_varid(ncid, "EVAPQZM", zmdqevap_id))
        call check( nf90_get_var(ncid, zmdqevap_id, zmdqevap_nc))

        call check( nf90_inq_varid(ncid, "RVMTEND_CLUBB", clubbdq_id))
        call check( nf90_get_var(ncid, clubbdq_id, clubbdq))
        call check( nf90_inq_varid(ncid, "STEND_CLUBB", clubbds_id))
        call check( nf90_get_var(ncid, clubbds_id, clubbds))
        call check( nf90_inq_varid(ncid, "RCMTEND_CLUBB", clubbdqc_id))
        call check( nf90_get_var(ncid, clubbdqc_id, clubbdqc))
        call check( nf90_inq_varid(ncid, "RIMTEND_CLUBB", clubbdqi_id))
        call check( nf90_get_var(ncid, clubbdqi_id, clubbdqi))
        call check( nf90_inq_varid(ncid, "TTENDICE", clubbdtadj_id))
        call check( nf90_get_var(ncid, clubbdtadj_id, clubbdtadj))
        call check( nf90_inq_varid(ncid, "QVTENDICE", clubbdqadj_id))
        call check( nf90_get_var(ncid, clubbdqadj_id, clubbdqadj))
        call check( nf90_inq_varid(ncid, "QCTENDICE", clubbdqiadj_id))
        call check( nf90_get_var(ncid, clubbdqiadj_id, clubbdqiadj))

        call check( nf90_inq_varid(ncid, "MPDT", mpdt_id))
        call check( nf90_get_var(ncid, mpdt_id, mpdt))
        call check( nf90_inq_varid(ncid, "MPDQ", mpdq_id))
        call check( nf90_get_var(ncid, mpdq_id, mpdq))
        call check( nf90_inq_varid(ncid, "MPDLIQ", mpdqc_id))
        call check( nf90_get_var(ncid, mpdqc_id, mpdqc))
        call check( nf90_inq_varid(ncid, "MPDICE", mpdqi_id))
        call check( nf90_get_var(ncid, mpdqi_id, mpdqi))

        ! Close the nc file
        call check( nf90_close(ncid))

        write(*, *) 'Finished reading outputs file.'

    end subroutine read_cam_outputs

    subroutine check(err_status)
        !! Check error status after netcdf call and print message for
        !! error codes.

        integer, intent(in) :: err_status
            !! error status from nf90 function

        if(err_status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(err_status))
        end if

    end subroutine check

end module cam_profile

