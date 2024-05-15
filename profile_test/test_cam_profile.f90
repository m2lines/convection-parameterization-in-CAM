program cam_profile_tests

  use cam_profile, only: read_cam_profile, read_cam_outputs
  use nn_interface_CAM, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, &
                              nn_convection_flux_CAM_finalize, fetch_sam_data
  use nf_out, only: nf_setup, nf_close, nf_set_t, nf_write_sam, nf_write_cam, nf_write_scalar

  implicit none

  character(len=136) :: nn_file = "NN_weights_YOG_convection.nc"
  character(len=136) :: sounding_file = "resources/SAM_sounding.nc"
  character(len=136) :: cam_profile_file = "scam_gateIII_215_YOG_testing.cam.h2.1974-08-30-00000.nc"
  character(len=136) :: cam_out_file = "scam_gateIII_215_YOG_testing.cam.h1.1974-08-30-00000.nc"

  !> SAM grid data
  real(8), dimension(48) :: lev_sam, int_sam, gamazsam, rhosam, zsam

  !> SAM input and output data
  real(8), dimension(1,32) :: t, qv, qc, qi
  real(8), dimension(1,32) :: plev
  real(8), dimension(1,33) :: pint
  real(8), dimension(1) :: ps, precsfc
  real(8), dimension(1,32) :: dqi, dqv, dqc, ds

  !> CAM input data
  real(8), dimension(1,1,32,101) :: t_nc, qv_nc, qc_nc, qi_nc
  real(8), dimension(32) :: plev_nc
  real(8), dimension(33) :: pint_nc
  real(8), dimension(1,1,101) :: ps_nc

  !> Tendency outputs from parameterisation on CAM
  real(8), dimension(1,1,32,101) :: yogdt_nc, yogdq_nc, yogdqi_nc, yogdqc_nc
  real(8), dimension(1,1,32,101) :: zmdt_nc, zmdq_nc, zmdqi_nc, zmdqc_nc, zmdtevap_nc, zmdqevap_nc
  real(8), dimension(1,1,32,101) :: clubbdq, clubbds, clubbdqc, clubbdqi, clubbdtadj, clubbdqadj, clubbdqiadj
  real(8), dimension(1,1,32,101) :: mpdt, mpdq, mpdqc, mpdqi
  real(8), dimension(1,1,101) :: zmprec
  real(8), dimension(1,32) :: yogdt, yogdq, zmdt, zmdq, zmdqi, zmdqc

  integer :: i

  ! Initialise the NN module
  call nn_convection_flux_CAM_init(nn_file, sounding_file)
  call fetch_sam_data(lev_sam, int_sam, gamazsam, rhosam, zsam)

  ! Read in the data that was output from CAM h2 and h1 fincl files
  call read_cam_profile(cam_profile_file, t_nc, plev_nc, qv_nc, qc_nc, qi_nc, pint_nc, ps_nc)
  call read_cam_outputs(cam_out_file, yogdt_nc, yogdq_nc, yogdqi_nc, yogdqc_nc, &
                        zmdt_nc, zmdq_nc, zmdqi_nc, zmdqc_nc, zmdtevap_nc, zmdqevap_nc, &
                        clubbdq, clubbds, clubbdqc, clubbdqi, clubbdtadj, clubbdqadj, clubbdqiadj, &
                        mpdt, mpdq, mpdqc, mpdqi, &
                        zmprec)

  ! Prepare the netcdf file for writing output
  call nf_setup(lev_sam(1:30), int_sam(1:31), plev(1,:), pint(1,:))

  ! Loop over timesteps
  do i = 1, 101
      call nf_set_t(i)

      t = t_nc(1,:,:,i)
      qv = qv_nc(1,:,:,i)
      qc = qc_nc(1,:,:,i)
      qi = qi_nc(1,:,:,i)
      ! P saved as hPa so convert back to Pa
      pint(1,:) = 100.0 * pint_nc
      plev(1,:) = 100.0 * plev_nc
      ps = ps_nc(1,:,i)

      yogdt = yogdt_nc(1,:,:,i)
      yogdq = yogdq_nc(1,:,:,i)
      zmdt = zmdt_nc(1,:,:,i)
      zmdq = zmdq_nc(1,:,:,i)
      zmdqi = zmdqi_nc(1,:,:,i)
      zmdqc = zmdqc_nc(1,:,:,i)

      ! Flip arrays in output to match definition of pressure coordinate from low-high
      ! used in parameterisation (instead of high-low from CAM).
      call nf_write_cam(t(1,size(t):1:-1), "YOG_T_IN")
      call nf_write_cam(qi(1,size(qi):1:-1), "YOG_QV_IN")
      call nf_write_cam(qc(1,size(qc):1:-1), "YOG_QC_IN")
      call nf_write_cam(qv(1,size(qv):1:-1), "YOG_QI_IN")

      ! Run parameterisation
      call nn_convection_flux_CAM(plev(:,32:1:-1), pint(:,33:1:-1), ps, &
                                  t(:,32:1:-1), qv(:,32:1:-1), qc(:,32:1:-1), qi(:,32:1:-1), &
                                  1004D0, &
                                  1200D0, &
                                  1, 32, &
                                  1, 1, 1, &
                                  precsfc, &
                                  dqi(:,32:1:-1), dqv(:,32:1:-1), dqc(:,32:1:-1), ds(:,32:1:-1))

      ! Flip arrays in output to match definition of pressure coordinate from low-high
      ! used in parameterisation (instead of high-low from CAM).
      call nf_write_cam(dqv(1,size(dqv):1:-1), "YOGDQ")
      call nf_write_cam(dqi(1,size(dqi):1:-1), "YOGDQICE")
      call nf_write_cam(dqc(1,size(dqc):1:-1), "YOGDQCLD")
      call nf_write_cam(ds(1,size(ds):1:-1)/1004D0, "YOGDT")  ! Scale by cp=1004 to match ZMDT
      call nf_write_scalar(precsfc(1), "YOGPREC")

      ! Write out YOG tendencies from CAM for comparison
      call nf_write_cam(yogdq_nc(1,1,size(zmdq):1:-1,i), "YOGDQ_NC")
      call nf_write_cam(yogdqi_nc(1,1,size(zmdqi):1:-1,i), "YOGDQICE_NC")
      call nf_write_cam(yogdqc_nc(1,1,size(zmdqc):1:-1,i), "YOGDQCLD_NC")
      call nf_write_cam(yogdt_nc(1,1,size(zmdt):1:-1,i), "YOGDT_NC")

      ! Write out ZM tendencies for comparison
      call nf_write_cam(zmdq_nc(1,1,size(zmdq):1:-1,i), "ZMDQ")
      call nf_write_cam(zmdqi_nc(1,1,size(zmdqi):1:-1,i), "ZMDQICE")
      call nf_write_cam(zmdqc_nc(1,1,size(zmdqc):1:-1,i), "ZMDQCLD")
      call nf_write_cam(zmdt_nc(1,1,size(zmdt):1:-1,i), "ZMDT")
      call nf_write_cam(zmdtevap_nc(1,1,size(zmdt):1:-1,i), "ZMDTEVAP")
      call nf_write_cam(zmdqevap_nc(1,1,size(zmdt):1:-1,i), "ZMDQEVAP")
      call nf_write_scalar(zmprec(1,1,i), "ZMPREC")

      call nf_write_cam(qv(1,:) + dqv(1,:) * 1200D0, "YOG_QV_OUT")
      call nf_write_cam(t(1,:) + ds(1,:) * 1200D0 / 1004D0, "YOG_T_OUT")
      call nf_write_cam(qv(1,:) + zmdq(1,:) * 1200D0, "ZM_QV_OUT")
      call nf_write_cam(t(1,:) + zmdt(1,:) * 1200D0, "ZM_T_OUT")

      ! Write out CLUBB tendencies for comparison
      call nf_write_cam(clubbdq(1,1,size(zmdq):1:-1,i), "CLUBBDQ")
      call nf_write_cam(clubbds(1,1,size(zmdq):1:-1,i) / 1004D0, "CLUBBDT")
      call nf_write_cam(clubbdqc(1,1,size(zmdq):1:-1,i), "CLUBBDQCLD")
      call nf_write_cam(clubbdqi(1,1,size(zmdq):1:-1,i), "CLUBBDQICE")
      call nf_write_cam(clubbdtadj(1,1,size(zmdq):1:-1,i), "CLUBBDT_ADJ")
      call nf_write_cam(clubbdqadj(1,1,size(zmdq):1:-1,i), "CLUBBDQ_ADJ")
      call nf_write_cam(clubbdqiadj(1,1,size(zmdq):1:-1,i), "CLUBBDQICE_ADJ")

      ! Write out MPHYS tendencies for comparison
      call nf_write_cam(mpdt(1,1,size(zmdq):1:-1,i) / 1004D0, "MPDT")
      call nf_write_cam(mpdq(1,1,size(zmdq):1:-1,i), "MPDQ")
      call nf_write_cam(mpdqc(1,1,size(zmdq):1:-1,i), "MPDQCLD")
      call nf_write_cam(mpdqi(1,1,size(zmdq):1:-1,i), "MPDQICE")

  end do

  ! Clean up
  call nn_convection_flux_CAM_finalize()
  call nf_close()

end program cam_profile_tests
