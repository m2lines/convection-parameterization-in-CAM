program cam_profile_tests

  use cam_profile, only: read_cam_profile, read_cam_outputs
  use nn_interface_CAM, only: nn_convection_flux_CAM, nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize

  implicit none

  character(len=136) :: nn_file = "NN_weights_YOG_convection.nc"
  character(len=136) :: sounding_file = "resources/SAM_sounding.nc"
  character(len=136) :: cam_profile_file = "scam_gateIII_215_ZMmod.cam.h2.1974-08-30-00000.nc"
  character(len=136) :: cam_out_file = "scam_gateIII_215_ZMmod.cam.h1.1974-08-30-00000.nc"

  real(8), dimension(1,1,32,101) :: t_nc, qv_nc, qc_nc, qi_nc
  real(8), dimension(32) :: plev_nc
  real(8), dimension(33) :: pint_nc
  real(8), dimension(1,1,101) :: ps_nc
  
  real(8), dimension(1,32) :: t, qv, qc, qi
  real(8), dimension(1,32) :: plev
  real(8), dimension(1,33) :: pint
  real(8), dimension(1) :: ps, precsfc
  real(8), dimension(1,32) :: dqi, dqv, dqc, ds

  real(8), dimension(1,1,32,101) :: yogdt_nc, yogdq_nc, zmdt_nc, zmdq_nc
  real(8), dimension(1,32) :: yogdt, yogdq, zmdt, zmdq

  integer :: i

  ! Initialise the NN module
  call nn_convection_flux_CAM_init(nn_file, sounding_file)

  call read_cam_profile(cam_profile_file, t_nc, plev_nc, qv_nc, qc_nc, qi_nc, pint_nc, ps_nc)
  call read_cam_outputs(cam_out_file, yogdt_nc, yogdq_nc, zmdt_nc, zmdq_nc)

  i = 5

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
  
  ! Run tests
  call nn_convection_flux_CAM(plev(:,32:1:-1), pint(:,33:1:-1), ps, &
                              t(:,32:1:-1), qv(:,32:1:-1), qc(:,32:1:-1), qi(:,32:1:-1), &
                              1004D0, &
                              1200D0, &
                              1, 32, &
                              1, 1, 1, &
                              precsfc, &
                              dqi(:,32:1:-1), dqv(:,32:1:-1), dqc(:,32:1:-1), ds(:,32:1:-1))

  write(*,*) "YOGDT:"
  write(*,*) ds / 1004.0
  write(*,*) "YOGDT CAM:"
  write(*,*) yogdt
  write(*,*) "ZMDT:"
  write(*,*) zmdt

  ! Clean up
  call nn_convection_flux_CAM_finalize()

end program cam_profile_tests
