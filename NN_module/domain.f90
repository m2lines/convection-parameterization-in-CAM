! Set the domain dimensionality, size and number of subdomains.

module domain

       integer, parameter :: YES3D = 1  ! Domain dimensionality: 1 - 3D, 0 - 2D
! wrb: settings for wide domain 5 km res run:
!       integer, parameter :: nx_gl = 576 ! Number of grid points in X
!       integer, parameter :: ny_gl = 1440 ! Number of grid points in Y
!       integer, parameter :: nz_gl = 48 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 12 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 24 ! No of subdomains in y

!       integer, parameter :: nx_gl = 64 ! Number of grid points in X
!       integer, parameter :: ny_gl = 160 ! Number of grid points in Y
!       integer, parameter :: nz_gl = 48 ! Number of pressure (scalar) levels
!       integer, parameter :: nsubdomains_x  = 2 ! No of subdomains in x
!       integer, parameter :: nsubdomains_y  = 4 ! No of subdomains in y

       integer, parameter :: nx_gl = 72 ! Number of grid points in X - Yani changed to 36 from 32
       integer, parameter :: ny_gl = 180 ! Number of grid points in Y
       integer, parameter :: nz_gl = 48 ! Number of pressure (scalar) levels
       integer, parameter :: nsubdomains_x  = 3 ! No of subdomains in x
       integer, parameter :: nsubdomains_y  = 12 ! No of subdomains in y

! Note:
!  * nx_gl and ny_gl should be a factor of 2,3, or 5 (see User's Guide)
!  * if 2D case, ny_gl = nsubdomains_y = 1 ;
!  * nsubdomains_x*nsubdomains_y = total number of processors
!  * if one processor is used, than  nsubdomains_x = nsubdomains_y = 1;

end module domain
