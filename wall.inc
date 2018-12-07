      integer(4), parameter :: nnphi = 90, nstepp=1
      integer(4), parameter :: nntheta = 45, nstept=1
      real(8), parameter :: rmajour = 1.20d0, zmajour = 0.0d0

      real(8) :: rmnw(0:5,-10:10), zmnw(0:5,-10:10)
      real(8) ::  aawall(0:nnphi,0:nntheta),thwall(0:nnphi,0:nntheta)
      integer(4) :: mmm2,nnn2

      common /comwmd/ rmnw,zmnw,aawall,thwall

