!#####################################################################
      subroutine readwmode
!#####################################################################
!
!   readwmode read fourier spectrum of the wall data
!   read Fourier mode (rmn,zmn) for wall
!   and scoreing wall data
!   rwall,zwall,rmajour,zmajour should be in MKS units

      implicit none

      include 'wall.inc' 

      real(8) :: rwall, zwall, xwall, ywall
      real(8) :: rr2, zz2, tmpth
      real(8) :: pi
      real(8) :: phi,theta, rread, zread

      integer(4) :: iphi,jtheta
      integer(4) :: nnn, i

! init rwall,zwall
      rwall = 0.0d0
      zwall = 0.0d0
      pi = atan(1.0d0)*4.0d0

      open(32,file='./wmode.dat')
      open(28,file='./wall3d.txt')
! read wmode
      do i=1,116
        read(32,*)mmm2,nnn2,rread,zread
        rmnw(mmm2,nnn2)=rread
        zmnw(mmm2,nnn2)=zread
        write(6,'(2(1x,1pe10.3))') rread,zread
      enddo

      mmm2=0
      nnn2=-10
      do nnn=-10,-1,1
        rmnw(mmm2,nnn)=0.0d0
        zmnw(mmm2,nnn)=0.0d0
      enddo

! calculate (rwall,zwall) 
      do iphi = 0, nnphi
        phi = real(iphi)/real(nnphi)*2.0d0*pi
        do jtheta = 0, nntheta
          theta = real(jtheta)/real(nntheta)*2.0d0*pi
          rwall = 0.0d0
          zwall = 0.0d0
          do mmm2=0,5
            do nnn2=-10,10
              rwall=rwall+ 
     &          rmnw(mmm2,nnn2)*cos(real(mmm2)*theta-4.0*real(nnn2)*phi)
              zwall=zwall+ 
     &          zmnw(mmm2,nnn2)*sin(real(mmm2)*theta-4.0*real(nnn2)*phi)
            end do
          end do

! calculate and score radius and thwall at (phiv,thetav) in torus coordinate

          rr2 = rwall - rmajour
          zz2 = zwall - zmajour

          aawall(iphi,jtheta) = sqrt(rr2*rr2 + zz2*zz2)
!          call findatan2(rr2,zz2,tmpth)
          thwall(iphi,jtheta) = datan2(zz2,rr2)

! output wall position
          if(mod(iphi,nstepp).eq.0.and.mod(jtheta,nstept).eq.0) then
            xwall = rwall * cos(phi)
            ywall = rwall * sin(phi)
            rr2 = rwall - rmajour
            zz2 = zwall - zmajour
            write(28,'(3(1pe12.4))') xwall, ywall, zwall
!            write(28,'(10(1pe12.4),2i4)') xwall, ywall, zwall, rwall, 
!     &      phi, aawall(iphi,jtheta), thwall(iphi,jtheta), rr2, zz2, 
!     &       datan2(zz2,rr2), iphi, jtheta
          endif
          if(mod(iphi,nstepp).eq.0.and.jtheta.eq.nntheta) write(28,*)

        enddo ! jtheta
      enddo ! iphi

      return
      end

!#####################################################################
      subroutine findwall(rpos,zpos,ppos,row,lout)
!#####################################################################
!
!   findwall is calculate wall position of Heliotron J device
!     taking into account its 3D shape
!       using fourier spectrum of the wall data
!
      implicit none

      include 'wall.inc'
! input
      real(8) :: rpos, zpos, ppos
! outout
      real(8) :: row
      integer(4) :: lout

! local
      real(8) :: rr3, zz3, ppos2, aa, aa2
      real(8) :: pi
      real(8) :: phi,theta

      integer(4) :: iphi2,mtheta

      real(8) :: dtheta2(0:nntheta)        ! local

! init lout
      lout = 0
      pi = atan(1.0d0)*4.0d0

      rr3 = rpos - rmajour
      zz3 = zpos - zmajour
      aa = sqrt(rr3*rr3 + zz3*zz3)
!      call findatan2(rr3,zz3,theta)
      theta = datan2(zz3,rr3)

! find aawall,thwall
      ppos2 = mod(ppos,2.0d0*pi)
      if(ppos2.lt.0.0d0) ppos2 = ppos2 + 2.0d0*pi

      iphi2 = int(ppos2/2.0d0/pi*nnphi)

!       do jtheta2=0,nntheta
!        dtheta2(jtheta2)=abs(theta - thwall(iphi2,jtheta2))
!      enddo
      dtheta2(0:nntheta)=abs(theta - thwall(iphi2,0:nntheta))

      mtheta = minloc(dtheta2,1) - 1

      aa2 = aawall(iphi2,mtheta)
      row = aa/aa2

      if(aa >= aa2) then
         lout = 1
         return
      end if

      return
      end

