        implicit none

        include 'wall.inc'

        real(8) :: x1(3),x2(3),x3(3),x4(3),x5(3),x6(3)
        real(8) :: v12(3),v13(3),v45(3),v46(3),vn123(3),vn456(3)
        real(8) :: p1,p2,p3,p4,p5,p6
        real(8) :: r1,r2,r3,r4,r5,r6
        real(8) :: z1,z2,z3,z4,z5,z6
        real(8) :: a1,a2,a3,a4,a5,a6
        real(8) :: t1,t2,t3,t4,t5,t6
        real(8) :: pi
        integer(4) :: i,j,k,l,m
        integer(4) :: lchk

        pi=4.0d0*atan(1.0d0)

        call readwmode('./wmode.dat')

        open(16,file='walltest.stl')
        write(16,*)'solid wallp1440t720'

! notice
        write(6,*) ' coordinate data (phi=0 for straight section)'
        write(6,*) ' output data unit in mm'
        write(6,*) ' OK? (Y for 0)'
        read(5,*) lchk
        if(lchk.ne.0) stop

        do i=0,nnphi-1
          do j=0,nntheta-1
! pi/4 is added to phi
            p1=real(i)/real(nnphi)*2.0d0*pi + pi/4.0d0
            p2=real(i+1)/real(nnphi)*2.0d0*pi + pi/4.0d0
            p3=real(i)/real(nnphi)*2.0d0*pi + pi/4.0d0
            p4=real(i+1)/real(nnphi)*2.0d0*pi + pi/4.0d0
            p5=real(i+1)/real(nnphi)*2.0d0*pi + pi/4.0d0
            p6=real(i)/real(nnphi)*2.0d0*pi + pi/4.0d0
            a1=aawall(i,j)
            a2=aawall(i+1,j)
            a3=aawall(i,j+1)
            a4=aawall(i+1,j)
            a5=aawall(i+1,j+1)
            a6=aawall(i,j+1)
            t1=thwall(i,j)
            t2=thwall(i+1,j)
            t3=thwall(i,j+1)
            t4=thwall(i+1,j)
            t5=thwall(i+1,j+1)
            t6=thwall(i,j+1)
            r1=a1*cos(t1)+rmajour  ;  z1=a1*sin(t1)+zmajour
            r2=a2*cos(t2)+rmajour  ;  z2=a2*sin(t2)+zmajour
            r3=a3*cos(t3)+rmajour  ;  z3=a3*sin(t3)+zmajour
            r4=a4*cos(t4)+rmajour  ;  z4=a4*sin(t4)+zmajour
            r5=a5*cos(t5)+rmajour  ;  z5=a5*sin(t5)+zmajour
            r6=a6*cos(t6)+rmajour  ;  z6=a6*sin(t6)+zmajour
            x1(1)=r1*cos(p1) ; x1(2)=r1*sin(p1) ; x1(3)=z1
            x2(1)=r2*cos(p2) ; x2(2)=r2*sin(p2) ; x2(3)=z2
            x3(1)=r3*cos(p3) ; x3(2)=r3*sin(p3) ; x3(3)=z3
            x4(1)=r4*cos(p4) ; x4(2)=r4*sin(p4) ; x4(3)=z4
            x5(1)=r5*cos(p5) ; x5(2)=r5*sin(p5) ; x5(3)=z5
            x6(1)=r6*cos(p6) ; x6(2)=r6*sin(p6) ; x6(3)=z6
! m -> mm
            x1=x1*1.0e3
            x2=x2*1.0e3
            x3=x3*1.0e3
            x4=x4*1.0e3
            x5=x5*1.0e3
            x6=x6*1.0e3

! set vector
            v12=x2-x1
            v13=x3-x1
            v45=x5-x4
            v46=x6-x4

! set normal vector
            call crossprod(v12,v13,vn123)
            call crossprod(v45,v46,vn456)

! output
            write(16,'(a14,3(1x,1pe14.6))')'  facet normal',vn123
            write(16,'(a14)')'    outer loop'
            write(16,'(a12,3(1x,1pe14.6))')'      vertex',x1
            write(16,'(a12,3(1x,1pe14.6))')'      vertex',x2
            write(16,'(a12,3(1x,1pe14.6))')'      vertex',x3
            write(16,'(a11)')'    endloop'
            write(16,'(a10)')'  endfacet'
            write(16,'(a14,3(1x,1pe14.6))')'  facet normal',vn456
            write(16,'(a14)')'    outer loop'
            write(16,'(a12,3(1x,1pe14.6))')'      vertex',x4
            write(16,'(a12,3(1x,1pe14.6))')'      vertex',x5
            write(16,'(a12,3(1x,1pe14.6))')'      vertex',x6
            write(16,'(a11)')'    endloop'
            write(16,'(a10)')'  endfacet'

          enddo         !j
        enddo           !i

        write(16,*)'endsolid wallp1440t720'
        close(16)

        stop
        end
