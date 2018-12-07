        subroutine ScalarProd(a,b,c)

        implicit none

        real(8) :: a(3),b(3),c

        c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

        return
        end

        subroutine CrossProd(a,b,c)

        implicit none

        real(8) :: a(3),b(3),c(3)

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

        return
        end

        function abs3d(a)

        implicit none

        real(8) ::  a(3), abs3d

        abs3d = sqrt(a(1)*a(1) + a(2)*a(2) + a(3)*a(3))

        return
        end

