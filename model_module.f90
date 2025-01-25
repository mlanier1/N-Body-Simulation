module model
    implicit none 
    double precision :: G = 6.67430E-20 ! gravitational constant
    integer, parameter :: num_bodies = 65, num_steps = 500000

contains

! Kinetic energy, potential energy, angular momentum components
! Input: m, x, y, z, vx, vy, vz, bodies
! Output: lx, ly, lz, kinetic, tot_potential, tot_kinetic, total_lx, total_ly, total_lz
subroutine momentum_calc(m, x, y, z, vx, vy, vz, lx, ly, lz, kinetic, &
    tot_potential, tot_kinetic, total_lx, total_ly, total_lz, bodies)
    integer, intent(in) :: bodies
    double precision, dimension(:), intent(in) :: m, x, y, z, vx, vy, vz
    double precision, dimension(:) :: lx, ly, lz, kinetic  ! Corrected array size
    double precision, intent(out) :: tot_potential, tot_kinetic, total_lx, total_ly, total_lz
    double precision :: dx, dy, dz, Euclidean
    integer :: i, j

    tot_potential = 0.0d0 ! total potential energy
    tot_kinetic = 0.0d0 ! total kinetic
    total_lx = 0.0 ! total lx
    total_ly = 0.0 ! total ly
    total_lz = 0.0 ! total lz
    lx = 0.0
    ly = 0.0
    lz = 0.0

    do i = 1, bodies ! Angular momentum and kinetic energy
        lx(i) = m(i) * (y(i) * vz(i) - z(i) * vy(i)) ! Comes from cross product
        ly(i) = m(i) * (z(i) * vx(i) - x(i) * vz(i))
        lz(i) = m(i) * (x(i) * vy(i) - y(i) * vx(i))
        total_lx = total_lx + lx(i)
        total_ly = total_ly + ly(i)
        total_lz = total_lz + lz(i)

        kinetic(i) = 0.5d0 * m(i) * (vx(i)**2 + vy(i)**2 + vz(i)**2) ! Kinetic formula including all three components
        tot_kinetic = tot_kinetic + kinetic(i)

        do j = i + 1, bodies ! Potential energy
            dx = x(j) - x(i) ! Distance between x coordinates
            dy = y(j) - y(i) ! Distance between y coords
            dz = z(j) - z(i) ! Distance between z coords
            Euclidean = (dx**2 + dy**2 + dz**2)**0.5 ! Euclidean norm
            tot_potential = tot_potential + ((-G * m(i) * m(j))) / Euclidean ! Total potential energy
        end do
    end do
end subroutine


    ! Acceleration due to gravity
    ! Input: bodies, x, y, z, m
    ! Output: ax, ay, az
    subroutine acceleration(bodies, x, y, z, m, ax, ay, az)
        integer, intent(in) :: bodies
        double precision, dimension(:), intent(in) :: x, y, z, m
        double precision, dimension(:) :: ax, ay, az
        double precision :: dx, dy, dz, Euclidean
        integer :: i, j

        ax = 0 ! Acceleration in x dir
        ay = 0 ! Acceleration in y dir
        az = 0 ! Acceleration in z dir

        do i = 1, bodies ! Iterates over all bodies in the system
            do j = 1, i - 1
                dx = x(j) - x(i) ! Distance between x coordinates
                dy = y(j) - y(i) ! Distance between y coords
                dz = z(j) - z(i) ! Distance between z coords

                Euclidean = (dx**2 + dy**2 + dz**2)**0.5 ! Euclidean norm
                ! Acceleration for j -> j pulled toward i
                ax(j) = -((G * m(i) * dx) / (Euclidean**3)) + ax(j)
                ay(j) = -((G * m(i) * dy) / (Euclidean**3)) + ay(j)
                az(j) = -((G * m(i) * dz) / (Euclidean**3)) + az(j)
                ! Acceleration for i -> i pulled toward j
                ax(i) = ((G * m(j) * dx) / (Euclidean**3)) + ax(i)
                ay(i) = ((G * m(j) * dy) / (Euclidean**3)) + ay(i)
                az(i) = ((G * m(j) * dz) / (Euclidean**3)) + az(i)
            end do
        end do
    end subroutine

    ! Verlet integration
    ! Input: bodies, steps, t, m, x, y, z, vx, vy, vz, ax, ay, az
    ! Output: updated x, y, z, vx, vy, vz
    subroutine Verlet(bodies, steps, t, m, x, y, z, vx, vy, vz, ax, ay, az)
        integer, intent(in) :: bodies, steps
        double precision, intent(in) :: t
        double precision, dimension(:), intent(inout) :: m
        double precision, dimension(:,:), intent(inout) :: x, y, z, vx, vy, vz, ax, ay, az
        integer :: i, j

        call acceleration(bodies, x(1,:), y(1,:), z(1,:), m,  ax(1,:), ay(1,:), az(1,:))

        ! Euler initial step
        do i = 1, bodies
            vx(2, i) = vx(1,i) + t * ax(1,i)
            x(2,i) = x(1,i) + t * (vx(1,i) + vx(2,i)) / 2.0d0
            vy(2, i) = vy(1,i) + t * ay(1,i)
            y(2,i) = y(1,i) + t * (vy(1,i) + vy(2,i)) / 2.0
            vz(2, i) = vz(1,i) + t * az(1,i)
            z(2,i) = z(1,i) + t * (vz(1,i) + vz(2,i)) / 2.0d0
        end do    

        ! Verlet i > 1
        do i = 2, steps - 1
            call acceleration(bodies, x(i,:), y(i,:), z(i,:), m,  ax(i,:), ay(i,:), az(i,:))
            do j = 1, bodies
                x(i + 1,j) = 2 * x(i,j) - x(i - 1,j) + ax(i,j) * t**2
                y(i + 1,j) = 2 * y(i,j) - y(i - 1,j) + ay(i,j) * t**2
                z(i + 1,j) = 2 * z(i,j) - z(i - 1,j) + az(i,j) * t**2
    
                vx(i + 1,j) = (x(i + 1,j) - x(i,j)) / t
                vy(i + 1,j) = (y(i + 1,j) - y(i,j)) / t
                vz(i + 1,j) = (z(i + 1,j) - z(i,j)) / t
            end do
        end do
    end subroutine
end module
