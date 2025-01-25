program lanier_n_body
  use model
  integer :: i, id
  real :: time_Start, time_End
  double precision, dimension(:,:), allocatable :: x, y, z, ax, ay, az, vx, vy, vz, lx, ly, lz
  double precision, dimension(:), allocatable :: m, kinetic
  double precision :: t, tot_potential, tot_kinetic, total_lx, total_ly, total_lz
  character(len=40) :: name

  ! Allocate arrays
  allocate(x(num_steps, num_bodies))
  allocate(y(num_steps, num_bodies))
  allocate(z(num_steps, num_bodies))
  allocate(ax(num_steps, num_bodies))
  allocate(ay(num_steps, num_bodies))
  allocate(az(num_steps, num_bodies))
  allocate(vx(num_steps, num_bodies))
  allocate(vy(num_steps, num_bodies))
  allocate(vz(num_steps, num_bodies))
  allocate(lx(num_steps, num_bodies))
  allocate(ly(num_steps, num_bodies))
  allocate(lz(num_steps, num_bodies))
  allocate(kinetic(num_steps))
  allocate(m(num_bodies))

  ! Open and read initial conditions from file
  open(unit=10, file='1980-11-12.txt', action='read')
  do i = 1, num_bodies
      read(10, *) name, id, m(i), x(1, i), y(1, i), z(1, i), vx(1, i), vy(1, i), vz(1, i)
      print *, "Read body ", i, " - Name: ", name, " ID: ", id, " Mass: ", m(i)
  end do
  close(10)

  days = 365 * 10 ! 50 years
  seconds = 60 * 60 * 24 * days ! seconds in a day
  t = seconds / num_steps

  ! Start Verlet timer
  call cpu_time(time_Start)
  ! Call Verlet subroutine
  call Verlet(num_bodies, num_steps, t, m(:), x(:,:), y(:,:), z(:,:), vx(:,:), vy(:,:), vz(:,:), ax(:,:), ay(:,:), az(:,:))
  ! End Verlet timer
  call cpu_time(time_End)

  ! Open output files
  open(unit=100, file='Lx.txt')
  open(unit=200, file='Ly.txt')
  open(unit=300, file='Lz.txt')
  open(unit=400, file='kinetic.txt')
  open(unit=500, file='potential.txt')
  open(unit=600, file='tot_energy.txt')

  ! Output results in a loop
  do i = 2, num_steps, 5
      call momentum_calc(m(:), x(i,:), y(i,:), z(i,:), vx(i,:), vy(i,:), vz(i,:),&
           lx(i,:), ly(i,:), lz(i,:), kinetic(:), tot_potential, tot_kinetic, &
           total_lx, total_ly, total_lz, num_bodies)
      write(100, *) t * i, total_lx
      write(200, *) t * i, total_ly
      write(300, *) t * i, total_lz
      write(400, *) t * i, tot_kinetic
      write(500, *) t * i, tot_potential
      write(600, *) t * i, tot_potential + tot_kinetic
  end do

  ! Close output files
  close(100)
  close(200)
  close(300)
  close(400)
  close(500)
  close(600)

  ! Open output files for individual bodies
  open(unit=1000, file = 'jupiter.txt')
  open(unit=2000, file ='io.txt')
  open(unit=3000, file ='europa.txt')
  open(unit=4000, file ='ganymede.txt')
  open(unit=5000, file ='callisto.txt')

!     open(unit=1000, file='sun.txt')
!   open(unit=2000, file='mercury.txt')
!   open(unit=3000, file='venus.txt')
!   open(unit=4000, file='earth.txt')
!   open(unit=5000, file='mars.txt')
!   open(unit=6000, file='jupiter.txt')
!   open(unit=7000, file='io.txt')
!   open(unit=8000, file='europa.txt')
!   open(unit=9000, file='ganymede.txt')
!   open(unit=1010, file='callisto.txt')
!   open(unit=1020, file='moon.txt')
!   open(unit=1030, file='smith.txt')

  ! Output individual body positions
  do i = 2, num_steps, 100
      write(1000, *) x(i, 1), y(i, 1), z(i, 1)
      write(2000, *) x(i, 2), y(i, 2), z(i, 2)
      write(3000, *) x(i, 3), y(i, 3), z(i, 3)
      write(4000, *) x(i, 4), y(i, 4), z(i, 4)
      write(5000, *) x(i, 5), y(i, 5), z(i, 5)
    !   write(6000, *) x(i, 6), y(i, 6), z(i, 6)
    !   write(7000, *) x(i, 7), y(i, 7), z(i, 7)
    !   write(8000, *) x(i, 8), y(i, 8), z(i, 8)
    !   write(9000, *) x(i, 9), y(i, 9), z(i, 9)
    !   write(1010, *) x(i, 10), y(i, 10), z(i, 10)
    !   write(1020, *) x(i, 11), y(i, 11), z(i, 11)
    ! !   write(1030, *) x(i, 12), y(i, 12), z(i, 12)
  end do

  ! Close output files for individual bodies
  close(1000)
  close(2000)
  close(3000)
  close(4000)
  close(5000)
!   close(6000)
!   close(7000)
!   close(8000)
!   close(9000)
!   close(1010)
!   close(1020)
!   close(1030)
  write(*,*) 'Number of bodies: ', num_bodies
  write(*,*) 'Time for Verlet: ', (time_End - time_Start) / 60, 'minutes'
end program lanier_n_body