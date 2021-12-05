!*****************************************************************************************
!>
!  Units test for 1d-6d tensor product b-spline interpolation.

    program bspline_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer(ip),parameter :: nx = 50     !! number of points in x
    integer(ip),parameter :: ny = 50     !! number of points in y
    integer(ip),parameter :: nx_new = 512     !! nbr of interpolated x-points
    integer(ip),parameter :: ny_new = 512     !! nbr of interpolated y-points
    integer(ip),parameter :: kx = 4     !! order in x
    integer(ip),parameter :: ky = 4     !! order in y

    integer(ip),parameter :: iknot = 0  !! automatically select the knots

    real(wp), parameter   :: x_min=0._wp, x_max=10._wp !! range og x-grid
    real(wp), parameter   :: y_min=0._wp, y_max=10._wp !! range og y-grid

    real(wp) :: x(nx),y(ny),tx(nx+kx),ty(ny+ky)
    real(wp) :: fcn_2d(nx,ny)

    real(wp),dimension(3*kx)                     :: w1_1d
    real(wp),dimension(ky)                       :: w1_2d
    real(wp),dimension(3*max(kx,ky))             :: w2_2d

    real(wp) :: tol,dx,dy,dx_new,dy_new,int_2d,int_2d_new, &
                sum_x
    real(wp),dimension(6) :: val,tru,err,errmax
    logical :: fail
    integer(ip) :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids
    integer(ip),dimension(6) :: iflag
    integer(ip) :: inbvx,inbvy,inbvz,inbvq,inbvr,inbvs
    integer(ip) :: iloy,iloz,iloq,ilor,ilos
    integer ::     count_rate, count_max,count,t_start,t_final
    real    ::     time_s

    real(wp), allocatable :: xout(:),yout(:),fout(:,:)

    call system_clock(count_max=count_max, count_rate=count_rate)

    call system_clock(t_start)

    print*, ""
    print*, "-----------2D-Bspline Interpolation--------------"
    print*, ""

    fail = .false.
    tol = 100 * epsilon(1.0_wp)
    idx = 0
    idy = 0

    allocate(xout(nx_new)); allocate(yout(ny_new));allocate(fout(nx_new,ny_new))

     dx = (x_max-x_min)/real(nx-1,wp); dy = (y_max-y_min)/real(ny-1,wp)
!x-grid
     do i=1,nx
        x(i) = x_min + real(i-1,wp)*dx
     end do
!y-grid
     do j=1,ny
        y(j) = y_min + real(j-1,wp)*dy
     enddo
!2d function to be interpolated
     do i=1,nx
       do j=1,ny
          fcn_2d(i,j) = f2(x(i),y(j))
       enddo
     end do

!defining the interpolated 2D-grid
     dx_new = (x_max-x_min)/real(nx_new-1,wp); dy_new = (y_max-y_min)/real(ny_new-1,wp)
     do i=1,nx_new
        xout(i) = x_min + real(i-1,wp)*dx_new
     enddo
     do j=1,ny_new
        yout(j) = y_min + real(j-1,wp)*dy_new
     enddo

     print*,"--Grid parameter----"
     print*, "--nbr of point:nx, nx_new",nx, nx_new
     print*, "--nbr of point:ny, ny_new",ny, ny_new
     print*, "--spatial step:dx, dx_new", dx, dx_new
     print*, "--spatial step:dy, dy_new", dy, dy_new
     print*,"--xmax,ymax",x_max,y_max
     print*,""
 
    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1

    ! initialize
    call db2ink(x,nx,y,ny,fcn_2d,kx,ky,iknot,tx,ty,fcn_2d,iflag(2))

    if (any(iflag/=0)) then
        do i=2,2
            if (iflag(i)/=0) then
                write(*,*) 'Error initializing ',i,'D spline: '//get_status_message(iflag(i))
            end if
        end do
    end if

    ! compute max error at interpolation points

     errmax = 0.0_wp
     do i=1,size(xout)
       do j=1,size(yout)
          call db2val(xout(i),yout(j),idx,idy,&
                     tx,ty,nx,ny,kx,ky,fcn_2d,val(2),iflag(2),&
                     inbvx,inbvy,iloy,&
                     w1_2d,w2_2d)

                     tru(2)    = f2(x(i),y(j))
                     err(2)    = abs(tru(2)-val(2))
                     errmax(2) = max(err(2),errmax(2))
                    
          fout(i,j) = val(2)
       enddo
     end do

    ! check max error against tolerance
    do i=2,2
        write(*,*) i,'D: max error:', errmax(i)
!        if (errmax(i) >= tol) then
!            write(*,*)  ' ** test failed ** '
!        else
!            write(*,*)  ' ** test passed ** '
!        end if
        write(*,*) ''
    end do

!Evaluate the integral
!without interpolation
   sum_x=0._wp
   do j=1,ny
      sum_x = sum_x + sum(fcn_2d(1:nx,j))
   enddo
   int_2d = sum_x*dx*dy

!with interpolation
   sum_x=0._wp
   do j=1,ny_new
      sum_x = sum_x + sum(fout(1:nx_new,j))
   enddo
   int_2d_new = sum_x*dx_new*dy_new

   print*,"--2D-Integral without interpolation", int_2d
   print*,"--2D-Integral with interpolation", int_2d_new

   call system_clock(t_final)

   time_s = real(t_final - t_start)/real(count_rate)
   print*,""
   print*, '--Time it takes (s)', time_s

!print out the output file
!    open(100,file='2D-TestBspline-Without-Interpolation.dat')
!    open(101,file='2D-TestBspline-With-Interpolation.dat')

!prinout data
!data without interpolation
!   do i=1,nx
!      do j=1,ny
!         write(100,*)x(i),y(j),fcn_2d(i,j)
!      enddo
!      write(100,*)
!   enddo
!data with interpolation
!   do i=1,size(xout)
!     do j=1,size(yout)
!      write(101,*)xout(i),yout(j),fout(i,j)
!     enddo
!     write(101,*)
!   enddo

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp) :: x
        !f1 = 0.5_wp * (x*exp(-x) + sin(x) )
         f1 = sin(x)
        end function f1

        real(wp) function f2(x,y) !! 2d test function
        implicit none
        real(wp) x,y,piov2
        !piov2 = 2.0_wp * atan(1.0_wp)
        !f2 = 0.5_wp * (y*exp(-x) + sin(piov2*y) )
         f2 = sin(sqrt(x*x + y*y))
        end function f2

        real(wp) function f3 (x,y,z) !! 3d test function
        implicit none
        real(wp) x,y,z,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f3 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) )
        end function f3

        real(wp) function f4 (x,y,z,q) !! 4d test function
        implicit none
        real(wp) x,y,z,q,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f4 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) + q )
        end function f4

        real(wp) function f5 (x,y,z,q,r) !! 5d test function
        implicit none
        real(wp) x,y,z,q,r,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f5 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) + q*r )
        end function f5

        real(wp) function f6 (x,y,z,q,r,s) !! 6d test function
        implicit none
        real(wp) x,y,z,q,r,s,piov2
        piov2 = 2.0_wp*atan(1.0_wp)
        f6 = 0.5_wp*( y*exp(-x) + z*sin(piov2*y) + q*r + 2.0_wp*s )
        end function f6

    end program bspline_test
