!*****************************************************************************************
!>
!  Units test for 1d-6d tensor product b-spline interpolation.

    program bspline_test

    use bspline_module
    use bspline_kinds_module, only: wp, ip

    implicit none

    integer(ip),parameter :: nx = 50     !! number of points in x
    integer(ip),parameter :: nx_new = 512     !! nbr of interpolated points
    integer(ip),parameter :: kx = 4     !! order in x

    integer(ip),parameter :: iknot = 0  !! automatically select the knots

    real(wp), parameter   :: x_min=0._wp, x_max=10._wp !! range og x-grid

    real(wp) :: x(nx),tx(nx+kx)
    real(wp) :: fcn_1d(nx)

    real(wp),dimension(3*kx)                     :: w1_1d

    real(wp) :: tol,dx,dx_new,int_1d,int_1d_new
    real(wp),dimension(6) :: val,tru,err,errmax
    logical :: fail
    integer(ip) :: i,j,k,l,m,n,idx,idy,idz,idq,idr,ids
    integer(ip),dimension(6) :: iflag
    integer(ip) :: inbvx,inbvy,inbvz,inbvq,inbvr,inbvs
    integer(ip) :: iloy,iloz,iloq,ilor,ilos
    integer ::     count_rate, count_max,count,t_start,t_final
    real    ::     time_s

    real(wp), allocatable :: xout(:),fout(:)

    call system_clock(count_max=count_max, count_rate=count_rate)

    call system_clock(t_start)

    print*, ""
    print*, "-----------1D-Bspline Interpolation--------------"
    print*, ""

!print out the output file
    open(100,file='1D-TestBspline-Without-Interpolation.dat')
    open(101,file='1D-TestBspline-With-Interpolation.dat')

    fail = .false.
    tol = 100 * epsilon(1.0_wp)
    idx = 0

    allocate(xout(nx_new)); allocate(fout(nx_new))

     dx = (x_max-x_min)/real(nx-1,wp)
     do i=1,nx
        x(i) = x_min + real(i-1,wp)*dx
     end do
     do i=1,nx
        fcn_1d(i) = f1(x(i))
     end do

!defining the interpolated grid
     dx_new = (x_max-x_min)/real(nx_new-1,wp)
     do i=1,nx_new
        xout(i) = x_min + real(i-1,wp)*dx_new
     enddo
     print*,"--Grid parameter----"
     print*, "--nbr of point:nx, nx_new",nx, nx_new
     print*, "--spatial step:dx, dx_new", dx, dx_new
     print*,"--xout(1),xout(nx_new)",xout(1),xout(nx_new)
     print*,"--x(a),x(nx)",x(1),x(nx)
     
    !have to set these before the first evaluate call:
    inbvx = 1

    ! initialize
    call db1ink(x,nx,fcn_1d,kx,iknot,tx,fcn_1d,iflag(1))

    if (any(iflag/=0)) then
        do i=1,1
            if (iflag(i)/=0) then
                write(*,*) 'Error initializing ',i,'D spline: '//get_status_message(iflag(i))
            end if
        end do
    end if

    ! compute max error at interpolation points

     errmax = 0.0_wp
     do i=1,size(xout)
       call db1val(xout(i),idx,&
           tx,nx,kx,fcn_1d,val(1),iflag(1),inbvx,&
           w1_1d)
                        tru(1)    = f1(x(i))
                        err(1)    = abs(tru(1)-val(1))
                        errmax(1) = max(err(1),errmax(1))
           fout(i) = val(1)
     end do

    ! check max error against tolerance
    do i=1,1
        write(*,*) i,'D: max error:', errmax(i)
!        if (errmax(i) >= tol) then
!            write(*,*)  ' ** test failed ** '
!        else
!            write(*,*)  ' ** test passed ** '
!        end if
        write(*,*) ''
    end do

!Evaluate the integral
   int_1d = sum(fcn_1d(1:nx))*dx
   int_1d_new = sum(fout(1:nx_new))*dx_new

   print*,"--1D-Integral without interpolation", int_1d
   print*,"--1D-Integral with interpolation", int_1d_new

   call system_clock(t_final)

   time_s = real(t_final - t_start)/real(count_rate)
   print*,""
   print*, '--Time it takes (s)', time_s

!prinout data
!data without interpolation
    do i=1,nx
      write(100,*)x(i),fcn_1d(i)
   enddo
!data with interpolation
   do i=1,size(xout)
      write(101,*)xout(i),fout(i)
   enddo
 
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
        piov2 = 2.0_wp * atan(1.0_wp)
        f2 = 0.5_wp * (y*exp(-x) + sin(piov2*y) )
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
