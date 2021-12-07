       subroutine read_input(nx,ny,nx_new,ny_new,kx,ky,iknot,&
                            x_min,x_max,y_min,y_max)

         use bspline_kinds_module, only: wp, ip

         integer(ip)   :: nx,ny,nz,nx_new,ny_new,nz_new,kx,ky,kz,iknot
         real(wp)      :: x_min,x_max,y_min,y_max,z_min,z_max
         character(80) :: grid_parameter,Bsp_parameter

         read(*,"(A)")grid_parameter
         read(*,*)
         read(*,*)nx,ny,nz
         read(*,*)nx_new,ny_new,nz_new
         read(*,*)
         read(*,*)x_min,x_max
         read(*,*)y_min,y_max
         read(*,*)z_min,z_max
         read(*,*)
         read(*,"(A)")Bsp_parameter
         read(*,*)kx,ky,kz
         read(*,*)iknot

       end subroutine read_input

    
