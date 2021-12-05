B-Spline-based Interpolation on a Regular Grid

# Brief description

This Bspline fortran code is taken and adapted
[from](http://jacobwilliams.github.io/bspline-fortran/)

The code is based on Bspline interpolation on a regular grid. The dependencies are eliminated from the original version. The code is organized according to the dimension of the data to be interpolated (i.e. 1D-3D). The code is originally written in a serial form and the goal here is to accelerate it with the use of OpenACC offloading.                                                 

## Subroutines

The core routines for the the Bspline code are:

```Fortran

!f(x)
subroutine db1ink(x,nx,fcn,kx,iknot,tx,bcoef,iflag)
subroutine db1val(xval,idx,tx,nx,kx,bcoef,f,iflag,inbvx,w0,extrap)

!f(x,y)
subroutine db2ink(x,nx,y,ny,fcn,kx,ky,iknot,tx,ty,bcoef,iflag)
subroutine db2val(xval,yval,idx,idy,tx,ty,nx,ny,kx,ky,bcoef,f,iflag,inbvx,inbvy,iloy,w1,w0,extrap)

!f(x,y,z)
subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,iknot,tx,ty,tz,bcoef,iflag)
subroutine db3val(xval,yval,zval,idx,idy,idz,tx,ty,tz,nx,ny,nz,kx,ky,kz,bcoef,f,iflag,inbvx,inbvy,inbvz,iloy,iloz,w2,w1,w0,extrap)
```

The ```ink``` routines compute the interpolant coefficients, and the ```val``` routines evalute the interpolant at the specified value of each coordinate.

# Compiling

A Makefile is provided for the compilation process. It can be compiled with NVHPC compiler, and also with GNU and Intel Fortran.

# Documentation

Further details about the documentation can be found [here](http://jacobwilliams.github.io/bspline-fortran/).

# License

The bspline-fortran source code and related files are distributed under a permissive free software [license](https://github.com/jacobwilliams/bspline-fortran/blob/master/LICENSE) (BSD-style).
