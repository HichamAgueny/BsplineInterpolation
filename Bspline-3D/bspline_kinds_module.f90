!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!### Description
!  Numeric kind definitions for BSpline-Fortran.

    module bspline_kinds_module

    use,intrinsic :: iso_fortran_env

    implicit none

    private

     integer,parameter,public :: wp = selected_real_kind(13,300)
     integer,parameter,public :: ip = selected_int_kind(9)

!    integer,parameter,public :: wp = selected_real_kind(27,2400)   !! Real working precision [8 bytes]

!    integer,parameter,public :: ip = selected_int_kind(18)    !! Integer working precision if not specified [4 bytes]

!*****************************************************************************************
    end module bspline_kinds_module
!*****************************************************************************************
