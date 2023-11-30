!-------------------------------------------------------------------------------
! "DEM2xyz v.3.0" (DEM manager tool)
! Copyright 2016-2022 (RSE SpA)
! "DEM2xyz v.3.0" authors and email contact are provided on the documentation 
! file.
! This file is part of DEM2xyz v.3.0 .
! DEM2xyz v.3.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DEM2xyz is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DEM2xyz. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! This file is copied and pasted from SPHERA v.8.0 (RSE SpA).
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: distance_point_line_2D  
! Description: Computation of the distance between a point and a plane.     
!-------------------------------------------------------------------------------
subroutine distance_point_line_2D(P0,P1_line,P2_line,dis,normal)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(IN) :: P0(2),P1_line(2),P2_line(2)
double precision,intent(INOUT) :: dis
double precision,intent(INOUT) :: normal(2)
double precision :: a,b,c
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
a = (P2_line(2) - P1_line(2))
b = -(P2_line(1) - P1_line(1)) 
c = -(a * P1_line(1) + b * P1_line(2))
normal(1) = a / dsqrt(a ** 2 + b ** 2)
normal(2) = b / dsqrt(a ** 2 + b ** 2)
dis = (a * P0(1) + b * P0(2) + c) / dsqrt(a ** 2 + b ** 2)
!------------------------
! Deallocations
!------------------------
return
end subroutine distance_point_line_2D

