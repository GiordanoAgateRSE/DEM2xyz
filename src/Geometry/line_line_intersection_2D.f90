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
! Program unit: line_line_intersection_2D 
! Description: Computation of the intersection between two lines in 2D.     
!-------------------------------------------------------------------------------
subroutine line_line_intersection_2D(P1_line1,P2_line1,P1_line2,P2_line2,Pint, &
                                     test)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: P1_line1(2),P2_line1(2),P1_line2(2),P2_line2(2)
logical,intent(out) :: test
double precision,intent(out) :: Pint(2)
double precision :: aa,bb,cc,dd,ee,ff,denom
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
test = .false.
Pint(:) = -9.d8
!------------------------
! Statements
!------------------------
aa = (P1_line1(1) * P2_line1(2)) - (P1_line1(2) * P2_line1(1))
bb = (P1_line2(1) * P2_line2(2)) - (P1_line2(2) * P2_line2(1))
cc = P1_line2(1) - P2_line2(1)
dd = P1_line1(1) - P2_line1(1)
ee = P1_line2(2) - P2_line2(2)
ff = P1_line1(2) - P2_line1(2)
denom = (dd * ee) - (ff * cc)
if ((denom>1.d-9).or.(denom<-1.d-9)) then
   test = .true.
   Pint(1) = ((aa * cc) - (dd * bb)) / denom
   Pint(2) = ((aa * ee) - (ff * bb)) / denom
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine line_line_intersection_2D

