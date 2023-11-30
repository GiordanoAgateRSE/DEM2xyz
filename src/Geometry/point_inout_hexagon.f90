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
! Program unit: point_inout_hexagon
! Description: Test to evaluate if a point lies inside or strictly outside a
!              generic hexagon. The hexagon is partitioned into 4 triangles
!              (P1P2P6,P2P5P6,P2P3P5,P3P4P5). A point is internal to the hexagon
!              if it is internal to one of its triangles.   
!-------------------------------------------------------------------------------
subroutine point_inout_hexagon(point,point_pol_1,point_pol_2,point_pol_3,      &
                               point_pol_4,point_pol_5,point_pol_6,test)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
double precision,intent(in) :: point_pol_3(2),point_pol_4(2),point_pol_5(2)
double precision,intent(in) :: point_pol_6(2)
integer(4),intent(inout) :: test
integer(4) :: test1,test2,test3,test4
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,         &
                                                        point_pol_1,           &
                                                        point_pol_2,           &
                                                        point_pol_3,           &
                                                        point_pol_4,           &
                                                        point_pol_5,           &
                                                        point_pol_6,test)
      implicit none
      integer(4),intent(in) :: n_sides
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
      double precision :: dis1,dis2
      double precision :: normal(2)
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
test1 = 0
test2 = 0
test3 = 0
test4 = 0
!------------------------
! Statements
!------------------------
call point_inout_convex_non_degenerate_polygon(point,3,point_pol_1,            &
                                               point_pol_2,point_pol_6,        &
                                               point_pol_6,point_pol_6,        &
                                               point_pol_6,test1)
test = test1
if (test/=1) then 
   call point_inout_convex_non_degenerate_polygon(point,3,point_pol_2,         &
                                                  point_pol_5,point_pol_6,     &
                                                  point_pol_6,point_pol_6,     &
                                                  point_pol_6,test2)
   test = test2
endif
if (test/=1) then
   call point_inout_convex_non_degenerate_polygon(point,3,point_pol_2,         &
                                                  point_pol_3,point_pol_5,     &
                                                  point_pol_5,point_pol_5,     &
                                                  point_pol_5,test3)
   test = test3
endif
if (test/=1) then
   call point_inout_convex_non_degenerate_polygon(point,3,point_pol_3,         &
                                                  point_pol_4,point_pol_5,     &
                                                  point_pol_5,point_pol_5,     &
                                                  point_pol_5,test4)
   test = test4
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine point_inout_hexagon

