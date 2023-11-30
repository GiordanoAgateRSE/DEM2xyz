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
! Description. “DEM2xyz v.3.0” (RSE SpA) reads a “.asc” input grid file and 
!              writes the 
!              associated DEM in a corresponding “xyz” file, possibly 
!              changing the spatial resolution (as requested by the user). 
!              In case the absolute value of the mean latitude is provided with
!              a non-negative value, the following conversion takes place 
!              "(lon,lat) in (°) to (X,Y) in (m)". In this case, an 
!              interpolation (weighted on the square of the distance) is 
!              carried out to provide a regular Cartesian output grid in (X,Y). 
!              The height of the DEM points which belong to the digging/filling 
!              regions (provided in input) is modified. After this treatment, 
!              each digging/filling region has null slope. 
!              Bathymetry is possibly extruded from the heights of the 
!              most upstream and downstream coastline points.
!              The bathymetry/reservoir extrusion is corrected in case the  
!              volume reservoir is provided as an input parameter.
!              Multiple reservoirs are admitted.
!              Digging regions cannot overlap each other.
!              Reservoir/bathymetry regions cannot overlap each other.
!              In case a digging region overlaps a reservoir region, the latter 
!              holds the priority.
!              In the presence of a volume correction, two reference shapes are 
!              available: "reservoir" and "volcanic lake".
!              The output rectangular DEM can be a sub-domain of the input DEM. 
!              The edge coordinates to cut the output DEM are provided in input.
!              The cut procedure is a post-processing task which does not 
!              affect the other tasks.
!              A translation vector is added to the output grid point positions 
!              so that the origin of the reference system in output is chosen 
!              by the user by means of the input file.
!              DEM2xyz v.3.0 is compatible with SPHERA v.9.0.0 (RSE SpA).
!              Variables:
!              input variables (ref. template of the main input file)
!              coastline(n_bathymetries,n_rows,n_col_out): logical flag to 
!                 detect the reservoir coastline.
!              dis_Pdown_Pint: distance between the most downstream point (m) 
!                 and the current intersection point
!              dis_down_up(n_bathymetries): distance between the most upstream 
!                 and downstream points (m)
!              dis_Pint_Pcoast: distance between "point_coast" and "Pint" 
!              dx,dy: spatial resolution -final values in (m)-
!              mat_z_in(n_row_n_col_in): input DEM
!              mat_z_out(n_row_n_col_out): output DEM
!              n_col_in: number of columns in the input DEM
!              n_col_out: number of columns in the output DEM
!              n_points_in: number of input vertices
!              n_points_out: number of output vertices
!              n_row: number of rows in the input/output DEM
!              Pint(2): intersection between the line r_down_up (passing for  
!                 both the upstream and the downstream points) and the line r_iC 
!                 (passing for the current reservoir "point" and the associated 
!                 coastline point).
!              point_coast(2): provided an inner reservoir point, to detect the 
!                 associated coast point which belongs to the same coast side 
!                 and is the closest to the line r_iC. This line passes for the 
!                 inner reservoir point and is parallel to the above normal.
!              reservoir(n_bathymetries,n_rows,n_col_out): logical flag to 
!                 detect the reservoir.
!              volume_res_est(n_bathymetries): first estimations of the 
!                 reservoir volumes
!              volume_res_corr(n_bathymetries): corrected estimations of the 
!                 reservoir volumes
!              weight(n_bathymetries,n_row,n_col_out): point weights for 
!                 reservoir volume corrections
!              weight_sum(n_bathymetries): weight sums 
!              x_in,y_in: horizontal coordinates in input -(m) or (°)- 
!              x_out,y_out: horizontal coordinates in output (m)
!              z_Pint: hieght of the point Pint (m)
!-------------------------------------------------------------------------------
PROGRAM DEM2xyz
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
logical :: test_logical
integer :: i_in,j_in,i_out,j_out,n_col_in,n_col_out,n_row,res_fact,n_points_in
integer :: n_points_out,i_aux,j_aux,n_digging_regions,i_reg,test_integer
integer :: i_bath,aux_integer,i_close,j_close,j2_out,i2_out,n_bathymetries
integer :: alloc_stat,i_out_cut_min,i_out_cut_max,j_out_cut_min,j_out_cut_max
double precision :: dx,dy,abs_mean_latitude,denom,distance,x_in,x_out,y_in,y_out
double precision :: dis,dis2,min_dis2,z_Pint,dis_Pint_Pcoast,aux_scalar
double precision :: dis_Pdown_Pint,aux_scalar_2,aux_scalar_3,x_inp_min,y_inp_min
double precision :: dis3,x_inp_cut_min,y_inp_cut_min,x_inp_cut_max,y_inp_cut_max
double precision :: dy_cut_min,dy_cut_max,dx_cut_min,x_trans_out,y_trans_out
double precision :: dx_cut_max
double precision :: point(2),point2(2),point_plus_normal(2),normal(2),normal2(2)
double precision :: point_coast(2),Pint(2)
logical,dimension(:),allocatable :: volume_flag
integer,dimension(:),allocatable :: n_digging_vertices,n_vertices_around_res
integer,dimension(:),allocatable :: weight_type(:)
double precision,dimension(:),allocatable :: z_digging_regions,z_FS,z_eps
double precision,dimension(:),allocatable :: z_downstream,volume_res_est
double precision,dimension(:),allocatable :: volume_res_inp,weight_sum
double precision,dimension(:),allocatable :: dis_down_up,volume_res_corr
double precision,dimension(:,:),allocatable :: mat_z_in,mat_z_out
double precision,dimension(:,:),allocatable :: pos_res_downstream
double precision,dimension(:,:),allocatable :: pos_res_upstream
logical,dimension(:,:,:),allocatable :: reservoir,coastline
double precision,dimension(:,:,:),allocatable :: digging_vertices,weight
double precision,dimension(:,:,:),allocatable :: vertices_around_res
character(len=100) :: char_aux,input_grid_file_name
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
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
interface
   subroutine point_inout_hexagon(point,point_pol_1,point_pol_2,point_pol_3,   &
                                  point_pol_4,point_pol_5,point_pol_6,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_hexagon
end interface
interface
   subroutine point_inout_pentagon(point,point_pol_1,point_pol_2,point_pol_3,  &
                                   point_pol_4,point_pol_5,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_pentagon
end interface
interface
   subroutine point_inout_quadrilateral(point,point_pol_1,point_pol_2,         &
                                        point_pol_3,point_pol_4,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_quadrilateral
end interface
interface
   subroutine distance_point_line_2D(P0,P1_line,P2_line,dis,normal)
      implicit none
      double precision,intent(in) :: P0(2),P1_line(2),P2_line(2)
      double precision,intent(inout) :: dis
      double precision,intent(inout) :: normal(2)
   end subroutine distance_point_line_2D
end interface
interface
   subroutine line_line_intersection_2D(P1_line1,P2_line1,P1_line2,P2_line2,   &
                                        Pint,test)
      implicit none
      double precision,intent(in) :: P1_line1(2),P2_line1(2),P1_line2(2)
      double precision,intent(in) :: P2_line2(2)
      logical,intent(out) :: test
      double precision,intent(out) :: Pint(2)
   end subroutine line_line_intersection_2D
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
n_digging_regions = 0
n_bathymetries = 0
!------------------------
! Statements
!------------------------
write(*,*) "DEM2xyz v.3.0 (RSE SpA) is running. DEM2xyz is a DEM manager tool. "
write(*,*) "Reading DEM file, DEM2xyz main input file and pre-processing. "
open(12,file='DEM2xyz.inp')
read(12,*) input_grid_file_name
read(12,*) res_fact,abs_mean_latitude,n_digging_regions,n_bathymetries
read(12,*) x_inp_cut_min,y_inp_cut_min
read(12,*) x_inp_cut_max,y_inp_cut_max
read(12,*) x_trans_out,y_trans_out
if (n_digging_regions>0) then
   if (.not.allocated(n_digging_vertices)) then
      allocate(n_digging_vertices(n_digging_regions),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "n_digging_vertices" failed; the execution',&
            ' terminates here.'
         stop
         else
            write(*,*) 'Allocation of "n_digging_vertices" is successfully ',  &
               'completed.'
      endif
   endif
   if (.not.allocated(z_digging_regions)) then
      allocate(z_digging_regions(n_digging_regions),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "z_digging_regions" failed; the execution', &
            ' terminates here.'
         stop
         else
            write(*,*) 'Allocation of "z_digging_regions" is successfully ',   &
               'completed.'
      endif
   endif
   if (.not.allocated(digging_vertices)) then
      allocate(digging_vertices(n_digging_regions,6,2),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "digging_vertices" failed; the execution',  &
            ' terminates here.'
         stop
         else
            write(*,*) 'Allocation of "digging_vertices" is successfully ',    &
               'completed.'
      endif
   endif
   n_digging_vertices(:) = 0
   z_digging_regions(:) = 0.d0
   digging_vertices(:,:,:) = 0.d0
   do i_reg=1,n_digging_regions
      read(12,*) z_digging_regions(i_reg),n_digging_vertices(i_reg)
      do j_aux=1,n_digging_vertices(i_reg)
         read (12,*) digging_vertices(i_reg,j_aux,1:2)
      enddo
   enddo
endif
if (n_bathymetries>0) then
   if (.not.allocated(volume_flag)) then
      allocate(volume_flag(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "volume_flag" failed; the execution',       &
            ' terminates here.'
         stop
         else
            write(*,*) 'Allocation of "volume_flag" is successfully ',         &
               'completed.'
      endif
   endif
   if (.not.allocated(n_vertices_around_res)) then
      allocate(n_vertices_around_res(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "n_vertices_around_res" failed; the ',      &
            'execution terminates here.'
         stop
         else
            write(*,*) 'Allocation of "n_vertices_around_res" is ',            &
               'successfully completed.'
      endif
   endif
   if (.not.allocated(z_downstream)) then
      allocate(z_downstream(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "z_downstream" failed; the execution ',     &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "z_downstream" is successfully completed.'
      endif
   endif
   if (.not.allocated(z_FS)) then
      allocate(z_FS(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "z_FS" failed; the execution terminates ',  &
            'here.'
         stop
         else
            write(*,*) 'Allocation of "z_FS" is successfully completed.'
      endif
   endif
   if (.not.allocated(z_eps)) then
      allocate(z_eps(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "z_eps" failed; the execution terminates ', &
            'here.'
         stop
         else
            write(*,*) 'Allocation of "z_eps" is successfully completed.'
      endif
   endif
   if (.not.allocated(volume_res_est)) then
      allocate(volume_res_est(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "volume_res_est" failed; the execution ',   &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "volume_res_est" is successfully ',      &
               'completed.'
      endif
   endif
   if (.not.allocated(volume_res_inp)) then
      allocate(volume_res_inp(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "volume_res_inp" failed; the execution ',   &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "volume_res_inp" is successfully ',      &
               'completed.'
      endif
   endif
   if (.not.allocated(volume_res_corr)) then
      allocate(volume_res_corr(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "volume_res_corr" failed; the execution ',  &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "volume_res_corr" is successfully ',     &
               'completed.'
      endif
   endif
   if (.not.allocated(weight_sum)) then
      allocate(weight_sum(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "weight_sum" failed; the execution ',       &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "weight_sum" is successfully completed.'
      endif
   endif
   if (.not.allocated(dis_down_up)) then
      allocate(dis_down_up(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "dis_down_up" failed; the execution ',      &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "dis_down_up" is successfully completed.'
      endif
   endif
   if (.not.allocated(weight_type)) then
      allocate(weight_type(n_bathymetries),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "weight_type" failed; the execution ',      &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "weight_type" is successfully completed.'
      endif
   endif
   if (.not.allocated(pos_res_downstream)) then
      allocate(pos_res_downstream(n_bathymetries,2),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "pos_res_downstream" failed; the ',         &
            'execution terminates here.'
         stop
         else
            write(*,*) 'Allocation of "pos_res_downstream" is successfully ',  &
               'completed.'
      endif
   endif
   if (.not.allocated(pos_res_upstream)) then
      allocate(pos_res_upstream(n_bathymetries,2),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "pos_res_upstream" failed; the ',           &
            'execution terminates here.'
         stop
         else
            write(*,*) 'Allocation of "pos_res_upstream" is successfully ',    &
               'completed.'
      endif
   endif
   if (.not.allocated(vertices_around_res)) then
      allocate(vertices_around_res(n_bathymetries,6,2),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "vertices_around_res" failed; the ',        &
            'execution terminates here.'
         stop
         else
            write(*,*) 'Allocation of "vertices_around_res" is successfully ', &
               'completed.'
      endif
   endif
   volume_flag(:) = .false.
   n_vertices_around_res(:) = 0
   z_downstream(:) = 0.d0
   z_FS(:) = 0.d0
   z_eps(:) = 0.d0
   volume_res_est(:) = 0.d0
   volume_res_inp(:) = 0.d0
   volume_res_corr(:) = 0.d0
   weight_sum(:) = 0.d0
   dis_down_up(:) = 0.d0
   weight_type(:) = 0
   pos_res_downstream(:,:) = 0.d0
   pos_res_upstream(:,:) = 0.d0
   vertices_around_res(:,:,:) = 0.d0
   do i_bath=1,n_bathymetries
      read(12,*) pos_res_downstream(i_bath,1),pos_res_downstream(i_bath,2),    &
         z_downstream(i_bath),pos_res_upstream(i_bath,1),                      &
         pos_res_upstream(i_bath,2),z_FS(i_bath),z_eps(i_bath),                &
         volume_flag(i_bath),weight_type(i_bath),n_vertices_around_res(i_bath)
      do j_aux=1,n_vertices_around_res(i_bath)
         read(12,*) vertices_around_res(i_bath,j_aux,1),                       &
            vertices_around_res(i_bath,j_aux,2)
      enddo
      if (volume_flag(i_bath).eqv..true.) then
         read(12,*) volume_res_inp(i_bath) 
      endif
! Distance between the most upstream and downstream points
      dis_down_up(i_bath) = dsqrt((pos_res_upstream(i_bath,1) -                &
                 pos_res_downstream(i_bath,1)) ** 2 +                          &
                 (pos_res_upstream(i_bath,2) - pos_res_downstream(i_bath,2))   &
                 ** 2)
   enddo
endif
close(12)
open(11,file=trim(input_grid_file_name))
read(11,*) char_aux,n_col_in
read(11,*) char_aux,n_row
read(11,*) char_aux,x_inp_min
read(11,*) char_aux,y_inp_min
read(11,*) char_aux,dy
if (.not.allocated(mat_z_in)) then
   allocate(mat_z_in(n_row,n_col_in),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(0,*) 'Allocation of "mat_z_in" failed; the execution terminates ', &
         'here.'
      stop
      else
         write(*,*) 'Allocation of "mat_z_in" is successfully completed.'
   endif
endif
mat_z_in = 0.d0
read(11,*)
do i_in=1,n_row
   read(11,*) mat_z_in(i_in,:)
enddo
close(11)
! Conversion from geograghic to cartographic coordinates and first assessment 
! of the number of output columns
dx_cut_min = x_inp_cut_min - x_inp_min
dy_cut_min = y_inp_cut_min - y_inp_min
dx_cut_max = x_inp_cut_max - x_inp_min
dy_cut_max = y_inp_cut_max - y_inp_min
if (abs_mean_latitude>=0.d0) then
! At this stage dx is not initialized and should be equal to dy
   call delta_lon_lat_to_delta_x_y(dy,dy,abs_mean_latitude,dx,dy)
   call delta_lon_lat_to_delta_x_y(dx_cut_min,dy_cut_min,abs_mean_latitude,    &
      dx_cut_min,dy_cut_min)
   call delta_lon_lat_to_delta_x_y(dx_cut_max,dy_cut_max,abs_mean_latitude,    &
      dx_cut_max,dy_cut_max)
   n_col_out = floor(n_col_in * dx / dy)
   else
      dx = dy
      n_col_out = n_col_in
endif
! Assessment of the indices to cut the output DEM. Pay attention to the format 
!    of the ".asc" files: see below.
i_out_cut_min = n_row - int(dy_cut_max / dy)
i_out_cut_max = n_row - int(dy_cut_min / dy)
j_out_cut_min = int(dx_cut_min / dx) + 1
j_out_cut_max = int(dx_cut_max / dx) + 1
if (.not.allocated(mat_z_out)) then
   allocate(mat_z_out(n_row,n_col_out),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(0,*) 'Allocation of "mat_z_out" failed; the execution terminates ',&
         'here.'
      stop
      else
         write(*,*) 'Allocation of "mat_z_out" is successfully completed.'
   endif
endif
mat_z_out = 0.d0
if (n_bathymetries>0) then
   if (.not.allocated(reservoir)) then
      allocate(reservoir(n_bathymetries,n_row,n_col_out),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "reservoir" failed; the execution ',        &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "reservoir" is successfully completed.'
      endif
   endif
   reservoir(:,:,:) = .false.
   if (.not.allocated(coastline)) then
      allocate(coastline(n_bathymetries,n_row,n_col_out),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "coastline" failed; the execution ',        &
            'terminates here.'
         stop
         else
            write(*,*) 'Allocation of "coastline" is successfully completed.'
      endif
   endif
   coastline(:,:,:) = .false.
   if (.not.allocated(weight)) then
      allocate(weight(n_bathymetries,n_row,n_col_out),STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Allocation of "weight" failed; the execution terminates ',&
            'here.'
         stop
         else
            write(*,*) 'Allocation of "weight" is successfully completed.'
      endif
   endif
   weight(:,:,:) = 0.d0
endif
write(*,*) "Possible grid interpolation, possible digging/filling DEM ",       &
   "regions, possible reservoir detection and writing xyz file (version ",     &
   "before the possible reservoir/batymetry extrusions). "
open(13,file="xyz_no_extrusion.txt")
write(13,'(a)') '           x(m)          y(m)           z(m)           z   '
n_points_in = n_row * n_col_in / res_fact / res_fact
n_points_out = (i_out_cut_max - i_out_cut_min + 1) * (j_out_cut_max -          &
               j_out_cut_min + 1) / res_fact / res_fact
write(*,'(a,i15)') 'Number of vertices in the input ".asc" file: ',n_points_in
write(*,'(a,i15)') 'Number of vertices in the output "xyz" file: ',n_points_out
! Pay attention: 
!    rows represent y, columns x
!    the ".asc" points are ordered from top-left to bottom-right: provided a 
!       given row / column, x increases / is constant and y is constant / 
!       decreases;  
do j_out=1,n_col_out,res_fact
   x_out = (j_out - 1) * dy + dy / 2.d0
   do i_out=1,n_row,res_fact
      y_out = (n_row + 1 - i_out) * dy - dy / 2.d0
      if (abs_mean_latitude>=0.d0) then
! Interpolation: inverse of the distance**2
         denom = 0.d0
         j_aux = nint((j_out - 0.5) * real(dy / dx) + 0.5)
         i_aux = i_out
         do j_in=(j_aux-1),(j_aux+1)
            x_in = (j_in - 1) * dx + dx / 2.d0
            do i_in=(i_aux-1),(i_aux+1)
               if ((i_in<1).or.(i_in>n_row).or.(j_in<1).or.(j_in>n_col_in))    &
                  cycle
               y_in = (n_row + 1 - i_in) * dy - dy / 2.d0
               distance = dsqrt((x_in - x_out) ** 2 + (y_in - y_out) ** 2)
               if (distance<=dy) then
                  mat_z_out(i_out,j_out) = mat_z_out(i_out,j_out) +            &
                                           mat_z_in(i_in,j_in) / distance ** 2
                  denom = denom + 1.d0 / distance ** 2
               endif
            enddo
         enddo
         if (denom/=0.d0) mat_z_out(i_out,j_out) = mat_z_out(i_out,j_out) /    &
                                                   denom
         else
! No interpolation (dx=dy)
            mat_z_out(i_out,j_out) = mat_z_in(i_out,j_out)
      endif
! Reservoir detection
      do i_bath=1,n_bathymetries
         test_integer = 0
         point(1) = x_out
         point(2) = y_out
         select case(n_vertices_around_res(i_bath))
            case(3)
               call point_inout_convex_non_degenerate_polygon(point,3,         &
                  vertices_around_res(i_bath,1,1:2),                           &
                  vertices_around_res(i_bath,2,1:2),                           &
                  vertices_around_res(i_bath,3,1:2),point,point,point,         &
                  test_integer)
            case(4)
               call point_inout_quadrilateral(point,                           &
                  vertices_around_res(i_bath,1,1:2),                           &
                  vertices_around_res(i_bath,2,1:2),                           &
                  vertices_around_res(i_bath,3,1:2),                           &
                  vertices_around_res(i_bath,4,1:2),test_integer)
            case(5)
               call point_inout_pentagon(point,                                &
                  vertices_around_res(i_bath,1,1:2),                           &
                  vertices_around_res(i_bath,2,1:2),                           &
                  vertices_around_res(i_bath,3,1:2),                           &
                  vertices_around_res(i_bath,4,1:2),                           &
                  vertices_around_res(i_bath,5,1:2),test_integer)
            case(6)
               call point_inout_hexagon(point,                                 &
                  vertices_around_res(i_bath,1,1:2),                           &
                  vertices_around_res(i_bath,2,1:2),                           &
                  vertices_around_res(i_bath,3,1:2),                           &
                  vertices_around_res(i_bath,4,1:2),                           &
                  vertices_around_res(i_bath,5,1:2),                           &
                  vertices_around_res(i_bath,6,1:2),test_integer)
         endselect
         if (test_integer>0) then
            if ((mat_z_out(i_out,j_out)<=(z_FS(i_bath)+z_eps(i_bath))).and.    &
               (mat_z_out(i_out,j_out)>=(z_FS(i_bath)-z_eps(i_bath)))) then
               reservoir(i_bath,i_out,j_out) = .true.
            endif
         endif
      enddo
! Detection of the digging/filling regions      
      do i_reg=1,n_digging_regions
         test_integer = 0
         point(1) = x_out
         point(2) = y_out
         select case(n_digging_vertices(i_reg))
            case(3)
               call point_inout_convex_non_degenerate_polygon(point,3,         &
                  digging_vertices(i_reg,1,1:2),digging_vertices(i_reg,2,1:2), &
                  digging_vertices(i_reg,3,1:2),point,point,point,test_integer)
            case(4)
               call point_inout_quadrilateral(point,                           &
                  digging_vertices(i_reg,1,1:2),digging_vertices(i_reg,2,1:2), &
                  digging_vertices(i_reg,3,1:2),digging_vertices(i_reg,4,1:2), &
                  test_integer)
            case(5)
               call point_inout_pentagon(point,                                &
                  digging_vertices(i_reg,1,1:2),digging_vertices(i_reg,2,1:2), &
                  digging_vertices(i_reg,3,1:2),digging_vertices(i_reg,4,1:2), &
                  digging_vertices(i_reg,5,1:2),test_integer)
            case(6)
               call point_inout_hexagon(point,digging_vertices(i_reg,1,1:2),   &
                  digging_vertices(i_reg,2,1:2),digging_vertices(i_reg,3,1:2), &
                  digging_vertices(i_reg,4,1:2),digging_vertices(i_reg,5,1:2), &
                  digging_vertices(i_reg,6,1:2),test_integer)
         endselect
         if (test_integer>0) then
            mat_z_out(i_out,j_out) = z_digging_regions(i_reg)
         endif
      enddo
      if ((i_out>=i_out_cut_min).and.(i_out<=i_out_cut_max).and.               &
         (j_out>=j_out_cut_min).and.(j_out<=j_out_cut_max)) then
         x_out = x_out + x_trans_out
         y_out = y_out + y_trans_out
         write(13,'(4(F15.4))') x_out,y_out,mat_z_out(i_out,j_out),            &
            mat_z_out(i_out,j_out)
         x_out = x_out - x_trans_out
         y_out = y_out - y_trans_out
      endif
   enddo
enddo
close(13)
if (n_bathymetries>0) then
   write(*,*) "Coastline detections. "
! Coastline detections
   do j_out=1,n_col_out,res_fact
      do i_out=1,n_row,res_fact
         do i_bath=1,n_bathymetries
            if (reservoir(i_bath,i_out,j_out).eqv..true.) then
               aux_integer = 0
               close_points: do j_close=j_out-1,j_out+1
                  do i_close=i_out-1,i_out+1
                     if ((i_close<1).or.(i_close>n_row).or.(j_close<1).or.     &
                        (j_close>n_col_out)) exit close_points
                     if (reservoir(i_bath,i_close,j_close).eqv..true.) then
                        aux_integer = aux_integer + 1
                        else
                           exit close_points
                     endif
                  enddo
               enddo close_points
               if (aux_integer<8) then
                  coastline(i_bath,i_out,j_out) = .true.
               endif
            endif
         enddo
      enddo
   enddo
! Bathymetry extrusions
   write(*,*) "Bathymetry extrusions. "
! Loop over the DEM output points
   do j_out=1,n_col_out,res_fact
      do i_out=1,n_row,res_fact
         do_extrusion: do i_bath=1,n_bathymetries
            if (reservoir(i_bath,i_out,j_out).eqv..true.) then
! To treat the inner reservoir points
               if (coastline(i_bath,i_out,j_out).eqv..true.) then
                  exit do_extrusion
               endif
! Position of current the inner reservoir point 
               point(1) = (j_out - 1) * dy + dy / 2.d0
               point(2) = (n_row + 1 - i_out) * dy - dy / 2.d0
! Distance (project on the horizontal) between a reservoir point and the line 
! passing for the upstream and the downstream points. 
! Unit vector perpendicular to the line. 
               call distance_point_line_2D(point,                              &
                  pos_res_downstream(i_bath,1:2),pos_res_upstream(i_bath,1:2), &
                  dis,normal)
! Position of a second point belonging to the line r_iC (beyond the centreline 
! "point")
               point_plus_normal(1) = point(1) + 1000.d0 * normal(1)
               point_plus_normal(2) = point(2) + 1000.d0 * normal(2)
! point_coast
               point_coast(:) = -9.d8
               min_dis2 = 9.d8
! Loop over the DEM output points
               do j2_out=1,n_col_out,res_fact
                  do i2_out=1,n_row,res_fact
                     if (coastline(i_bath,i2_out,j2_out).eqv..true.) then
! To treat the coast points
! Position of the generic coast point
                        point2(1) = (j2_out - 1) * dy + dy / 2.d0
                        point2(2) = (n_row + 1 - i2_out) * dy - dy / 2.d0
! Distance between the coast point and the line r_iC
                        call distance_point_line_2D(point2,point,              &
                           point_plus_normal,dis2,normal2)
! Distance (project on the horizontal) between the coastline point and the line 
! passing for the upstream and the downstream points. 
                        call distance_point_line_2D(point2,                    &
                           pos_res_downstream(i_bath,1:2),                     &
                           pos_res_upstream(i_bath,1:2),dis3,normal2)                        
! To update the position of the associated coast point
                        aux_scalar = dabs(dis2)
                        if ((aux_scalar<min_dis2).and.((dis*dis3)>=0.d0)) then
                           min_dis2 = aux_scalar
                           point_coast(:) = point2(:)
                        endif
                     endif
                  enddo
               enddo
               if (point_coast(1)<-8.9d8) then
                  write(*,*) 'The following inner reservoir point cannot be ', &
                     'associated to any coast point: ',point(1:2),'The ',      &
                     'program stops here. '
                  stop
               endif
! Pint
               call line_line_intersection_2D(pos_res_downstream(i_bath,1:2),  &
                  pos_res_upstream(i_bath,1:2),point,point_coast,Pint,         &
                  test_logical)
               if (test_logical.eqv..false.) then
                  write(*,*) 'Error. The intersection between the present ',   &
                     'two lines does not provide an unique point. The ',       &
                     'program stops here. '
                  write(*,*) 'pos_res_downstream(i_bath,1:2): ',               &
                     pos_res_downstream(i_bath,1:2)
                  write(*,*) 'pos_res_upstream(i_bath,1:2): ',                 &
                     pos_res_upstream(i_bath,1:2)
                  write(*,*) 'point(1:2): ',point(:)
                  write(*,*) 'point_coast(1:2): ',point_coast(:)
                  stop
               endif
! Bathymetry height at the intersection point Pint (linear interpolation)
               dis_Pdown_Pint = dsqrt((pos_res_downstream(i_bath,1) - Pint(1)) &
                                ** 2 + (pos_res_downstream(i_bath,2) -         &
                                Pint(2)) ** 2)
               z_Pint = z_downstream(i_bath) + dabs(dis_Pdown_Pint /           &
                        dis_down_up(i_bath)) * (z_FS(i_bath) -                 &
                        z_downstream(i_bath))
! Distance (projected along the horizontal) between the coast point and the 
! line r_down_up (or the point Pint)
               dis_Pint_Pcoast = dsqrt((point_coast(1) - Pint(1)) ** 2 +       &
                                 (point_coast(2) - Pint(2)) ** 2)
! Bathymetry height at the inner reservoir point (linear interpolation)
               mat_z_out(i_out,j_out) = z_Pint + min(dabs(dis /                &
                                        dis_Pint_Pcoast),1.d0) * (z_FS(i_bath) &
                                        - z_Pint)
! To update the estimated reservoir volume
               volume_res_est(i_bath) = volume_res_est(i_bath) + (z_FS(i_bath) &
                                        - mat_z_out(i_out,j_out)) * (dy ** 2)
               if (volume_flag(i_bath).eqv..true.) then
! Volume correction is active
! To compute the weight for the volume correction
                  if (weight_type(i_bath)==1) then
                     weight(i_bath,i_out,j_out) = max((mat_z_out(i_out,j_out)  &
                                                  - z_downstream(i_bath)),0.d0)
                     elseif (weight_type(i_bath)==2) then
                        aux_scalar_3 = dabs(dis_Pdown_Pint /                   &
                                       dis_down_up(i_bath))
                        if (aux_scalar_3<=0.5d0) then
                           aux_scalar = 0.d0
                           aux_scalar_2 = 1.d0
                           else
                              aux_scalar = 1.d0
                              aux_scalar_2 = -1.d0
                        endif
                        weight(i_bath,i_out,j_out) = 2.d0 * (aux_scalar +      &
                                                     aux_scalar_2 *            &
                                                     aux_scalar_3) * (1.d0 -   &
                                                     min(dabs(dis /            &
                                                     dis_Pint_Pcoast),1.d0))
                        write(*,*) "weight(i_bath,i_out,j_out): ",             &
                           weight(i_bath,i_out,j_out)
                        write(*,*) "aux_scalar: ",aux_scalar
                        write(*,*) "aux_scalar_2",aux_scalar_2
                        write(*,*) "aux_scalar_3",aux_scalar_3
                        write(*,*) "dis",dis
                        write(*,*) "dis_Pint_Pcoast",dis_Pint_Pcoast
                        else
                           write(*,*) "A reservoir volume correction is ",     &
                              "requested, but no admissible weight type is ",  &
                              "selected. The program stops here. "
                           stop 
                  endif
! To update the weight sum
                  weight_sum(i_bath) = weight_sum(i_bath) +                    &
                                       weight(i_bath,i_out,j_out)
               endif
            endif
         enddo do_extrusion
      enddo
   enddo
! Bathymetry correction on the actual reservoir volumes (for inner reservoir 
! points)
   write(*,*) "Possible bathymetry corrections on the actual reservoir volumes."
   do i_bath=1,n_bathymetries
      write(*,*) "Reservoir ",i_bath,": first volume estimation is ",          &
         volume_res_est(i_bath),"m**3 . "
      if (volume_flag(i_bath).eqv..true.) then
! Volume correction is active
         write(*,*) "Reservoir ",i_bath,": input/corrected volume is ",        &
            volume_res_inp(i_bath),"m**3 . "
         else
            write(*,*) "Reservoir ",i_bath,": no input/corrected volume is ",  &
            "provided . "
      endif
   enddo
   do j_out=1,n_col_out,res_fact
      do i_out=1,n_row,res_fact
         do_correction: do i_bath=1,n_bathymetries
            if (volume_flag(i_bath).eqv..true.) then
! Volume correction is active
! Bathymetry correction is applied in case a reservoir volume is provided in 
! input
               if (reservoir(i_bath,i_out,j_out).eqv..true.) then
! To treat inner reservoir points
                  if (coastline(i_bath,i_out,j_out).eqv..true.) then
                     exit do_correction
                  endif
                  mat_z_out(i_out,j_out) = mat_z_out(i_out,j_out) +            &
                                           (volume_res_est(i_bath) -           &
                                           volume_res_inp(i_bath)) *           &
                                           (weight(i_bath,i_out,j_out) /       &
                                           weight_sum(i_bath)) / (dy ** 2)
! To update the corrected reservoir volume
                  volume_res_corr(i_bath) = volume_res_corr(i_bath) +          &
                                           (z_FS(i_bath) -                     &
                                           mat_z_out(i_out,j_out)) * (dy ** 2)
               endif
            endif
         enddo do_correction
      enddo
   enddo
! Possible writing of the corrected reservoir volumes
   do i_bath=1,n_bathymetries
      if (volume_flag(i_bath).eqv..true.) then
! Volume correction is active
         write(*,*) "Reservoir ",i_bath,": corrected volume estimation is ",   &
            volume_res_corr(i_bath),"m**3 . "
      endif
   enddo
   write(*,*) "Possible writing of the xyz file (version after the ",          &
      "reservoir/batymetry extrusions). "
   open(14,file="xyz_with_extrusions.txt")
   write(14,'(a)') '           x(m)          y(m)           z(m)           z   '
   do j_out=1,n_col_out,res_fact
      do i_out=1,n_row,res_fact
         if ((i_out>=i_out_cut_min).and.(i_out<=i_out_cut_max).and.            &
            (j_out>=j_out_cut_min).and.(j_out<=j_out_cut_max)) then
            x_out = (j_out - 1) * dy + dy / 2.d0 + x_trans_out
            y_out = (n_row + 1 - i_out) * dy - dy / 2.d0 + y_trans_out
            write(14,'(4(F15.4))') x_out,y_out,mat_z_out(i_out,j_out),         &
               mat_z_out(i_out,j_out)
         endif
      enddo
   enddo
   close(14)
   open(15,file="weight_bathymetry_1.txt")
   write(15,'(a)') '           x(m)          y(m)           z(m)           z   '
   do j_out=1,n_col_out,res_fact
      do i_out=1,n_row,res_fact
         if ((i_out>=i_out_cut_min).and.(i_out<=i_out_cut_max).and.            &
            (j_out>=j_out_cut_min).and.(j_out<=j_out_cut_max)) then
            x_out = (j_out - 1) * dy + dy / 2.d0 + x_trans_out
            y_out = (n_row + 1 - i_out) * dy - dy / 2.d0 + y_trans_out
            write(15,'(4(F15.4))') x_out,y_out,weight(1,i_out,j_out),          &
               weight(1,i_out,j_out)
         endif
      enddo
   enddo
   close(15)
endif
!------------------------
! Deallocations
!------------------------
if(allocated(mat_z_in)) then
   deallocate(mat_z_in,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(0,*) 'Deallocation of "mat_z_in" failed; the execution terminates',&
         ' here. '
      stop
      else
         write(*,'(1x,a)') 'Deallocation of "mat_z_in" is successfully ',      &
            'completed. '
   endif
endif
if(allocated(mat_z_out)) then
   deallocate(mat_z_out,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(0,*) 'Deallocation of "mat_z_out" failed; the execution ',         &
         'terminates here. '
      stop
      else
         write(*,'(1x,a)') 'Deallocation of "mat_z_out" is successfully ',     &
            'completed. '
   endif
endif
if (n_digging_regions>0) then
   if(allocated(n_digging_vertices)) then
      deallocate(n_digging_vertices,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "n_digging_vertices" failed; the ',       &
            'execution terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "n_digging_vertices" is ',      &
               'successfully completed. '
      endif
   endif
   if(allocated(z_digging_regions)) then
      deallocate(z_digging_regions,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "z_digging_regions" failed; the ',        &
            'execution terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "z_digging_regions" is ',       &
               'successfully completed. '
      endif
   endif
   if(allocated(digging_vertices)) then
      deallocate(digging_vertices,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "digging_vertices" failed; the ',         &
            'execution terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "digging_vertices" is ',        &
               'successfully completed. '
      endif
   endif
endif
if (n_bathymetries>0) then
   if(allocated(n_vertices_around_res)) then
      deallocate(n_vertices_around_res,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "n_vertices_around_res" failed; the ',    &
            'execution terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "n_vertices_around_res" is ',   &
               'successfully completed. '
      endif
   endif
   if(allocated(z_downstream)) then
      deallocate(z_downstream,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "z_downstream" failed; the execution ',   &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "z_downstream" is ',            &
               'successfully completed. '
      endif
   endif
   if(allocated(z_FS)) then
      deallocate(z_FS,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "z_FS" failed; the execution terminates ',&
            'here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "z_FS" is successfully ',       &
               'completed. '
      endif
   endif
   if(allocated(z_eps)) then
      deallocate(z_eps,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "z_eps" failed; the execution terminates',&
            ' here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "z_eps" is successfully ',      &
               'completed. '
      endif
   endif
   if(allocated(volume_res_est)) then
      deallocate(volume_res_est,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "volume_res_est" failed; the execution ', &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "volume_res_est" is ',          &
               'successfully completed. '
      endif
   endif
   if(allocated(volume_res_inp)) then
      deallocate(volume_res_inp,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "volume_res_inp" failed; the execution ', &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "volume_res_inp" is ',          &
               'successfully completed. '
      endif
   endif
   if(allocated(volume_res_corr)) then
      deallocate(volume_res_corr,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "volume_res_corr" failed; the execution ',&
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "volume_res_corr" is ',         &
               'successfully completed. '
      endif
   endif
   if(allocated(weight_sum)) then
      deallocate(weight_sum,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "weight_sum" failed; the execution ',     &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "weight_sum" is successfully ', &
               'completed. '
      endif
   endif
   if(allocated(dis_down_up)) then
      deallocate(dis_down_up,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "dis_down_up" failed; the execution ',    &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "dis_down_up" is successfully ',&
               'completed. '
      endif
   endif
   if(allocated(pos_res_downstream)) then
      deallocate(pos_res_downstream,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "pos_res_downstream" failed; the ',       &
            'execution terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "pos_res_downstream" is ',      &
               'successfully completed. '
      endif
   endif
   if(allocated(pos_res_upstream)) then
      deallocate(pos_res_upstream,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "pos_res_upstream" failed; the execution',&
            ' terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "pos_res_upstream" is ',        &
               'successfully completed. '
      endif
   endif
   if(allocated(vertices_around_res)) then
      deallocate(vertices_around_res,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "vertices_around_res" failed; the ',      &
            'execution terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "vertices_around_res" is ',     &
               'successfully completed. '
      endif
   endif
   if(allocated(reservoir)) then
      deallocate(reservoir,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "reservoir" failed; the execution ',      &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "reservoir" is successfully ',  &
               'completed. '
      endif
   endif
   if(allocated(coastline)) then
      deallocate(coastline,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "coastline" failed; the execution ',      &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "coastline" is successfully ',  &
               'completed. '
      endif
   endif
   if(allocated(volume_flag)) then
      deallocate(volume_flag,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "volume_flag" failed; the execution ',    &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "volume_flag" is successfully ',&
               'completed. '
      endif
   endif
   if(allocated(weight)) then
      deallocate(weight,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "weight" failed; the execution ',         &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "weight" is successfully ',     &
               'completed. '
      endif
   endif
   if(allocated(weight_type)) then
      deallocate(weight_type,STAT=alloc_stat)
      if (alloc_stat/=0) then
         write(0,*) 'Deallocation of "weight_type" failed; the execution ',    &
            'terminates here. '
         stop
         else
            write(*,'(1x,a)') 'Deallocation of "weight_type" is successfully ',&
               'completed. '
      endif
   endif
endif
write(*,*) "DEM2xyz has terminated. "
end program DEM2xyz
