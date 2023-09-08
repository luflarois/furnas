module filesMod
   !! ## Manipulate files
   !!
   !! ![](https://i.ibb.co/LNqGy3S/logo-Monan-Color-75x75.png)
   !! ## MONAN
   !!
   !! Author: luiz F. Rodrigues [LFR]
   !!
   !! E-mail: <mailto:luizfrodrigues@protonmail.com>
   !!
   !! Date: 16Agosto2023 14:31
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !!
   !! Manipulate files
   !!
   !! ** History**:
   !!
   !! ---
   !! ** Licence **:
   !!
   !!  <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   !!  This program is free software: you can redistribute it and/or modify
   !!  it under the terms of the GNU General Public License as published by
   !!  the  Free  Software  Foundation, either version 3 of the License, or
   !!  (at your option) any later version.
   !!
   !!  This program is distributed in the hope that it  will be useful, but
   !!  ** WITHOUT  ANY  WARRANTY **;  without  even  the   implied   warranty  of
   !!  **MERCHANTABILITY** or **FITNESS FOR A  PARTICULAR PURPOSE**.  See  the, GNU
   !!  GNU General Public License for more details.
   !!
   !!  You should have received a copy  of the GNU General  Public  License
   !!  along with this program.  If not, see [GNU Public License](https://www.gnu.org/licenses/gpl-3.0.html).
   !!
   use netcdf
   use utilsMod, only: &
        date_secs_ymdt &
      , hypsometric

   implicit none
   include "netcdf.inc"
   !use modConstants, only: only_list
   character(len=*), parameter :: p_source_name = 'filesMod.F90'
   !! Source code name
   character(len=*), parameter :: p_module_name = 'filesMod'
   !! module name

   type t_mp
      character(len=256) :: prefix
      character(len=256) :: suffix
      integer :: year_i
      integer :: month_i
      integer :: day_i
      integer :: year_f
      integer :: month_f
      integer :: day_f
   end type t_mp
   type(t_mp) :: mp

   real, allocatable, dimension(:, :, :)   :: rshort !2d
   real, allocatable, dimension(:, :, :)   :: totpcp !2d
   real, allocatable, dimension(:, :, :)   :: spress !2d
   real, allocatable, dimension(:, :, :, :) :: rh     !3d
   real, allocatable, dimension(:, :, :, :) :: temp   !3d
   real, allocatable, dimension(:, :, :, :) :: u_wind !3d
   real, allocatable, dimension(:, :, :, :) :: v_wind !3d
   real, allocatable, dimension(:, :, :, :) :: geo    !3d

   real, allocatable, dimension(:, :, :)   :: rshort_p !2d
   real, allocatable, dimension(:, :, :)   :: totpcp_p !2d
   real, allocatable, dimension(:, :, :)   :: spress_p !2d
   real, allocatable, dimension(:, :, :, :) :: rh_p     !3d
   real, allocatable, dimension(:, :, :, :) :: temp_p   !3d
   real, allocatable, dimension(:, :, :, :) :: u_wind_p !3d
   real, allocatable, dimension(:, :, :, :) :: v_wind_p !3d
   real, allocatable, dimension(:, :, :, :) :: geo_p    !3d

   integer :: info3(4), info2(4), nhours
   real, allocatable :: levels(:), lat(:), lon(:), time_count(:)
   real, allocatable :: tp(:, :, :), ssrd(:, :, :), z(:, :, :, :)
   real, allocatable :: maguv(:, :, :, :), maguv_p(:, :, :, :)
   real, allocatable :: diruv(:, :, :, :), diruv_p(:, :, :, :)

   integer :: time_init
   integer :: time_end
   integer :: start_pos
   integer :: end_pos
   integer :: h_ini
   integer :: h_end
   integer :: h_count

   type t_est
      character(len=50) :: name
      real :: lat
      real :: lon
      integer :: plat
      integer :: plon
   end type t_est
   integer, parameter :: num_station = 50
   type(t_est) :: data_station(num_station)

   private
   public :: t_mp, OpenNetCdfFile, readNamelist, readInfoNetCDF_fase1, getDimsNetCDF, closeNetCDFFile
   public :: readNetCDFS, time_init, time_end, info2, rshort, totpcp, spress, rh, temp, u_wind, v_wind, geo
   public :: rshort_p, totpcp_p, spress_p, rh_p, temp_p, u_wind_p, v_wind_p, geo_p, readStations
   public :: getStationPos, getDimsVars, maguv, maguv_p, writeStationsCSV, mp, writePercNetcdfFile
   public :: writeStationsNetCDF, nhours, diruv, diruv_p

contains

#define check_nf90(x) check_nf_err(x, __LINE__, __FILE__)

   subroutine check_nf_err(iErr, line, filename)

      implicit none

      INTEGER, INTENT(IN) :: iErr
      INTEGER, INTENT(IN) :: line
      CHARACTER(LEN=*), INTENT(IN) :: filename

      if (iErr /= nf90_noerr) then
         print *, "NF90 ERROR: ", trim(nf90_strerror(iErr))
         print *, "At line ", line, " in file ", filename
         call exit(iErr)
      end if

   end subroutine check_nf_err

   function writeSeriesNetcdfFile(sn,lo,lt,rs,tp,mv,dv,r,t,sp, hv) result(ok)
      !! ## Write Netcdf file to percentis computed
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 25Agosto2023 14:40
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Write Netcdf file to percentis computed
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'writeSeriesNetcdfFile'

      !! Function name

      ! Variables (input):
      character(len=*), intent(in) :: sn
      ! Station name
      real, intent(in) :: lo
      !! Station longitude
      real, intent(in) :: lt
      !! Station latitude
      real, intent(in) :: rs(:)
      !! Short Radiation
      real, intent(in) :: sp(:)
      !! Surface pressure
      real, intent(in) :: tp(:)
      !! Total precipitation
      real, intent(in) :: mv(:,:)
      !! Wind speed
      real, intent(in) :: dv(:,:)
      !! Wind direction
      real, intent(in) :: r(:,:)
      !! Moisture
      real, intent(in) :: t(:,:)
      !! Temperature
      integer, intent(in) :: hv(:)

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
     
      integer :: ncid, lon_dimid, lat_dimid, lev_dimid, time_dimid, diruv_dimid
      integer :: temp_varid, maguv_varid, rh_varid, geo_varid, diruv_varid
      integer :: rshort_varid, tp_varid, sp_varid
      integer :: lon_varid, lat_varid, lev_varid,tim_varid
      integer :: dimids_2d(2), dimids_1d(1), dimids_lev(1), dimids_0d(1)
      integer :: iErr, perc
      character(len=10) :: perc_str
 
      ! Create a new netCDF file
      print *, "Creating "//sn//".nc file"
      call check_nf90(nf90_create(trim(mp%suffix)//"/"//sn//".nc", NF90_CLOBBER, ncid))
      ! Define dimensions
      !print *, "Defining dimensions"
      call check_nf90(nf90_def_dim(ncid, "longitude", 1, lon_dimid))
      call check_nf90(nf90_def_dim(ncid, "latitude", 1, lat_dimid))
      call check_nf90(nf90_def_dim(ncid, "level", size(levels), lev_dimid))
      call check_nf90(nf90_def_dim(ncid, "time", nhours+1, time_dimid))
      dimids_2d = (/ lev_dimid, time_dimid /)
      dimids_lev = (/ lev_dimid /)
      dimids_1d = (/ time_dimid /)
      dimids_0d = (/ 1 /)

         ! rshort_p !2d
         ! totpcp_p !2d
         ! spress_p !2d
         ! rh_p     !3d
         ! temp_p   !3d
         ! maguv_p  !3d
         ! geo_p    !3d

         !print *, "Defining variables"

         call check_nf90(nf90_def_var(ncid, "longitude", NF90_REAL, dimids_0d, lon_varid))
         call check_nf90(nf90_def_var(ncid, "latitude", NF90_REAL, dimids_0d, lat_varid))
         call check_nf90(nf90_def_var(ncid, "level", NF90_REAL, dimids_lev, lev_varid))
         call check_nf90(nf90_def_var(ncid, "time", NF90_REAL, dimids_1d, tim_varid))
         
         call check_nf90(nf90_def_var(ncid, "temp", NF90_REAL, dimids_2d, temp_varid))
         call check_nf90(nf90_def_var(ncid, "rh", NF90_REAL, dimids_2d, rh_varid))
         call check_nf90(nf90_def_var(ncid, "maguv", NF90_REAL, dimids_2d, maguv_varid))
         call check_nf90(nf90_def_var(ncid, "dirguv", NF90_REAL, dimids_2d, diruv_varid))
         call check_nf90(nf90_def_var(ncid, "rshort", NF90_REAL, dimids_1d, rshort_varid))
         call check_nf90(nf90_def_var(ncid, "precip", NF90_REAL, dimids_1d, tp_varid))
         call check_nf90(nf90_def_var(ncid, "spress", NF90_REAL, dimids_1d, sp_varid))

         call check_nf90(nf90_put_att(ncid, tim_varid, "standard_name", "time"))
         call check_nf90(nf90_put_att(ncid, tim_varid, "long_name", "time"))
         call check_nf90(nf90_put_att(ncid, tim_varid, "units", "hours since 1900-01-01 00:00:00.0"))
         call check_nf90(nf90_put_att(ncid, tim_varid, "calendar", "gregorian"))
         call check_nf90(nf90_put_att(ncid, tim_varid, "axis", "T"))

         ! End define mode
         call check_nf90(nf90_enddef(ncid))         

         ! Write data to file
         call check_nf90(nf90_put_var(ncid, lon_varid, lo))
         call check_nf90(nf90_put_var(ncid, lat_varid, lt))
         call check_nf90(nf90_put_var(ncid, lev_varid, levels(:)))
         call check_nf90(nf90_put_var(ncid, tim_varid, hv(:)))

         call check_nf90(nf90_put_var(ncid, temp_varid, t))
         call check_nf90(nf90_put_var(ncid, rh_varid, r))
         call check_nf90(nf90_put_var(ncid, maguv_varid, mv))
         call check_nf90(nf90_put_var(ncid, diruv_varid, dv))
         call check_nf90(nf90_put_var(ncid, rshort_varid, rs))
         call check_nf90(nf90_put_var(ncid, tp_varid, tp))
         call check_nf90(nf90_put_var(ncid, sp_varid, sp))

         ! Close the netCDF file
         call check_nf90(nf90_close(ncid))

      !print *, "Done series NC File"
      ok = 0

   end function writeSeriesNetcdfFile



   function writePercNetcdfFile(outdate_in, outdate_fi) result(ok)
      !! ## Write Netcdf file to percentis computed
      !!
      !! Author: Rodrigues, L.F. [LFR], Freitas, M.F. [MFF]
      !!
      !! E-mail: luiz.rodrigues@inpe.br, mateus.ffreitas@hotmail.com
      !!
      !! Date: 25Agosto2023 14:40
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Write Netcdf file to percentis computed
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'writePercNetcdfFile'

      !! Function name

      ! Variables (input):
      character(len=*), intent(in) :: outdate_in, outdate_fi

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
     
      integer :: ncid, lon_dimid, lat_dimid, lev_dimid, time_dimid
      integer :: temp_varid, maguv_varid, rh_varid, geo_varid, diruv_varid
      integer :: rshort_varid, tp_varid, sp_varid
      integer :: lon_varid, lat_varid, lev_varid
      integer :: dimids_3d(4), dimids_2d(3), dimids_lat(1), dimids_lon(1), dimids_lev(1), dimids_tim(1)
      integer :: iErr, perc
      character(len=10) :: perc_str
   
      do perc=1,9
         write(perc_str, '(I10)') perc*10
 
         ! Create a new netCDF file
         print *, "Creating P"//trim(adjustl(perc_str))//" NC file"
         call check_nf90(nf90_create(trim(mp%suffix)//"/P"//trim(adjustl(perc_str))//".nc", NF90_CLOBBER, ncid))

         ! Define dimensions
         print *, "Defining dimensions"
         call check_nf90(nf90_def_dim(ncid, "longitude", size(lon), lon_dimid))
         call check_nf90(nf90_def_dim(ncid, "latitude", size(lat), lat_dimid))
         call check_nf90(nf90_def_dim(ncid, "level", size(levels), lev_dimid))
         call check_nf90(nf90_def_dim(ncid, "time", 1, time_dimid))

         dimids_3d = (/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)
         dimids_2d = (/ lon_dimid, lat_dimid, time_dimid /)
         dimids_lat = (/ lat_dimid /)
         dimids_lon = (/ lon_dimid /)
         dimids_lev = (/ lev_dimid /)
         dimids_tim = (/ time_dimid /)

         ! rshort_p !2d
         ! totpcp_p !2d
         ! spress_p !2d
         ! rh_p     !3d
         ! temp_p   !3d
         ! maguv_p  !3d
         ! geo_p    !3d

         print *, "Defining variables"

         call check_nf90(nf90_def_var(ncid, "longitude", NF90_REAL, dimids_lon, lon_varid))
         call check_nf90(nf90_def_var(ncid, "latitude", NF90_REAL, dimids_lat, lat_varid))
         call check_nf90(nf90_def_var(ncid, "level", NF90_REAL, dimids_lev, lev_varid))
         
         call check_nf90(nf90_def_var(ncid, "temp", NF90_REAL, dimids_3d, temp_varid))
         call check_nf90(nf90_def_var(ncid, "geo", NF90_REAL, dimids_3d, geo_varid))
         call check_nf90(nf90_def_var(ncid, "rh", NF90_REAL, dimids_3d, rh_varid))
         call check_nf90(nf90_def_var(ncid, "maguv", NF90_REAL, dimids_3d, maguv_varid))
         call check_nf90(nf90_def_var(ncid, "diruv", NF90_REAL, dimids_3d, diruv_varid))
         call check_nf90(nf90_def_var(ncid, "rshort", NF90_REAL, dimids_2d, rshort_varid))
         call check_nf90(nf90_def_var(ncid, "precip", NF90_REAL, dimids_2d, tp_varid))
         call check_nf90(nf90_def_var(ncid, "spress", NF90_REAL, dimids_2d, sp_varid))

         ! End define mode
         call check_nf90(nf90_enddef(ncid))         

         ! Write data to file
         call check_nf90(nf90_put_var(ncid, lon_varid, lon(:)))
         call check_nf90(nf90_put_var(ncid, lat_varid, lat(:)))
         call check_nf90(nf90_put_var(ncid, lev_varid, levels(:)))

         call check_nf90(nf90_put_var(ncid, temp_varid, temp_p(perc,:,:,:)))
         call check_nf90(nf90_put_var(ncid, geo_varid, geo_p(perc,:,:,:)))
         call check_nf90(nf90_put_var(ncid, rh_varid, rh_p(perc,:,:,:)))
         call check_nf90(nf90_put_var(ncid, maguv_varid, maguv_p(perc,:,:,:)))
         call check_nf90(nf90_put_var(ncid, diruv_varid, diruv_p(perc,:,:,:)))
         call check_nf90(nf90_put_var(ncid, rshort_varid, rshort_p(perc,:,:)))
         call check_nf90(nf90_put_var(ncid, tp_varid, totpcp_p(perc,:,:)))
         call check_nf90(nf90_put_var(ncid, sp_varid, spress_p(perc,:,:)))

         ! Close the netCDF file
         call check_nf90(nf90_close(ncid))
      end do

      print *, "Done NC File"
      ok = 0

   end function writePercNetcdfFile


   function writeNetcdfFileTest() result(ok)
      !! ## Write Netcdf file to percentis computed
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 25Agosto2023 14:40
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Write Netcdf file to percentis computed
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'writeNetcdfFileTest'

      !! Function name

      ! Variables (input):

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
     
      integer :: ncid, lon_dimid, lat_dimid, lev_dimid, time_dimid
      integer :: temp_varid, maguv_varid, rh_varid, time_varid
      integer :: dimids(4)
      integer :: status
      real, dimension(14,3,7,12) :: temp_data, maguv_data, rh_data

      ! Initialize data arrays with some values
      temp_data = 0.0
      maguv_data = 1.0
      rh_data = 2.0

      ! Create a new netCDF file
      status = nf90_create("output.nc", NF90_CLOBBER, ncid)

      ! Define dimensions
      status = nf90_def_dim(ncid, "longitude", 14, lon_dimid)
      status = nf90_def_dim(ncid, "latitude", 3, lat_dimid)
      status = nf90_def_dim(ncid, "level", 7, lev_dimid)
      status = nf90_def_dim(ncid, "time", 12, time_dimid)

      status = nf90_def_var(ncid, "time", NF90_REAL, time_dimid, time_varid)
      ! Define variables
      dimids = (/ lon_dimid, lat_dimid, lev_dimid, time_dimid /)
      status = nf90_def_var(ncid, "temp", NF90_REAL, dimids, temp_varid)
      status = nf90_def_var(ncid, "maguv", NF90_REAL, dimids, maguv_varid)
      status = nf90_def_var(ncid, "rh", NF90_REAL, dimids, rh_varid)



      ! End define mode
      status = nf90_enddef(ncid)

      ! Write data to file
      status = nf90_put_var(ncid, temp_varid, temp_data)
      status = nf90_put_var(ncid, maguv_varid, maguv_data)
      status = nf90_put_var(ncid, rh_varid, rh_data)

      ! Close the netCDF file
      status = nf90_close(ncid)

   end function writeNetcdfFileTest

   function writeGradsFile() result(ok)
      !! ## Write grads file to percentis computed
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 25Agosto2023 14:40
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Write grads file to percentis computed
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'writeGradsFile'

      !! Function name

      ! Variables (input):

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      integer :: number, irec, recordLen, nvar, k
      integer :: output_byte_size, npercs, nvars, v, p
      real :: bytes_in_float
      !# bytes_in_float
      character(len=8), allocatable :: varn(:, :)
      character(len=2) :: pC
      character(len=17) :: comp_data

      inquire (iolength=output_byte_size) bytes_in_float
      recordLen = output_byte_size*size(lat)*size(lon)

      !write(comp_data,fmt='(I4.4,"-",I2.2,"-",I2.2,"-",I6.6)') dateTime%year &
      !,dateTime%month,dateTime%day,dateTime%hour

   end function writeGradsFile

   function writeStationsNetCDF() result(ok)
      !! ## Write the stations in a NetCDF File
      !!
      !! Author: Rodrigues, L. F. [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 08Setembro2023 08:54
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Write the stations in a NetCDF File
      !!
      !! ** History**:
      !!
      !! --- 
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!
   
      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'writeStationsNetCDF' 
   
      !! Function name
   
      ! Variables (input):
   
   
      ! Local variables:
      integer :: ok
      !! function Output: ok
      integer :: i, t, l
      integer :: hour_val(nhours+1)
      
      i = 0
      do t=time_init,time_end
         i = i+1
         hour_val(i) = t
      end do

      write(*,*) "=== Write NetCDF files temporal series. Total=",num_station
      do i = 1, num_station
         ok = writeSeriesNetcdfFile( trim(data_station(i)%name) &
         , data_station(i)%lon, data_station(i)%lat &
         , rshort(data_station(i)%plon, data_station(i)%plat,:) &
         , totpcp(data_station(i)%plon, data_station(i)%plat,:) &
         , maguv(data_station(i)%plon, data_station(i)%plat,:,:) &
         , diruv(data_station(i)%plon, data_station(i)%plat,:,:) &
         , rh(data_station(i)%plon, data_station(i)%plat,:,:) &
         , temp(data_station(i)%plon, data_station(i)%plat,:,:) &
         , spress(data_station(i)%plon, data_station(i)%plat,:) &
         , hour_val)
      end do
   
   end function writeStationsNetCDF


   function writeStationsCSV(outdate_in, outdate_fi) result(ok)
      !! ## Write the csv files with station info
      !!
      !! Author: Rodrigues, L. F.
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 18Agosto2023 16:04
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Write the csv files with station info
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'writeStationsCSV'

      !! Function name

      ! Variables (input):
      character(len=*), intent(in) :: outdate_in, outdate_fi

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      integer :: i, t, l
      character(len=256) :: filename, fileName_pf
      integer :: iyear1, imonth1, idate1, ihour1
      character(len=19) :: dates

      write(*,*) "=== Write CSV files. Total=",num_station
      do i = 1, num_station
         write (*, fmt='(A,A,2(F9.4,1X),2(I3.3,1X))') "Writing data for the station " &
            , trim(data_station(i)%name) &
            , data_station(i)%lat, data_station(i)%lon, data_station(i)%plat, data_station(i)%plon
         filename = trim(mp%suffix)//"/"//trim(data_station(i)%name)//"_s_"//outdate_in//"_"//outdate_fi//"DT1.0h.csv"
         fileName_pf = trim(mp%suffix)//"/"//trim(data_station(i)%name)//"_v_"//outdate_in//"_"//outdate_fi//"DT1.0h.csv"

         !Abre o arquivo de superfície
         !trim(mp%suffix)//"/"//
         open (unit=33, file=trim(filename), status="replace", action="write")
         write (33, fmt='("# Local: ",A,", LAT=",F9.5,", LON=",F9.5)') trim(data_station(i)%name) &
            , data_station(i)%lat, data_station(i)%lon
         write (unit=33, fmt='(A)') "# Sinal/Unidade: rshort [w/m2],precip [mm],v0m [m/s],v10m [m/s],rh [%],tempc [C],press [mb])"
         write(unit=33,fmt='(A)') "               Tempo,      rshort,      precip,         v0m,        v10m,          rh,       tempc,       press,"

         ! Abre o arquivo de ventos
         !trim(mp%suffix)//"/"//
         open (unit=34, file=trim(filename_pf), status="replace", action="write")
         write (34, fmt='("# Local: ",A,", LAT=",F9.5,", LON=",F9.5)') trim(data_station(i)%name) &
            , data_station(i)%lat, data_station(i)%lon
         write (unit=34, fmt='(A)') "# Sinal/Unidade: MAGUV, unidade 'm/s', Desc: 'Magnitude do Vento'"
         write (unit=34, fmt='(A)') "              Tempo,   0.00m,  10.00m,  50.00m,  75.00m, 100.00m, 150.00m, 200.00m, 300.00m,"

         !Laço no tempo
         do t = time_init, time_end
            call date_secs_ymdt(real(t, 8)*3600.0, iyear1, imonth1, idate1, ihour1)
            write (dates, fmt='(I4.4,"-",I2.2,"-",I2.2,"T",I2.2,":00:00")') iyear1, imonth1, idate1, ihour1/10000
            !Escreve superfície
            write (unit=33, fmt='(A19,",",7(F12.6,","))') dates &
               , rshort(data_station(i)%plon, data_station(i)%plat,t) &
               , totpcp(data_station(i)%plon, data_station(i)%plat,t) &
               , maguv(data_station(i)%plon, data_station(i)%plat,8,t) &
               , maguv(data_station(i)%plon, data_station(i)%plat,7,t) &
               , rh(data_station(i)%plon, data_station(i)%plat,8,t) &
               , temp(data_station(i)%plon, data_station(i)%plat,8,t) &
               , spress(data_station(i)%plon, data_station(i)%plat,t)
               ! do l=1,8
               !    write (88,*) data_station(i)%lat, data_station(i)%lon,l,levels(l) &
               !    ,spress(data_station(i)%plat,data_station(i)%plon,t) &
               !    ,temp(data_station(i)%plat,data_station(i)%plon,l,t) &
               !    ,hypsometric(spress(data_station(i)%plat,data_station(i)%plon,t) &
               !    ,levels(l) &
               !    ,temp(data_station(i)%plat, data_station(i)%plon,l,t))
               ! end do
            ! Escreve vento
            write (unit=34, fmt='(A19,",",8(F12.6,","))') dates, maguv(data_station(i)%plon, data_station(i)%plat,8:1:-1, t)
         end do
         close (unit=33)
         close (unit=34)
      end do

   end function writeStationsCSV

   function getStationPos() result(ok)
      !! ## Verify d=the x and y position os station in lat/lon
      !!
      !! Author: Luiz Flávio Rodrigues
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 18Agosto2023 15:35
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Verify d=the x and y position os station in lat/lon
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'getStationPos'

      !! Function name

      ! Variables (input):

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      integer :: i, j

      do i = 2, size(lat)
         do j = 1, num_station
            if (data_station(j)%lat <= lat(i - 1) .and. data_station(j)%lat >= lat(i)) then
               data_station(j)%plat = i
            end if
         end do
      end do
      do i = 2, size(lon)
         !if (i == 248) write (78, *) i, lon(i - 1), lon(i)
         do j = 1, num_station
            if (data_station(j)%lon >= lon(i - 1) .and. data_station(j)%lon <= lon(i)) then
               data_station(j)%plon = i
            end if
         end do
      end do
      write (*, fmt='(A)') "* Stations inventory:"
      do i = 1, num_station
         write (*, fmt='(A50,1X,2(I3,1X),2(F13.5,1X))') data_station(i)%name, data_station(i)%plat, data_station(i)%plon &
            , data_station(i)%lat, data_station(i)%lon
      end do
      write (*, fmt='(A)') "------------------------------------"

   end function getStationPos

   function readStations(stationFile) result(ok)
      !! ## Read the station file
      !!
      !! Author: Luiz F Rodrigues
      !!
      !! E-mail: <mailto:luisfrodrigues@protonmail.com>
      !!
      !! Date: 18Agosto2023 15:14
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Read the station file
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'readStations'

      !! Function name

      ! Variables (input):
      character(len=*) :: stationFile

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      character(len=50) :: dummy, line
      real :: lat, lon
      integer :: i

      open (unit=33, file=trim(stationFile), status="old", action="read")
      read (33, *) dummy
      do i = 1, num_station
         read (33, *) data_station(i)%name, data_station(i)%lat, data_station(i)%lon
      end do
      close (unit=33)

   end function readStations

   subroutine readNetCDFS(year, mp)
      !! ## read the NetCDF files for a especific year
      !!
      !! Author: Luiz Flavio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 16Agosto2023 21:32
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! read the NetCDF files for a especific year
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_name = 'readNetCDFS'
      !! subroutine name

      ! Variables (input, output, inout)
      integer, intent(in) :: year
      !! Year to be readed
      type(t_mp) :: mp

      ! Local variables:
      integer :: id, i
      integer :: info1(2)
      character(len=4) :: yn

      write (yn, fmt='(I4)') year

      write (*, *) "-------------------------------------"
      write (*, fmt='(A,I4)') "Processing year ", year
      print *, " "
      print *, " "

      id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_superficie-"//yn//".nc")
      !Getting info1 from NetCDF
      ! 1 = ndims
      ! 2 = nvars
      info1 = readInfoNetCDF_fase1(id)
      print *, "info1 = readInfoNetCDF_fase1"  ,info1(1), info1(2)
      info2 = getDimsNetCDF(id, info1(1))
      i = getDimsVars(id, info2, "s") !S is surface
      i = getVarNamesAndValues(id, info1(2), year)
      i = closeNetCDFFile(id)

      id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-geopotential_"//yn//".nc")
      info1 = readInfoNetCDF_fase1(id)
      info2 = getDimsNetCDF(id, info1(1))
      i = getDimsVars(id, info2, "p") !p is pressure
      i = getVarNamesAndValues(id, info1(2), year)
      i = closeNetCDFFile(id)

      id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-temperature_"//yn//".nc")
      info1 = readInfoNetCDF_fase1(id)
      info2 = getDimsNetCDF(id, info1(1))
      i = getDimsVars(id, info2, "p") !p is pressure
      i = getVarNamesAndValues(id, info1(2), year)
      i = closeNetCDFFile(id)

      id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-relative_humidity_"//yn//".nc")
      info1 = readInfoNetCDF_fase1(id)
      info2 = getDimsNetCDF(id, info1(1))
      i = getDimsVars(id, info2, "p") !p is pressure
      i = getVarNamesAndValues(id, info1(2), year)
      i = closeNetCDFFile(id)

      ! id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-temperature_"//yn//".nc")
      ! info1 = readInfoNetCDF_fase1(id)
      ! info2 = getDimsNetCDF(id, info1(1))
      ! i = getDimsVars(id, info2, "p") !p is pressure
      ! i = getVarNamesAndValues(id, info1(2), year)
      ! i = closeNetCDFFile(id)

      id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-u_component_of_wind_"//yn//".nc")
      info1 = readInfoNetCDF_fase1(id)
      info2 = getDimsNetCDF(id, info1(1))
      i = getDimsVars(id, info2, "p") !p is pressure
      i = getVarNamesAndValues(id, info1(2), year)
      i = closeNetCDFFile(id)

      id = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-v_component_of_wind_"//yn//".nc")
      info1 = readInfoNetCDF_fase1(id)
      info2 = getDimsNetCDF(id, info1(1))
      i = getDimsVars(id, info2, "p") !p is pressure
      i = getVarNamesAndValues(id, info1(2), year)
      i = closeNetCDFFile(id)

      print *, ""
      print *, ""

   end subroutine readNetCDFS

   function getVarNamesAndValues(id, nvars, year) result(ok)
      !! ## Get the names an positions of vars in NetCDF
      !!
      !! Author: Luiz Fláfio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 17Agosto2023 12:24
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Get the names an positions of vars in NetCDF
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'getVarNamesAndValues'

      !! Function name

      ! Variables (input):
      integer, intent(in) :: nvars
      integer, intent(in) :: id
      integer, intent(in) :: year

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      integer :: i, iErr, nat, an, t, result, ii, jj
      character(len=32) :: varName, atName
      real :: scaleFactor, addOfset
      !real, allocatable :: vteste(:, :, :)
      !integer :: iyear1,imonth1,idate1,ihour1
      !character(len=13) :: cd

      do i = 1, nvars
         iErr = nf90_Inquire_Variable(id, i, name=varName, nAtts=nat)
         call check_nf90(iErr)
         select case (trim(varName))
         case ("time")
            print *, '===================> passei no time', info2
            allocate (time_count(info2(4)))
            iErr = nf90_get_var(id, i, time_count)
            call check_nf90(iErr)
            start_pos = findloc(time_count, real(time_init), 1)
            end_pos = findloc(time_count, real(time_end), 1)
            if (start_pos == 0) start_pos = 1
            if (end_pos == 0) end_pos = size(time_count)
            print *,"start end " , start_pos, end_pos
            h_ini = int(time_count(start_pos))
            h_end = int(time_count(end_pos))
            h_count = h_end - h_ini + 1
         end select
      end do
      !allocate (vteste(size(lon), size(lat), time_init:time_end))
      !print *,h_ini,h_end
      !call date_secs_ymdt(real(h_ini+15,8)*3600.0,iyear1,imonth1,idate1,ihour1)
      !write(cd,fmt='(I4.4,"-",I2.2,"-",I2.2,"-",I2.2)') iyear1,imonth1,idate1,ihour1/10000
      !Get names of vars in netCDF file and total of atributes
      do i = 1, nvars
         iErr = nf90_Inquire_Variable(id, i, name=varName, nAtts=nat)
         call check_nf90(iErr)
         print *, "Variable ", i, "name ", varName, nat

         select case (trim(varName))
         case ("ssrd")
            write (*, fmt='(4(A,I7))') "reading ssrd/rshort from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
                  print *, 'scaleFactor=', scaleFactor
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
                  print *, 'addOfset=', addOfset
               end if
            end do
            !print *, info2(1), info2(2), start_pos
            iErr = nf90_get_var(id, i, rshort(:,:,h_ini:h_end), start=(/1, 1, start_pos/) &
                                , count=(/info2(1), info2(2),h_count/))
            call check_nf90(iErr)
            
            rshort(:, :, h_ini:h_end) = (rshort(:, :, h_ini:h_end)*scaleFactor + addOfset)/3600.0
             
            print *,'h_ini:h_end:',h_ini,h_end
             do ii=1, 1 !size(rshort,1)
                do jj=h_ini,h_end
                   print *, rshort(ii,1,jj)
                   !write(77,fmt='(A,1X,2(I5,1X,F8.2,1X),2(F16.5,1X))') &
                   !print*,cd,ii,lat(ii),jj,lon(jj),vteste(h_ini+15,ii,jj),rshort(h_ini+15,ii,jj)
                end do
             end do
             !stop "depois de ssrd"

         case ("tp")
            write (*, fmt='(4(A,I7))') "reading tp/totpcp from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
                  print *, 'scaleFactor=', scaleFactor
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
                  print *, 'addOfset=', addOfset
               end if
            end do
            iErr = nf90_get_var(id, i, totpcp(:, :, h_ini:h_end), start=(/1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), h_count/))
            call check_nf90(iErr)
            totpcp(:, :, h_ini:h_end) = totpcp(:, :, h_ini:h_end)*scaleFactor + addOfset

         case ("sp")
            write (*, fmt='(4(A,I7))') "reading sp/spress from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
                  print *, 'scaleFactor=', scaleFactor
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
                  print *, 'addOfset=', addOfset
               end if
            end do
            iErr = nf90_get_var(id, i, spress(:, :, h_ini:h_end), start=(/1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), h_count/))
            call check_nf90(iErr)
            spress(:, :, h_ini:h_end) = (spress(:, :, h_ini:h_end)*scaleFactor + addOfset)/100.0
   
         case ("z")
            write (*, fmt='(4(A,I7))') "reading z/geo from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
               end if
            end do
            
            iErr = nf90_get_var(id, i, geo(:, :, :, h_ini:h_end), start=(/1, 1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), info2(3), h_count/))
            call check_nf90(iErr)

            geo(:, :, :, h_ini:h_end) = geo(:, :, :, h_ini:h_end)*scaleFactor + addOfset

         case ("r")
            write (*, fmt='(4(A,I7))') "reading r/rh from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
               end if
            end do
            iErr = nf90_get_var(id, i, rh(:, :, :, h_ini:h_end), start=(/1, 1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), info2(3), h_count/))
            call check_nf90(iErr)

            rh(:, :, :, h_ini:h_end) = rh(:, :, :, h_ini:h_end)*scaleFactor + addOfset

         case ("t")
            write (*, fmt='(4(A,I7))') "reading t/temp from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
               end if
            end do
            iErr = nf90_get_var(id, i, temp(:, :, :, h_ini:h_end), start=(/1, 1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), info2(3), h_count/))
            call check_nf90(iErr)

            temp(:, :, :, h_ini:h_end) = (temp(:, :, :, h_ini:h_end)*scaleFactor + addOfset) - 273.15

         case ("u")
            write (*, fmt='(4(A,I7))') "reading u/u_wind from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
               end if
            end do
            iErr = nf90_get_var(id, i, u_wind(:, :, :, h_ini:h_end), start=(/1, 1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), info2(3), h_count/))
            call check_nf90(iErr)

            u_wind(:, :, :, h_ini:h_end) = u_wind(:, :, :, h_ini:h_end)*scaleFactor + addOfset

         case ("v")
            write (*, fmt='(4(A,I7))') "reading v/v_wind from ", int(time_count(start_pos)), " to ", int(time_count(end_pos)) &
               , " of data from ", time_init, " to ", time_end
            do an = 1, nat
               iErr = nf90_inq_attname(id, i, an, atName)
               call check_nf90(iErr)
               if (trim(atName) == 'scale_factor') then
                  iErr = nf90_get_att(id, i, atName, scaleFactor)
                  call check_nf90(iErr)
               elseif (trim(atName) == 'add_offset') then
                  iErr = nf90_get_att(id, i, atName, addOfset)
                  call check_nf90(iErr)
               end if
            end do
            iErr = nf90_get_var(id, i, v_wind(:, :, :, h_ini:h_end), start=(/1, 1, 1, start_pos/) &
                                , count=(/info2(1), info2(2), info2(3), h_count/))
            call check_nf90(iErr)

            v_wind(:, :, :, h_ini:h_end) = v_wind(:, :, :, h_ini:h_end)*scaleFactor + addOfset

         end select
      end do
      deallocate (time_count)

      ok = 1

   end function getVarNamesAndValues

   function readNamelist() result(mp)
      !! ## Read the namelist file
      !!
      !! Author: Luiz Flávio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 16Agosto2023 15:08
      !!
      !! #####Version: version
      !!
      !! ---
      !! **Full description**:
      !!
      !! Read the namelist file
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'readNamelist'
      !! Function name

      ! Variables (input):

      ! Local variables:
      type(t_mp) :: mp      !! function Output: mp
      integer :: yeari, monthi, dayi
      integer :: yearf, monthf, dayf
      character(len=256) :: prefix, suffix

      NAMELIST /MODEL/ prefix, suffix, yeari, monthi, dayi, yearf, monthf, dayf
      ! Code:

      open (unit=33, file="model_param.nml", status="old")
      read (33, MODEL)
      close (33)

      mp%prefix = prefix
      mp%suffix = suffix
      mp%year_i = yeari
      mp%month_i = monthi
      mp%day_i = dayi
      mp%year_f = yearf
      mp%month_f = monthf
      mp%day_f = dayf

   end function readNamelist

   function OpenNetCdfFile(fileName) result(id)
      !! ## Open a NetCDF File and return the ID
      !!
      !! Author: Luiz Fláfio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 16Agosto2023 14:41
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Open a NetCDF File and return the ID
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'OpenNetCdfFile'

      !! Function name

      ! Variables (input):
      character(len=*), intent(in) :: fileName

      ! Local variables:
      integer :: id
      !! function Output: id
      integer :: iErr, ncid

      ! Code:
      iErr = nf90_open(path=trim(fileName), mode=nf90_nowrite, ncid=ncid)
      call check_nf90(iErr)    
      id = ncid

   end function OpenNetCdfFile

   function closeNetCDFFile(id) result(ok)
      !! ## Close de file NetCDF
      !!
      !! Author: Luiz Fláfio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 16Agosto2023 21:19
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Close de file NetCDF
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'closeNetCDFFile'

      !! Function name

      ! Variables (input):
      integer, intent(in) :: id

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      integer :: iErr

      iErr = nf90_close(id)
      call check_nf90(iErr)
      ok = 0

   end function closeNetCDFFile

   function readInfoNetCDF_fase1(id) result(info1)
      !! ## read nvars and ndims of NetCDF file
      !!
      !! Author: Luiz Fláfio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 16Agosto2023 20:45
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! read nvars and ndims of NetCDF file
      !!
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'readInfoNetCDF_fase1'

      !! Function name

      ! Variables (input):
      integer, intent(in) :: id
      integer :: iErr, ndims, nvars

      ! Local variables:
      integer :: info1(2)
      !! function Output: info1

      ! get info about netCDF file
      iErr = nf90_inquire(id, ndims, nvars)
      call check_nf90(iErr)
      !print *, "ndims nvars ",ndims, nvars

      info1(1) = ndims
      info1(2) = nvars

   end function readInfoNetCDF_fase1

   function getDimsNetCDF(id, ndims) result(info2)
      !! ## Get the NetCDF dimensions
      !!
      !! Author: Luiz Fláfio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 16Agosto2023 20:57
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Get the NetCDF dimensions
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'getDimsNetCDF'

      !! Function name

      ! Variables (input):
      integer, intent(in) :: id
      integer, intent(in) :: ndims

      ! Local variables:
      integer :: info2(4)
      !! function Output: variable_out
      integer :: i, iErr
      character(len=32) :: name
      integer :: lenDim

      do i = 1, ndims
         iErr = nf90_Inquire_Dimension(id, i, name, lenDim)
         call check_nf90(iErr)
         !print *,'Dimension ',i,'name ',trim(name)
         select case (trim(name))
         case ('level')
            !print *, "level:",i,lenDim
            info3(3) = i
            info2(3) = lenDim
         case ('longitude')
            info3(1) = i
            info2(1) = lenDIm
         case ('latitude')
            info3(2) = i
            info2(2) = lenDIm
         case ('time')
            info3(4) = i
            info2(4) = lenDIm
         end select
      end do

   end function getDimsNetCDF

   function getDimsVars(id, info2, tipo) result(ok)
      !! ## Get the dimensional variables lon, lat, lev and time
      !!
      !! Author: Luiz Flávio Rodrigues [LFR]
      !!
      !! E-mail: <mailto:luizfrodrigues@protonmail.com>
      !!
      !! Date: 17Agosto2023 11:35
      !!
      !! #####Version: 0.1.0
      !!
      !! ---
      !! **Full description**:
      !!
      !! Get the dimensional variables lon, lat, lev and time
      !!
      !! ** History**:
      !!
      !! ---
      !! ** Licence **: Under the terms of the GNU General Public version 3
      !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
      !!

      implicit none
      ! Parameters:
      character(len=*), parameter :: p_procedure_ame = 'getDimsVars'

      !! Function name

      ! Variables (input):
      integer, intent(in) :: id
      integer, intent(in) :: info2(:)
      character, intent(in) :: tipo

      ! Local variables:
      integer :: ok
      !! function Output: variable_out
      integer :: iErr

      !print *, "GetDimsVars ", info2
      if (.not. allocated(levels)) then
          allocate (levels(info2(3)), lat(info2(2)), lon(info2(1)))
          !print *, "Allocated levels"
      end if
      
      if (tipo .ne. "s") then
         !Get pressure levels
         !print *, info3(3)
         iErr = nf90_get_var(id, info3(3), levels)
         call check_nf90(iErr)
         !print *, levels
      end if

      iErr = nf90_get_var(id, info3(1), lon)
      call check_nf90(iErr)
      iErr = nf90_get_var(id, info3(2), lat)
      call check_nf90(iErr)

      !print *, "!Info3:", info3

      ok = 0
   end function getDimsVars

end module filesMod
