program main
   !! ## Estatística ERA5 para Furnas
   !!
   !! Author: Rodrigues, L.F. [LFR], Freitas, M.F. [MFF]
   !!
   !! E-mail: luiz.rodrigues@inpe.br, mateus.ffreitas@hotmail.com
   !!
   !! Date: 21Junho2023 11:59
   !!
   !! #####Version: 0.1.0
   !!
   !! ---
   !! **Full description**:
   !! Esse programa realiza a estatística necessária por Furnas para os dados obtidos
   !! de rodadas do modelo BRAMS
   !!
   !! ** History**:
   !!
   !! --- 
   !! ** Licence **: Under the terms of the GNU General Public version 3
   !!   <img src="https://www.gnu.org/graphics/gplv3-127x51.png width="63">
   !!
   use filesMod, only: t_mp  & !Type model param
                     , OpenNetCdfFile &
                     , readNamelist &
                     , readInfoNetCDF_fase1 &
                     , getDimsNetCDF &
                     , closeNetCDFFile &
                     , readNetCDFS &
                     , time_init &
                     , time_end &
                     , rshort &
                     , totpcp &
                     , spress &
                     , rh    &
                     , temp  &
                     , u_wind &
                     , v_wind &
                     , geo    &
                     , info2 &
                     , rshort_p &
                     , totpcp_p &
                     , spress_p &
                     , geo_p &
                     , temp_p &
                     , rh_p &
                     , u_wind_p &
                     , v_wind_p &
                     , readStations &
                     , getStationPos &
                     , getDimsVars &
                     , maguv &
                     , maguv_p &
                     , writeStationsCSV &
                     , mp &                     
                     , writePercNetcdfFile &
                     , writeStationsNetCDF &
                     , nhours &
                     , diruv &
                     , diruv_p

   use utilsMod, only: date_make_big &
                     , date_abs_secs

   use sortsMod, only: double_insertion_sort

   use modStatistic, only: quantil

   implicit none
   
   integer :: id1,ny,nYears, stat, i, y, j, perc, k, t
   character(len=4) :: yn, yni, ynf
   character(len=17) :: fSuffix
   character(len=14) :: outdate_in, outdate_fi
   real(kind=8) :: seconds1,seconds2
   integer :: info1(2)
   !! String of Year with 4 chars
   real(8), allocatable :: dados(:),dados_out(:)
   character(10) :: hora

   real :: startTime, stopTime

   call cpu_time(startTime)

   mp = readNamelist() !read the namelist with info about dates and prefix/suffix

   write(fSuffix,fmt='(I4,I2,I2"-",I4,I2,I2)') mp%year_i,mp%month_i,mp%day_i,mp%year_f,mp%month_f,mp%day_f
   write(yn,fmt='(I4)') mp%year_i

   call date_make_big (mp%year_i,mp%month_i,mp%day_i,0,outdate_in)
   call date_make_big (mp%year_f,mp%month_f,mp%day_f,230000,outdate_fi)
   write(*,fmt='(A,A," to ",A)') "Dates to compute - from ",outdate_in,outdate_fi
   call date_abs_secs(outdate_in,seconds1)
   call date_abs_secs(outdate_fi,seconds2)
   time_init = seconds1/3600
   time_end = seconds2/3600
   nhours = int(ceiling((seconds2-seconds1)/3600)) ! +1
   print *, time_init, " ", seconds1, " ", time_end, " ", seconds2 
   write(*,fmt='(A,I20)') "Number of hours: ",nhours

   nYears = mp%year_f-mp%year_i+1 !Verify number of years
   write(*,fmt='(A,I2)') "Number fo Years between datas: ",nYears

   !Reading the stations
   i = readStations("./estacoes.csv")

   id1 = OpenNetCdfFile(trim(mp%prefix)//"/era5-niveis_de_pressao-geopotential_"//yn//".nc")
   !Getting info1 from NetCDF
   ! 1 = ndims
   ! 2 = nvars
   info1 = readInfoNetCDF_fase1(id1)
   !Getting info2 from NetCDF
   ! 1 = nx
   ! 2 = ny
   ! 3 = nz
   ! 4 = ntimes
   !print *, "readInfoNetCDF", info1
   
   info2 = getDimsNetCDF(id1,info1(1))
   !print *, "!Info2:", info2

   allocate(rshort(info2(1),info2(2),time_init:time_end))
   allocate(totpcp(info2(1),info2(2),time_init:time_end))
   allocate(spress(info2(1),info2(2),time_init:time_end))
   allocate(geo   (info2(1),info2(2),info2(3),time_init:time_end))
   allocate(rh    (info2(1),info2(2),info2(3),time_init:time_end))
   allocate(temp  (info2(1),info2(2),info2(3),time_init:time_end))
   allocate(u_wind(info2(1),info2(2),info2(3),time_init:time_end))
   allocate(v_wind(info2(1),info2(2),info2(3),time_init:time_end))    

   i = getDimsVars(id1, info2, "s") 
   i = closeNetCDFFile(id1)

   i = getStationPos()

   i = 0
   do y=mp%year_i, mp%year_f
      i = i+1
      call readNetCDFS(y,mp)
   end do

   allocate(dados(size(geo,4)),dados_out(size(geo,4)))
   allocate(maguv(info2(1),info2(2),info2(3),time_init:time_end))
   allocate(diruv(info2(1),info2(2),info2(3),time_init:time_end))

   allocate(rshort_p(9,info2(1),info2(2)))
   allocate(totpcp_p(9,info2(1),info2(2)))
   allocate(spress_p(9,info2(1),info2(2)))
   allocate(maguv_p (9,info2(1),info2(2),info2(3)))
   allocate(diruv_p (9,info2(1),info2(2),info2(3)))
   allocate(geo_p   (9,info2(1),info2(2),info2(3)))
   allocate(temp_p  (9,info2(1),info2(2),info2(3)))
   allocate(u_wind_p(9,info2(1),info2(2),info2(3)))
   allocate(v_wind_p(9,info2(1),info2(2),info2(3)))
   allocate(rh_p    (9,info2(1),info2(2),info2(3)))

   !Making estatistics
   print *,""
   print *,"-----------------------------------"
   write(*,*) "Making Rshort Percentil"
   do i=1,info2(1)
      do j=1,info2(2)
         dados_out = rshort(i,j,:)
         !sorting data
         call double_insertion_sort(dados_out)
         !Making the percentil
         do perc=1,9
            rshort_p(perc,i,j) = quantil(dados_out, real(perc)/10.0)
         end do
      end do
   end do
   
   write(*,*) "Making totpcp Percentil"
   do i=1,info2(1)
      do j=1,info2(2)

         dados_out = totpcp(i,j,:)
         !sorting data
         call double_insertion_sort(dados_out)
         !Making the percentil
         do perc=1,9
            totpcp_p(perc,i,j) = quantil(dados_out, real(perc)/10.0)
         end do
      end do
   end do 

   write(*,*) "Making spress Percentil"
   do i=1,info2(1)
      do j=1,info2(2)

         dados_out = spress(i,j,:)
         !sorting data
         call double_insertion_sort(dados_out)
         !Making the percentil
         do perc=1,9
            spress_p(perc,i,j) = quantil(dados_out, real(perc)/10.0)
         end do
      end do
   end do 

   write(*,*) "Making geo Percentil"
   do i=1,info2(1)
      do j=1,info2(2)
         do k=1,info2(3)
            dados_out = geo(i,j,k,:)
            !sorting data
            call double_insertion_sort(dados_out)
            !Making the percentil
            do perc=1,9
               geo_p(perc,i,j,k) = quantil(dados_out, real(perc)/10.0)
            end do
         end do
      end do
   end do      

   write(*,*) "Making temp Percentil"
   do i=1,info2(1)
      do j=1,info2(2)
         do k=1,info2(3)
            dados_out = temp(i,j,k,:)
            !sorting data
            call double_insertion_sort(dados_out)
            !Making the percentil
            do perc=1,9
               temp_p(perc,i,j,k) = quantil(dados_out, real(perc)/10.0)
            end do
         end do
      end do
   end do

   write(*,*) "Making rh Percentil"
   do i=1,info2(1)
      do j=1,info2(2)
         do k=1,info2(3)
            dados_out = rh(i,j,k,:)
            !sorting data
            call double_insertion_sort(dados_out)
            !Making the percentil
            do perc=1,9
               rh_p(perc,i,j,k) = quantil(dados_out, real(perc)/10.0)
            end do
         end do
      end do
   end do

   write(*,*) "Computing maguv from winds"
   maguv = sqrt(u_wind**2+v_wind**2)

   write(*,*) "Computing diruv from winds"
   diruv = ((atan2(u_wind,v_wind)*180.0)/3.14159)

   write(*,*) "Making maguv Percentil"
   do i=1,info2(1)
      do j=1,info2(2)
         do k=1,info2(3)
            dados_out = maguv(i,j,k,:)
            !sorting data
            call double_insertion_sort(dados_out)
            !Making the percentil
            do perc=1,9
               maguv_p(perc,i,j,k) = quantil(dados_out, real(perc)/10.0)
            end do
         end do
      end do
   end do

   write(*,*) "Making diruv Percentil"
   do i=1,info2(1)
      do j=1,info2(2)
         do k=1,info2(3)
            dados_out = diruv(i,j,k,:)
            !sorting data
            call double_insertion_sort(dados_out)
            !Making the percentil
            do perc=1,9
               diruv_p(perc,i,j,k) = quantil(dados_out, real(perc)/10.0)
            end do
         end do
      end do
   end do


   i = writePercNetcdfFile(outdate_in,outdate_fi)
   i = writeStationsNetCDF()

   call cpu_time(stopTime)
   write(*, '(A, F18.6)') 'Elapsed time, s : ',  (stopTime - startTime)

end program main
