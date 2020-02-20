!=================================================================== 
! Program reads aerosol profiles from MATCH-hourly netCDF files
! 
! Name: read_hourly_match.f90
!
! Modules used: netcdf
!
! Author: Ryan Scott, SSAI, ryan.c.scott@nasa.gov
!
! February 18, 2020
!===================================================================
! To do:
! 
! (0) Figure out how to compile this code on AMI...
! 
! (1) Instead of reading MATCH file directly from local dir, need 
!     to read the data the from appropriate ASDC_archive dir on AMI.
!     To do this, need to use PCF, etc. .......?
! 
! (2) Combine constituents & extract what's needed for Fu-Liou 
!
! (3) DONE: Match the aerosol data to specific CERES FOVs ?
!     This include HOUR time step and LAT/LON
!
!===================================================================
program read_hourly_match 
 
! use netcdf module - see README for compilation instructions
use netcdf

! no implicit variables
implicit none

! name of file to be read
character (len = *), parameter :: FILE_NAME = "CER_MATCH-hourly_Terra-Aqua-MODIS_Edition4_402402.20190101.nc"

! dimensions of data
integer, parameter :: nlon = 192, nlat = 94, ntime = 24, nlev = 28
real :: plev(nlev)
real :: lon(nlon), lat(nlat)

! ID for netcdf file and select variables
integer :: pvarid                ! pressure variable ID
integer :: latvarid, lonvarid    ! lat/lon variable ID
integer :: ncid                  ! nc file ID
integer :: i, j, p, k, t         ! loop indices

! Artificial data
! Extract MATCH data at the time and location of each SSF FOV
! CERES fov lat, lon, hr - these will need to be read in
real :: fov_lon = -145
real :: fov_lat = -20
integer :: fov_hr = 0

!=====================================================

! MATCH aerosol vertical profile (mavp) data structure
type mavptype
real :: aero_prof(11,nlon,nlat,nlev,ntime)  ! aerosol profile array, types 1-11
real :: aero_prof7(7,nlon,nlat,nlev,ntime)  ! aerosol profile array, types 1-11 combined into 7 categories
character*7 aero_type11(11)                 ! aerosol profile type label, 11 constituents
character*7 aero_type7(7)                   ! aerosol profile type label, 7 constituents obtained by combining categories
end type mavptype

type (mavptype) mpro                        ! MATCH profile data structure
 
! MATCH aerosol optical depth (AOD) data structure
type aodtype
real :: aod_array(12,nlon,nlat,ntime)       ! AOD array, total + 11 aerosol types
real :: aod_array8(8,nlon,nlat,ntime)       ! AOD array, total + 7 combined aerosol types
character*8 aod_type12(12)                  ! AOD type label, total + 11 aerosol types
character*8 aod_type8(8)                    ! AOD type label, total + 7 combined aerosol types
end type aodtype

type (aodtype) aod                          ! AOD data structure

integer :: iprofvarid(1:11)                 ! integer aerosol profile variable id, 11 aerosol types
integer :: iaodvarid(1:12)                  ! integer aod variable id, total + 11 aerosol types
       
!=====================================================

! aerosol type strings
data mpro%aero_type11 / &
'DSTQ01',        & !1                                                                                                                 
'DSTQ02',        & !2                                           
'DSTQ03',        & !3                                                                                                                  
'DSTQ04',        & !4                                                                                                               
'SO4'  ,         & !5                                                                                                                
'SSLT',          & !6                                                                                                              
'BCPHI',         & !7                                                                                                              
'BCPHO',         & !8                                                                                             
'OCPHI',         & !9
'OCPHO',         & !10
'VOLC' /           !11

data aod%aod_type8 / &
'Total',   &    
'DustSm',  &
'DustLg',  &
'SO4',     &
'SSLT',    &
'Soot',    &
'Solub',   &
'Insol' /     

data aod%aod_type12 / &
'AEROD',     &
'DSTODX01',  &
'DSTODX02',  &
'DSTODX03',  &
'DSTODX04',  &
'SO4OD',     &
'SSLTOD',    &
'BCPHIOD',   &
'BCPHOOD',   &
'OCPHIOD',   &
'OCPHOOD',   &
'VOLCOD' /

!==============================================

! open the netcdf file in read mode
call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

! get the id for each variable based on its name
call check( nf90_inq_varid(ncid, "lon", lonvarid) )
call check( nf90_inq_varid(ncid, "lat", latvarid) )
call check( nf90_inq_varid(ncid, "lev", pvarid) )

! print nc file and variable ids
print*, "===================================="
print*, "  NetCDF file and variable IDs...   "
print*, "===================================="

! read the data
call check( nf90_get_var(ncid,lonvarid, lon) )
call check( nf90_get_var(ncid,latvarid, lat) )
call check( nf90_get_var(ncid,pvarid, plev) )

!********************************************************************************************************
! dimensions are reversed in the netCDF file - this line is old but left here FYI
!call check( nf90_get_var(ncid, varid1, data3d, start=(/1,1,1/), count=(/nlon,nlat,ntime/) ) )
!call check( nf90_get_var(ncid, varid2, data4d, start=(/1,1,1,1/), count=(/nlon,nlat,nlev,ntime/) ) )
!********************************************************************************************************

! read aerosol vertical profiles
do i = 1,11
  ! get variable ids for each aerosol-type profile
  call check( nf90_inq_varid(ncid, mpro%aero_type11(i), iprofvarid(i) ) )
  print*, "Getting...", mpro%aero_type11(i)," variable id...", iprofvarid(i)
  ! get vertical aerosol profiles for each type and store in single 5-d array
  call check( nf90_get_var(ncid, iprofvarid(i), mpro%aero_prof(i,:,:,:,:), start=(/1,1,1,1/), count=(/nlon,nlat,nlev,ntime/) ) )
end do

! read AODs
do i = 1,12
  ! get variable ids for each aerosol type OD
  call check( nf90_inq_varid(ncid, aod%aod_type12(i), iaodvarid(i) ) )
  print*, "Getting...", aod%aod_type12(i),"variable id...", iaodvarid(i) 
  ! get aerosols optical depths for each type and store in single 4-d array
  call check( nf90_get_var(ncid, iaodvarid(i), aod%aod_array(i,:,:,:), start=(/1,1,1/), count=(/nlon,nlat,ntime/) ) )
end do

print*, "===================================="
print*, " Successfully read 4D & 5D data..."
print*, " AOD(type,nlon,nlat,ntime)          "
print*, " Profile(type,nlon,nlat,nlev,ntime) "
print*, "===================================="

print*, "===================================="
print*, " Printing aerosol optical depth...  "
print*, "===================================="

print*, " Suppressed Output..."
!do t = 1,12
!  do i = 1,nlon
!    do j = 1,
!      print*,"data3d(",i,j,",1) = ", data3d(i,j,1)
!    end do
!  end do
!end do 

print*, "===================================="
print*, "Printing aerosol vertical profiles..."
print*, "===================================="


do t = 1,11            ! aerosol type
  do i = 31,31         ! longitude
     do j = 31,31      ! latitude
        do k = 1,1     ! time
          print*, "============================================================================"
          print*, "Hour:",k, "Type:  ",mpro%aero_type11(t), "VarID:  ",iaodvarid(t)
          print*, "============================================================================"
          do p=1,nlev  ! pressure levels
            print*, "Lon:", lon(i),"Lat:", lat(j), "Pres:",plev(p), "Aero:", mpro%aero_prof(t,i,j,p,k)
          end do
        end do
     end do
  end do
end do


print*, "======================================="
print*, "Combining AODs into 7 types + total AOD"
print*, "======================================="

aod%aod_array8(1,:,:,:) = SUM(aod%aod_array(2:12,:,:,:),1) !Total - should equal AEROOD
aod%aod_array8(2,:,:,:) = aod%aod_array(2,:,:,:)           !DustSm
aod%aod_array8(3,:,:,:) = SUM(aod%aod_array(3:5,:,:,:),1)  !DustLg
aod%aod_array8(4,:,:,:) = aod%aod_array(6,:,:,:)           !SO4
aod%aod_array8(5,:,:,:) = aod%aod_array(7,:,:,:)           !SSLT - Sea Salt
aod%aod_array8(6,:,:,:) = SUM(aod%aod_array(8:9,:,:,:),1)  !SOOT - Black Carbon
aod%aod_array8(7,:,:,:) = aod%aod_array(10,:,:,:)          !Solub
aod%aod_array8(8,:,:,:) = aod%aod_array(11,:,:,:)          !Insolub

!print*, "Computed total AOD"
!print*, aod%aod_array8(1,1:10,1,1)
!print*, "Read from file"
!print*, aod%aod_array(1,1:10,1,1)

print*, "======================================"

print*, "MATCH Latitude centers"
print*, lat

print*, "MATCH Longitude centers"
print*, lon

print*, "MATCH Pressure levels"
print*, plev

print*, "======================================"
print*, "Matching MATCH to CERES footprints... "
print*, "======================================"


! CERES FOV hr is offset from MATCH by 1
fov_hr = fov_hr + 1

print*, "Artificial CERES FOV lon:", fov_lon
print*, "Artificial CERES FOV lon:", fov_lat

print*, "MATCH longitude index, MATCH longitude"
print*, lon_match_ceresfov(fov_lon), lon(lon_match_ceresfov(fov_lon))

print*, "MATCH latitude index, MATCH latitude "
print*, lat_match_ceresfov(fov_lat), lat(lat_match_ceresfov(fov_lat))

print*, "======================================"

! Get the aerosol profiles aods matched to the CERES FOV

print*, "Getting aerosol profile at CERES FOV location..."

do t = 1,11              ! aerosol type                                                                                              
  do k = fov_hr,fov_hr   ! time
    print*, "============================================================================"
    print*, "Hour:",k, "Type:  ",mpro%aero_type11(t), "VarID:  ",iaodvarid(t)
    print*, "============================================================================"
    do p=1,nlev  ! pressure levels                                                                                              
      print*, "Lon:", lon(lon_match_ceresfov(fov_lon)),"Lat:", lat(lat_match_ceresfov(fov_lat)), "Pres:",plev(p), "Aero:", &
                      mpro%aero_prof(t,lon_match_ceresfov(fov_lon),lat_match_ceresfov(fov_lat),p,k)
    end do
  end do
end do

! close the file, freeing all resources
call check( nf90_close(ncid) )

!================================================================

contains

!**********************************************************                                                                     
! SUBROUTINE:                                                                                                                          
! Check netcdf-fortran functions for errors
!**********************************************************
subroutine check(status)
integer, intent (in) :: status

if(status /= nf90_noerr) then
    print *, trim(nf90_strerror(status))
    stop "Stopped"
    end if
end subroutine check

!**********************************************************
! FUNCTION:
! Calculate CERES FOV & MATCH long distance, match by index
!********************************************************** 
function lon_match_ceresfov(fovlon) result(ilon)
integer :: i            !loop index
integer :: ilon         !min dist index
real :: dist_lon(nlon)  !distance array
real :: fovlon          !fov lon

! if CERES FOV lon < 0 convert range to 0-360
if (fovlon<0) fovlon = fovlon+360
! compute distance between CERES FOV & MATCH lon centers
dist_lon = abs(lon-fovlon)
! find the minimum distance and extract index
do i=1,nlon
  if (dist_lon(i) == minval(dist_lon)) ilon = i 
end do
end function lon_match_ceresfov

!********************************************************** 
! FUNCTION                                                                        
! Calculate CERES FOV & MATCH lat distance, match by index
!**********************************************************                                                                           
function lat_match_ceresfov(fovlat) result(ilat)
integer :: i            !loop index      
integer :: ilat         !index of profile closest to fov                                                                          
real :: dist_lat(nlat)  !distance array                                                                                             
real :: fovlat          !fov lat                                                                                                      
! compute distance between CERES FOV & MATCH lat centers                                                                        
dist_lat = abs(lat-fovlat)
! find the minimum distance and extract index                                                                                         
do i=1,nlat
  if (dist_lat(i) == minval(dist_lat)) ilat = i
end do
end function lat_match_ceresfov


!********************************************************** 
end program read_hourly_match









