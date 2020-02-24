!=================================================================== 
!
! Name: read_hourly_match.f90
!
! Purpose: Program reads aerosol vertical profiles (kg/kg) and AOD
!          from hourly MATCH netCDF files. Profiles are reported on
!          pressure layers, not pressure levels.
!
! Modules Used: netcdf
!
! Author:  Ryan Scott, SSAI, ryan.c.scott@nasa.gov
! Updated: February 24, 2020
!
!===================================================================
! To do:
! 
! (0) Figure out how to compile this script on AMI...
! 
! (1) Instead of reading MATCH file directly from local dir, need 
!     to read the data the from appropriate ASDC_archive dir on AMI.
!     To do this, need to use pcf module, etc.
!     Fred knows these ins and outs... 
!
! (2) DONE: Match the aerosol data to specific CERES FOVs...
!     This include HOUR time step and LAT/LON
!     This will need to be modified to read this info from SSF data
!
! (3) IN PROGRESS: Process vertical profiles for input to Fu-Liou...
!                  Combine constituents                 
!                  Convert profile into AOT profile
!                  Re-grid onto Fu-Liou model profile 
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
integer, parameter :: nlon = 192, nlat = 94, ntime = 24, nlay = 28, nlev = nlay+1
real :: play(nlay), plev(nlev)
real :: lon(nlon), lat(nlat)

! ID for netcdf file and select variables
integer :: iplayvarid, iplevvarid  ! integer pressure lay/lev variable ID
integer :: ilatvarid, ilonvarid    ! integer lat/lon variable ID

integer :: iprofvarid(1:11)      ! integer aerosol profile variable id, 11 aerosol types                                 
integer :: iaodvarid(1:12)       ! integer aod variable id, total + 11 aerosol types 

integer :: ncid                  ! nc file ID
integer :: i, j, p, k, t         ! loop indices


!===========================================  
! MATCH aerosol vertical profile (mavp) data structure definition
type mavptype
real :: aero_prof(11,nlon,nlat,nlay,ntime)  ! aerosol profile array, types 1-11
real :: aero_prof7(7,nlon,nlat,nlay,ntime)  ! aerosol profile array, types 1-11 combined into 7 categories
character*7 aero_type11(11)                 ! aerosol profile type label, 11 constituents
character*10 aero_type7(7)                  ! aerosol profile type label, 7 constituents obtained by combining categories
end type mavptype

type (mavptype) mavp                        ! MATCH profile data structure variable

!=========================================== 
! MATCH aerosol optical depth (AOD) data structure definition
type aodtype
real :: aod_array(12,nlon,nlat,ntime)       ! AOD array, total + 11 aerosol types
real :: aod_array8(8,nlon,nlat,ntime)       ! AOD array, total + 7 combined aerosol types
character*8 aod_type12(12)                  ! AOD type label, total + 11 aerosol types
character*8 aod_type8(8)                    ! AOD type label, total + 7 combined aerosol types
end type aodtype

type (aodtype) aod                          ! AOD data structure variable

!===========================================  

! Artificial data (temporary)

! CERES fov lat, lon, hr - these will need to be read in
real :: fov_lon
real :: fov_lat
integer :: fov_hr
integer ::  mfov_hr     ! match fov hr = ceres fov hr + 1 

! Prompt user for artificial CERES FOV and hr
! This will need to be replaced by actually reading this info

write(*,*) "Enter CERES FOV longitude:"
read(*,*) fov_lon
write(*,*) "Enter CERES FOV latitude:"
read(*,*) fov_lat
write(*,*) "Enter CERES hour:"
read(*,*) fov_hr

!=====================================================

! aerosol type strings
data mavp%aero_type11 / &
'DSTQ01',        & !1  Dust small
'DSTQ02',        & !2  Dust med-small                                          
'DSTQ03',        & !3  Dust med-large                                                                                          
'DSTQ04',        & !4  Dust large                                                                                                     
'SO4'  ,         & !5  Sulfate                    
'SSLT',          & !6  Sea salt, d'Almedia maritime                                                                           
'BCPHI',         & !7  Hydrophilic black carbon                                                                               
'BCPHO',         & !8  Hydrophobic black carbon
'OCPHI',         & !9  Hydrophilic organic carbon
'OCPHO',         & !10 Hydrophobic organic carbon
'VOLC' /           !11 Volcanic

data mavp%aero_type7 / &
'DustSm',     &    ! DSTQ01
'DustLg',     &    ! DSTQ02+DSTQ03+DSTQ04
'OPAC SUSO',  &    ! VOLC + SO4(strato)
'SSLT',       &    ! d'Almedia maritime
'OPAC SOOT',  &    ! BCPHI+BCPHO
'OPAC WASO',  &    ! OCPHI + SO4(tropo)
'OPAC INSO' /      ! OCPHO

data aod%aod_type8 / &
'Total',   &       ! Total AOD
'DustSm',  &       ! Dust small 
'DustLg',  &       ! Dust large (med+large)
'SO4',     &       ! Sulfate
'SSLT',    &       ! Sea Salt
'Soot',    &       ! Black carbon
'Solub',   &       ! Organic carbon 
'Insol' /          ! Organic 

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
call check( nf90_inq_varid(ncid, "lon", ilonvarid) )   ! lon var id
call check( nf90_inq_varid(ncid, "lat", ilatvarid) )   ! lat var id
call check( nf90_inq_varid(ncid, "ilev", iplevvarid) ) ! lev var id - nc files poorly named
call check( nf90_inq_varid(ncid, "lev" , iplayvarid) ) ! lay var id


! print nc file and variable ids
print*, "===================================="
print*, "  NetCDF file and variable IDs...   "
print*, "===================================="

! read the data
call check( nf90_get_var(ncid,ilonvarid, lon) )
call check( nf90_get_var(ncid,ilatvarid, lat) )
call check( nf90_get_var(ncid,iplevvarid, plev) )
call check( nf90_get_var(ncid,iplayvarid, play) )


!********************************************************************************************************
! dimensions are reversed in the netCDF file - this code is old but left here as an FYI
!call check( nf90_get_var(ncid, varid1, data3d, start=(/1,1,1/), count=(/nlon,nlat,ntime/) ) )
!call check( nf90_get_var(ncid, varid2, data4d, start=(/1,1,1,1/), count=(/nlon,nlat,nlev,ntime/) ) )
!********************************************************************************************************

! read aerosol vertical profiles
do i = 1,11
  ! get variable ids for each aerosol-type profile
  call check( nf90_inq_varid(ncid, mavp%aero_type11(i), iprofvarid(i) ) )
  print*, "Getting...", mavp%aero_type11(i)," variable id...", iprofvarid(i)
  ! get vertical aerosol profiles for each type and store in single 5-d array
  call check( nf90_get_var(ncid, iprofvarid(i), mavp%aero_prof(i,:,:,:,:), start=(/1,1,1,1/), count=(/nlon,nlat,nlay,ntime/) ) )
end do

! read AODs
do i = 1,12
  ! get variable ids for each aerosol type OD
  call check( nf90_inq_varid(ncid, aod%aod_type12(i), iaodvarid(i) ) )
  print*, "Getting...", aod%aod_type12(i),"variable id...", iaodvarid(i) 
  ! get aerosol optical depths for each type and store in single 4-d array
  call check( nf90_get_var(ncid, iaodvarid(i), aod%aod_array(i,:,:,:), start=(/1,1,1/), count=(/nlon,nlat,ntime/) ) )
end do

print*, "===================================="
print*, " Successfully read 4D & 5D data..."
print*, " AerosOD(type,nlon,nlat,ntime)          "
print*, " Profile(type,nlon,nlat,nlev,ntime) "
print*, "===================================="

print*, "===================================="
print*, " Printing aerosol optical depth...  "
print*, "===================================="

print*, " Suppressed output..."
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

print*, "Suppressed output..."
!do t = 1,11            ! aerosol type
!  do i = 31,31         ! longitude
!     do j = 31,31      ! latitude
!        do k = 1,1     ! time
!          print*, "========================================================================================="
!          print*, "Hour:",k, "Type:  ",mavp%aero_type11(t), "VarID:  ",iaodvarid(t)
!          print*, "========================================================================================="
!          do p=1,nlay  ! pressure layers
!            print*, "Lon:", lon(i),"Lat:", lat(j), "Pres:",play(p), "Aero:", mavp%aero_prof(t,i,j,p,k)
!          end do
!        end do
!     end do
!  end do
!end do


print*, "======================================="
print*, "Combining Ed4 aerosol profiles...      "
print*, "======================================="

! 7 species
mavp%aero_prof7(1,:,:,:,:) = mavp%aero_prof(1,:,:,:,:)                                  ! 1. DustSm = DSTQ01
mavp%aero_prof7(2,:,:,:,:) = sum(mavp%aero_prof(2:4,:,:,:,:),1)                         ! 2. DustLg = DSTQ02+DSTQ03+DSTO4
mavp%aero_prof7(4,:,:,:,:) = mavp%aero_prof(6,:,:,:,:)                                  ! 4. SSLT
mavp%aero_prof7(5,:,:,:,:) = sum(mavp%aero_prof(7:8,:,:,:,:),1)                         ! 5. OPAC SOOT = BCPHI+BCPHO
mavp%aero_prof7(7,:,:,:,:) = mavp%aero_prof(10,:,:,:,:)                                 ! 7. OPAC INSO

! assumes 200 mb tropopause
do p=1,nlay
  ! stratosphere
  if (p <= 10) then
    mavp%aero_prof7(3,:,:,p,:) = mavp%aero_prof(5,:,:,p,:) + mavp%aero_prof(11,:,:,p,:) ! 3. OPAC SUSO = SO4(strato)+VOLC
    mavp%aero_prof7(6,:,:,p,:) = mavp%aero_prof(9,:,:,p,:)                              ! 6. OPAC WASO = OCPHI
  ! troposphere
  else if (p > 10) then
    mavp%aero_prof7(3,:,:,p,:) = mavp%aero_prof(11,:,:,p,:)                             ! 3. OPAC SUSO = VOLC only - 0 in troposphere
    mavp%aero_prof7(6,:,:,p,:) = mavp%aero_prof(9,:,:,p,:) + mavp%aero_prof(5,:,:,p,:)  ! 6. OPAC WASO = OCPHI + SO4(tropo)
  end if
end do

print*, "Done combining profiles..."

print*, "======================================="
print*, "Combining AODs into 7 types + total AOD"
print*, "======================================="

! total + 7 species
aod%aod_array8(1,:,:,:) = sum(aod%aod_array(2:12,:,:,:),1) ! Total - should equal AEROOD
aod%aod_array8(2,:,:,:) = aod%aod_array(2,:,:,:)           ! DustSm
aod%aod_array8(3,:,:,:) = sum(aod%aod_array(3:5,:,:,:),1)  ! DustLg
aod%aod_array8(4,:,:,:) = aod%aod_array(6,:,:,:)           ! SO4
aod%aod_array8(5,:,:,:) = aod%aod_array(7,:,:,:)           ! SSLT - Sea Salt
aod%aod_array8(6,:,:,:) = sum(aod%aod_array(8:9,:,:,:),1)  ! SOOT - Black Carbon
aod%aod_array8(7,:,:,:) = aod%aod_array(10,:,:,:)          ! Solub
aod%aod_array8(8,:,:,:) = aod%aod_array(11,:,:,:)          ! Insolub

! This checks out to ~4 decimal places
! Printing out only 10 numbers...
!print*, "Computed total AOD"
!print*, aod%aod_array8(1,1:10,1,1)
!print*, "Read from file"
!print*, aod%aod_array(1,1:10,1,1)

print*, "======================================="

print*, "MATCH Latitude Centers"
print*, lat

print*, "MATCH Longitude Centers"
print*, lon

print*, "MATCH Pressure Levels"
print*, plev

print*, "MATCH Pressure Layers"
print*, play

print*, "======================================"
print*, "Matching MATCH data to CERES FOV...   "
print*, "======================================"

! CERES FOV hr is offset from MATCH by 1
mfov_hr = fov_hr + 1

print*, "Artificial CERES Hour   :", fov_hr
print*, "Artificial CERES FOV lon:", fov_lon
print*, "Artificial CERES FOV lon:", fov_lat

print*, "MATCH longitude index, MATCH longitude"
print*, lon_match_ceresfov(fov_lon), lon(lon_match_ceresfov(fov_lon))

print*, "MATCH latitude index, MATCH latitude "
print*, lat_match_ceresfov(fov_lat), lat(lat_match_ceresfov(fov_lat))

print*, "======================================"

! Get the aerosol profiles / aods? matched to the CERES FOV

!print*, "Getting aerosol profiles at CERES FOV location..."

!do t = 1,11                ! aerosol type 
!   do k = mfov_hr,mfov_hr  ! time
!      print*, "========================================================================================="
!      print*, "Hour:",k, "Type:  ",mavp%aero_type11(t), "VarID:  ",iaodvarid(t)
!      print*, "========================================================================================="
!      do p=1,nlay          ! pressure layers  
!         print*, "Lon:", lon(lon_match_ceresfov(fov_lon)),"Lat:", lat(lat_match_ceresfov(fov_lat)), "Pres:",play(p), "Aero:", &
!                         mavp%aero_prof(t,lon_match_ceresfov(fov_lon),lat_match_ceresfov(fov_lat),p,k)
!      end do
!   end do
!end do


print*, "Getting combined aerosol profiles at CERES FOV location..."

do t = 1,7                 ! aerosol type combined
   do k = mfov_hr,mfov_hr  ! time
      print*, "========================================================================================="
      print*, "Hour:",k, "Type:  ",mavp%aero_type7(t)
      print*, "========================================================================================="
      do p = 1,nlay        ! pressure layers
         print*, "Lon:", lon(lon_match_ceresfov(fov_lon)),"Lat:", lat(lat_match_ceresfov(fov_lat)), "Pres:",play(p), "Aero:", &
                         mavp%aero_prof7(t,lon_match_ceresfov(fov_lon),lat_match_ceresfov(fov_lat),p,k)
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









