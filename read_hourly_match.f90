!=================================================================== 
!
! Name:    read_hourly_match.f90
!
! Purpose: This program reads aerosol vertical profiles (kg/kg) and
!          AOD information from hourly MATCH netCDF files for 11 
!          different aerosol constituents.
!
!          Artificial CERES FOV lat, lon & swath hour are used
!          (input by user) to extract necessary aerosol data for 
!          Fu-Liou computation of CERES footprint-level sfc fluxes.
!          This will soon be replaced by reading in SSF footprints.
!
!          Aerosol profiles are combined for CERES Edition 4 
!          according to OPAC standards. Aerosol profiles need to be 
!          regridded from the MATCH to the Fu-Liou vertical grid, 
!          which may change...
!         
!          MATCH profiles are reported for pressure layers, not levels.
!
! Modules: netcdf
!
! Author:  Ryan Scott, SSAI, ryan.c.scott@nasa.gov
! Updated: February 25, 2020
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
! (2) IN PROGRESS: Match the aerosol data to specific CERES FOVs...
!                  This includes the HOUR time step and LAT/LON
!                      > Read SSF data directly...
!
! (3) IN PROGRESS: Construct Fu-Liou model vertical profile...
!                  Combine constituents and get into Fu-Liou form
!                      > Convert profile into AOT profile
!                      > Re-grid onto Fu-Liou model profile 
!
!===================================================================
program read_hourly_match 
 
! use netcdf module - see README for compilation instructions
use netcdf

! no implicit variables
implicit none

! name of MATCH file to be read
character (len = *), parameter :: FILE_NAME = "CER_MATCH-hourly_Terra-Aqua-MODIS_Edition4_402402.20190101.nc"

! dimensions of MATCH data
integer, parameter :: nlon = 192   ! number of longitude boxes
integer, parameter :: nlat = 94    ! number of latitude boxes
integer, parameter :: ntime = 24   ! number of hourly time steps
integer, parameter :: nlay = 28    ! number of pressure layers
integer, parameter :: nlev = nlay+1! number of pressure levels
integer, parameter :: np = 11      ! number of particle constituents
integer, parameter :: npc = 7      ! number of particle constituents combined
real :: lon(nlon), lat(nlat)       ! MATCH grid box lon, lat
real :: play(nlay), plev(nlev)     ! MATCH pressure layers, levels
real :: hybi(nlev), hybm(nlay)     ! MATCH sigma coordinates for levels & layer midpoints

! IDs for file & select variables
integer :: iplayvarid, iplevvarid  ! integer pressure lay/lev variable IDs
integer :: ihybivarid, ihybmvarid  ! integer sigma coord variable IDs
integer :: ilatvarid, ilonvarid    ! integer lat/lon variable IDs
integer :: iprofvarid(1:np)        ! integer aerosol profile variable ID, 11 aerosol types                                 
integer :: iaodvarid(1:np+1)       ! integer aerosol OD variable ID

integer :: ncid                    ! netCDF (nc) file ID
integer :: i, j, p, k, t           ! array loop indices - lon, lat, lay, hr, type

!===========================================  
! MATCH aerosol vertical profile (mavp) data structure definition
type mavptype
real :: aero_array(np,nlon,nlat,nlay,ntime)  ! aerosol array, types 1-11
real :: aero_array7(npc,nlon,nlat,nlay,ntime)! aerosol array, types 1-11 combined into 7 categories
character*7 aero_type11(np)                  ! aerosol type label, 11 constituents
character*10 aero_type7(npc)                 ! aerosol type label, 7 constituents obtained by combining categories
real :: aero_profs(np,nlay)                  ! aerosol profiles extracted at CERES FOV LAT/LON/HR
real :: aero_profs7(npc,nlay)                ! aerosol profiles extracted at CERES FOV LAT/LON/HR, combined into 7 categories
end type mavptype

type (mavptype) mavp                         ! MATCH aerosol profile data structure variable

!=========================================== 
! MATCH aerosol optical depth (AOD) data structure definition
type aodtype
real :: aod_array(np+1,nlon,nlat,ntime)      ! AOD array, total + 11 aerosol types
real :: aod_array8(npc+1,nlon,nlat,ntime)    ! AOD array, total + 7 combined aerosol types
character*8 aod_type12(np+1)                 ! AOD type label, total + 11 aerosol types
character*8 aod_type8(npc+1)                 ! AOD type label, total + 7 combined aerosol types
real :: aod_fov(np)                          ! AOD extracted at CERES FOV LAT/LON/HR
end type aodtype

type (aodtype) aod                           ! MATCH AOD data structure variable

!===============================================
!Parameters required to regrid to Fu-Liou

integer, parameter :: mlev = 30 ! max # of Fu-Liou model levels
real :: dpp(mlev)               ! pressure layer thickness profile (size nv = nv1 - 1)
real :: p1, p2
real :: dp          
real :: psfc, px, pt, pb
real :: asm1(np), nsm1
real :: plast
real :: ttt_comb
real :: sdpp       

! Fu-Liou model profile (mpro) data structure definition

type mvptype
real :: fovlon            ! CERES FOV lon
real :: fovlat            ! CERES FOV lat
integer :: nv1            ! vertical levels = nv + 1
integer :: mlev           ! max # model levels
real :: pp(mlev)          ! pressure profile levels
real :: dpp(mlev-1)       ! pressure layer thickness profile
real :: mvp(np,mlev)      !
real :: mvp_comb(npc,mlev)! 
end type mvptype

type (mvptype) mpro

!================================================  
! Artificial CERES footprint data (temporary)

real :: fov_lon       ! CERES FOV longitude
real :: fov_lat       ! CERES FOV latitude
integer :: fov_hr     ! CERES swath hour
integer ::  mfov_hr   ! match fov hr = ceres fov hr + 1 

! Prompt user for artificial CERES geolocation & time
! This will need to be replaced eventually...

write(*,*) "Enter artificial CERES FOV longitude:"
read(*,*) fov_lon
write(*,*) "Enter artificial CERES FOV latitude:"
read(*,*) fov_lat
write(*,*) "Enter artificial CERES hour:"
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


! Set parameters of Fu-Liou model profile...
! location of CERES FOV
mpro%fovlon = fov_lon
mpro%fovlat = fov_lat

! specify model # vertical levels
mpro%nv1 = 30
!mpro%nv1 = 10

! specify model pressure levels                                                                                                                      
mpro%pp(1:mpro%nv1) =                                &  ! Fu-Liou equivalent ( fi%pp(1:fi%nv +1 )
(/0.01, 1., 5., 10., 30. ,50., 70.,100.,125.,150.,   &
175.,200.,250.,300.,350.,400.,450.,500.,550.,600.,   &
650.,700.,750.,800.,850.,900.,950.,975.,995.,1005./)

!mpro%pp(1:mpro%nv1) = (/1., 10., 50., 100., 200., 400., 600., 700., 900., 1000./)

psfc = mpro%pp(mpro%nv1)

!==============================================

! open the netcdf file in read mode
call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )

! get the id for each variable based on its name
call check( nf90_inq_varid(ncid, "lon", ilonvarid) )   ! lon var id
call check( nf90_inq_varid(ncid, "lat", ilatvarid) )   ! lat var id
call check( nf90_inq_varid(ncid, "ilev", iplevvarid) ) ! lev var id - nc files poorly named...
call check( nf90_inq_varid(ncid, "lev" , iplayvarid) ) ! lay var id
call check( nf90_inq_varid(ncid, "hybi", ihybivarid) ) ! hybi var id - lev
call check( nf90_inq_varid(ncid, "hybm", ihybmvarid) ) ! hybm var id - lay

! print nc file and variable ids
print*, "===================================="
print*, "  NetCDF file and variable IDs...   "
print*, "===================================="

! read the data
call check( nf90_get_var(ncid,ilonvarid, lon) )
call check( nf90_get_var(ncid,ilatvarid, lat) )
call check( nf90_get_var(ncid,iplevvarid, plev) )
call check( nf90_get_var(ncid,iplayvarid, play) )
call check( nf90_get_var(ncid,ihybivarid, hybi) )
call check( nf90_get_var(ncid,ihybmvarid, hybm) )

!********************************************************************************************************
! dimensions are reversed in the netCDF file - this code is old but left here as an FYI
!call check( nf90_get_var(ncid, varid1, data3d, start=(/1,1,1/), count=(/nlon,nlat,ntime/) ) )
!call check( nf90_get_var(ncid, varid2, data4d, start=(/1,1,1,1/), count=(/nlon,nlat,nlev,ntime/) ) )
!********************************************************************************************************

! read aerosol arrays
do i = 1,11
  ! get variable ids for each aerosol-type profile
  call check( nf90_inq_varid(ncid, mavp%aero_type11(i), iprofvarid(i) ) )
  print*, "Getting...", mavp%aero_type11(i)," variable id...", iprofvarid(i)
  ! get vertical aerosol profiles for each type and store in single 5-d array
  call check( nf90_get_var(ncid, iprofvarid(i), mavp%aero_array(i,:,:,:,:), start=(/1,1,1,1/), count=(/nlon,nlat,nlay,ntime/) ) )
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
print*, " Successfully read 4D & 5D data...  "
print*, " AerosOD(type,nlon,nlat,ntime)      "
print*, " Profile(type,nlon,nlat,nlev,ntime) "



!print*, "Getting aerosol profiles at CERES FOV location..."  
!do t = 1,11                ! aerosol type                                                                                                          
!   do k = mfov_hr,mfov_hr  ! time 
!      print*, "========================================================================================="                                            
!      print*, "MATCH Hour:",k, "Type:  ",mavp%aero_type11(t), "VarID:  ",iaodvarid(t)                                                                 
!      print*, "========================================================================================="                                 
!      do p=1,nlay          ! pressure layers                                                                                                      
!         print*, "Lon:", lon(lon_match_ceresfov(fov_lon)),"Lat:", lat(lat_match_ceresfov(fov_lat)), "Pres:",play(p), "Aero:", &                        
!                         mavp%aero_prof(t,lon_match_ceresfov(fov_lon),lat_match_ceresfov(fov_lat),p,k)                                              
!      end do                                                                                                                                   
!   end do                                                                                                                                           
!end do

!print*, "Printing combined aerosol profiles at CERES FOV location..."                                                                                 
!do t = 1,7                 ! aerosol type combined                                                                                                    
!   do k = mfov_hr,mfov_hr  ! time                                                                                                                    
!      print*, "========================================================================================="                                             
!      print*, "MATCH Hour:",k, "Type:  ",mavp%aero_type7(t)                                                                                           
!      print*, "========================================================================================="                                            
!      do p = 1,nlay        ! pressure layers                                                                                                          
!         print*, "Lon:", lon(lon_match_ceresfov(fov_lon)),"Lat:", lat(lat_match_ceresfov(fov_lat)), "Pres:",play(p), "Aero:", &                  
!                         mavp%aero_array7(t,lon_match_ceresfov(fov_lon),lat_match_ceresfov(fov_lat),p,k)                                        
!      end do                                                                                                                                          
!   end do                                                                                                                                      
!end do                                                                                                                                         


!print*, "===================================="
!print*, "Printing aerosol vertical profiles..."
!print*, "===================================="

!do t = 1,11            ! aerosol type
!  do i = 31,31         ! longitude
!     do j = 31,31      ! latitude
!        do k = 1,1     ! time
!          print*, "========================================================================================="
!          print*, "Hour:",k, "Type:  ",mavp%aero_type11(t), "VarID:  ",iaodvarid(t)
!          print*, "========================================================================================="
!          do p=1,nlay  ! pressure layers
!            print*, "Lon:", lon(i),"Lat:", lat(j), "Pres:",play(p), "Aero:", mavp%aero_array(t,i,j,p,k)
!          end do
!        end do
!     end do
!  end do
!end do

!print*, "======================================="
!print*, "Combining Ed4 aerosol arrays into 7 types...      "
!print*, "======================================="

! 7 species
!mavp%aero_array7(1,:,:,:,:) = mavp%aero_array(1,:,:,:,:)                                  ! 1. DustSm = DSTQ01
!mavp%aero_array7(2,:,:,:,:) = sum(mavp%aero_array(2:4,:,:,:,:),1)                         ! 2. DustLg = DSTQ02+DSTQ03+DSTO4
!mavp%aero_array7(4,:,:,:,:) = mavp%aero_array(6,:,:,:,:)                                  ! 4. SSLT
!mavp%aero_array7(5,:,:,:,:) = sum(mavp%aero_array(7:8,:,:,:,:),1)                         ! 5. OPAC SOOT = BCPHI+BCPHO
!mavp%aero_array7(7,:,:,:,:) = mavp%aero_array(10,:,:,:,:)                                 ! 7. OPAC INSO

! assumes 200 mb tropopause
!do p=1,nlay
!  ! stratosphere
!  if (p <= 10) then
!    mavp%aero_array7(3,:,:,p,:) = mavp%aero_array(5,:,:,p,:) + mavp%aero_array(11,:,:,p,:) ! 3. OPAC SUSO = SO4(strato)+VOLC
!    mavp%aero_array7(6,:,:,p,:) = mavp%aero_array(9,:,:,p,:)                              ! 6. OPAC WASO = OCPHI
!  ! troposphere
!  else if (p > 10) then
!    mavp%aero_array7(3,:,:,p,:) = mavp%aero_array(11,:,:,p,:)                             ! 3. OPAC SUSO = VOLC only - 0 in troposphere
!    mavp%aero_array7(6,:,:,p,:) = mavp%aero_array(9,:,:,p,:) + mavp%aero_array(5,:,:,p,:)  ! 6. OPAC WASO = OCPHI + SO4(tropo)
!  end if
!end do

!print*, "Done combining aerosol arrays into new OPAC types..."

!print*, "======================================="
!print*, "Combining AOD arrays into 7 types + total AOD"
!print*, "======================================="

! total + 7 species
!aod%aod_array8(1,:,:,:) = sum(aod%aod_array(2:12,:,:,:),1) ! Total - should equal AEROOD
!aod%aod_array8(2,:,:,:) = aod%aod_array(2,:,:,:)           ! DustSm
!aod%aod_array8(3,:,:,:) = sum(aod%aod_array(3:5,:,:,:),1)  ! DustLg
!aod%aod_array8(4,:,:,:) = aod%aod_array(6,:,:,:)           ! SO4
!aod%aod_array8(5,:,:,:) = aod%aod_array(7,:,:,:)           ! SSLT - Sea Salt
!aod%aod_array8(6,:,:,:) = sum(aod%aod_array(8:9,:,:,:),1)  ! SOOT - Black Carbon
!aod%aod_array8(7,:,:,:) = aod%aod_array(10,:,:,:)          ! Solub
!aod%aod_array8(8,:,:,:) = aod%aod_array(11,:,:,:)          ! Insolub

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

print*, "MATCH Sigma Levels"
print*, hybi

print*, "MATCH Pressure Layers (Midpoints)"
print*, play

print*, "MATCH Sigma Layers (Midpoints)"
print*, hybm

print*, "Fu-Liou Pressure Levels"
print*, mpro%pp(1:mpro%nv1)

print*, "======================================"
print*, "Matching MATCH data to CERES FOV...   "
print*, "======================================"

! CERES FOV hr is offset from MATCH by 1 hr
mfov_hr = fov_hr + 1

print*, "Artificial CERES Hour   :", fov_hr
print*, "Artificial CERES FOV lon:", fov_lon
print*, "Artificial CERES FOV lon:", fov_lat

print*, "MATCH longitude index, MATCH longitude"
print*, ilon_match_ceresfov(fov_lon), lon(ilon_match_ceresfov(fov_lon))

print*, "MATCH latitude index, MATCH latitude "
print*, jlat_match_ceresfov(fov_lat), lat(jlat_match_ceresfov(fov_lat))

print*, "======================================"
print*, "Getting aerosol profiles at CERES FOV "
print*, "======================================"

do t = 1,np       ! aerosol types
   do p = 1,nlay  ! pressure layers          ! call functions to match aerosol profiles to CERES FOV
      mavp%aero_profs(t,p) = mavp%aero_array(t,ilon_match_ceresfov(fov_lon),jlat_match_ceresfov(fov_lat),p,mfov_hr)
   end do
end do

print*, "======================================="
print*, "Combining aerosol profiles at FOV...   "
print*, "======================================="

! 7 species                                                                                                                                            
mavp%aero_profs7(1,:) = mavp%aero_profs(1,:)                                  ! 1. DustSm = DSTQ01                                       
mavp%aero_profs7(2,:) = sum(mavp%aero_profs(2:4,:),1)                         ! 2. DustLg = DSTQ02+DSTQ03+DSTO4                           
mavp%aero_profs7(4,:) = mavp%aero_profs(6,:)                                  ! 4. SSLT                                                   
mavp%aero_profs7(5,:) = sum(mavp%aero_profs(7:8,:),1)                         ! 5. OPAC SOOT = BCPHI+BCPHO                              
mavp%aero_profs7(7,:) = mavp%aero_profs(10,:)                                 ! 7. OPAC INSO

! assumes 200 mb tropopause                                                                                                                              
do p=1,nlay                                                                                                                                             
  ! stratosphere                                                                                                                                        
  if (p <= 10) then                                                                                                                                     
    mavp%aero_profs7(3,p) = mavp%aero_profs(5,p) + mavp%aero_profs(11,p) ! 3. OPAC SUSO = SO4(strato)+VOLC                            
    mavp%aero_profs7(6,p) = mavp%aero_profs(9,p)                         ! 6. OPAC WASO = OCPHI                                        
  ! troposphere                                                                                                                                         
  else if (p > 10) then                                                                                                                                 
    mavp%aero_profs7(3,p) = mavp%aero_profs(11,p)                        ! 3. OPAC SUSO = VOLC only - 0 in troposphere                 
    mavp%aero_profs7(6,p) = mavp%aero_profs(9,p) + mavp%aero_profs(5,p)  ! 6. OPAC WASO = OCPHI + SO4(tropo)                          
  end if                                                                                                                                                
end do 

print*, "======================================================================================================================"
print*, "      DustSm           DustLg         OPAC SUSO           SSLT          OPAC SOOT        OPAC WASO        OPAC INSO"
print*, "======================================================================================================================"
do p = 1,nlay
   write(*,*) ( mavp%aero_profs7(t,p), t=1,7 )
end do
print*, "======================================================================================================================"

!===============================================================
! The following code has been modified from MATCH_ed4_hourly_SYNI.f90...

mpro%mvp = 0.0
plast = 0.0                         ! mpro%pp(1) <- TOA

OUTLEV : do i = 1,mpro%nv1-1        ! loop over Fu-Liou model layers starting at TOA
            p1 = mpro%pp(i)         ! p1 < p2 
            p2 = mpro%pp(i+1)       ! p2
            dpp(i) = p2-p1          ! calculate pressure layer thickness profile

            asm1=0    
            nsm1=0   

            hybi(1) = 0.00          ! TOA at zero pressure

INLEV  : do k = 1,nlay              ! loop over MATCH layers starting at TOA
            px = hybm(k)*psfc       ! multiply MATCH sigma layer midpoints by Fu-Liou psfc
            pt = hybi(k)*psfc       ! multiply MATCH sigma level by Fu-Liou psfc
            pb = hybi(k+1)*psfc     ! set p bottom equal to Fu-Liou surface pressure

            if ( k .eq. nlay)  pb = psfc  
            dp = 0

            
            if( (p1 >= pt .and. p2 <= pb ) .or. &
                (pt >= p1 .and. pb <= p2 ) .or. &
                (p1 >= pt .or. p2 >= pb )  .and. ( p1<=pb .or. p2<= pb) ) then

            dp = pb - plast 
             
            !print*, pt, pb, dp, plast 
                                                                                                                   
            if ( dp == 0 ) cycle  ! exit inner do loop
            
            if ( p2 < pb ) dp = p2 - plast  

            plast = plast + dp

            ! change aero_profs7 to aero_profs and npc to np to do raw constituents...
            asm1(1:npc) = asm1(1:npc) + mavp%aero_profs7(1:npc,k)*dp
            nsm1 = nsm1+dp

            !print*, nsm1

            end if

end do INLEV
            
            if( nsm1 > 0 ) mpro%mvp(i,1:npc) =  asm1(1:npc)/ nsm1

            ! print out the resulting profiles on Fu-Liou grid
            print*, mpro%mvp(i,1:npc)

enddo OUTLEV










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
function ilon_match_ceresfov(fovlon) result(ilon)
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
end function ilon_match_ceresfov

!********************************************************** 
! FUNCTION                                                                        
! Calculate CERES FOV & MATCH lat distance, match by index
!**********************************************************                                                                           
function jlat_match_ceresfov(fovlat) result(ilat)
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
end function jlat_match_ceresfov


!********************************************************** 
end program read_hourly_match









