!=================================================================== 
!
! Name:    read_hourly_match.f90
!
! Purpose: This program reads vertical aerosol profiles (kg/kg) and
!          aerosol optical depth (AOD) information from hourly MATCH
!          netCDF files for 11 aerosol constituents. First raw data
!          data are matched in time and space and extracted at the
!          location of CERES footprints. The profiles are then re-
!          gridded from MATCH to Fu-Liou pressure layers. They are 
!          then converted into AOD profiles as required for input to 
!          the Fu-Liou RT model for CRS flux calculations.
!
! Modules: netcdf
!
! Author:  Ryan Scott, SSAI, ryan.c.scott@nasa.gov
! Updated: February 28, 2020
!
!===================================================================
! To do:
! 
! (0) Figure out how to compile this script on AMI...
! 
! (1) Instead of reading MATCH file from local dir, read the data 
!     from appropriate ASDC_archive dir on AMI. To do this, need to
!     use pcf module, etc. Fred can help w/ these ins and outs... 
!
! (2) Read SSF data directly instead of artificially.
!
! (3) Modularize this code.
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
!character*8 aod_type8(npc+1)                 ! AOD type label, total + 7 combined aerosol types
real :: aod_fov(np)                          ! AOD extracted at CERES FOV LAT/LON/HR
end type aodtype

type (aodtype) aod                           ! MATCH AOD data structure variable

!===============================================
!Parameters required to regrid to Fu-Liou

integer, parameter :: mlev = 30 ! max # of Fu-Liou model levels
real :: dpp(mlev)               ! pressure layer thickness profile (size nv = nv1 - 1)
real :: p1, p2                  ! used in regrid loop - Fu-Liou levels
real :: dp                      ! used in regrid loop - pressure thickness
real :: psfc                    ! used in regrid loop - Fu-Liou surface pressure
real :: pt, pb                  ! used in regrid loop - MATCH p layers - top, bottom
real :: asm1(np), nsm1          ! used in regrid loop 
real :: plast                   ! used in regrid loop - last/previous pressure
real :: sdpp                    ! sum of (Fu-Liou p-layer thicnkess profile, dpp)
real :: tta(np), tta_comb(npc)  ! used to combine profiles
real :: ttt(np)                 ! used to combine profiles
real :: tts(np)                 ! used to combine profiles
real :: amr(np)                 ! X = total AOD / total mass
integer :: icxA(np)             ! pointer to combine constituents
integer :: icxT(np)             ! pointer to combine constituents   
integer :: icxS(np)             ! pointer to combine constituents
real :: ptrop = 200.            ! tropopause pressure, assume 200 mb
integer :: ilev_trop            ! Fu-Liou pressure profile tropospause location index

!=============================================== 
! Fu-Liou model profile (mpro) data structure definition
type mvptype
real :: fovlon            ! CERES FOV lon
real :: fovlat            ! CERES FOV lat
integer :: nv1            ! vertical levels = nv + 1
integer :: mlev           ! max # model levels
real :: pp(mlev)          ! pressure profile levels
real :: dpp(mlev-1)       ! pressure layer thickness profile
real :: mvp(np,mlev)      ! MATCH 11 aerosol prof on Fu-Liou levels at CERES footprint
real :: mvp_comb(npc,mlev)! MATCH 7  aerosol prof on Fu-Liou levels at CERES footprint
real :: aods_all(0:np)    ! AODs 11 + tot at CERES footprint
real :: aods_comb(0:npc)  ! AODs  7 + tot at CERES footprint 
integer :: ityp_comb(npc) ! Fu-Liou aerosol types
end type mvptype

type (mvptype) mpro

!================================================  
! Artificial CERES footprint data
real :: fov_lon       ! CERES FOV longitude
real :: fov_lat       ! CERES FOV latitude
integer :: fov_hr     ! CERES swath hour
integer :: mfov_hr    ! match fov hr = ceres fov hr + 1 

! Prompt user for artificial CERES FOV geolocation & time
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

! pointers for combining profiles
!          1  2  3  4  5  6  7  8  9 10 11          
data icxA /1, 2, 2, 2, 0, 4, 5, 5, 6, 7, 3/ ! full atmosphere
data icxT /0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0/ ! troposphere
data icxS /0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0/ ! stratosphere  

! set parameters of Fu-Liou model profile...
! location of CERES FOV
mpro%fovlon = fov_lon
mpro%fovlat = fov_lat

! specify # of model vertical levels
mpro%nv1 = 30

! specify model pressure levels
mpro%pp(1:mpro%nv1) =                                &  ! Fu-Liou equivalent ( fi%pp(1:fi%nv +1 )
(/0.01, 1., 5., 10., 30. ,50., 70.,100.,125.,150.,   &
175.,200.,250.,300.,350.,400.,450.,500.,550.,600.,   &
650.,700.,750.,800.,850.,900.,950.,975.,995.,1005./)

! Test code on multiple Fu-Liou model profiles...
! second Fu-Liou test pressure level profile...
!mpro%nv1 = 10
!mpro%pp(1:mpro%nv1) = (/1., 10., 50., 100., 200., 400., 600., 700., 900., 1000./)

! sfc pressure from lowest Fu-Liou level
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

! close the file, freeing all resources
call check( nf90_close(ncid) )

print*, "==================================================="
print*, " Successfully read 4D & 5D data...                 "
print*, " Aerosol optical depth (type,nlon,nlat,ntime)      "
print*, " Aerosol profiles      (type,nlon,nlat,nlev,ntime) "
print*, "==================================================="

!print*, "Aerosol profiles at CERES FOV location..."  
!do t = 1,np                ! aerosol type                                                                                                          
!   do k = 1,1  ! time 
!      print*, "========================================================================================="                                            
!      print*, "MATCH Hour:",k, "Type:  ",mavp%aero_type11(t), "VarID:  ",iaodvarid(t)                                                                 
!      print*, "========================================================================================="                                 
!      do p=1,nlay          ! pressure layers                                                                                                      
!         print*, "Lon:", lon(ilon_match_ceresfov(fov_lon)),"Lat:", lat(jlat_match_ceresfov(fov_lat)), "Pres:",play(p), "Aero:", &                        
!                         mavp%aero_array(t,ilon_match_ceresfov(fov_lon),jlat_match_ceresfov(fov_lat),p,k)                                              
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

!print*, "MATCH Latitude Centers"
!print*, lat

!print*, "MATCH Longitude Centers"
!print*, lon

!print*, "MATCH Pressure Levels"
!print*, plev

!print*, "MATCH Sigma Levels"
!print*, hybi

!print*, "MATCH Pressure Layers (Midpoints)"
!print*, play

!print*, "MATCH Sigma Layers (Midpoints)"
!print*, hybm

!print*, "Fu-Liou Pressure Levels"
!print*, mpro%pp(1:mpro%nv1)

!print*, "MATCH sigma layers times psfc"
!print*, psfc*hybm

!print*, "MATCH sigma levels times psfc"
!print*, psfc*hybi

print*, "======================================================="
print*, "Getting MATCH aerosol profiles/AOD at CERES FOV...     "
print*, "======================================================="

! CERES FOV hr is offset from MATCH by 1 hr
mfov_hr = fov_hr + 1

print*, "Artificial CERES Hour   :", fov_hr
print*, "Artificial CERES FOV Lon:", fov_lon
print*, "Artificial CERES FOV Lat:", fov_lat
print*, "----------------------------------------"

print*, "MATCH Longitude index, MATCH Longitude"
print*, ilon_match_ceresfov(fov_lon), lon(ilon_match_ceresfov(fov_lon))

print*, "MATCH Latitude index, MATCH Latitude "
print*, jlat_match_ceresfov(fov_lat), lat(jlat_match_ceresfov(fov_lat))

print*, "======================================================="

! aerosol profiles
do t = 1,np       ! aerosol (t)ypes
   do p = 1,nlay  ! (p)ressure layers         ! call functions that extract MATCH data at the CERES FOV
      mavp%aero_profs(t,p) = mavp%aero_array(t,ilon_match_ceresfov(fov_lon),jlat_match_ceresfov(fov_lat),p,mfov_hr)
   end do
end do

! aerosol optical depth (index 0 refers to total AOD)
do t = 0,np          !indicies... aod%aod_array(1:12) but want mpro%aods_all(0:11)
   mpro%aods_all(t) = aod%aod_array(t+1,ilon_match_ceresfov(fov_lon),jlat_match_ceresfov(fov_lat),mfov_hr)
end do
 
print*, "================================================================="
print*, "Aerosol profiles on native MATCH pressure layers (showing 1-9)..."
print*, "================================================================="

do p = 1,nlay
   write(*,*) (mavp%aero_profs(t,p), t=1,9) ! just showing 9 since it fits on my small screen
end do

print*, "======================================================================"
print*, "Aerosol profiles regridded to Fu-Liou pressure levels (showing 1-9)..."
print*, "======================================================================"

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
            pt = hybi(k)*psfc       ! p-top = MATCH sigma level * Fu-Liou psfc
            pb = hybi(k+1)*psfc     ! p-bottom = next MATCH sigma level * Fu-Liou psfc

            if ( k .eq. nlay)  pb = psfc  ! p-bottom = psfc at end of profile
            dp = 0

            ! compare Fu-Liou & MATCH pressure level positions 
            ! to determine dp at each step
            if( (p1 >= pt .and. p2 <= pb ) .or. &
                (pt >= p1 .and. pb <= p2 ) .or. &
                (p1 >= pt .or. p2 >= pb )  .and. ( p1<=pb .or. p2<= pb) ) then
               
            dp = pb - plast 
                                                                                                                   
            if ( dp == 0 ) cycle  ! exit INLEV do loop
            
            ! if Fu-Liou < MATCH pressure
            if ( p2 < pb ) dp = p2 - plast  

            plast = plast + dp ! increment plast by dp

            ! compute weighted avg
            asm1(1:np) = asm1(1:np) + mavp%aero_profs(1:np,k)*dp ! for all constituents, multiply by dp thickness
            nsm1 = nsm1+dp

            end if

end do INLEV

            ! assign aerosol value to Fu-Liou model layer
            if( nsm1 > 0 ) mpro%mvp(1:np,i) =  asm1(1:np)/ nsm1

end do OUTLEV


! print the resulting profiles on Fu-Liou model layers
do p = 1,mpro%nv1-1
   write(*,*) ( mpro%mvp(t,p), t=1,9)
end do


print*, "================================================================="
print*, "Aerosol profiles converted to AOT (showing 1-9)..."
print*, "================================================================="

! normalize Fu-Liou dpp profile by TOA-surface p difference
dpp(1:mpro%nv1-1) = dpp(1:mpro%nv1-1)/sum(dpp(1:mpro%nv1-1))  ! can comment out with no effect... why?

! do for all aerosol (t)ypes...
do t=1,np

   ! multiply constituent profiles by normalized dpp profile
   mpro%mvp(t, 1:mpro%nv1-1 )= mpro%mvp(t, 1:mpro%nv1-1)* dpp(1:mpro%nv1-1)

   ! 'X' = Total AOD / Total Mass
   amr(t) = mpro%aods_all(t)/sum( mpro%mvp(t, 1:mpro%nv1-1 ) )

   ! multiply profile by 'X' to get AOD profile
   if (amr(t) > 0) mpro%mvp(t, 1:mpro%nv1-1 ) = mpro%mvp(t, 1:mpro%nv1-1 )* amr(t)

end do

! print the resulting AOD profiles
do p = 1,mpro%nv1-1
   write(*,*) (mpro%mvp(t,p), t=1,9) ! just showing 9 since it fits on my small screen                                                                
end do

print*, "=============================================================="
print*, "Comparing 11 AODs (i) sum AOD profiles (ii) MATCH... "
print*, "=============================================================="

! retrieve index of Fu-Liou tropopause (assuming p = 200mb)
do i =1,mpro%nv1-1     
   if ( mpro%pp(i) >= ptrop ) exit
      ilev_trop = i
end do

tta_comb = 0
tta = 0          ! total thickness - atmosphere
ttt = 0          ! total thickness - troposphere
tts = 0          ! total thickness - stratosphere

mpro%mvp_comb  = 0
mpro%aods_comb = 0
mpro%aods_comb(0) = mpro%aods_all(0)  ! total AOD


do t = 1,np

   ! AOD for each atmospheric layer
   tta(t) = sum( mpro%mvp(t, 1:mpro%nv1-1 )  )          ! sum over full atmospheric profile
   tts(t) = sum( mpro%mvp(t, 1:ilev_trop )  )           ! sum over stratospheric portion                            
   ttt(t) = sum( mpro%mvp(t, ilev_trop+1:mpro%nv1-1 ) ) ! sum over tropospheric portion                                                                 
   
   ! these values should be equal
   print*, tta(t), mpro%aods_all(t)

   if ( icxA(t) > 0 .and. tta(t) > 0 ) then
      tta_comb( icxA(t) ) = tta_comb( icxA(t) ) + tta(t)
      mpro%mvp_comb( icxA(t), 1:mpro%nv1-1 ) =  &
      mpro%mvp_comb( icxA(t), 1:mpro%nv1-1 ) + mpro%mvp(t, 1:mpro%nv1-1  )
      mpro%aods_comb ( icxA(t)) = mpro%aods_comb ( icxA(t)) + mpro%aods_all(t) * tta(t)/tta(t)   ! total profile
   endif

   if ( icxS(t) > 0 .and. tts(t) > 0) then
      tta_comb( icxS(t) ) = tta_comb( icxS(t) ) + tts(t)
      mpro%mvp_comb( icxS(t), 1:ilev_trop  ) =  &
      mpro%mvp_comb( icxS(t), 1:ilev_trop  ) + mpro%mvp(t, 1:ilev_trop )
      mpro%aods_comb ( icxS(t)) = mpro%aods_comb ( icxS(t)) + mpro%aods_all(t) * tts(t)/tta(t)   ! stratosphere
   endif

   if ( icxT(t) > 0 .and. ttt(t) > 0) then
      tta_comb( icxT(t) ) = tta_comb( icxT(t) ) + ttt(t)
      mpro%mvp_comb( icxT(t), ilev_trop+1:mpro%nv1-1  ) =  &
      mpro%mvp_comb( icxT(t), ilev_trop+1:mpro%nv1-1  ) + mpro%mvp(t, ilev_trop+1:mpro%nv1-1  )
      mpro%aods_comb ( icxT(t)) = mpro%aods_comb ( icxT(t)) + mpro%aods_all(t) * ttt(t)/tta(t)   ! troposphere
   endif

end do


! renormalize profiles by total AOD and convert to percentages
! this allows division by zero, which yields NaN values...
!do k = 1, mpro%nv1-1 
!   mpro%mvp(1:np,k)       = 100*mpro%mvp     (1:np,k ) / tta(1:np)
!   mpro%mvp_comb(1:npc,k) = 100*mpro%mvp_comb(1:npc,k) / tta_comb(1:npc)
!end do

! renormalize profiles by total AOD and convert to percentages
! ...but don't divide by zero...
do t = 1,np
   do k = 1,mpro%nv1-1
      if (tta(t)>0) mpro%mvp(t,k) = 100*mpro%mvp(t,k) / tta(t)
   end do
end do

do t = 1,npc
   do k = 1,mpro%nv1-1
      if (tta_comb(t)>0) mpro%mvp_comb(t,k) = 100*mpro%mvp_comb(t,k) / tta_comb(t)
   end do
end do


! get combined AODs
mpro%aods_comb (1:npc) = mpro%aods_comb (1:npc)
mpro%aods_comb (0)     = sum(mpro%aods_comb(1:npc))

mpro%ityp_comb(1:npc) =  (/ 24, 25, 18, 1, 11, 10, 9 /) ! Fu-Liou aerosol types

print*, "======================================================================================================================"
print*, "      DustSm           DustLg         OPAC SUSO           SSLT          OPAC SOOT        OPAC WASO        OPAC INSO"
print*, "======================================================================================================================"

do p = 1,mpro%nv1-1
   write(*,*) ( mpro%mvp_comb(t,p), t=1,npc)
end do

print*, "================================================="
print*, "Checking to ensure each profile sums to 100 % ..."
print*, "================================================="
do t = 1,npc
   write(*,*) sum(mpro%mvp_comb(t,1:mpro%nv1-1))
enddo



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









