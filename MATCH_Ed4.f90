 module MATCH_Ed4
!=================================================================== 
! 
!  Name:    MATCH_Ed4.f90
!
!  Purpose: This module contains several subroutines and functions
!           to read MATCH aerosol vertical profiles (kg/kg) & aerosol
!           optical depth information from daily, hourly resolution
!           netCDF files for 11 aerosol constituents; to match the
!           data in time & space with CERES footprints; to regrid
!           the profiles from MATCH pressure layers to Fu-Liou
!           pressure layers; and to combine profiles and AODs into
!           7 OPAC aerosol types.
!
!  Subroutines: ingest_match_aerosol_arrays <- reads MATCH data from ASDC netCDF file
!               match_with_ceres_fov        <- extracts MATCH data at the CERES FOV
!               regrid_from_match_to_fu     <- places MATCH aerosols on Fu-Liou model profile
!               convert_profiles_to_aot     <- converts aerosol mixing ratio profile to AOT profile
!               combine_aerosol_profiles    <- returns AOT profiles and total AOT for input to Fu-Liou
!  
!  Functions:   lat_ind_match_ceresfov      <- returns latitude index of MATCH gridbox
!               lon_ind_match_ceresfov      <- returns longitude index of MATCH gridbox
!
!  Author:  Ryan Scott, SSAI
!           ryan.c.scott@nasa.gov
!
!  Updated: March 6, 2020
!
!===================================================================
   
! no implicit variables
implicit none

! dimensions of MATCH data
integer, parameter :: nlon = 192   ! number of longitude boxes
integer, parameter :: nlat = 94    ! number of latitude boxes
integer, parameter :: ntime = 24   ! number of hourly time steps
integer, parameter :: nlay = 28    ! number of pressure layers
integer, parameter :: nlev = nlay+1! number of pressure levels
integer, parameter :: np = 11      ! number of particlate constituents
integer, parameter :: npc = 7      ! number of particlate constituents - combined
real :: lon(nlon)                  ! MATCH grid box lon  
real :: lat(nlat)                  ! MATCH grid box lat   
real :: play(nlay), plev(nlev)     ! MATCH pressure layers, levels
real :: hybi(nlev), hybm(nlay)     ! MATCH sigma coordinates for levels & layer midpoints

! IDs for file & select variables
integer :: ncid                    ! netCDF (nc) file ID
integer :: iplayvarid, iplevvarid  ! integer pressure layer, level variable IDs
integer :: ihybivarid, ihybmvarid  ! integer sigma coordinate variable IDs
integer :: ilatvarid, ilonvarid    ! integer latitude, longitude variable IDs
integer :: iprofvarid(1:np)        ! integer aerosol profile variable ID, 11 aerosol types                                 
integer :: iaodvarid(1:np+1)       ! integer aerosol OD variable ID

! For inquiring dimension name / lengths from file 
!integer :: ilatdimid, ilondimid, iplaydimid, iplevdimid, itimedimid
!integer :: latlen, lonlen, playlen, plevlen, timelen
!character*10 :: latname
!character*10 :: lonname
!character*10 :: playname
!character*10 :: plevname
!character*10 :: timename

integer :: i, j, p, k, t           ! array loop indices - lon, lat, lay, hr, type

!===========================================  
! MATCH aerosol vertical profile (mavp) data structure definition
type mavptype
   real :: aero_array(np,nlon,nlat,nlay,ntime)  ! aerosol array, types 1-11
   character*7 aero_type11(np)                  ! aerosol type label, 11 constituents
   character*10 aero_type7(npc)                 ! aerosol type label, 7 constituents obtained by combining categories
   real :: aero_profs(np,nlay)                  ! aerosol profiles extracted at CERES FOV LAT/LON/HR
end type mavptype

type (mavptype) mavp                            ! MATCH aerosol vertical profile data structure variable

!=========================================== 
! MATCH aerosol optical depth (AOD) data structure definition
type aodtype
   real :: aod_array(np+1,nlon,nlat,ntime)      ! AOD array, total + 11 aerosol types
   real :: aod_array8(npc+1,nlon,nlat,ntime)    ! AOD array, total + 7 combined aerosol types
   character*8 aod_type12(np+1)                 ! AOD type label, total + 11 aerosol types
   real :: aod_fov(np)                          ! AOD extracted at CERES FOV LAT/LON/HR
end type aodtype

type (aodtype) aod                              ! MATCH AOD data structure variable ...of course AOD=AOT

!===============================================
!Parameters required to regrid to Fu-Liou

integer, parameter :: mlev = 200! max # of Fu-Liou model levels
real :: dpp(mlev)               ! pressure layer thickness profile (size nv = nv1 - 1)
real :: p1, p2                  ! used in regrid loop - Fu-Liou levels
real :: dp                      ! used in regrid loop - pressure thickness
real :: psfc                    ! used in regrid loop - Fu-Liou surface pressure
real :: pt, pb                  ! used in regrid loop - MATCH p layers - top, bottom
real :: asm1(np), nsm1          ! used in regrid loop to calculate Fu-Liou value
real :: plast                   ! used in regrid loop - last/previous pressure
real :: tta(np), tta_comb(npc)  ! total thickness - full atmosphere: individual species, combined
real :: ttt(np)                 ! total thickness - troposphere
real :: tts(np)                 ! total thickness - stratosphere
real :: amr(np)                 ! X = total AOD / total mass
integer :: icxA(np)             ! pointer to combine constituents - full atmos
integer :: icxT(np)             ! pointer to combine constituents - troposphere
integer :: icxS(np)             ! pointer to combine constituents - stratosphere
real :: ptrop = 200.            ! tropopause pressure, assume = 200 mb
integer :: ilev_trop            ! index of Fu-Liou pressure profile tropospause location

!=============================================== 
! Fu-Liou model profile (mpro) data structure definition
type mvptype
   real :: rlon                 ! CERES FOV lon
   real :: rlat                 ! CERES FOV lat
   real :: fovhr                ! CERES FOV hour
   integer :: nv1               ! vertical levels = nv + 1
   integer :: mlev              ! max # model levels
   real :: pp(mlev)             ! pressure profile levels
   real :: dpp(mlev-1)          ! pressure layer thickness profile
   real :: mvp(np,mlev)         ! MATCH 11 aerosol prof on Fu-Liou levels at CERES footprint
   real :: mvp_comb(npc,mlev)   ! MATCH 7  aerosol prof on Fu-Liou levels at CERES footprint
   real :: mvp_comb1(mlev,npc)  ! switch dimensions at the end - necessary!
   real :: aods_all(0:np)       ! AODs 11 + tot at CERES footprint
   real :: aods_comb(0:npc)     ! AODs  7 + tot at CERES footprint 
   integer :: ityp_comb(npc)    ! Fu-Liou aerosol types
end type mvptype

type (mvptype) mpro             ! model profile variable


integer :: match_hr             ! match hr = ceres fov hr + 1

!================================================

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
'DustSm',     &    ! DSTQ01               (small dust)
'DustLg',     &    ! DSTQ02+DSTQ03+DSTQ04 (large dust)
'OPAC SUSO',  &    ! VOLC + SO4(strato)   (OPAC SUSO)
'SSLT',       &    ! d'Almedia maritime   (sea salt)
'OPAC SOOT',  &    ! BCPHI+BCPHO          (OPAC SOOT)
'OPAC WASO',  &    ! OCPHI + SO4(tropo)   (OPAC WASO)
'OPAC INSO' /      ! OCPHO                (OPAC INSO)

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

! pointers for combining aerosol data
! types    1  2  3  4  5  6  7  8  9 10 11
data icxA /1, 2, 2, 2, 0, 4, 5, 5, 6, 7, 3/ ! full atmosphere
data icxT /0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0/ ! troposphere
data icxS /0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0/ ! stratosphere

! MATCH latitude centers
!real, parameter :: lat(nlat) =                                   &
!(/-88.54195,-86.65317,-84.75323, -82.85078,-80.94736,-79.0434,   &
!-77.13935, -75.23505, -73.33066, -71.42619,-69.52167,-67.617,    &
!-65.71251, -63.8079, -61.90326, -59.99861,-58.09395,-56.1892,    &
!-54.2846, -52.37991, -50.47522, -48.57052, -46.66582,-44.7611,   &
!-42.8564, -40.95169, -39.04697, -37.14225, -35.23753,-33.3328,   &
!-31.42808, -29.52336, -27.61863, -25.7139, -23.80917,-21.9044,   &
!-19.99971,-18.09498,-16.19024,-14.28551, -12.38078,-10.4760,     &
!-8.571308, -6.666573, -4.761838, -2.857103, -0.952368,0.952368,  &
!2.857103, 4.761838, 6.666573,8.571308,10.47604,12.38078,14.285,  &
!16.19024, 18.09498, 19.99971,21.90444,23.80917,25.7139, 27.618,  &
!29.52336, 31.42808, 33.33281,35.23753,37.14225,39.04697,40.951,  &
!42.8564, 44.76111, 46.66582,48.57052,50.47522,52.37991, 54.284,  &
!56.18928, 58.09395, 59.99861,61.90326,63.8079,65.71251, 67.617,  &
!69.52167, 71.42619, 73.33066,75.23505,77.13935,79.04349,80.947,  &
!82.85078, 84.75323, 86.65317,88.54195 /)

! MATCH longitude centers
!real, parameter :: lon(nlon) =                                   &
!(/0.,1.875,3.75,5.625,7.5,9.375,11.25,13.125,15.,16.875,18.75,   &
!20.625,22.5,24.375,26.25,28.125,30.,31.875,33.75,35.625,37.5,    &
!39.375,41.25,43.125,45.,46.875,48.75,50.625,52.5,54.375,56.25,   &
!58.125,60.,61.875,63.75,65.625,67.5,69.375,71.25,73.125,75.,     &
!76.875,78.75,80.625,82.5,84.375,86.25,88.125,90.,91.875,93.75,   &
!95.625,97.5,99.375,101.25,103.125,105.,106.875,108.75,110.625,   &
!112.5,114.375,116.25,118.125,120.,121.875,123.75,125.625,127.5,  &
!129.375,131.25,133.125,135.,136.875,138.75,140.625,142.5,144.375,&
!146.25,148.125,150.,151.875,153.75,155.625,157.5,159.375,161.25, &
!163.125,165.,166.875,168.75,170.625,172.5,174.375,176.25,178.125,&
!180.,181.875,183.75,185.625,187.5,189.375,191.25,193.125,195.,   &
!196.875,198.75,200.625,202.5,204.375,206.25,208.125,210.,211.875,&
!213.75,215.625,217.5,219.375,221.25,223.125,225.,226.875,228.75, &
!230.625,232.5,234.375,236.25,238.125,240.,241.875,243.75,245.625,&
!247.5,249.375,251.25,253.125,255.,256.875,258.75,260.625,262.5,  &
!264.375,266.25,268.125,270.,271.875,273.75,275.625,277.5,279.375,&
!281.25,283.125,285.,286.875,288.75,290.625,292.5,294.375,296.25, &
!298.125,300.,301.875,303.75,305.625,307.5,309.375,311.25,313.125,&
!315.,316.875,318.75,320.625,322.5,324.375,326.25,328.125,330.,   &
!331.875,333.75,335.625,337.5,339.375,341.25,343.125,345.,346.875,&
!348.75,350.625,352.5,354.375,356.25,358.125/)

! BELOW : DEVELOPMENT PURPOSES ONLY

! set parameters of Fu-Liou model profile...
! location of CERES FOV
!mpro%fovlon = fov_lon
!mpro%fovlat = fov_lat

! specify # of model vertical levels
!mpro%nv1 = 30

! specify model pressure profile levels
!mpro%pp(1:mpro%nv1) =                                &  ! Fu-Liou equivalent ( fi%pp(1:fi%nv +1 )
!(/0.01, 1., 5., 10., 30. ,50., 70.,100.,125.,150.,   &
!175.,200.,250.,300.,350.,400.,450.,500.,550.,600.,   &
!650.,700.,750.,800.,850.,900.,950.,975.,995.,1005./)

! Another Fu-Liou pressure level profile for testing...
! Uncomment the next two lines to overwrite the above definitions

!mpro%nv1 = 10
!mpro%pp(1:mpro%nv1) = (/1., 10., 50., 100., 200., 400., 600., 700., 900., 1000./)

! sfc pressure from lowest Fu-Liou level
!psfc = mpro%pp(mpro%nv1)

!==========================================================
!==========================================================

contains

!**********************************************************
!**********************************************************
! SUBROUTINE:
! Check netcdf-fortran functions for errors
!**********************************************************
!**********************************************************


subroutine check(status)

use netcdf
  
integer, intent (in) :: status

if(status /= nf90_noerr) then
   print*, trim(nf90_strerror(status))
   print*, status
   stop "Stopped"
end if


end subroutine check


!**********************************************************
!**********************************************************
! SUBROUTINE ingest_match_aerosol_arrays:
!
! Read daily file of hourly resolution MATCH aerosol
! data from appropriate file from the ASDC. Loads data
! into 5-d (profiles) and 4-d (aods) arrays.
!
!**********************************************************
!**********************************************************

subroutine ingest_match_aerosol_arrays

use netcdf
use pcf, only    : GetFileName
!use ceres_status, only: DIRECT_FILE, FMT_FILE, READ_FILE, SEQ_FILE, UNFMT_FILE, WRITE_FILE
!use io, only: CloseFile, OpenFile

integer, parameter :: incpath_yes = 1             ! path appended to file name? 
integer, parameter :: lid_file_match = 4          ! PCF MATCH file logic ID
integer            :: ios                         ! input output status
integer            :: iver                        ! 
character*300 match_filename                      ! to be retrieved by GetFileName 

iver = 1

call GetFileName (logicID = lid_file_match,      & ! I
                  filename= match_filename,      & ! O
                  status  = ios,                 & 
                  incpath = incpath_yes,         &
                  version = iver)  

print*, ios, match_filename 

if ( ios .ne. 0 ) stop ' ERROR OPENING MATCH FILE... '

! open the netcdf file in read mode
call check( nf90_open(match_filename, NF90_NOWRITE, ncid) )

! get the id for each variable based on its name
call check( nf90_inq_varid(ncid, "lon", ilonvarid) )   ! lon var id
call check( nf90_inq_varid(ncid, "lat", ilatvarid) )   ! lat var id
call check( nf90_inq_varid(ncid, "ilev", iplevvarid) ) ! lev var id - nc files poorly named...
call check( nf90_inq_varid(ncid, "lev" , iplayvarid) ) ! lay var id
call check( nf90_inq_varid(ncid, "hybi", ihybivarid) ) ! hybi var id - lev
call check( nf90_inq_varid(ncid, "hybm", ihybmvarid) ) ! hybm var id - lay

! dimension inquiry
!call check( nf90_inq_dimid(ncid, "lon",ilondimid) )
!call check( nf90_inq_dimid(ncid, "lat",ilatdimid) )
!call check( nf90_inq_dimid(ncid, "lev",iplaydimid))
!call check( nf90_inq_dimid(ncid,"ilev",iplevdimid))
!call check( nf90_inq_dimid(ncid,"time",itimedimid))

!print*,"lat dim id", ilatdimid
!print*,"lon dim id", ilondimid
!print*,"lev dim id", iplaydimid
!print*,"ilevdim id", iplevdimid
!print*,"timedim id", itimedimid

!call check( nf90_inquire_dimension(ncid, ilondimid, lonname, lonlen) )
!call check( nf90_inquire_dimension(ncid, ilatdimid, latname, latlen) )
!call check( nf90_inquire_dimension(ncid,iplaydimid,playname,playlen) )
!call check( nf90_inquire_dimension(ncid,iplevdimid,plevname,plevlen) )
!call check( nf90_inquire_dimension(ncid,itimedimid,timename,timelen) )

!print*, "Latitude Name:", latname, "Latitude Length:", latlen
!print*, "Longitude Name:", lonname, "Longitude Length:", lonlen
!print*, "P Layer Name: ", playname, "P Layer Length:", playlen
!print*, "P Level Name: ", plevname, "P Level Length:", plevlen
!print*, "Time/Hr Name: ", timename, "Time Length: " , timelen

! print nc file and variable ids
print*, "========================================"
print*, " MATCH NetCDF file and variable IDs...  "
print*, "========================================"

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

! read AODs (total + 11 types)
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
print*, " Successfully read MATCH aerosol arrays...         "
print*, " Aerosol optical depth (type,nlon,nlat,ntime)      "
print*, " Aerosol profiles      (type,nlon,nlat,nlev,ntime) "
print*, "==================================================="


end subroutine ingest_match_aerosol_arrays


!**********************************************************
!**********************************************************
! SUBROUTINE match_with_ceres_fov:
!
! Extracts aerosol data coincident with the CERES FOV by
! matching them via lat/lon indices
!
! Calls functions: lon_ind_match_ceresfov
!                  lat_ind_match_ceresfov
!
! which accept CERES FOV lon and lat as inputs
!**********************************************************
!**********************************************************


subroutine match_with_ceres_fov


print*, "======================================================="
print*, "Getting MATCH aerosol profiles/AOD at CERES FOV...     "
print*, "======================================================="

! CERES FOV hr is offset from MATCH by 1 hr
 match_hr = mpro%fovhr + 1

print*, "MATCH Hour   :", match_hr
print*, "CERES SwathHr:", mpro%fovhr
print*, "CERES FOV Lon:", mpro%rlon
print*, "CERES FOV Lat:", mpro%rlat
print*, "----------------------------------------"

print*, "MATCH Longitude index, MATCH Longitude"
print*, lon_ind_match_ceresfov(mpro%rlon), lon(lon_ind_match_ceresfov(mpro%rlon))

print*, "MATCH Latitude index, MATCH Latitude "
print*, lat_ind_match_ceresfov(mpro%rlat), lat(lat_ind_match_ceresfov(mpro%rlat))

print*, "======================================================="

! aerosol profiles
do t = 1,np       ! aerosol (t)ypes
   do p = 1,nlay  ! (p)ressure layers         ! call functions that extract MATCH data at the CERES FOV
      mavp%aero_profs(t,p) = mavp%aero_array(t,lon_ind_match_ceresfov(mpro%rlon),lat_ind_match_ceresfov(mpro%rlat),p,match_hr)
   end do
end do

! aerosol optical depth for 11 species (index 0 refers to total AOD)
do t = 0,np   !note on indicies: aod%aod_array(1:12) but want mpro%aods_all(0:11)
   mpro%aods_all(t) = aod%aod_array(t+1,lon_ind_match_ceresfov(mpro%rlon),lat_ind_match_ceresfov(mpro%rlat),match_hr)
end do

print*, "================================================================="
print*, "Aerosol profiles on native MATCH pressure layers (showing 1-9)..."
print*, "================================================================="

do p = 1,nlay
   write(*,*) (mavp%aero_profs(t,p), t=1,9) ! just showing 9 since it fits on my small screen
end do

end subroutine match_with_ceres_fov


!**********************************************************
!**********************************************************
! SUBROUTINE regrid_from_match_to_fu:
!
! Regrids the MATCH aerosol vertical profiles (reported on
! 28 pressure layers) to Fu-Liou RT model pressure layers
! (variable number of Fu layers)
!
!**********************************************************
!**********************************************************


subroutine regrid_from_match_to_fu


print*, "======================================================================"
print*, "Aerosol profiles regridded to Fu-Liou pressure layers (showing 1-9)..."
print*, "======================================================================"

psfc = mpro%pp(mpro%nv1)

!print*, "psfc:", psfc

mpro%mvp = 0.0
plast = 0.0                       ! mpro%pp(1) <- TOA

OUTLEV : do i = 1,mpro%nv1-1      ! loop over Fu-Liou model layers starting at TOA
            p1 = mpro%pp(i)       ! p1 < p2
            p2 = mpro%pp(i+1)     ! p2
            dpp(i) = p2-p1        ! calculate pressure layer thickness profile

            asm1=0
            nsm1=0

            hybi(1) = 0.00        ! TOA at zero pressure

INLEV  : do k = 1,nlay            ! loop over MATCH layers starting at TOA
            pt = hybi(k)*psfc     ! p-top = MATCH sigma level * Fu-Liou psfc
            pb = hybi(k+1)*psfc   ! p-bottom = next MATCH sigma level * Fu-Liou psfc

            if ( k .eq. nlay)  pb = psfc  ! p-bottom = psfc at end of profile
            dp = 0

            ! compare Fu-Liou & MATCH pressure level positions
            ! to determine dp at each step
            if( (p1 >= pt .and. p2 <= pb ) .or. &
                (pt >= p1 .and. pb <= p2 ) .or. &
                (p1 >= pt .or. p2 >= pb )  .and. ( p1<=pb .or. p2<= pb) ) then
               
            dp = pb - plast
                                                                                                                   
            if ( dp == 0 ) cycle   ! exit INLEV do loop
            
            ! if Fu-Liou < MATCH pressure
            if ( p2 < pb ) dp = p2 - plast

            plast = plast + dp     ! increment plast by dp

            ! compute weighted avg
            asm1(1:np) = asm1(1:np) + mavp%aero_profs(1:np,k)*dp ! for all constituents, multiply by dp thickness
            nsm1 = nsm1 + dp

            end if

end do INLEV

            ! assign aerosol value to Fu-Liou model layer
            if (nsm1 > 0) mpro%mvp(1:np,i) =  asm1(1:np) / nsm1

end do OUTLEV

!do ip=1,np
!      if ( sum( mpro%mvp(1:mpro%nv1-1,ip)) < 1.0E-30 ) mpro%mvp(ilev_trop,ip) =1.0E-30
!end do


! print the resulting profiles on Fu-Liou model layers
do p = 1,mpro%nv1-1
   write(*,*) ( mpro%mvp(t,p), t=1,9)
end do


end subroutine regrid_from_match_to_fu


!**********************************************************
!**********************************************************
! SUBROUTINE:
! Convert aerosol mass mixing ratio profiles to AOT profile
!**********************************************************
!**********************************************************


subroutine convert_profiles_to_aot


print*, "================================================================="
print*, "Aerosol profiles converted to AOT (showing 1-9)...               "
print*, "================================================================="

! normalize Fu-Liou dpp profile by TOA-surface p difference
dpp(1:mpro%nv1-1) = dpp(1:mpro%nv1-1)/sum(dpp(1:mpro%nv1-1))

! do for all aerosol (t)ypes...
do t=1,np

   ! multiply constituent profiles by normalized dpp profile
   mpro%mvp(t, 1:mpro%nv1-1 )= mpro%mvp(t, 1:mpro%nv1-1) * dpp(1:mpro%nv1-1)

   ! 'X' = Total AOD / Total Mass
   amr(t) = mpro%aods_all(t)/sum( mpro%mvp(t, 1:mpro%nv1-1 ) )

   ! multiply profile by 'X' to get AOD profile
   if (amr(t) > 0) mpro%mvp(t, 1:mpro%nv1-1 ) = mpro%mvp(t, 1:mpro%nv1-1 )* amr(t)

end do

! print the resulting AOD profiles
do p = 1,mpro%nv1-1
   write(*,*) (mpro%mvp(t,p), t=1,9) ! just showing 9 since it fits on my small 13" screen
end do


end subroutine convert_profiles_to_aot


!**********************************************************
!**********************************************************
! SUBROUTINE:
! Combining aerosol profiles into OPAC types for Ed4
!**********************************************************
!**********************************************************


subroutine combine_aero_profiles


print*, "=============================================================="
print*, "Comparing 11 AODs: (i) sum over AOD profiles, (ii) MATCH...   "
print*, "=============================================================="

!do t=1,np
!   if( sum( mpro%mvp(1:mpro%nv1-1,t)) < 1.0E-30 ) mpro%mvp(ilev_trop,t) =1.0E-30
!end do


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

   ! AOD of each atmospheric layer
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

! renormalize profiles by total AOD and convert to %
! but don't divide by zero... to avoid NaN output
do t = 1,np
   do k = 1,mpro%nv1-1
  if (tta(t)>0)  mpro%mvp(t,k) = 100*mpro%mvp(t,k) / tta(t)
   end do
end do

do t = 1,npc
   do k = 1,mpro%nv1-1
   if (tta_comb(t)>0)  mpro%mvp_comb(t,k) = 100*mpro%mvp_comb(t,k) / tta_comb(t)
   end do
end do


! if profile sums to zero replace with extremely small non-zero values...
do t=1,npc
      if ( sum( mpro%mvp_comb(t,1:mpro%nv1-1)) < 1.0E-30 ) mpro%mvp_comb(t,ilev_trop) =1.0E-30
end do



! get combined AODs
mpro%aods_comb(1:npc) = mpro%aods_comb(1:npc)
mpro%aods_comb(0)     = sum(mpro%aods_comb(1:npc))

! Fu-Liou aerosol types
mpro%ityp_comb(1:npc) =  (/ 24, 25, 18, 1, 11, 10, 9 /)


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


! change dimensions to those required by Fu-Liou fi%
do k = mpro%nv1-1,1,-1
   do t = 1,npc
      mpro%mvp_comb1(k,t) = mpro%mvp_comb(t,k)
   end do
end do

!print*, "Vertical profiles for 7 species for input to Fu-Liou..."
!print*, "Profile order is from TOA to surface"
!do t = 1,npc  !dimensions now: model pressure layer, aerosol type
!   write(*,*) (mpro%mvp_comb1(k,t), k = 1,mpro%nv1-1)
!end do

print*, "AODs combined:", mpro%aods_comb

!where (mpro%aods_comb == 0.0) mpro%aods_comb = 1E-9


end subroutine combine_aero_profiles


!**********************************************************
!**********************************************************
!
! FUNCTION lon_ind_match_ceresfov:
!
! Input:  CERES FOV longitude
! Output: Lon index of closest MATCH gridbox
!
! Calculate the distance from the CERES FOV longitude to
! MATCH longitude centers and return the index of the
! MATCH gridbox nearest to the CERES FOV.
!
!**********************************************************
!**********************************************************


function lon_ind_match_ceresfov(fovlon) result(ilon)
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
end function lon_ind_match_ceresfov


!**********************************************************
!********************************************************** 
!FUNCTION lat_ind_match_ceresfov:
!
! Input:  CERES FOV latitude
! Output: Lat index of closest MATCH gridbox
!
! Calculate the distance from the CERES FOV latitude to
! MATCH latitude centers and return the index of the
! MATCH gridbox nearest to the CERES FOV. 
!
!**********************************************************
!**********************************************************


function lat_ind_match_ceresfov(fovlat) result(ilat)
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
end function lat_ind_match_ceresfov


!********************************************************** 

end module MATCH_Ed4










