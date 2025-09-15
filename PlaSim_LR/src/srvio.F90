  module srvio
! module srvio
!
! author: FastOpt GmbH
!
! purpose: provide routines for direct access of service-formatted
!          data files in connection with PlaSim
!
! contains:
!      subroutine initsrvio(kunit,yname)
!      subroutine postsrvio
!      subroutine writesrv(krec,kcode,klev,fwt)
!      subroutine readsrv(krec,kcode,klev,frd)
!      subroutine code2field(kcode,klev,frd)
!      subroutine diag_geop(pgeop,klev)
! public: initsrvio,postsrvio,writesrv,readsrv,code2field
! private: everything else
!
    use pumamod
    implicit none
    integer,parameter :: ioprec=8            ! precision of data fields
    integer(kind=4)   :: ircw1,ircw2,ircw3,ircw4 ! fortran record control words
!    integer           :: ihead(8)            ! SERVICE header always kind=4
    integer (kind=4)          :: ihead(8)            ! SERVICE header 
	real(kind=ioprec) :: fio(NLON*NLAT)      ! data field 
    integer :: iunit
    real :: doro(NHOR) ! gridpoint orography for calculation of geopot. height.
!!FL
!   add variable from simba:
    
    real :: dveg(NHOR)
    real :: dlai(NHOR)
    real :: dcveg(NHOR)
    real :: dcsoil(NHOR)
!!

! constants from burn5.cpp
    real,parameter :: mars_rd =      (189.0 )
    real,parameter :: rkbol =  (1.380658e-23)
    real,parameter :: rnavo = (6.0221367e+23)
    real,parameter :: r =     (rkbol * rnavo)
    real,parameter :: rmd =   (28.9644)
    real,parameter :: rmv =   (18.0153)
    real,parameter :: earth_rd = (1000. * r / rmd)
    real,parameter :: earth_grav = 9.80665
    real,parameter :: mars_grav =  3.7
    real:: vtmp,rga
    real :: rd ! should be plasimmod:gascon, but postprocessor uses slightly different value
    private
    public initsrvio,postsrvio,writesrv,readsrv,code2field
  contains

    subroutine initsrvio(kunit,yname,ystatus)
! Initialisation of module: opening file, computation of internal constants
! input:
!       integer :: kunit            ! unit number of data file
!       character(len=*) :: yname   ! name of data file
!       character(len=*) :: ystatus ! status of data file
!
      integer :: kunit
      character(len=*) :: yname
      character(len=*) :: ystatus
      integer :: recl
      iunit = kunit
! open file: fortran blocked SERVICE-format file in direct access mode
!      if(mypid==NROOT) then
       inquire(iolength=recl) ircw1,ihead(:),ircw2,ircw3,fio(:),ircw4
       open (iunit,file=yname,ACCESS='DIRECT',RECL=recl,STATUS=ystatus, convert='BIG_ENDIAN')
!      endif

! set up physical constants as defined in burn5.cpp
      if ( mars == 1 ) then
         rd =  mars_rd
         rga = 1. / mars_grav
      else
         rd = earth_rd
         rga = 1. / earth_grav
      end if
      vtmp = (rv / rd) - 1.0
! provide gridpoint orography
      doro(:) = 0.
      call sp2fc(so(:),doro)
      call fc2gp(doro(:),NLON,NLPP)
      doro(:) = doro(:) * CV * CV
      return
    end subroutine initsrvio

    subroutine postsrvio
! Finalisation of module, closing of data file
     if(mypid==NROOT) then
      close (iunit)
     endif
    end subroutine postsrvio

    subroutine writesrv(krec,kcode,klev,fwt)
! write service-header and data to file
! input:
!       integer :: krec      ! number of direct access record
!       integer :: kcode     ! code number of data field
!       integer :: klev      ! level identifier of data field
!       real    :: fwt(NHOR) ! data to be written
      integer, intent(in) :: krec,kcode,klev
      integer :: nmin,nhour,nday,nmonth,nyear
      integer :: istep,irec
      real, intent(in) :: fwt(NLON*NLAT)
      character(len=256) :: yname

      if(mypid /= NROOT) then
       print*,'ERROR: writesrv not called by root... stopped'
       stop
      endif

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear) ! timestep to date conversion

      ihead(1) = kcode
      ihead(2) = klev
      ihead(3) = nday + 100 * nmonth + 10000 * nyear
!      ihead(4) = nmin + 100 * nhour
      ihead(4) = 100 * (nmin + 100 * nhour)  ! alternative time format including seconds: hhmmss
      ihead(5) = NLON
      ihead(6) = NLAT
      ihead(7) = 1
      ihead(8) = n_days_per_year

      ircw1 = kind(ihead(1))*8
      ircw2 = ircw1
      ircw3 = ioprec*NLON*NLAT
      ircw4 = ircw3
! fortran blocked SERVICE-format file in direct access mode
      inquire(unit=iunit,name=yname)
      print '(1x,a,i8,a,a)','writing obs record # ',krec,' to ',trim(yname)
!      write (iunit,rec=krec) ircw1,ihead(:),ircw2,ircw3,real(fwt(:),kind=ioprec),ircw4 ! gives problems with lf95 --dbl
      fio(:) = fwt(:) ! precison conversion
      write (iunit,rec=krec) ircw1,ihead(:),ircw2,ircw3,fio(:),ircw4
!
!      inquire(unit=iunit,name=yname)
       print '(1x,a,i8,a,a)','wrote obs record # ',krec,' to ',trim(yname)
!
      return
    end subroutine writesrv


    subroutine readsrv(krec,kcode,klev,frd)
! read service-header and data from file
! input:
!       integer :: krec  ! number of direct access record
!       integer :: kcode ! code number of data field
!       integer :: klev  ! level identifier of data field
! output:
!       real :: frd(NLON*NLAT) ! data read
!
      integer, intent(in) :: krec,kcode,klev
      integer :: nmin,nhour,nday,nmonth,nyear
      integer :: istep
      integer :: jhead(6)
      real, intent(out) :: frd(NLON*NLAT)
      character(len=256) :: yname

      ! if(mypid /= NROOT) then
      !  print*,'ERROR: readsrv not called by root... stopped'
      !  stop
      ! endif

      istep = nstep
      call ntomin(istep,nmin,nhour,nday,nmonth,nyear) ! timestep to date conversion

      jhead(1) = kcode
      jhead(2) = klev
      jhead(3) = nday + 100 * nmonth + 10000 * nyear
!      jhead(4) = nmin + 100 * nhour
      jhead(4) = 100 * (nmin + 100 * nhour)  ! alternative time format including seconds: hhmmss
      jhead(5) = NLON
      jhead(6) = NLAT

! fortran blocked SERVICE-format file in direct access mode
      !print*, "iunit in readsrv", iunit
      ! print*, "ihead in readsrv", ihead(:) 
      read (iunit,rec=krec) ircw1,ihead(:),ircw2,ircw3,fio(:),ircw4
      !print*, "record in readsrv", krec
          
      if (any(ihead((/1,2,3,5,6/)).ne.jhead((/1,2,3,5,6/))).or.(ihead(4).ne.jhead(4).and.100*ihead(4).ne.jhead(4))) then
        inquire(unit=iunit,name=yname)
        print'(1x,a,a)',"header mismatch while reading from ",trim(yname)
        print*,'record =',krec
        print'(a,6i9)','wanted: ',jhead
        print'(a,6i9)','found:  ',ihead(1:6)
      end if
      if (ircw1.ne.ircw2) then
        print *,'RCW mismatch 1',ircw1,ircw2
        stop
      end if
      if (ircw3.ne.ircw4) then
        print *,'RCW mismatch 2',ircw3,ircw4
        stop
      end if
      
      frd(:) = fio(:) ! precision conversion

      return
    end subroutine readsrv


    subroutine code2field(kcode,klev,frd)
! get data corresponding to a specified code number and level
! input:
!     integer :: kcode  ! code number of requested data field
!     integer :: klev   ! level identifier of requested data field
! output:
!     real :: frd(NHOR) ! requested data field
!
      integer, intent(in) :: kcode,klev
      real :: zsave
      real, intent(out) :: frd(NHOR)

      frd = 0.

      select case (kcode)

!     **********************************
!     * mixed-layer depth (from ocean) *
!     **********************************
      case (110)
         frd(:) = dmld

!     *******************  
!     * air temperature *
!     *******************
      case (130)
         frd = dt(:,klev)

!     **********  
!     * u-wind *
!     **********
      case (131)
         frd = du(:,klev)

!     **********  
!     * v-wind *
!     **********
      case (132)
         frd = dv(:,klev)

!     ******************  
!     * spec. humidity *
!     ******************
      case (133)
         frd = dq(:,klev)

!     ********************  
!     * surface pressure *
!     ********************
      case (134)
         frd = dp(:)

!     ******************  
!     * atm. vorticity *
!     ******************
      case (138)
!$taf incompl sz
         call mpgasp(sz(:,klev),szp(:,klev),1)
         zsave = sz(3,klev)
         sz(3,klev) = sz(3,klev) - EZ ! convert absolute to relative vorticity
         call sp2fc(sz(:,klev),frd(:))
         call fc2gp(frd(:),NLON,NLPP)
         sz(3,klev) = zsave ! restore absolute vorticity
         frd = frd * ww

!     ***********************
!     * surface temperature *
!     ***********************
      case (139)
         frd(:)=dt(:,NLEP)

!     ****************
!     * soil wetness *
!     ****************
      case (140)
         frd(:) = dwatc

!     **************
!     * snow depth *
!     **************
      case (141)
         frd(:) = dsnow

!! Armin's code ==============================
      case(164)
         frd(:) = dcc(:,NLEP)   !acc(:)/real(naccuout)

      case (176)
         frd(:) = dswfl(:,NLEP) ! assol(:)/real(naccuout) !
! conflicts with nafter 
! conflicts with nafter !     *****************************
! conflicts with nafter !     * surface thermal radiation *
! conflicts with nafter !     *****************************
      case (177)
         frd(:) =  dlwfl(:,NLEP) ! asthr(:)/real(naccuout) !
!         frd(:) =  dlwfl(:,klev) ! asthr(:)/real(naccuout) !
! conflicts with nafter 
! conflicts with nafter !     ***********************
! conflicts with nafter !     * top solar radiation *
! conflicts with nafter !     ***********************
       case (178)
         frd(:) = dswfl(:,1) ! atsol(:)/real(naccuout) !
! conflicts with nafter 
! conflicts with nafter !     *************************
! conflicts with nafter !     * top thermal radiation *
! conflicts with nafter !     *************************
       case (179)
          frd(:) = dlwfl(:,1) ! atthr(:)/real(naccuout) !
 
!!===============================

! conflicts with nafter !     **********************
! conflicts with nafter !     * large scale precip *
! conflicts with nafter !     **********************
! conflicts with nafter       case (142)
! conflicts with nafter          frd(:) = aprl(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *********************
! conflicts with nafter !     * convective precip *
! conflicts with nafter !     *********************
! conflicts with nafter       case (143)
! conflicts with nafter          frd(:) = aprc(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *************
! conflicts with nafter !     * snow fall *
! conflicts with nafter !     *************
! conflicts with nafter       case (144)
! conflicts with nafter          frd(:) = aprs(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     **********************
! conflicts with nafter !     * sensible heat flux *
! conflicts with nafter !     **********************
! conflicts with nafter       case (146)
! conflicts with nafter          frd(:) = ashfl(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ********************
! conflicts with nafter !     * latent heat flux *
! conflicts with nafter !     ********************
! conflicts with nafter       case (147)
! conflicts with nafter          frd(:) = alhfl(:)/real(naccuout)

!     ***********************  
!     * ln surface pressure *
!     ***********************
      case (152)
         frd = log(dp(:))

!     *******************  
!     * atm. divergence *
!     *******************  
      case (155)
!$taf incompl sd
         call mpgasp(sd(:,klev),sdp(:,klev),1)
         call sp2fc(sd(:,klev),frd(:))
         call fc2gp(frd(:),NLON,NLPP)
         frd = frd * ww

!     ***********************
!     * geopotential height *
!     ***********************

      case(156)

         call diag_geop(frd,klev)

!     ************************
!     * liquid water content *
!     ************************
      case (161)
         frd(:) = dql(:,klev)

!     *************
!     * u-star**3 *
!     *************
      case (159)
         frd(:) = dust3

! conflicts with nafter !     **********
! conflicts with nafter !     * runoff *
! conflicts with nafter !     **********
! conflicts with nafter       case (160)
! conflicts with nafter          frd(:) = aroff(:)/real(naccuout)

!     ***************
!     * cloud cover *
!     ***************
      case (162)
         frd(:) = dcc(:,klev)

! conflicts with nafter !     *********************
! conflicts with nafter !     * total cloud cover *
! conflicts with nafter !     *********************
! conflicts with nafter       case(164)
! conflicts with nafter          frd(:) = acc(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ***************************
! conflicts with nafter !     * surface air temperature *
! conflicts with nafter !     ***************************
! conflicts with nafter       case (167)
! conflicts with nafter          frd(:) = atsa(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ******************************
! conflicts with nafter !     * surface temperature (accu) *
! conflicts with nafter !     ******************************
! conflicts with nafter       case (169)
! conflicts with nafter          frd(:) = ats0(:)/real(naccuout)

!     *************************
!     * deep soil temperature *
!     *************************
      case (170)
         frd(:) = dtd5

!     *****************
!     * land sea mask *
!     *****************
      case (172)
         frd(:) = dls

!     *********************
!     * surface roughness *
!     *********************
      case (173)
         frd(:) = dz0

!     **********
!     * albedo *
!     **********
      case (175)
         frd(:) = dalb

! conflicts with nafter !     ***************************
! conflicts with nafter !     * surface solar radiation *
! conflicts with nafter !     ***************************
! conflicts with nafter       case (176)
! conflicts with nafter          frd(:) = assol(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *****************************
! conflicts with nafter !     * surface thermal radiation *
! conflicts with nafter !     *****************************
! conflicts with nafter       case (177)
! conflicts with nafter          frd(:) = asthr(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ***********************
! conflicts with nafter !     * top solar radiation *
! conflicts with nafter !     ***********************
! conflicts with nafter       case (178)
! conflicts with nafter          frd(:) = atsol(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *************************
! conflicts with nafter !     * top thermal radiation *
! conflicts with nafter !     *************************
! conflicts with nafter       case (179)
! conflicts with nafter          frd(:) = atthr(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ************
! conflicts with nafter !     * u-stress *
! conflicts with nafter !     ************
! conflicts with nafter       case (180)
! conflicts with nafter          frd(:) = ataux(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *************
! conflicts with nafter !     * v- stress *
! conflicts with nafter !     *************
! conflicts with nafter       case (181)
! conflicts with nafter          frd(:) = atauy(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ***************
! conflicts with nafter !     * evaporation *
! conflicts with nafter !     ***************
! conflicts with nafter       case (182)
! conflicts with nafter          frd(:) = aevap(:)/real(naccuout)

!     *********************
!     * soil temperature *
!     *********************
      case (183)
         frd(:) = dtsoil

!     ********************
!     * vegetation cover *
!     ********************
      case (199)
         frd(:) = dveg

!     *******************
!     * leaf area index *
!     *******************
      case (200)
         frd(:) = dlai

! conflicts with nafter !     ********************
! conflicts with nafter !     * top solar upward *
! conflicts with nafter !     ********************
! conflicts with nafter       case (203)
! conflicts with nafter          frd(:) = atsolu(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     ************************
! conflicts with nafter !     * surface solar upward *
! conflicts with nafter !     ************************
! conflicts with nafter       case (204)
! conflicts with nafter          frd(:) = assolu(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     **************************
! conflicts with nafter !     * surface thermal upward *
! conflicts with nafter !     **************************
! conflicts with nafter       case (205)
! conflicts with nafter          frd(:) = asthru(:)/real(naccuout)

!     *******************************
!     * soil temperatures level 2-4 *
!     *******************************
      case (207)
         frd(:) = dtd2
      case (208)
         frd(:) = dtd3
      case (209)
         frd(:) = dtd4

!     *****************
!     * sea ice cover *
!     *****************
      case (210)
         frd(:) = dicec

!     *********************
!     * sea ice thickness *
!     *********************
      case (211)
         frd(:) = diced

!     ****************
!     * forest cover *
!     ****************
      case (212)
         frd(:) = dforest

! conflicts with nafter !     *************
! conflicts with nafter !     * snow melt *
! conflicts with nafter !     *************
! conflicts with nafter       case (218)
! conflicts with nafter          frd(:) = asmelt(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *********************
! conflicts with nafter !     * snow depth change *
! conflicts with nafter !     *********************
! conflicts with nafter       case (221)
! conflicts with nafter          frd(:) = asndch(:)/real(naccuout)

!     ******************
!     * field capacity *
!     ******************
      case (229)
         frd(:) = dwmax

! conflicts with nafter !     *****************************************
! conflicts with nafter !     * vertical integrated specific humidity *
! conflicts with nafter !     *****************************************
! conflicts with nafter       case (230)
! conflicts with nafter          frd(:) = aqvi(:)/real(naccuout)

!     ****************
!     * glacier mask *
!     ****************
      case (232)
         frd(:) = dglac

!     *********************
!     ***   S I M B A   ***
!     *********************

! conflicts with nafter !     ****************************
! conflicts with nafter !     * gross primary production *
! conflicts with nafter !     ****************************
! conflicts with nafter       case (300)
! conflicts with nafter          frd(:) = agpp(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     **************************
! conflicts with nafter !     * net primary production *
! conflicts with nafter !     **************************
! conflicts with nafter       case (301)
! conflicts with nafter          frd(:) = anpp(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *********************
! conflicts with nafter !     * light limited GPP *
! conflicts with nafter !     *********************
! conflicts with nafter       case (302)
! conflicts with nafter          frd(:) = agppl(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *********************
! conflicts with nafter !     * water limited GPP *
! conflicts with nafter !     *********************
! conflicts with nafter       case (303)
! conflicts with nafter          frd(:) = agppw(:)/real(naccuout)

!     *********************
!     * vegetation carbon *
!     *********************
      case (304)
         frd(:) = dcveg

!     ***************
!     * soil carbon *
!     ***************
      case (305)
         frd(:) = dcsoil

! conflicts with nafter !     ************************
! conflicts with nafter !     * no growth allocation *
! conflicts with nafter !     ************************
! conflicts with nafter       case (306)
! conflicts with nafter          frd(:) = anogrow(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *****************************
! conflicts with nafter !     * heterotrophic respiration *
! conflicts with nafter !     *****************************
! conflicts with nafter       case (307)
! conflicts with nafter          frd(:) = aresh(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     *********************
! conflicts with nafter !     * litter production *
! conflicts with nafter !     *********************
! conflicts with nafter       case (308)
! conflicts with nafter          frd(:) = alitter(:)/real(naccuout)
! conflicts with nafter 
! conflicts with nafter !     **************
! conflicts with nafter !     * water loss *
! conflicts with nafter !     **************
! conflicts with nafter       case (309)
! conflicts with nafter          frd(:) = awloss(:)/real(naccuout)

      case (278)             ! top of atmosphere energy balance
         frd(:) = dswfl(:,1)+ dlwfl(:,1) ! atsol(:)/real(naccuout) !

      case default
         print *,"This code not implemented in code2field: ",kcode
         stop
      end select
      return
    end subroutine code2field


    subroutine diag_geop(pgeop,klev)
!
! purpose: compute the geopotential on level klev (model level if klev<=NLEV, pressure level if otherwise)
!          in the same way as the plasim postprocessor burn5 does it.
!
! input:
!      integer :: klev     ! level identifier OR pressure level of requested geopotential
! output:
!      real :: pgeop(NHOR) ! requested geopotential
!
      implicit none
      real :: pgeop(NHOR)
      integer :: klev ! pressure or model level
      integer :: ilev
      integer :: jhor, jlev
      real :: zgeop(NHOR,NLEV+1)
      real :: alpha,tstar,tmsl,ZALP,ZALPAL
! dt is present
! dq is present
! dp is present
      if (klev.le.NLEV) then
         ilev = max(2,klev)
      else
         ilev = 2 
      end if
! surface 
      zgeop(:,NLEV+1) = doro(:) ! surface geopotential
! atm. half levels
      do jlev=NLEV,ilev,-1
         zgeop(:,jlev)=zgeop(:,jlev+1)+ &
              RD*dt(:,jlev)*(1.+VTMP*dq(:,jlev))*log(sigmah(jlev)/sigmah(jlev-1))
      end do
! top level (if pressure level or level 1 is requested)
      if (klev.gt.NLEV.or.klev==1) zgeop(:,1)=zgeop(:,2)+ RD*dt(:,1)*(1.+VTMP*dq(:,1))* 2. * log(2.)
!
! model or pressure level requested?   
!
      if (klev.le.NLEV) then ! assume model level requested
         pgeop(:) = zgeop(:,klev) * rga ! m**2/s**2 -> m
      else
!
! assume pressure level is requested and klev is in Pa
!
         do jhor=1,NHOR
! now find the neighbouring levels jlev and jlev+1
            jlev = 1
            do while (dp(jhor)*sigmah(jlev+1) < klev .and. jlev < NLEV-2)
               jlev = jlev + 1
            end do
            if (dp(jhor)*sigmah(jlev+1) < klev) then
!
! extrapolate below surface level
!
               alpha = RD * alr * rga
               tstar = (1.0 + alpha * (1./sigma(NLEV) - 1.0)) * dt(jhor,NLEV)
               if (tstar < 255.0) tstar = 0.5 * (255.0 + tstar)
               tmsl = tstar + alr * rga * doro(jhor);
               if (tmsl > 290.5 .and. tstar > 290.5) then
                  tstar = 0.5 * (290.5 + tstar)
                  tmsl  = tstar
               end if
               if (tmsl > 290.5 .and. tstar <= 290.5) tmsl = 290.5;
               if (tmsl-tstar < 0.000001 .and. tstar-tmsl < 0.000001) then
                  alpha = 0.0
               else 
                  if (doro(jhor) > 0.0001 .or. doro(jhor) < -0.0001) alpha = RD * (tmsl-tstar) / doro(jhor)
               end if
               ZALP   = log(klev/dp(jhor))
               ZALPAL = ZALP * alpha
               pgeop(jhor) = (doro(jhor) - RD * tstar * ZALP * (1.0 + ZALPAL * (0.5 + ZALPAL/6.0))) * rga
            else
! interpolate linearly in ln(p) coordinates
!            pgeop(jhor)= rga * &
!               (  (zgeop(jhor,jlev+1)-zgeop(jhor,jlev)) *   &
!                 log(klev/sigma(jlev)/dp(jhor)) /           &
!                 log(sigma(jlev+1)/sigma(jlev))+zgeop(jhor,jlev) )
!
! interpolate linearly in p coordinates
!
               pgeop(jhor) = rga * &
                    (  zgeop(jhor,jlev+1) + ( real(klev)/dp(jhor) - sigmah(jlev) ) &
                    * (zgeop(jhor,jlev+2) - zgeop(jhor,jlev+1)) &
                    / (    sigmah(jlev+1) -     sigmah(jlev)) )
            end if
         end do
      end if
    end subroutine diag_geop
      
  end module srvio
