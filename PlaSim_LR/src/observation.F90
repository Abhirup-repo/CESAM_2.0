!-*-f90-*-

! This file is cpp-preprocessed. Be sure to edit observation.F90
      module observation !FastOpt
! module observation
!
! author: FastOpt GmbH
! modified by larry 03.01.2013
!
! purpose: provide routines for computing the cost function based on observations
!
! contains:
!        subroutine initobs
!        subroutine reinitobs
!        subroutine genav
!        subroutine resetav
!        subroutine writeobs(step)
!        subroutine addcost(step)
!        subroutine postobs
! private: obscode,obslev,lobsav,iobsav,irec,iobs,miss,fobserr
! public: everything else
!
        use pumamod, only : NHOR,NLON,NLAT,NLPP,mypid,NROOT,nud
        implicit none
        integer, parameter :: nobs = 41 ! number of variables in single observation
        integer :: obscode(nobs) ! list of codes in a single observation
        integer :: obslev(nobs)  ! list of levels corresponding to codes
        logical :: lobsav(nobs)  ! is this field a time average (T/F)?
        integer :: iobsav(nobs)  ! index of time averaged quantity
        real,allocatable :: fobsav(:,:) ! time averaged field, partial
        real :: fobserr(NHOR,nobs) ! observational uncertainties, time independent
        real :: fobserrfac(nobs) ! factor on observational errors (1/k**2) for experimenting
        integer :: countav       ! counter for averaging period
        integer :: uflux = 60 !201
        real :: gausswsum        ! sum of gaussw
        real :: cost,costm
        integer :: irec
        integer :: iobs
        integer :: nobsav
        real,parameter :: miss = -9.d+33 ! cdo default missing value for srv_32
        character(len=*), parameter :: yname = 'observations.srv'
        
        private obscode,obslev,lobsav,iobsav,irec,iobs,miss,fobserr

      contains 

        subroutine initobs
!
! set up structure of observation - this defines most of the cost function!
!
! IMPORTANT: The number and sequence of records must reflect exactly the structure of
!            the data in "observations.srv".
!
! files read (optional): cost-function-definition-atm.dat
! files written:         cost-function-definition-atm.dat
!
          use pumamod, only : gaussw,ncost,nstepe,nstep0,nafter,mypid,NROOT
          use srvio, only : initsrvio
          integer :: jdim
          integer :: recl
          integer :: iav
          integer :: junit
          logical :: lexist
          character(len=15) :: ynam
          character(len=32),parameter :: yfile = 'cost-function-definition-atm.dat'
          integer, external :: get_free_unit
!
          integer :: initialize 
          real :: zfobserr(nobs)
!
! set up structure of observation - this defines most of the cost function!
!
! IMPORTANT: The number and sequence of records must reflect exactly the structure of 
!            the data in "observations.srv".
!
          if(mypid==NROOT) then
           initialize=1 
           inquire(file=yfile,exist=lexist)
           if (lexist) initialize=0
           junit=1993
           open(junit,file=yfile,form='formatted')
          endif
          call mpbci(initialize)
          if(initialize==0) then
           if(mypid==NROOT) then
             print*,'Reading cost function definition from: ',trim(yfile)
             read(junit,*) jdim
             if (jdim.ne.nobs) stop 'incompatible dimension!'
             read(junit,*) ynam
             if(trim(ynam).ne.'obscode') stop 'wrong tag!'
             read(junit,1000) obscode(:)
             read(junit,*) ynam
             if(trim(ynam).ne.'obslev') stop 'wrong tag!'
             read(junit,1000) obslev(:)
             read(junit,*) ynam
             if(trim(ynam).ne.'lobsav') stop 'wrong tag!'
             read(junit,1010) lobsav(:)
             read(junit,*) ynam
             if (trim(ynam)=='fobserrfac') then
                read(junit,*) fobserrfac(:)
                print*,'Distorting per observation errors with the following factors:'
                print 1020,fobserrfac(:)
                read(junit,*) ynam
             else
                fobserrfac(:) = 1.
             end if
             if(trim(ynam).ne.'fobserr') stop 'wrong tag!'
             read(junit,1020) zfobserr(:)
           endif
           call mpbcrn(zfobserr,nobs)            
           call mpbcrn(fobserrfac,nobs)
           call mpbci(jdim) 
           call mpbcin(obscode,nobs)
           call mpbcin(obslev,nobs)
           call mpbcln(lobsav,nobs)
           fobserr(:,:) = spread(zfobserr(:),dim=1,ncopies=NHOR)

          else
           obscode(:) = -1 ! mark as uninitialised
           obslev(:)  = -1 ! mark as uninitialised
           fobserr(:,:) = -1.  ! mark as uninitialised
           lobsav(:) = .true. ! eventually mark quantity as time average (to be computed!)
           fobserrfac(:) = 1. ! factor on observational errors for experimenting (1/k**2)
! order in dats:  138 155 130 133 134 #210 141
! vorticity
           obscode(1:10) = 138   ! code
           obslev(1) = 1     ! level
           obslev(2) = 2     ! level
           obslev(3) = 3     ! level
           obslev(4) = 4     ! level
           obslev(5) = 5     ! level
           obslev(6) = 6     ! level
           obslev(7) = 7     ! level
           obslev(8) = 8     ! level
           obslev(9) = 9     ! level
           obslev(10) = 10     ! level
           fobserr(:,1:10) = 1./(2.e-5**2)  ! set up first guess inverse error covariance 'matrix'
! test, FastOpt ! divergence
! test, FastOpt           obscode(11:20) = 155   ! code
! test, FastOpt           obslev(11:20) = obslev(1:10)     ! level
! test, FastOpt           fobserr(:,11:20) = 1./(3.e-6**2)  ! set up first guess inverse error covariance 'matrix'
! temperature
!          obscode(21:30) = 130   ! code
!          obslev(21:30) = obslev(1:10)  ! level
!          fobserr(:,21:30) =  1./(3.**2)  ! set up first guess inverse error covariance 'matrix'
           obscode(11:19) = 130   ! code
           obslev(11:19) = obslev(2:10)  ! level
           fobserr(:,11:19) =  1./(3.**2)  ! set up first guess inverse error covariance 'matrix'
!!! surface pressure
           obscode(20) = 134   ! code
           obslev(20) = 0  ! level
!          fobserr(:,41) = NLEV/(5.e+2**2)  ! set up first guess inverse error covariance 'matrix'
           fobserr(:,20) = 1./(5.e+2**2)  ! set up first guess inverse error covariance 'matrix'
           if(mypid==NROOT) then
            print*,"Writing cost function definition to: ",trim(yfile)
            write(junit,*)nobs
            write(junit,*)"obscode"
            write(junit,1000) obscode
            write(junit,*)"obslev"
            write(junit,1000) obslev
            write(junit,*)"lobsav"
            write(junit,1010) lobsav
            if (any(fobserrfac.ne.1.)) then
             write(junit,*)"fobserrfac"
             write(junit,1020) fobserrfac(:)
            end if
            write(junit,*)"fobserr"
            write(junit,1020) fobserr(1,:)
           endif
          end if
          if(mypid==NROOT) then
           close(junit)
          endif
!
! digest the above
!
! check for completeness
           if (any(obscode(:).lt.0).or.any(obslev(:).lt.0).or.any(fobserr(:,:).lt.0.)) then
              call mpabort ('definition of observations contains errors!')
           end if
! allocate field for time-averaged observations
          nobsav=count(lobsav(:))
          if (nobsav.gt.0) then
             allocate(fobsav(NHOR,nobsav))
             fobsav(:,:) = 0.
          else                     ! workaround to avoid error on store
             allocate(fobsav(1,1)) ! workaround to avoid error on store
          end if
          countav = 0. ! initialise step counter for averages
! set up index array for time-averaged observations
          iav=0
          do iobs=1,nobs
             if (lobsav(iobs)) then
                iav = iav + 1
                iobsav(iav) = iobs
             end if
          end do
! open file for direct access SERVICE i/o
          uflux = 1994

          call initsrvio(uflux,yname,'OLD')
! set gausswsum
          gausswsum = sum(gaussw(1:nlpp))
         !  call mpsumr(gausswsum,1)
! reset cost function
          cost = 0.
1000      format (10(i6,1x))
1010      format (10(L6,1x))
1020      format (10(e16.10,1x))

        end subroutine initobs

        subroutine reinitobs
! reset counters and cost function for repeated cost function evaluation
          if (allocated(fobsav)) fobsav(:,:) = 0.
          countav = 0.
! reset cost function
          cost = 0.
          costm =0.
        end subroutine reinitobs

        subroutine genav
! compute average for model fields needed for cost function
          use pumamod, only : NHOR,NLON,NLAT,mypid,NROOT
          use srvio, only : code2field
          real :: frd(NHOR)
          do iobs=1,nobs
             if (lobsav(iobs)) then ! accumulate for averaging
                call code2field(obscode(iobs),obslev(iobs),frd(:)) ! copy model field specified to frd(:)
                fobsav(:,iobsav(iobs)) = fobsav(:,iobsav(iobs)) + frd(:)
             end if
          end do
          countav = countav + 1
        end subroutine genav

        subroutine resetav
! reset averages and counter for computation of averages
          if (allocated(fobsav)) fobsav(:,:) = 0.
          countav = 0
        end subroutine resetav

        subroutine writeobs(step)
! write pseudo-observational data
! input:
!          integer :: step ! number of cost function contribution

          use pumamod, only : NHOR,NUGP,mypid,NROOT
          use srvio, only : writesrv,code2field
          integer :: step,iobs,ilev
          real :: fwtf(NUGP)
          real :: fwt(NHOR)
          do iobs=1,nobs
             if (lobsav(iobs)) then
                fwt(:) = fobsav(:,iobsav(iobs)) / countav ! compute time average and store in fwt()
             else
                call code2field(obscode(iobs),obslev(iobs),fwt(:)) ! copy specified model field to fwt(:)
             end if
             call mpgagp(fwtf,fwt,1)
             irec  = nobs*(step-1)+iobs
             if(mypid==NROOT) then
              call writesrv(irec,obscode(iobs),obslev(iobs),fwtf(:)) ! write fwt(:)
             endif
          end do
          call resetav
        end subroutine writeobs

        subroutine addcost(step)
! compute contribution to cost function
! input:
!          integer :: step ! number of cost function contribution
!
          use pumamod, only : NUGP,NLON,NLAT,NHOR,gaussw,lunstable,nstep0,nstep,ncost,NROOT,mypid
          use srvio, only : readsrv,code2field
          implicit none
          integer :: step
          real :: fm(NHOR)
          real :: frdf(NUGP)
          real :: frd(NHOR)
          real, dimension(NLON) :: fh
          real :: ch(1),ch1, gmm, gmo
          integer :: jlat,nval
          logical :: lval(nlon) ! value mask
!
          if(mypid==NROOT) then
           write (nud,'(1x,a,i5,a,i5)') 'Cost function contribution #',step,' step #',nstep-nstep0+1
          endif
!! $TAF STORE countav = obs1, REC= (nstep-nstep0+1)/ncost
!!$TAF STORE countav = obs1
          do iobs=1,nobs
             irec  = nobs*(step-1)+iobs
             if (mypid==NROOT) then
                call readsrv(irec,obscode(iobs),obslev(iobs),frdf(:)) ! read obs to frd(:)
             else
                frdf = 0.
             endif
             call mpscgp(frdf,frd,1)
             if (lobsav(iobs)) then
                fm(:) = fobsav(:,iobsav(iobs)) !/ countav ! compute requested model average
             else
                call code2field(obscode(iobs),obslev(iobs),fm(:)) ! copy requested model field  to fm(:)
             end if
!!$TAF STORE fm,frd= obsn, REC=nobs*(nstep-nstep0+1)/ncost
!$TAF STORE fm,frd= obsn
             ch = 0.
             nval = 0.
             do jlat=1,NLPP
                lval(:)  = (abs(frd((jlat-1)*NLON+1:jlat*NLON)/miss-1.).gt.1.e-3 .and. &
                     fobserr((jlat-1)*NLON+1:jlat*NLON,iobs) .gt. 0.) ! missing value mask on current parallel - data
                nval     = nval + count(lval(:)) ! number of valid values on current parallel
                where (lval)
                   fh(:) = fm((jlat-1)*NLON+1:jlat*NLON) - frd((jlat-1)*NLON+1:jlat*NLON)
                   fh(:) = fh(:) ** 2 * fobserr((jlat-1)*NLON+1:jlat*NLON,iobs)
                end where
                ch = ch + gaussw(jlat) * sum(fh(:),mask=lval)
             end do
!             call mpsumi(nval,1)
!             call mpsumr(ch,1)
             if(mypid==NROOT) then
!!$TAF STORE nval,ch = obsn, REC=nobs*(nstep-nstep0+1)
              write (nud,'(1x,a,i3,a,i4,a,i6,a,e20.12,a,i5)') 'Cost function contribution of obs ',iobs, &
                  ' code ',obscode(iobs),' level ',obslev(iobs),' = ', &
                  ch(1) * 0.5 * fobserrfac(iobs) * NLAT /gausswsum / nval,' missing:',nlon*nlat-nval
              cost = cost + ch(1) * 0.5 * fobserrfac(iobs) * NLAT / gausswsum / nval ! mean square deviation per grid cell
             endif
!             call mpbcr(cost)
          end do
          if(mypid==NROOT) then
           if(.not.cost.lt.huge(cost).or..not.(-cost).lt.huge(cost)) then
             print*,"detecting NaN in cost function!"
             lunstable = .true.
           end if
          endif
          call resetav
        end subroutine addcost

        subroutine postobs
          use srvio, only : postsrvio
          if (allocated(fobsav)) deallocate(fobsav)
          call postsrvio
        end subroutine postobs

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                          Nuddging part 
! ----------------------------------------------------------------------------------------------
      subroutine nudgini
      use pumamod
      implicit none
      !  for data type  
      integer,parameter :: ioprec=4 ! prcision of data files
      integer(kind=4) :: ircw1,ircw2,ircw3,ircw4 ! fortran record control words
      integer(kind=4)   :: ihead(8) ,jhead(8)         ! SERVICE header
      integer :: recl, junit=50
      character(len=50),parameter :: nfile = 'observations_nudg.srv'
      real(kind=ioprec):: fio1(NLON*NLAT), fio2(NLON*NLAT) ! single data field in single precision
         inquire(iolength=recl) ircw1,ihead(:),ircw2,ircw3,fio1(:),ircw4
         ! open (junit,file=nfile, FORM='UNFORMATTED',ACCESS='DIRECT',RECL=recl,STATUS='old',action='read',convert='SMALL_ENDIAN')
         open (unit=junit, file=nfile, form='UNFORMATTED', access='DIRECT', recl=recl, status='OLD', action='READ', convert='LITTLE_ENDIAN')
!         read(junit, rec=1) ircw1, ihead, ircw2, ircw3, fio1, ircw4
!         print *, 'Header values (ihead): ', ihead
!         print *, 'RECL =', recl
      end subroutine nudgini

      subroutine getnudg(kstep,kcode,frd)
         use pumamod
!         use srvio , only : readsrv,code2field
         implicit none
         integer, intent(in) :: kcode ,kstep !( 138 155 130 133)
         real, dimension(NLON*NLAT,NLEV), intent(out) :: frd
         real :: itr1,itr2  ! coefficient for temporal interploation
      
         integer :: nmin,nhour,nday,nmonth,nyear,kyday,nrec
         integer :: junit=50, klev
      !  for data type  
         integer,parameter :: ioprec=4 ! prcision of data files
         integer(kind=4) :: ircw1,ircw2,ircw3,ircw4 ! fortran record control words
         integer(kind=4)   :: ihead(8) ,jhead(8)         ! SERVICE header
         real(kind=ioprec):: fio1(NLON*NLAT), fio2(NLON*NLAT) ! single data field in single precision
         integer :: recl
         character(len=50) :: ncname
   
         if(mypid/=NROOT) then
            print*,'getnudg can only be in ROOT process!'
            stop
         endif
!         print*, "KSTEP in getnudge", kstep
         call ntomin(kstep+1,nmin,nhour,nday,nmonth,nyear)  
         call mmdd2yday(kyday,nyear,nmonth,nday)
         
         ! 26.07.2023, yulia: read nudging fields for multiple years from the same file
         if(nyear>=1980)then
            kyday=kyday+360*(nyear-1980)
         endif
         ! 26.07.2023: edit end
         
         nrec=(kyday-1)*24+((nhour)*60+nmin)/60 +1 
         itr2=(mod(real(nhour),1.)+real(nmin)/60.)/1.
         itr1=1-itr2  ! for nrec record


!         print*, "NRECORD",nrec

         ! print*,'simulating time: ', nday,nhour,nmin,nrec, itr1,itr2
         !  print*,'data is in: ',nrec
         
            jhead(1) = kcode     
            jhead(3) = nday + 100 * nmonth + 10000 * nyear
      !      ihead(4) = nmin + 100 * nhour
            jhead(4) = 100 * (nmin + 100 * nhour)  ! alternative time format including seconds: hhmmss
            jhead(5) = NLON
            jhead(6) = NLAT
            jhead(7) = 1
            jhead(8) = n_days_per_year


         if(kcode==138)then   ! vorticity
         do klev=1,10
            jhead(2) = klev
            read (junit,rec=(nrec-1)*41+klev) ircw1,ihead(:),ircw2,ircw3,fio1(:),ircw4
            if ( any(ihead((/1,2/)).ne. jhead((/1,2/))).or.( (ihead(4).gt.jhead(4)).and.(ihead(3).gt.jhead(3)))) then
            
            print'(1x,a,a)',"header mismatch while reading from nudg file"
            print*,'record =',(nrec-1)*40+klev
            print'(a,6i9)','wanted: ',jhead(1:6)
            print'(a,6i9)','found:  ',ihead(1:6)
            stop
            end if
      !    frd(:,klev)=fio1(:)
            read (junit,rec=nrec*41+klev) ircw1,ihead(:),ircw2,ircw3,fio2(:),ircw4
            frd(:,klev)=itr1*fio1(:)+itr2*fio2(:) 
         enddo
         elseif(kcode==155)then  ! divergence
         do klev=1,10
            jhead(2) = klev
            read (junit,rec=(nrec-1)*41+klev+10) ircw1,ihead(:),ircw2,ircw3,fio1(:),ircw4
            if (any(ihead((/1,2/)).ne. jhead((/1,2/))).or.( (ihead(4).gt.jhead(4)).and.(ihead(3).gt.jhead(3)))) then
            
            print'(1x,a,a)',"header mismatch while reading from nudg file"
            print*,'record =',(nrec-1)*40+klev+10
            print'(a,6i9)','wanted: ',jhead(1:6)
            print'(a,6i9)','found:  ',ihead(1:6)
            stop
            end if
      !       frd(:,klev)=fio1(:)
            read (junit,rec=nrec*41+klev+10) ircw1,ihead(:),ircw2,ircw3,fio2(:),ircw4
            frd(:,klev)=itr1*fio1(:)+itr2*fio2(:) 
         enddo

         elseif(kcode==130)then              
            do klev=1,10
            jhead(2) = klev
            read (junit,rec=(nrec-1)*41+klev+20) ircw1,ihead(:),ircw2,ircw3,fio1(:),ircw4
            if (any(ihead((/1,2/)).ne. jhead((/1,2/))).or.( (ihead(4).gt.jhead(4)).and.(ihead(3).gt.jhead(3)))) then
            
            print'(1x,a,a)',"header mismatch while reading from nudg file"
            print*,'record =',(nrec-1)*41+klev+20
            print'(a,6i9)','wanted: ',jhead(1:6)
            print'(a,6i9)','found:  ',ihead(1:6)
            stop
            end if
      !     frd(:,klev)=fio1(:)
            read (junit,rec=nrec*41+klev+20) ircw1,ihead(:),ircw2,ircw3,fio2(:),ircw4
            frd(:,klev)=itr1*fio1(:)+itr2*fio2(:) 
         enddo
         elseif(kcode==133)then
            do klev=1,10
            jhead(2) = klev
            read (junit,rec=(nrec-1)*41+klev+30) ircw1,ihead(:),ircw2,ircw3,fio1(:),ircw4
            if (any(ihead((/1,2/)).ne. jhead((/1,2/))).or.( (ihead(4).gt.jhead(4)).and.(ihead(3).gt.jhead(3)))) then
            
            print'(1x,a,a)',"header mismatch while reading from nudg file"
            print*,'record =',(nrec-1)*41+klev+30
            print'(a,6i9)','wanted: ',jhead(1:6)
            print'(a,6i9)','found:  ',ihead(1:6)
            stop
            endif
      !      frd(:,klev)=fio1(:)
            read (junit,rec=nrec*41+klev+30) ircw1,ihead(:),ircw2,ircw3,fio2(:),ircw4
            frd(:,klev)=itr1*fio1(:)+itr2*fio2(:) 
         enddo
         else
            print*,'kcode not supported: ', kcode
            stop
         endif 
      !endif
      end subroutine getnudg
      
      subroutine getnudg2d(kstep,kcode,frd)
         use pumamod
!         use srvio , only : readsrv,code2field

         implicit none
         integer, intent(in) :: kcode ,kstep !( 134 )
         real, dimension(NLON*NLAT), intent(out) :: frd
         real :: itr1,itr2  ! coefficient for temporal interploation
         integer, external :: get_free_unit 
         integer :: nmin,nhour,nday,nmonth,nyear,kyday,nrec
         integer ::  klev
      !  for data type  
         integer,parameter :: ioprec=4 ! prcision of data files
         integer(kind=4) :: ircw1,ircw2,ircw3,ircw4 ! fortran record control words
         integer(kind=4)   :: ihead(8) ,jhead(8)         ! SERVICE header
         real(kind=ioprec):: fio1(NLON*NLAT), fio2(NLON*NLAT) ! single data field in single precision
         integer :: recl, junit=50

      ! observations start from 00:00-31-12-2000 to 00-00-02-01-2002
      ! just for coding convinence 
      ! open file
      if(mypid/=NROOT) then
         print*,'getnudg can only be in ROOT process!'
         stop
      endif

      ! if (makenudge==1) then
      ! get hhddmmyy for former step
         call ntomin(kstep+1,nmin,nhour,nday,nmonth,nyear)  
         call mmdd2yday(kyday,nyear,nmonth,nday)
         
         ! 26.07.2023, yulia: read nudging fields for multiple years from the same file
         if(nyear>=1980)then
            kyday=kyday+360*(nyear-1980)
         endif
         ! 26.07.2023: edit end
         nrec=(kyday-1)*24+((nhour)*60+nmin)/60 +1 

         itr2=(mod(real(nhour),1.)+real(nmin)/60.)/1.
         itr1=1-itr2  ! for nrec record

         
            jhead(1) = kcode     
            jhead(3) = nday + 100 * nmonth + 10000 * nyear
      !      ihead(4) = nmin + 100 * nhour
            jhead(4) = 100 * (nmin + 100 * nhour)  ! alternative time format including seconds: hhmmss
            jhead(5) = NLON
            jhead(6) = NLAT
            jhead(7) = 1
            jhead(8) = n_days_per_year


         if(kcode==134)then   ! surface pressure 
            read (junit,rec=(nrec-1)*41+41) ircw1,ihead(:),ircw2,ircw3,fio1(:),ircw4
            if ( any(ihead((/1/)).ne. jhead((/1/))).or.( (ihead(4).gt.jhead(4)).and.(ihead(3).gt.jhead(3)))) then
            
            print'(1x,a,a)',"header mismatch while reading from nudg file"
            print*,'record =',(nrec-1)*41+klev
            print'(a,6i9)','wanted: ',jhead(1:6)
            print'(a,6i9)','found:  ',ihead(1:6)
            stop
            end if
      !   frd(:)=fio1(:)
            read (junit,rec=nrec*41+41) ircw1,ihead(:),ircw2,ircw3,fio2(:),ircw4
            !print *, "nudding is working"
            frd(:)=itr1*fio1(:)+itr2*fio2(:) 

      endif
      !endif
      end subroutine getnudg2d


      end module observation


