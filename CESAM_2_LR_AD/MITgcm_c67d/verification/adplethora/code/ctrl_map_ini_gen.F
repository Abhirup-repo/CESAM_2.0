#include "CTRL_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine ctrl_map_ini_gen3D(xxFileCur, wFileCur, xxDummyCur,
     & boundsVec, paramFld3d, maskFld3d, paramSmooth, mythid )

      use mo_mapping_plethora, only : cont2mod_oce_xyz
      use mo_mapping, only : ix

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#ifdef ALLOW_CTRL
# include "CTRL_SIZE.h"
# include "ctrl.h"
# include "CTRL_GENARR.h"
# include "ctrl_dummy.h"
# include "optim.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
# include "tamc_keys.h"
#endif /* ALLOW_AUTODIFF */

c     == routine arguments ==

      integer mythid
      character*(*) wFileCur,xxFileCur
      _RL boundsVec(5),tmpMax,xxDummyCur

      _RL wFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RL xxFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RL paramFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RS maskFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      integer paramSmooth

c     == local variables ==

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer il

      logical doglobalread
      logical ladinit

      character*( 80) fnamegeneric
      _RS dummyRS(1)

c     == external ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

#ifdef ALLOW_AUTODIFF_TAMC
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = 0
          ikey = (act3 + 1) + act4*max3
#endif /* ALLOW_AUTODIFF_TAMC */

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
c--   only do interior, and exchange outside
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

      doglobalread = .false.
      ladinit      = .false.

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO k = 1,nr
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx
           xxFld3d(i,j,k,bi,bj)=0. _d 0
           wFld3d(i,j,k,bi,bj)=0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL MDS_READ_FIELD(wFileCur,ctrlprec,.FALSE.,'RL',
     &                    nR,1,nR,wFld3d,dummyRS,1,mythid)
      _EXCH_XYZ_RL( wFld3d, mythid )

      il=ilnblnk( xxFileCur )
      write(fnamegeneric(1:80),'(2a,i10.10)')
     &     xxFileCur(1:il),'.',optimcycle
c #ifdef ALLOW_AUTODIFF
c       call active_read_xyz( fnamegeneric, xxFld3d, 1,
c      & doglobalread, ladinit, optimcycle, mythid, xxDummyCur )
c #else
c       CALL READ_REC_XYZ_RL( fnamegeneric, xxFld3d, 1, 1, myThid )
c #endif
C
      if (wFileCur=="wthetaLev") then
         ix=28672 
      elseif (wFileCur=="wsaltLev") then
         ix=28672 + 44*90*15
      else
         ix=0
      endif
c      print*,"The value of wthetaLev", sum(wFld3d)
C???ION???    forward run: no perturbation in ocean
      call cont2mod_oce_xyz(xxFld3d,mythid) ! FastOpt
C???ION???
C
      IF ( .NOT.ctrlSmoothCorrel3D ) THEN

c avoid xx larger than boundsVec(5) X uncertainty
c       if ( boundsVec(5).GT.0. ) then
c       do bj = jtlo,jthi
c         do bi = itlo,ithi
c           do k = 1,nr
c             do j = jmin,jmax
c               do i = imin,imax
c       if ( (maskFld3d(i,j,k,bi,bj).NE.0.).AND.
c      & (wFld3d(i,j,k,bi,bj).GT.0.) ) then
c        tmpMax=boundsVec(5)/sqrt(wFld3d(i,j,k,bi,bj))
c       if ( abs(xxFld3d(i,j,k,bi,bj)).GT.tmpMax ) then
c        xxFld3d(i,j,k,bi,bj)=sign(tmpMax,xxFld3d(i,j,k,bi,bj))
c       else
c        xxFld3d(i,j,k,bi,bj)=xxFld3d(i,j,k,bi,bj)
c       endif
c       endif
c               enddo
c             enddo
c           enddo
c        enddo
c       enddo
c       endif
      if ( boundsVec(5).GT.0. ) then
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do k = 1,nr
            do j = jmin,jmax
              do i = imin,imax
      if ( (maskFld3d(i,j,k,bi,bj).NE.0.).AND.
     & (wFld3d(i,j,k,bi,bj).GT.0.) ) then
          xxFld3d(i,j,k,bi,bj)=xxFld3d(i,j,k,bi,bj)
     &                /sqrt( wFld3d(i,j,k,bi,bj))
      else
          xxFld3d(i,j,k,bi,bj)=0.
      endif
              enddo
            enddo
          enddo
       enddo
      enddo
      endif

#ifdef ALLOW_SMOOTH
# ifdef ALLOW_SMOOTH_CTRL3D
      if (useSMOOTH) call smooth3D(xxFld3d,paramSmooth,mythid)
      write(fnamegeneric(1:80),'(2a,i10.10)')
     & xxFileCur(1:il),'.smooth.',optimcycle
      CALL MDS_WRITE_FIELD(fnamegeneric,ctrlprec,.FALSE.,.FALSE.,
     & 'RL',nr,1,nr,xxFld3d,dummyRS,1,optimcycle,mythid)
# endif
#endif

      do bj = jtlo,jthi
        do bi = itlo,ithi
          do k = 1,nr
            do j = jmin,jmax
              do i = imin,imax
      paramFld3d(i,j,k,bi,bj) = paramFld3d(i,j,k,bi,bj)
     & + xxFld3d(i,j,k,bi,bj)
              enddo
            enddo
          enddo
       enddo
      enddo

      ELSE !IF ( .NOT.ctrlSmoothCorrel3D ) THEN

#ifdef ALLOW_SMOOTH
c apply Weaver And Courtier correlation operator
      if ( paramSmooth.NE.0 .AND. useSMOOTH ) then
       call smooth_correl3d(xxFld3d,paramSmooth,mythid)
      endif
#endif

      do bj = jtlo,jthi
        do bi = itlo,ithi
          do k = 1,nr
            do j = jmin,jmax
              do i = imin,imax
c scale param adjustment
      if ( (maskFld3d(i,j,k,bi,bj).NE.0.)
     & .AND. (wFld3d(i,j,k,bi,bj).GT.0.) ) then
      xxFld3d(i,j,k,bi,bj)=xxFld3d(i,j,k,bi,bj)
     & /sqrt( wFld3d(i,j,k,bi,bj) )
      else
      xxFld3d(i,j,k,bi,bj)=0.
      endif
      paramFld3d(i,j,k,bi,bj) = paramFld3d(i,j,k,bi,bj)
     & + xxFld3d(i,j,k,bi,bj)
              enddo
            enddo
          enddo
       enddo
      enddo

      ENDIF !IF ( .NOT.ctrlSmoothCorrel3D ) THEN

      end subroutine ctrl_map_ini_gen3D

      subroutine ctrl_map_ini_gen3D_log(xxFileCur, wFileCur, xxDummyCur,
     & boundsVec, paramFld3d, maskFld3d, paramSmooth, mythid )


      use mo_mapping, only : ix
      use mo_mapping_plethora, only : cont2mod_oce_xyz

      implicit none

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#ifdef ALLOW_CTRL
# include "CTRL_SIZE.h"
# include "ctrl.h"
# include "CTRL_GENARR.h"
# include "ctrl_dummy.h"
# include "optim.h"
#endif
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
# include "tamc_keys.h"
#endif /* ALLOW_AUTODIFF */

c     == routine arguments ==

      integer mythid
      character*(*) wFileCur,xxFileCur
      _RL boundsVec(5),tmpMax,xxDummyCur

      _RL wFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RL xxFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RL xxFld3d1(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RL paramFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      _RL maskFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      integer paramSmooth

c     == local variables ==

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer il

      logical doglobalread
      logical ladinit

      character*( 80) fnamegeneric

      _RS dummyRS(1)
c     == external ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

#ifdef ALLOW_AUTODIFF_TAMC
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = 0
          ikey = (act3 + 1) + act4*max3
#endif /* ALLOW_AUTODIFF_TAMC */

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
c--   only do interior, and exchange outside
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

      doglobalread = .false.
      ladinit      = .false.


c      call mdsreadfield(wFileCur,ctrlprec,'RL',nR,wFld3d,1,mythid)
c      _EXCH_XYZ_RL( wFld3d, mythid )

      CALL MDS_READ_FIELD(wFileCur,ctrlprec,.FALSE.,'RL',
     &                    nR,1,nR,wFld3d,dummyRS,1,mythid)
      _EXCH_XYZ_RL( wFld3d, mythid )

      il=ilnblnk( xxFileCur )
      write(fnamegeneric(1:80),'(2a,i10.10)')
     &     xxFileCur(1:il),'.',optimcycle
c      call active_read_xyz( fnamegeneric, xxFld3d, 1,
c     & doglobalread, ladinit, optimcycle, mythid, xxDummyCur )
c      print*,'initial 3d is ', sum(xxFld3d)
      if (wFileCur=="wkapgmFld")  then
         ix=28672  + 44*90*15*2
      elseif (wFileCur=="wkaprediFld") then
         ix=28672 + 44*90*15*3
      else
         ix=0
      endif
         
      call cont2mod_oce_xyz(xxFld3d,myThid) ! FastOpt

      do bj = jtlo,jthi
        do bi = itlo,ithi
          do k = 1,nr
            do j = jmin,jmax
              do i = imin,imax
      if ( (maskFld3d(i,j,k,bi,bj).NE.0.).AND.
     & (wFld3d(i,j,k,bi,bj).GT.0.).and.
     & (paramFld3d(i,j,k,bi,bj).gt.0.) ) then
        paramFld3d(i,j,k,bi,bj)=(paramFld3d(i,j,k,bi,bj))
     &  *exp(xxFld3d(i,j,k,bi,bj)*log(1.+1./(paramFld3d(i,j,k,bi,bj)
     &  *sqrt( wFld3d(i,j,k,bi,bj)))))
c      else
c          paramFld3d(i,j,k,bi,bj)=paramFld3d(i,j,k,bi,bj)
      endif 
              enddo
            enddo
          enddo
       enddo
      enddo



c avoid param out of [boundsVec(1) boundsVec(4)]
      CALL CTRL_BOUND_3D(paramFld3d,maskFld3d,boundsVec,myThid)

      IF ( ctrlSmoothCorrel3D ) THEN
        write(fnamegeneric(1:80),'(2a,i10.10)')
     &    xxFileCur(1:il),'.effective.',optimcycle
        CALL MDS_WRITE_FIELD(fnamegeneric,ctrlprec,.FALSE.,.FALSE.,
     &   'RL',nr,1,nr,paramFld3d,dummyRS,1,optimcycle,mythid)
      ENDIF

      RETURN
      END
      
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      subroutine ctrl_map_ini_gen2D(xxFileCur, wFileCur, xxDummyCur,
     & boundsVec, paramFld2d, maskFld3d, paramSmooth, mythid )

      use mo_mapping, only : ix
      use mo_mapping_plethora, only : cont2mod_oce_xy

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#ifdef ALLOW_CTRL
# include "CTRL_SIZE.h"
# include "ctrl.h"
# include "ctrl_dummy.h"
# include "optim.h"
#endif
#ifdef ALLOW_AUTODIFF
#include "tamc.h"
#include "tamc_keys.h"
#endif /* ALLOW_AUTODIFF */

c     == routine arguments ==

      integer mythid
      character*(*) wFileCur,xxFileCur
      _RL boundsVec(5),tmpMax,xxDummyCur

      _RL wFld2d(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)
      _RL xxFld2d(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)
      _RL paramFld2d(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)
      _RS maskFld3d(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
      integer paramSmooth

c     == local variables ==

      integer bi,bj
      integer i,j,k
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer il

      logical doglobalread
      logical ladinit

      character*( 80) fnamegeneric
      _RS dummyRS(1)

c     == external ==

      integer  ilnblnk
      external ilnblnk

c     == end of interface ==

#ifdef ALLOW_AUTODIFF_TAMC
          act3 = myThid - 1
          max3 = nTx*nTy
          act4 = 0
          ikey = (act3 + 1) + act4*max3
#endif /* ALLOW_AUTODIFF_TAMC */

      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
c--   only do interior, and exchange outside
      jmin = 1
      jmax = sny
      imin = 1
      imax = snx

      doglobalread = .false.
      ladinit      = .false.

      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
           xxFld2d(i,j,bi,bj)=0. _d 0
           wFld2d(i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL MDS_READ_FIELD(wFileCur,ctrlprec,.FALSE.,'RL',
     &                    1,1,1,wFld2d,dummyRS,1,mythid)
      _EXCH_XY_RL( wFld2d, mythid )

      il=ilnblnk( xxFileCur )
      write(fnamegeneric(1:80),'(2a,i10.10)')
     &     xxFileCur(1:il),'.',optimcycle
c #ifdef ALLOW_AUTODIFF
c       call active_read_xy( fnamegeneric, xxFld2d, 1,
c      & doglobalread, ladinit, optimcycle, mythid, xxDummyCur )
c #else
c       CALL READ_REC_XY_RL( fnamegeneric, xxFld2d, 1, 1, myThid )
c #endif
      ix=0
      call cont2mod_oce_xy(xxFld2d,myThid) ! FastOpt

      IF ( .NOT.ctrlSmoothCorrel2D ) THEN

c avoid xx larger than boundsVec(5) X uncertainty
      if ( boundsVec(5).GT.0.) then
      do bj = jtlo,jthi
        do bi = itlo,ithi
            do j = jmin,jmax
              do i = imin,imax
      if ( (maskFld3d(i,j,1,bi,bj).NE.0.).AND.
     & (wFld2d(i,j,bi,bj).GT.0.) ) then
       tmpMax=boundsVec(5)/sqrt(wFld2d(i,j,bi,bj))
      if ( abs(xxFld2d(i,j,bi,bj)).GT.tmpMax ) then
       xxFld2d(i,j,bi,bj)=sign(tmpMax,xxFld2d(i,j,bi,bj))
      else
       xxFld2d(i,j,bi,bj)=xxFld2d(i,j,bi,bj)
      endif
      endif
              enddo
            enddo
       enddo
      enddo
      endif

#ifdef ALLOW_SMOOTH
# ifdef ALLOW_SMOOTH_CTRL2D
      if (useSMOOTH) call smooth2d(xxFld2d,maskFld3d,paramSmooth,mythid)
      write(fnamegeneric(1:80),'(2a,i10.10)')
     & xxFileCur(1:il),'.smooth.',optimcycle
      CALL MDS_WRITE_FIELD(fnamegeneric,ctrlprec,.FALSE.,.FALSE.,
     & 'RL',1,1,1,xxFld2d,dummyRS,1,optimcycle,mythid)
# endif
#endif

      do bj = jtlo,jthi
        do bi = itlo,ithi
            do j = jmin,jmax
              do i = imin,imax
      paramFld2d(i,j,bi,bj) = paramFld2d(i,j,bi,bj)
     & + xxFld2d(i,j,bi,bj)
              enddo
            enddo
       enddo
      enddo

      ELSE !IF ( .NOT.ctrlSmoothCorrel2D ) THEN

#ifdef ALLOW_SMOOTH
c apply Weaver And Courtier correlation operator
      if ( paramSmooth.NE.0 .AND. useSMOOTH ) then
       call smooth_correl2d(xxFld2d,maskFld3d,paramSmooth,mythid)
      endif
#endif

      do bj = jtlo,jthi
        do bi = itlo,ithi
            do j = jmin,jmax
              do i = imin,imax
c scale param adjustment
      if ( (maskFld3d(i,j,1,bi,bj).NE.0.)
     & .AND. (wFld2d(i,j,bi,bj).GT.0.) ) then
      xxFld2d(i,j,bi,bj)=xxFld2d(i,j,bi,bj)
     & /sqrt( wFld2d(i,j,bi,bj) )
      else
      xxFld2d(i,j,bi,bj)=0.
      endif
      paramFld2d(i,j,bi,bj) = paramFld2d(i,j,bi,bj)
     & + xxFld2d(i,j,bi,bj)
              enddo
            enddo
       enddo
      enddo

      ENDIF !IF ( .NOT.ctrlSmoothCorrel2D ) THEN

      CALL CTRL_BOUND_2D(paramFld2d,maskFld3d,boundsVec,myThid)

      IF ( ctrlSmoothCorrel2D ) THEN
        write(fnamegeneric(1:80),'(2a,i10.10)')
     &   xxFileCur(1:il),'.effective.',optimcycle
        CALL MDS_WRITE_FIELD(fnamegeneric,ctrlprec,.FALSE.,.FALSE.,
     &   'RL',1,1,1,paramFld2d,dummyRS,1,optimcycle,mythid)
      endif

      RETURN
      END
