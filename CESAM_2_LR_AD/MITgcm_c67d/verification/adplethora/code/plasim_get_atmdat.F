C $Header: /u/gcmpack/MITgcm/pkg/exf/exf_mapfields.F,v -1.23 2009/04/28 18:15:33 jmc Exp $
C $Name: checkpoint62y $

#include "CPP_OPTIONS.h"
#include "CPP_EEOPTIONS.h"
c#include "ECCO_CPPOPTIONS.h"

c$taf SUBROUTINE MPI_CART_COORDS   INPUT = 1,2,3
c$taf SUBROUTINE MPI_CART_COORDS  OUTPUT = 4
c$taf SUBROUTINE MPI_CART_COORDS  ACTIVE = 
c$taf SUBROUTINE MPI_CART_COORDS  DEPEND = 

      subroutine plasim_get_atmdat(mythid,kstepo)

c     ==================================================================
c     SUBROUTINE plasim_get_atmdat
c     ==================================================================
c
c     o run the planet simulator
c     o map forcing fields from planet simulator onto ocean arrays.
c
c     started: Larry 1-Jul-2011
c
c     ==================================================================
c     SUBROUTINE plasim_get_atmdat
c     ==================================================================
c     

      use cplmod
      use cplmod, only : ndatim_cpl=>ndatim
      use pumamod, only: acc,aevap,alhfl,
     &     aprc,aprl,aprs,aqvi,ashfl,asmelt,
     &     asndch,assol,assolu,asthr,asthru,ats0,atsa,atsol,atsolu,
     &     atthr,bm1,dcc,delt,delt2,dflux,dforest,
     &     dftu,dfd,dfu,dicec,dlhfl,dlwfl,dlhdt,dshdt,gq,
     &     dp,dprc,dprl,dprs,dq,dql,dqt,dqvi,drhs,drunoff,
     &     dshfl,dsmelt,dsndch,dsnow,dswfl,
     &     dt,dtaux,dtauy,dtsa,du,dv,dwatc,dwmax,
     &     dz0,naccuout,ncost,ndivdamp,nesp,nhor,nkits,nlem,
     &     nlep,nlev,nstep,nspp,nveg,
     &     sd,sdm,sdp,sp,spm,spp,sq,sqm,sqp,st,stm,stp,sz,
     &     szm,szp,ndatim_puma=>ndatim
      use icemod, only: naccuo,xpme,xswfl,xtaux,xtauy,
     &     cheat_ice=>cheat,xmld,xsndch,xflxice2,xheat,xicec,xicecc,
     &     xiced,xoflux,xprs,xqmelt,xroff,xsnow,xsst,xtfreeze,xts,
     &     nstep_ice=>nstep
      use seamod, only : cheata,clhdta,clhfla,clwfla,cpmea,croffa,
     &     cshdta,cshfla,cswfla,ctauxa,ctauya,cust3a,
     &     cprsa,cicec,ciced,cmld,csmelt,csndch,csnow,csst,cts,dsst,
     &     naccua,dts_sea=>dts,dqs_sea=>dqs
      use landmod, only: devap,dalb,
     &     nlsoil,driver,dsnowz,dsoilt,dtclim,dtsm,dwater,dwclim,
     &     dts_land=>dts,dqs_land=>dqs
      use vegmod, only : valb_max,vz0_min,zvalb,zvrhs,dvsoil,
     &     aresh,anogrow,anpp,alitter,agpp,agppl,agppw,zvz0,
     &     dcsoil,dcveg,zgpp,zgppl,zgppw,znpp,zres,zlitter,zveg
      use radmod, only : gmu0,dqo3,gdist2
      use rainmod, only : icclev,icctot,dprcl,dprll,dprscl,dprsll
      use oceanmod, only : yifluxr,ysst,nlev_oce,nstep_ocean=>nstep,
     &     ndatim_ocean=>ndatim
      use observation, only : genav,addcost,writeobs,countav,nobs,cost,
     &     fobsav,nobsav
      use mo_mapping, only : ix
      use mo_mapping_plethora, only : cont2mod_oce_xy
      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "GRID.h"
c
#ifdef ALLOW_AUTODIFF
# include "ctrl.h"
# include "ctrl_dummy.h"
# include "ecco_cost.h"        
#endif

#ifdef ALLOW_AUTODIFF_TAMC
# include "tamc.h"
# include "tamc_keys.h"
#endif
# if (defined (ALLOW_USE_MPI) || defined (ALWAYS_USE_MPI))
#  include "EESUPPORT.h"
# else
      INTEGER nProcs
      PARAMETER (nProcs=1)
# endif
c     == routine arguments ==

c     mythid - thread number for this instance of the routine.

      integer mythid
      integer kstepo

c     == local variables ==

      integer bi,bj,i,j
      integer itlo,ithi
      integer jtlo,jthi
      integer jmin,jmax
      integer imin,imax
      integer ipx, ipy
      integer iglobe, jglobe
      integer iproc
      integer mpirc
      integer mpicrd(2)
c
      integer ktoplev
      real(kind=8) cheat(1:Nx,1:Ny)
      real(kind=8) cpme(1:Nx,1:Ny)
      real(kind=8) croff(1:Nx,1:Ny)
      real(kind=8) ctaux(1:Nx,1:Ny)
      real(kind=8) ctauy(1:Nx,1:Ny)
      real(kind=8) cswr(1:Nx,1:Ny)
      real(kind=8) cps(1:Nx,1:Ny)
      real(kind=8) cil(1:Nx,1:Ny)
      real(kind=8) xxFld2d(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)
      INTEGER nflxdim
      PARAMETER (nflxdim=Nx*Ny*360*39)
      integer ix_atm

c     == end of interface ==
c--   Initialise local variables.
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

      jmin = 1
      jmax = sny
      imin = 1
      imax = snx
c
c     A) RUN PLASIM
c

c???ION??? 23.06.2015  coupling time step: atm 10*48min = ocn 480 min = 8h
!SiS substituted 10 -> kcplstep in pumamod      call plasim_step(10,kstepo,myProcId)

      call plasim_step(10,kstepo,myProcId)

C$TAF STORE cost,countav,gdist2,naccua,naccuo,ndivdamp,nhfla,npmea, 
C$TAF&      nroffa,nssto,nstep,nstep,nswra,ntaua,zveg,ahfla,apmea,  
C$TAF&      aroffa,aswra,atauxa,atauya,cheat_ice,cheata,clhdta,clhfla,
C$TAF&      clwfla,cpmea,cprsa,croffa,cshdta,cshfla,cswfla,ctauxa,  
C$TAF&      ctauya,cust3a,dalb,dcc,dcsoil,dcveg,dfd,dflux,dforest,  
C$TAF&      dicec,dlhdt,dlwfl,dp,dprcl,dprll,dprscl,dprsll,dq,dql,  
C$TAF&      dqo3,dqt,drhs,driver,drunoff,dshdt,dsnow,dsnowz,  
C$TAF&      dsoilt,dsst,dswfl,dt,dtaux,dtauy,dtclim,dtsa, 
C$TAF&      dtsm,du,dv,dvsoil,dwatc,dwater,dwclim,dwmax,dz0,fobsav,
C$TAF&      gmu0,gq,icclev,icctot,ndatim,ossta,sd,sdm,sdp,sp,spm,  
C$TAF&      spp,sq,sqm,sqp,st,stm,stp,sz,szm,szp,xheat,xicec, 
C$TAF&      xicecc,xiced,xmld,xoflux,xpme,xprs,xqmelt,xroff,xsnow,
C$TAF&      xsst,xswfl,xtaux,xtauy,xtfreeze,xts,yifluxr,ysst,
C$TAF&      dqs_land,dqs_sea, dts_land,dts_sea
C$TAF&      = mastercp_lev1, REC = kstepo
c
c     B) GET FORCING FIELDS FROM COUPLER AND OVERWRITE
c        CLIMATOLOGICAL FORCING
c
!!      print*,'IN plasim_get_atmdat.F:'
!      print*,'VORHER: '
!      print*,'QNET: ',maxval(qnet),minval(qnet)
!      print*,'QSWR: ',maxval(qsw),minval(qsw)
!      print*,'PME: ',maxval(EmPmR),minval(EmPmR)
!      print*,'USTR:',maxval(fu),minval(fu)
!      print*,'VSTR:',maxval(fv),minval(fv)
!      print*,'PS:',maxval(pload),minval(pload)
!!      print*,'IL:',maxval(sIceLoad),minval(sIceLoad)
c
#if (defined (ALLOW_USE_MPI) || defined (ALWAYS_USE_MPI))
c--     get coordinates of processor myProcId
        call MPI_Cart_coords(
     &          mpiComm, myProcId, 2, mpicrd
     &        , mpirc
     &        )
        ipx = 1 + mpicrd(1)
        ipy = 1 + mpicrd(2)
#else
        ipx = 1
        ipy = 1
#endif
c$TAF loop = parallel
      DO bj = jtlo,jthi
c$TAF loop = parallel
        DO bi = itlo,ithi
         call cgethflo(cheat,Nx,Ny,bi,bj,myProcId)
         call cgetpmeo(cpme,Nx,Ny,bi,bj,myProcId)
         call cgetroffo(croff,Nx,Ny,bi,bj,myProcId)
         call cgettauo(ctaux,ctauy,Nx,Ny,bi,bj,myProcId)
#ifdef SHORTWAVE_HEATING
         call cgetswro(cswr,Nx,Ny,bi,bj,myProcId)
#endif
#ifdef ATMOSPHERIC_LOADING
         call cgetpso(cps,Nx,Ny,bi,bj,myProcId)
         call cgetimo(cil,Nx,Ny,bi,bj,myProcId)
#endif

            do j = jmin,jmax
              do i = imin,imax
                jglobe = j + sNy*(bj-1) + sNy*nSy*(ipy-1)
                iglobe = i + sNx*(bi-1) + sNx*nSx*(ipx-1)
                Qnet(i,j,bi,bj)=cheat(iglobe,jglobe) 
                EmPmR(i,j,bi,bj)=(cpme(iglobe,jglobe)
     &                           +croff(iglobe,jglobe))*rhoConstFresh
                fu(i,j,bi,bj)=ctaux(iglobe,jglobe) 
                fv(i,j,bi,bj)=ctauy(iglobe,jglobe) 
#ifdef SHORTWAVE_HEATING
                Qsw(i,j,bi,bj)=cswr(iglobe,jglobe) 
#endif
#ifdef ATMOSPHERIC_LOADING
                pload(i,j,bi,bj)=ccps(iglobe,jglobe) 
                sIceLoad(i,j,bi,bj)=cil(iglobe,jglobe) 
#endif

              ENDDO
            ENDDO

        ENDDO
      ENDDO


#ifdef ALLOW_HFLUX_CONTROL
C     print*,'in plasim_get_atmdat check-point 1'
C only root cpu share the data with the other cpu.
C     if (myProcId==nroot) then
          ix_atm=28672+237600           ! Offset from the  32*64*14 atmos + 44*90*15*4 ocean
c         ix_atm=28672+118800
          ix=ix_atm+floor(real(kstepo-3)/3.)*Nx*Ny
          call cont2mod_oce_xy(xxFld2d,myThid) ! FastOpt          
C        print*,'in plasim_get_atmdat check-point 2'
C	 endif
C         print*,'in plasim_get_atmdat check-point 3'
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do j = 1,sny
             do i = 1,snx
                if (whflux(i,j,bi,bj) .ne. 0.) then
                   Qnet(i,j,bi,bj) = Qnet(i,j,bi,bj)
     &              +xxFld2d(i,j,bi,bj)
     &                  /sqrt(whflux(i,j,bi,bj))
                endif
            enddo
          enddo
        enddo
      enddo
C      print*,'in plasim_get_atmdat check-point 4'
#endif
#ifdef ALLOW_SFLUX_CONTROL
         ix_atm=28672+237600+nflxdim         ! Offset from the  32*64*14 atmos ctrl vector
         ix=ix_atm+floor(real(kstepo-3)/3.)*Nx*Ny
         call cont2mod_oce_xy(xxFld2d,myThid) ! FastOpt          
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do j = 1,sny
             do i = 1,snx
                if (wsflux(i,j,bi,bj) .ne. 0.) then
                   EmPmR(i,j,bi,bj) = EmPmR(i,j,bi,bj)
     &              +xxFld2d(i,j,bi,bj)
     &                  /sqrt(wsflux(i,j,bi,bj))
               endif
            enddo
          enddo
        enddo
      enddo
#endif
#ifdef ALLOW_USTRESS_CONTROL
         ix_atm=28672+237600+nflxdim*2         ! Offset from the  32*64*14 atmos ctrl vector
         ix=ix_atm+floor(real(kstepo-3)/3.)*Nx*Ny
         call cont2mod_oce_xy(xxFld2d,myThid) ! FastOpt          
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do j = 1,sny
             do i = 1,snx
                if (wtauu(i,j,bi,bj) .ne. 0.) then
                    fu(i,j,bi,bj) = fu(i,j,bi,bj)
     &              +xxFld2d(i,j,bi,bj)
     &                  /sqrt(wtauu(i,j,bi,bj))
                endif
           enddo
          enddo
        enddo
      enddo
#endif
#ifdef ALLOW_VSTRESS_CONTROL
         ix_atm=28672+237600+nflxdim*3         ! Offset from the  32*64*14 atmos ctrl vector
         ix=ix_atm+floor(real(kstepo-3)/3.)*Nx*Ny
         call cont2mod_oce_xy(xxFld2d,myThid) ! FastOpt          
      do bj = jtlo,jthi
        do bi = itlo,ithi
          do j = 1,sny
             do i = 1,snx
                if (wtauv(i,j,bi,bj) .ne. 0.) then
                    fv(i,j,bi,bj) = fv(i,j,bi,bj)
     &              +xxFld2d(i,j,bi,bj)
     &                  /sqrt(wtauv(i,j,bi,bj))
               endif
            enddo
          enddo
        enddo
      enddo
#endif

c     Update the tile edges.

      _EXCH_XY_RS(  qnet, mythid )
      _EXCH_XY_RS( empmr, mythid )
       CALL EXCH_UV_XY_RS(fu, fv, .TRUE., myThid)
#ifdef SHORTWAVE_HEATING
      _EXCH_XY_RS(   Qsw, mythid )
#endif
#ifdef ATMOSPHERIC_LOADING
      _EXCH_XY_RS( pload, mythid )
      _EXCH_XY_RS(sIceLoad, mythid )
#endif
!!      print*,'NACHHER '
!      print*,'QNET: ',maxval(qnet),minval(qnet)
!      print*,'QSW: ',maxval(qsw),minval(qsw)
!      print*,'PME: ',maxval(EmPmR),minval(EmPmR)
!      print*,'USTR:',maxval(fu),minval(fu)
!      print*,'VSTR:',maxval(fv),minval(fv)
!      print*,'PS:',maxval(pload),minval(pload)
!!      print*,'IL:',maxval(sIceLoad),minval(sIceLoad)
c
      return
      end
