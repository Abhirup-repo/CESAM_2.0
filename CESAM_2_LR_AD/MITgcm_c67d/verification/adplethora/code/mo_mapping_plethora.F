#include "CPP_OPTIONS.h"

#undef UNIFORM_PERT
#undef DEBUG

      module mo_mapping_plethora
! author: FastOpt GmbH
! purpose: map the control vector onto the ocean variables
! contains:
!          subroutine cont2mod_oce_xy(active_var,myThid)
!          subroutine cont2mod_oce_xyz(active_var,myThid)
!          subroutine cont2mod_oce_xz(active_var,myThid)
!          subroutine cont2mod_oce_yz(active_var,myThid)
! private: jinc
! public: everything else
!
        use mo_mapping, only : ix, xsave, ncontrol=>nx, mypid, nroot
        implicit none
!     == global variables ==
#include "EEPARAMS.h"
#include "SIZE.h"
!     private variables
        integer, private :: jinc
        contains
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          subroutine cont2mod_oce_xy(active_var,myThid)
! purpose: map control vector on 'active_var' (2D horizontal)
! input:   myThid
! output:
!         _RL     active_var(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)
!
            integer :: myThid
            _RL     active_var(1-olx:snx+olx,1-oly:sny+oly,nsx,nsy)
#ifdef UNIFORM_PERT
            jinc = 1
            active_var = xsave(ix+1)
#else
            jinc = Nx * Ny
            call  SCATTER_VEC_R8(xsave(ix+1),
     &            active_var(1:snx,1:sny,:,:),
     &            snx*sny*nsx*nsy,
     &            myThid)    
#endif
!            ix = ix + jinc
#ifdef DEBUG
            if(mypid == NROOT)
     &           print*,'perturbing xy, ix=',ix,' with',xsave(ix)
#endif
            if (ix.gt.ncontrol) 
     &           stop 'fix dimension of control variable vector!'
          end subroutine cont2mod_oce_xy
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          subroutine cont2mod_oce_xyz(active_var,myThid)
! purpose: map control vector on 'active_var' (3D)
! input:   myThid
! output:
!         _RL active_var(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
!
            integer :: myThid,kid 
            _RL active_var(1-olx:snx+olx,1-oly:sny+oly,nr,nsx,nsy)
#ifdef UNIFORM_PERT
            jinc = 1
            active_var = xsave(ix+1)
#else 
          jinc = Nx * Ny
          print*, "SUM of xsave in cont2mod", sum(xsave(28672:88072))
          do kid=1,Nr
          call  SCATTER_VEC_R8(xsave((kid-1)*jinc+ix+1:kid*jinc+ix+1),
     &           active_var(1:snx,1:sny,kid,:,:),
     &           snx*sny*nsx*nsy,
     &            myThid) 
          enddo
             ix=ix+jinc*Nr

#endif
!            ix = ix + jinc
#ifdef DEBUG
            if(mypid == NROOT)
     &           print*,'perturbing xyz, ix=',ix,' with',xsave(ix)
#endif
            if (ix.gt.ncontrol)
     &           stop 'fix dimension of control variable vector!'
          end subroutine cont2mod_oce_xyz
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          subroutine cont2mod_oce_xz(active_var,myThid)
! purpose: map control vector on 'active_var' (2D zonal?)
! input:   myThid
! output:
!         _RL     active_var(1-olx:snx+olx,nsx,nsy)
!
            integer :: myThid
            _RL     active_var(1-olx:snx+olx,nsx,nsy)
#ifdef UNIFORM_PERT
            jinc = 1
            active_var = xsave(ix+1)
#else
            jinc = Nx
            call scatter_vector(snx*nsx*nsy,
     &           xsave(ix+1),
     &           active_var(1:snx,:,:),
     &           myThid)
#endif
!            ix = ix + jinc
#ifdef DEBUG
            if(mypid == NROOT) 
     &           print*,'perturbing xz, ix=',ix,' with',xsave(ix)
#endif
            if (ix.gt.ncontrol)
     &           stop 'fix dimension of control variable vector!'
          end subroutine cont2mod_oce_xz
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          subroutine cont2mod_oce_yz(active_var,myThid)
! purpose: map control vector on 'active_var' (2D meridional?)
! input:   myThid
! output:
!         _RL     active_var(1-oly:sny+oly,nsx,nsy)
!
            integer :: myThid
            _RL     active_var(1-oly:sny+oly,nsx,nsy)
#ifdef UNIFORM_PERT
            jinc = 1
            active_var = xsave(ix+1)
#else
            jinc = Ny
            call scatter_vector(sny*nsx*nsy,
     &           xsave(ix+1),
     &           active_var(1:sny,:,:),
     &           myThid)
#endif
!            ix = ix + jinc
#ifdef DEBUG
            if(mypid == NROOT)
     &           print*,'perturbing yz, ix=',ix,' with',xsave(ix)
#endif
            if (ix.gt.ncontrol)
     &           stop 'fix dimension of control variable vector!'
          end subroutine cont2mod_oce_yz

      end module mo_mapping_plethora
