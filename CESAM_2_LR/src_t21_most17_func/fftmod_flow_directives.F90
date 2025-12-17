!$taf SUBROUTINE gp2fc INPUT   = 1,2,3
!$taf SUBROUTINE gp2fc OUTPUT  = 1
!$taf SUBROUTINE gp2fc ACTIVE  = 1
!$taf SUBROUTINE gp2fc DEPEND  = 2,3
!$taf SUBROUTINE gp2fc ADNAME  = adgp2fc
!$taf SUBROUTINE gp2fc FTLNAME  = gp2fc

!$taf SUBROUTINE fc2gp INPUT   = 1,2,3
!$taf SUBROUTINE fc2gp OUTPUT  = 1
!$taf SUBROUTINE fc2gp ACTIVE  = 1
!$taf SUBROUTINE fc2gp DEPEND  = 2,3
!$taf SUBROUTINE fc2gp ADNAME  = adfc2gp
!$taf SUBROUTINE fc2gp FTLNAME  = fc2gp

      subroutine adgp2fc(ada,n,lot)
! interface routine for adjoint fft
! input:
!      integer n       ! leading dimension of ada
!      integer lot     ! second dimension of ada
!      REAL ada(n,lot) ! adjoint data field
! output:
!      REAL ada(n,lot) ! adjoint data field
!
      use pumamod
      implicit none
      integer n,lot
      REAL ada(n,lot)
      ada(1,:)=ada(1,:)/real(n)
      ada(2,:)=0.
      ada(2*NTP1+1:,:)=0.       ! these are supposed to be zero
      ada(3:2*NTP1,:)=ada(3:2*NTP1,:)/real(n)/2.
      call fc2gp(ada,n,lot) ! fc2gp is transpose of gp2fc
      return
      end

      
      subroutine adfc2gp(ada,n,lot)
! interface routine for adjoint fft
! input:
!      integer n       ! leading dimension of ada
!      integer lot     ! second dimension of ada
!      REAL ada(n,lot) ! adjoint data field
! output:
!      REAL ada(n,lot) ! adjoint data field
!
      use pumamod
      implicit none
      integer n,lot
      REAL ada(n,lot)
      call gp2fc(ada,n,lot) ! gp2fc is transpose of fc2gp
      ada(1,:)=ada(1,:)*real(n)
      ada(3:2*NTP1,:)=ada(3:2*NTP1,:)*2.*real(n)
      return 
      end

