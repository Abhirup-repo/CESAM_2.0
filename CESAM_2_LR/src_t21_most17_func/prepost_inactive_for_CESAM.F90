!-*-f90-*-

! This file is cpp-preprocessed. Be sure to edit prepost.F90

!================================================
! This subroutine sets the number
! of control variables
!================================================
SUBROUTINE NUMBMOD( N )
  use mo_mapping, only : nx
  implicit none
  INTEGER N
  N = nx
END SUBROUTINE NUMBMOD

!================================================
! The subroutine "INITMOD" is called before the
! optimization. It must set a first guess
! of the parameter vector.
! It may also contain the initialization of
! the model.
!================================================
SUBROUTINE INITMOD( N, X )
  use observation, only: initobs
  use pumamod, only : pumaversion,nlon
  use mo_mapping
  implicit none
  INTEGER N
  REAL X(N)

  call mpstart
  !
  !     set version identifier via FPP/CPP; macro definition must contain quotes
  !
  pumaversion=SVNVERSION
  !
  ! initialise model
  !
  call prolog
  !
  ! fix some values
  !
  !
  if (nlandt==0) then
     where (dtcl.lt.1)
        dtcl(:,0:13) = spread(dt(:,nlep),dim=2,ncopies=14)
     end where
     if (any(dtcl(:,0:13)<1.)) stop 'still too cold spots in dtcl!'
     if (any(dtcl(:,0:13)>374.)) stop 'way too hot spots in dtcl!'
  end if
  !
  ! initialise cost function
  !
  call initobs
  !
  ! save initial state
  !
  call write_internal_restart('internal_restart')
  !
  ! pack arguments for model
  !
  call init_global_mapping
  call mod2cont(n,x)
END SUBROUTINE INITMOD

!================================================
! The subroutine "POSTMOD" is called after the function is evalutated
! It should contain the output of the function value
!================================================
SUBROUTINE POSTMOD( N, X, FC)
  use observation, only: postobs
  
  INTEGER N
  REAL X(N), FC
  INTEGER I
  call epilog
  call mpstop
  !
  ! end cost function
  ! 
  call postobs
  PRINT *,'**************************************'
  PRINT *,'***     FUNCTION VALUE              **'
  PRINT *,'**************************************'
  PRINT 9010, FC
  PRINT 9020
  DO I = 1,N
     PRINT 9030, I, X(I)
  ENDDO
9010 FORMAT(1X,'The value of the cost function is : ',E15.9)
9020 FORMAT(1X,'The parameter values are :')
9030 FORMAT(1X,I5,1X,E15.9,1X,E15.9)
  
  !      call dealloc_unknown
END SUBROUTINE POSTMOD


SUBROUTINE setfunc( n, m )
  use mo_mapping, only : nx,ny
  IMPLICIT NONE
  INTEGER n, m
  n = nx
  m = ny
END SUBROUTINE setfunc


SUBROUTINE initfunc( n, x )
  IMPLICIT NONE
  INTEGER n
  REAL x(n)
  CALL initmod( n, x )
END SUBROUTINE initfunc

SUBROUTINE postfunc( n, x, m, y )
  use observation, only: postobs
  IMPLICIT NONE

  INTEGER :: n, m
  REAL :: x(n), y(m)
  INTEGER :: j
  
  call epilog
  call postobs
  call mpstop

  PRINT *,'**************************************'
  PRINT *,'***         Function VALUE          **'
  PRINT *,'**************************************'
  DO J = 1,M
     PRINT '(1X,I4,20(1X,E12.6))', J, Y(J)
  ENDDO
  
END SUBROUTINE postfunc

