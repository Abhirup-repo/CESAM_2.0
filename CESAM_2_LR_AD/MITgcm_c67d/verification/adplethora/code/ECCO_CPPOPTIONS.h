#ifndef ECCO_CPPOPTIONS_H
#define ECCO_CPPOPTIONS_H

C-- Collect here, in a single option-file, options to control which optional
C   features to compile in packages AUTODIFF, COST, CTRL, ECCO, CAL and EXF.
C   If used, this option-file needs to be directly included in CPP_OPTIONS.h
C   Although this method, inherited from ECCO setup, has been traditionally
C   used for all adjoint built, work is in progess to allow to use the
C   standard method (each of the above pkg get its own options from its
C   specific option-file) also for adjoint built.

C ********************************************************************
C ***                         ECCO Package                         ***
C ********************************************************************

C  see pkg/ecco/ECCO_OPTIONS.h file
C  Allow use of legacy ecco/ctrl codes
#define ECCO_CTRL_DEPRECATED

C  Cost function output format
#undef ALLOW_ECCO_OLD_FC_PRINT

C  Allow for generic cost function and integral terms
c#define ALLOW_GENCOST_CONTRIBUTION
C  Allow for 3 dimensional generic terms
c#define ALLOW_GENCOST3D

C ********************************************************************
C ***                  Adjoint Support Package                     ***
C ********************************************************************

C  Include/exclude code in order to be able to automatically
#define ALLOW_AUTODIFF_TAMC

C  Checkpointing as handled by TAMC
#define ALLOW_TAMC_CHECKPOINTING

C  Extend to 2-level checkpointing
#define AUTODIFF_2_LEVEL_CHECKPOINT

C  Extract adjoint state
#undef ALLOW_AUTODIFF_MONITOR

C  Use divided adjoint to split adjoint computations
#undef ALLOW_DIVIDED_ADJOINT
#undef ALLOW_DIVIDED_ADJOINT_MPI

C  Tape settings
C#define ALLOW_AUTODIFF_WHTAPEIO
C#define AUTODIFF_USE_OLDSTORE_2D
C#define AUTODIFF_USE_OLDSTORE_3D
C#define EXCLUDE_WHIO_GLOBUFF_2D
C#define ALLOW_INIT_WHTAPEIO


C ********************************************************************
C ***                     Calendar Package                         ***
C ********************************************************************

CPH >>>>>> THERE ARE NO MORE CAL OPTIONS TO BE SET <<<<<<

C ********************************************************************
C ***                Cost function Package                         ***
C ********************************************************************

#define ALLOW_COST
#define ALLOW_ECCO
#undef ALLOW_COST_ATLANTIC_HEAT
#undef ALLOW_COST_ATLANTIC_HEAT_DOMASS
c#define ALLOW_ECCO_EVOLUTION

#define ALLOW_THETA_COST_CONTRIBUTION
#define ALLOW_SALT_COST_CONTRIBUTION
c AB: 22.01.25


c# define ALLOW_SSH_MEAN_COST_CONTRIBUTION
c# define ALLOW_SSH_TPANOM_COST_CONTRIBUTION
c# define ALLOW_SSH_ERSANOM_COST_CONTRIBUTION
c# define ALLOW_SSH_GFOANOM_COST_CONTRIBUTION
c# if (defined (ALLOW_SSH_MEAN_COST_CONTRIBUTION) || \
c      defined (ALLOW_SSH_TPANOM_COST_CONTRIBUTION) || \
c      defined (ALLOW_SSH_ERSANOM_COST_CONTRIBUTION))
c#  define ALLOW_SSH_COST_CONTRIBUTION
c# endif

#define GENERIC_BAR_MONTH
c#define ALLOW_THETA0_COST_CONTRIBUTION

#define ALLOW_KAPREDI_COST_CONTRIBUTION
#define ALLOW_KAPGM_COST_CONTRIBUTION

C
#define ALLOW_SST_COST_CONTRIBUTION
c#define ALLOW_SSS_COST_CONTRIBUTION
#define ALLOW_ARGO_THETA_COST_CONTRIBUTION
#define ALLOW_ARGO_SALT_COST_CONTRIBUTION
c#define ARGO_PTEMP
 
C  
#undef ALLOW_HFLUX_COST_CONTRIBUTION
#undef ALLOW_SFLUX_COST_CONTRIBUTION
#undef ALLOW_USTRESS_COST_CONTRIBUTION
#undef ALLOW_VSTRESS_COST_CONTRIBUTION


C ********************************************************************
C ***               Control vector Package                         ***
C ********************************************************************

c AB: 22.11.24
c This is the flux correction part which is activated in plasim_get_atmdat.F

c AB: 22.01.24
c Additional controls are added in the code 

C  I/O and pack settings
#undef ALLOW_NONDIMENSIONAL_CONTROL_IO
#undef ALLOW_PACKUNPACK_METHOD2

C  Initial values.
#define ALLOW_THETA0_CONTROL
#define ALLOW_SALT0_CONTROL

C Sets of controls
#undef ALLOW_GENTIM2D_CONTROL
#undef ALLOW_GENARR2D_CONTROL
#undef ALLOW_GENARR3D_CONTROL

#undef ALLOW_HFLUX_CONTROL
#undef ALLOW_SFLUX_CONTROL
#undef ALLOW_USTRESS_CONTROL
#undef ALLOW_VSTRESS_CONTROL

#define ALLOW_KAPGM_CONTROL
#define ALLOW_KAPREDI_CONTROL

C ********************************************************************
C ***             External forcing Package                         ***
C ********************************************************************

C-- see pkg/exf/EXF_OPTIONS.h file

#endif /* ECCO_CPPOPTIONS_H */
