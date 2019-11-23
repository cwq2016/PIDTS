!=================================================================================================================================
! Copyright (c) 2019  Mr. Wenqian Chen
! This file is part of PIDTS, a parallel-in-time implementation of the dual time stepping (DTS) method.
!
! PIDTS is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! PIDTS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. IF not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
!=================================================================================================================================
!> The MOD_Global contains variables about
!> -- constant parameters
!> -- geometry and mesh
!> -- precompute variables for spatial discretization.
!> -- Input/output parameters
!> -- variables for tempral discretization.
!=================================================================================================================================
MODULE MOD_Global
IMPLICIT NONE
! parameters
INTEGER ,PARAMETER :: RP=8
INTEGER ,PARAMETER :: nVar=4
REAL(RP),PARAMETER :: PI=ACOS(-1.D0)
CHARACTER(100)     ::PATH_INPUT='INPUT'    !< The path of the input directory

! Geometry
REAL(RP)::Lenx,Leny   !< The length and width of the enclosure.
NAMELIST /Rectangle/                Lenx,Leny


! Fluctuating temperature on left sidewall
REAL(RP)::Amplitude,Period 
NAMELIST /Pulsating_Temperature/ Amplitude,Period

! Mesh
INTEGER :: Nx,Ny
NAMELIST /Grid_Fine/              Nx,Ny

! Physical preperty
REAL(RP)::Pr,Ra,RE,BetaP
NAMELIST /Physical_Property/        Pr,Ra,RE,BetaP

! Implicit solver
REAL(RP)::EPS
REAL(RP)::CFL
NAMELIST /PTMU_Paras/    EPS,CFL

!   Output parameters
INTEGER       :: IterPrint       !< Write intermidate result every pseudo time iterations 
CHARACTER(100):: PathOutput      !< Path of output directory
LOGICAL       :: WriteResult     !< Determine whether to write the results of each time step
NAMELIST /Output_Control/  IterPrint,PathOutput,WriteResult

! Unsteady definations
REAL(RP)::DT !< physical Time step
INTEGER ::NT !< the number of Time steps to be solved
NAMELIST /Unsteady_parameters/ DT,NT

! PIDTS parameters

INTEGER::Npt=100  ! The number of pseudo time iterations in one pseudo time marching unit
NAMELIST /Time_Parallel/ Npt

! Fields
REAL(RP),ALLOCATABLE,DIMENSION(:,:,:  ):: Uone,U0,Utmp,U_refine,Mid,Resi !< Solutions
REAL(RP),ALLOCATABLE,DIMENSION(:,:,:,:):: Uall !< copies of solutions at previous three Time steps

! Pre-built variables used in CHebyshev pseudospectral discretization 
REAL(RP) :: Kappa(nVar)
REAL(RP) :: xCoef,yCoef
REAL(RP) :: hx,hy    
REAL(RP),ALLOCATABLE,DIMENSION(:  ):: x,y
REAL(RP),ALLOCATABLE,DIMENSION(:,:):: DX,DY,DXX,DYY


! Simulation time
REAL(RP)      ::Time,TStart,TEnd           !< simulating Time  
CHARACTER(100)::TimeDir                    !< the directory for storing solutions 


END MODULE MOD_Global
