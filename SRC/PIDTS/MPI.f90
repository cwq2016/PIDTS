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
!> The module contains some MPI parameters and a MPI start routine .
!=================================================================================================================================
MODULE MOD_MPI
USE MOD_Global,ONLY:RP
USE MPI    
IMPLICIT NONE
PUBLIC
!---------------------------------------------------------------------------------------------------------------------------------
! PUBLIC VARIABLES
INTEGER :: nProc                   ! The number of CPU cores
INTEGER :: MyID                    ! CPU core ID
REAL(RP):: Walltime0               ! initial wall time
INTEGER :: STATUS(MPI_STATUS_SIZE) ! temporary variables
INTEGER :: IERR                    ! temporary variables
!=================================================================================================================================

CONTAINS
!=================================================================================================================================
!> Start MPI and get Walltime 
!=================================================================================================================================
SUBROUTINE Init_MPI
CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nProc,IERR)
CALL MPI_Comm_rank(MPI_COMM_WORLD,MyID,IERR)
Walltime0=MPI_WTIME()
END SUBROUTINE

!=================================================================================================================================
!> Stop MPI
!=================================================================================================================================
SUBROUTINE STOPMPI
CALL MPI_Finalize(IERR)
END SUBROUTINE

END MODULE