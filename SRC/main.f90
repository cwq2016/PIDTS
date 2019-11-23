!==================================================================================
!> Main program of the PIDTS method for natural convection problem
!==================================================================================
PROGRAM MAIN
USE MOD_Global
USE MOD_MPI
USE MOD_FileIO
USE MOD_PIDTS
USE MOD_Post
USE MOD_BC
IMPLICIT NONE
INTEGER::i

! Initlize MPI
CALL Init_MPI
    
! Initlize part of parameters in MOD_Global
CALL Init_Parameters

! Initlize boundary conditions
CALL INIT_BC

! Initlize postprocess 
CALL Init_Post

! each CPU core independently call InitFields to avoid IO error
DO i=0,nProc-1
    CALL MPI_Barrier(MPI_COMM_WORLD, IERR)
    IF(MyID==i) CALL InitFields
ENDDO
    
! Start and end time to be solved 	
TStart=Time
TEnd  =Time+DBLE(NT)*DT
    
! Root CPU records parameters	
IF (MyID==0)THEN
    CALL  Record_Parameters
    CALL  Record_Time
ENDIF

! start calculation
CALL PIDTS

! Close files of recording time
CALL Record_Time_Close

CALL MPI_Barrier(MPI_COMM_WORLD, IERR)

! Merge all files of recording time 
IF(MyID==0) CALL Merge_Time_File

CALL StopMPI
    
PRINT*,'MyID=',MyID,'has finish'
    
END PROGRAM 

    