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
!> Routines for implementing the parallel inverted dual time stepping method.
!> Each CPU core is employed to solve one time step with in the chunk.
!> In each iteration, all time step update its solution by one pseudo time marching unit(PTMU),
!> and then CPU cores communicate with others to request latest solutions at previous time steps.
!=================================================================================================================================
MODULE MOD_PIDTS  
USE MOD_Global
USE MOD_MPI
IMPLICIT NONE
PRIVATE
SAVE

!---------------------------------------------------------------------------------------------------------------------------------
! PUBLIC VARIABLES
PUBLIC::PIDTS
!---------------------------------------------------------------------------------------------------------------------------------
! PRIVATE VARIABLES
INTEGER::DataCount  ! The count of communication date(the solutions of one time step)
!=================================================================================================================================

CONTAINS

!=================================================================================================================================
!> The main program of the PIDTS method.
!=================================================================================================================================
SUBROUTINE PIDTS
USE MOD_Global
USE MOD_ErrorAnalysis,ONLY:ErrorCriterion
USE MOD_FileIO,    ONLY: Record_Time,Print_Result
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER ::MyRank       !< Used to denote the location of current CPU core in the chunk
REAL(RP)::MyTime       !< The physical time being solved
REAL(RP)::ChunkLen     !< The physical time duration of the chunk
LOGICAL ::HeadUpdate   !< Used to denote whether the chunk head has converged
LOGICAL ::Restart      !< Used to denote whether to start the solving of a new time step
INTEGER ::HeadShift    !< USed to denote the shift of CPU cores in the chunk
INTEGER ::i            !< Loop variable
!=================================================================================================================================

!-------------- initialization ----------------------

DataCount=(Nx+1)*(Ny+1)*nVar !The count of communication date

! Define the solutions of the required time steps for current time steps
IF (MyID<3)THEN
    DO i=-2,-MyID
        Uall(:,:,:,i)=Uall(:,:,:,i+MyID)
    ENDDO
ENDIF

! Record the start time by root CPU
IF (MyID==0)THEN
    Time=TStart
    CALL Record_Time()
ENDIF

ChunkLen=DBLE(nProc)*DT ! Chunk length 
    
!------------------ Startup -------------------------
! Serially update all time steps in the chunk by one PTMU

PRINT*,'START-UP FOR CHUNK...'
MyRank=MyID                     ! At the begining, MyRank=MyID
HeadShift=0
MyTime=TStart+DBLE(MyRank+1)*DT ! The time step handled by each 

DO i=0,nProc-1
    IF(MyID==i) Uone=Uall(:,:,:,0)
    Time=MyTime
    CALL PTMU
    CALL Commucation(HeadShift,MyRank)
    CALL MPI_Barrier(MPI_COMM_WORLD, IERR) 
ENDDO

!------------------ Iteration -------------------------    
HeadUpdate=.FALSE.
Restart=.FALSE.
ErrorCriterion=10.D0
!迭代求解过程
PRINT*,'ITERATION...'
DO WHILE(.TRUE.)
       
	! the CPU core change to solve a new time step
    IF(Restart)THEN
        Uone=Uall(:,:,:,0)
        Restart=.FALSE.
    ENDIF
	
	! Update by one PTMU
    Time=MyTime
    CALL PTMU
	
    ! The CPU with (MyRank==0) handling the chunk head reach convergence
    IF (MyRank==0 .AND. ErrorCriterion<EPS) Then
            
            WRITE(*,200)MyID,MyTime
200             format('MyID=',I2.1,': Time=',F5.1,',has convergenced')
                
            
			! Update status variables
            HeadUpdate=.TRUE.
            Restart=.TRUE.
                
            !record time 
            CALL Record_Time()
                
            ! Copy the convergent solution the temporary variables
			! in order to collectlively print results.
			
            Utmp=Uone
                
            ! Special treatment for nproc<4.
			! As each time step require the solutions at previous three time steps,
			! the CPU core handling the chunke head (t) will turn to hanle a new time step (t+nproc).
		    ! Therefore, current convergent solution will the the requied solution at the new time step.
            IF (nProc<=3)THEN
                DO i=-2,0
                    IF(i+nProc<=0) Uall(:,:,:,i)=Uall(:,:,:,i+nProc)
                ENDDO
                Uall(:,:,:,-nProc+1)=Uone
            ENDIF

            ! Update the time be saved
            MyTime=MyTime+ChunkLen
			
			! set ErrorCriterion as a relatively large value
            ErrorCriterion=10
    ENDIF

	! Share the new HeadUpdate with all CPU cores
    CALL MPI_BCAST(HeadUpdate, 1, MPI_LOGICAL, Rank_to_ID(HeadShift, 0), MPI_COMM_WORLD, IERR)
    ! Note that calling MPI_BCAST will implicitly call MPI_Barrier
        
	! As Chunk head has just reached convergence, each CPU core receieves this information from HeadUpdate.
	! Then, all CPU cores will update their status variables accordingly.
	! If the chunk have move forward for nProc time steps, all CPU cores will collectively write results.
    IF(HeadUpdate)THEN
        ! Check whether the chunk have move forward for nProc time steps
        IF(HeadShift==nProc-1)THEN
            
            Time=MyTime-ChunkLen ! The time of output solution
            U0=Utmp              ! Copy solution to U0, as Print_Result can only write solution stored in U0.
            
			! Check whether to write to write results
			IF(WriteResult)Then
				CALL Print_Result
			ENDIF
            
            !Check whether to stop according to the time handled by the CPU with (MyRank=0)
            Time=MyTime-DBLE(MyRank-0)*DT
            IF (MyRank==0) Time=MyTime-ChunkLen
            IF (Time>=TEnd) EXIT ! Stop PIDTS, and return to main program.
        ENDIF
            

		! Send the converegnt at the chunk head to other CPUs
        CALL Commucation(HeadShift,MyRank)
        CALL MPI_Barrier(MPI_COMM_WORLD, IERR)
            
		! Update MyRank and HeadUpdate
        MyRank=MyRank-1
        IF (MyRank<0) MyRank=MyRank+nProc
        HeadUpdate=.FALSE.
    ENDIF
	
    ! Update HeadShift
    HeadShift=MyID-MyRank
    IF (HeadShift<0) HeadShift=HeadShift+nProc
	
	! Communication
    CALL Commucation(HeadShift,MyRank)
ENDDO
END SUBROUTINE PIDTS

!=================================================================================================================================
!> Transform myrank to to CPU id, according to HeadShift
!=================================================================================================================================
FUNCTION Rank_to_ID(HeadShift,MyRank)
IMPLICIT NONE
INTEGER::Rank_to_ID,HeadShift,MyRank
IF (0<=MyRank .AND. MyRank<nProc)THEN
    Rank_to_ID=MOD(MyRank+HeadShift, nProc)
ELSE
    Rank_to_ID=MyRank
ENDIF    
END FUNCTION Rank_to_ID

!=================================================================================================================================
!> Communication among CPU cores.
!> Each CPU core will [send] current solution to its later three CPUs, if they exist.
!> Each CPU core will [receieve] solutions to its previous three CPUs, if they exist.
!=================================================================================================================================
SUBROUTINE Commucation(HeadShift, MyRank)
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)::HeadShift    !< The shift of CPU cores in the chunk 
INTEGER,INTENT(IN)::MyRank       !< The location of COU core in the chunk
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER::Dest0,Dest1,Dest2       !< the ID of later    three CPU cores
INTEGER::Source0,Source1,Source2 !< the ID of previous three CPU cores
!=================================================================================================================================

Dest0  =Rank_to_ID(HeadShift, MyRank+1)
Dest1  =Rank_to_ID(HeadShift, MyRank+2)
Dest2  =Rank_to_ID(HeadShift, MyRank+3)
Source0=Rank_to_ID(HeadShift, MyRank-1)
Source1=Rank_to_ID(HeadShift, MyRank-2)
Source2=Rank_to_ID(HeadShift, MyRank-3)
	
! COMMUNICATION
IF (0<=Dest0 .AND. Dest0<nProc)THEN
    CALL MPI_SEND( Uone(0,0,1), DataCount, MPI_DOUBLE_PRECISION, Dest0, MyID+1000, MPI_COMM_WORLD,IERR)
ENDIF
IF (0<=Dest1 .AND. Dest1<nProc)THEN
    CALL MPI_SEND( Uone(0,0,1), DataCount, MPI_DOUBLE_PRECISION, Dest1, MyID+2000, MPI_COMM_WORLD,IERR)
ENDIF
IF (0<=Dest2 .AND. Dest2<nProc)THEN
    CALL MPI_SEND( Uone(0,0,1), DataCount, MPI_DOUBLE_PRECISION, Dest2, MyID+3000, MPI_COMM_WORLD,IERR)
ENDIF
        
IF (0<=Source0 .AND. Source0<nProc)THEN
    CALL MPI_RECV( Uall(0,0,1,0), DataCount, MPI_DOUBLE_PRECISION, Source0, Source0+1000, MPI_COMM_WORLD,STATUS,IERR)
ENDIF
IF (0<=Source1 .AND. Source1<nProc)THEN
    CALL MPI_RECV( Uall(0,0,1,-1), DataCount, MPI_DOUBLE_PRECISION, Source1, Source1+2000, MPI_COMM_WORLD,STATUS,IERR)
ENDIF
IF (0<=Source2 .AND. Source2<nProc)THEN
    CALL MPI_RECV( Uall(0,0,1,-2), DataCount, MPI_DOUBLE_PRECISION, Source2, Source2+3000, MPI_COMM_WORLD,STATUS,IERR)
ENDIF 
END SUBROUTINE Commucation
    
!=================================================================================================================================
!> Pseudo time marching unit,
!> Containing the marching Npt of pseudo time steps
!=================================================================================================================================
SUBROUTINE PTMU
USE MOD_NSEqs ,    ONLY: PseudoTimeIntegrator
USE MOD_FileIO,    ONLY: Print_StatusLine
USE MOD_ErrorAnalysis 
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER::ipt
!=================================================================================================================================

ErrorCriterion=1.
DO ipt=1,Npt
    U0=Uone
    CALL PseudoTimeIntegrator
    CALL Calerr 
	! Print intermediate result, aiming to debug the code
    IF (MOD(ipt,IterPrint)==0)THEN
	    CALL Calallerr
	    CALL Print_StatusLine(ipt)  
    ENDIF
	
	! Reach convergence 
    IF (ErrorCriterion<eps) EXIT
ENDDO
END SUBROUTINE
    

END MODULE