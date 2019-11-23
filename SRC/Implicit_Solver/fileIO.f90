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
!> Routines to read and write solutions
!=================================================================================================================================
MODULE MOD_FileIO
USE MOD_Global
IMPLICIT NONE
PUBLIC


CONTAINS
!=================================================================================================================================
!> Write convergent solution U0 to file
!=================================================================================================================================
SUBROUTINE Print_Result
USE MOD_NSEqs
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER i,j
CHARACTER(100)::FileName
!=================================================================================================================================

! Compute residual
CALL Unsteady_Residual

!输出solution至各个时刻对应的文件夹中
TimeDir=''
WRITE(TimeDir,'(F10.3)')Time
TimeDir=TRIM(ADJUSTL(PathOutput))//'/'//'Time='//TRIM(ADJUSTL(TimeDir))
CALL system('mkdir "'//TRIM(ADJUSTL(TimeDir))//'"')

! Write to a tecplot file
FileName=''
WRITE(FileName,'("RESULT.plt")')
OPEN(UNIT=20,FILE=TRIM(ADJUSTL(TimeDir))//'/'//TRIM(ADJUSTL(FileName)))
WRITE(20,*) 'title="result"'
WRITE(20,*) 'variables="x","y","P","u","v","t","RES1","RES2","RES3","RES4"'
WRITE(20,"(a7,i3,a3,i3,a8)")&
	'zone,j=',Ny+1,',i=',Nx+1,',f=point'
DO j=0,Ny
	DO i=0,Nx
		WRITE(20,"(2(f8.4,1x),8(f21.16,1x))")x(i),y(j),U0(i,j,1),U0(i,j,2),U0(i,j,3),U0(i,j,4),&
		                                    & Resi(i,j,1),Resi(i,j,2),Resi(i,j,3),Resi(i,j,4)
	ENDDO
ENDDO
close(20)
END SUBROUTINE

!=================================================================================================================================
!> Write status to terminal
!=================================================================================================================================
SUBROUTINE Print_StatusLine(ipt)
USE MOD_ErrorAnalysis
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)::ipt
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP)::times
!==================================================================================

CALL CPU_TIME(times)
WRITE(*,100)               times,ipt,ErrorCriterion, ErrorL2(1), ErrorL2(2), ErrorL2(3),ErrorL2(4)
!WRITE(FILE_FID_ERROR,100) times,ipt,ErrorCriterion, ErrorL2(1), ErrorL2(2), ErrorL2(3),ErrorL2(4)
100 FORMAT(F10.2,I10,5(1X,F10.5))
END SUBROUTINE

!=================================================================================================================================
!> Record input parameters
!=================================================================================================================================
SUBROUTINE Record_Parameters
use MOD_Global
IMPLICIT NONE

IF(Time==0)THEN
	OPEN(UNIT=4,FILE=TRIM(ADJUSTL(PathOutput))//'/COMPU_INFO.dat',STATUS='REPLACE')
	WRITE(4,*)"----------------------------------------------------------------------"
	WRITE(4,*)TRIM(ADJUSTL(Date_Time_Str()))
	WRITE(4,*)"----------------------------------------------------------------------"
ELSE
	OPEN(UNIT=4,FILE=TRIM(ADJUSTL(PathOutput))//'/COMPU_INFO.dat',STATUS='OLD',POSITION='APPEND')
	WRITE(4,*)"----------------------------------------------------------------------"
	WRITE(4,*)TRIM(ADJUSTL(Date_Time_Str()))
	WRITE(4,*)"----------------------------------------------------------------------"
ENDIF
WRITE(4,NML=Rectangle)
WRITE(4,NML=Pulsating_Temperature)
WRITE(4,NML=Physical_Property)
WRITE(4,NML=Grid_Fine)
WRITE(4,NML=PTMU_Paras)
WRITE(4,NML=Unsteady_parameters)
WRITE(4,NML=Time_Parallel)
WRITE(4,NML=Output_Control)
CLOSE(4)
END SUBROUTINE

 
!=================================================================================================================================
!> Record the consumed wall time and other postprocess results of each time step.
!> As each CPU handle a time step, each CPU will record time step assigned to itself in a independent file.
!> At last, A routine [Merge_Time_File] will merge all the files generated by all CPUs to one file.
!=================================================================================================================================
SUBROUTINE Record_Time
USE MOD_Global
USE MOD_MPI
USE MOD_POST
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER,SAVE  ::STAT(0:23)=0     !< Stat of file IO, 24 is the allowable maixmum number of CPU cores
REAL(RP)      ::NuL              !< the Nusselt number at the right sidewall
INTEGER       ::FileFlag,IOFlag  !< Flags
CHARACTER(100)::FileName         !< File name
CHARACTER(100),external::Date_Time_Str !< A string to describe date and time
!=================================================================================================================================

! Each has an independent file and also file ID
FileFlag=MyID*100+7
FileName=''
! File name
WRITE(FileName,'(I3.1)')MyID
FileName=TRIM(ADJUSTL(PathOutput))//'/Time'//TRIM(ADJUSTL(FileName))//'.dat'
! Check whether the file has be opened
IF(STAT(MyID)==0)THEN
	OPEN(FileFlag,FILE=TRIM(ADJUSTL(FileName)),STATUS='REPLACE',IOSTAT=IOFlag)
END IF

!Compute NuL
NuL=ComputeNu()

! Write 
WRITE(FileFlag,100)Time,MPI_WTIME()-Walltime0, MyID, NuL
100 FORMAT(F11.3,2X,F11.3,2X,I2.1,2X,F31.24)

! Update stat
STAT(MyID)=STAT(MyID)+1
END SUBROUTINE Record_Time

!================================================================================================================================
!> Close the file generated by each core.
!================================================================================================================================
SUBROUTINE Record_Time_Close
USE MOD_Global
USE MOD_MPI
IMPLICIT NONE
CLOSE(MyID*100+7)
END SUBROUTINE Record_Time_Close

!================================================================================================================================
!> Merge all the files generated by all CPUs to one file.
!================================================================================================================================
SUBROUTINE Merge_Time_File
USE MOD_Global
USE MOD_MPI
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER::i
CHARACTER(100)::FileName,Line
INTEGER::STAT
LOGICAL::EndFlag(0:nProc-1),Alive
!================================================================================================================================

! Check whether a time file exists
INQUIRE(FILE=TRIM(ADJUSTL(PathOutput))//'/Time.dat',EXIST=Alive)
IF(Alive)THEN
	OPEN(777,FILE=TRIM(ADJUSTL(PathOutput))//'/Time.dat',ACTION='WRITE',STATUS='OLD',POSITION='APPEND')
ELSE
	OPEN(777,FILE=TRIM(ADJUSTL(PathOutput))//'/Time.dat',ACTION='WRITE',STATUS='NEW')
ENDIF
! Open file 
DO i=0,nProc-1
	FileName=''
	WRITE(FileName,'(I3.1)')i
	FileName=TRIM(ADJUSTL(PathOutput))//'/Time'//TRIM(ADJUSTL(FileName))//'.dat'
	OPEN(i*100+11,FILE=TRIM(ADJUSTL(FileName)),ACTION='READ',STATUS='OLD')
ENDDO

! Write timestamp
WRITE(777,*)"!----------------------------------------------------------------------"
WRITE(777,*)"!",TRIM(ADJUSTL(Date_Time_Str()))
WRITE(777,*)"!----------------------------------------------------------------------"

i=0
Line=''
READ(i*100+11,'(A)',IOSTAT=STAT)Line
WRITE(777,'(A)')Line
Line=''
READ(i*100+11,'(A)',IOSTAT=STAT)Line
WRITE(777,'(A)')Line

EndFlag=.FALSE.
DO WHILE(ANY(.NOT.EndFlag))
	DO i=0,nProc-1
		IF (EndFlag(i)) CYCLE
		Line=''
		READ(i*100+11,'(A)',IOSTAT=STAT)Line
		IF (STAT/=0)THEN
			EndFlag(i)=.TRUE.
			CLOSE(i*100+11,STATUS='DELETE')
			EXIT
		ENDIF
		WRITE(777,'(A)')Line
	ENDDO
ENDDO
CLOSE(777)

END SUBROUTINE Merge_Time_File

!================================================================================================================================
!> Read solution from a Tecplot file.
!================================================================================================================================
SUBROUTINE Read_File
USE MOD_Global
IMPLICIT NONE
!----------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER::i,j
!================================================================================================================================

OPEN(1,file=TRIM(ADJUSTL(TimeDir))//'/RESULT.plt',ACTION='READ')
read(1,*)
read(1,*)
read(1,*)
DO j=0,Ny
	DO i=0,Nx
		read(1,*)x(i),y(j),Utmp(i,j,1),Utmp(i,j,2),Utmp(i,j,3),Utmp(i,j,4),&
		         & Resi(i,j,1),Resi(i,j,2),Resi(i,j,3),Resi(i,j,4)
	ENDDO
ENDDO
close(1)

END SUBROUTINE

!================================================================================================================================
!> Generate a timestamp.
!================================================================================================================================
FUNCTION Date_Time_Str()
IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
CHARACTER(100)::Date_Time_Str
INTEGER       ::Datetime(8)
character*10  ::b(3)
!================================================================================================================================

CALL date_and_time(b(1), b(2), b(3), Datetime)

WRITE(Date_Time_Str,100)Datetime(1),Datetime(2),Datetime(3),Datetime(5),Datetime(6),Datetime(7)
100 FORMAT(I4.4,"/",I2.2,"/",I2.2,", ",I2.2,":",I2.2,":",I2.2)
END FUNCTION

END MODULE