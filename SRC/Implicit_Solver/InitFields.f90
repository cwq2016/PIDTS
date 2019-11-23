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
!> The subroutine is used to initialize soluionts:
!> Case 1. Directly set solutions to zero except for boundary
!> Case 2. Read from files of the solutions at adjacent three time steps
!=================================================================================================================================
SUBROUTINE InitFields
use MOD_Global
USE MOD_FileIO, ONLY: Read_File
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
character(100)::str
INTEGER       ::ISTAT,id
REAL(RP)      ::TIME_LATEST0,TIME_LATEST1,TIME_LATEST2 !< The adjacent three time steps
!=================================================================================================================================

! Write the list of the output directory to a temporary file
CALL system('dir '//TRIM(ADJUSTL(PathOutput))//' /ad /b'//' >'//TRIM(ADJUSTL(PathOutput))//'/info.log')
OPEN(111,file=TRIM(ADJUSTL(PathOutput))//'/info.log',action='read',IOSTAT=ISTAT)

! Look up the newest three time steps
TIME_LATEST0=0.D0
TIME_LATEST1=0.D0
TIME_LATEST2=0.D0
Time=0.D0
DO WHILE(ISTAT==0)
	read(111,*,IOSTAT=ISTAT)str
	id=index(str,'Time=')
	IF (id/=1) CYCLE
	id=index(str,'=')
	str(1:id)=''
	read(str,*)Time
	IF (TIME_LATEST0<Time)THEN
		TIME_LATEST0=Time
	ELSEIF(TIME_LATEST1<Time)THEN
		TIME_LATEST1=Time
	ELSEIF(TIME_LATEST2<Time)THEN
		TIME_LATEST2=Time
	ENDIF
END DO
Time=TIME_LATEST0
CLOSE(111,STATUS='DELETE')

! Initialize solutions
IF( abs(Time)<1E-10 )THEN !Directly set solutions to zero except for boundary
	CALL TIME0
ELSE !Read from files of the solutions at adjacent three time steps
	TimeDir=''
	WRITE(TimeDir,'(F10.3)')TIME_LATEST0
	TimeDir=TRIM(ADJUSTL(PathOutput))//'/'//'Time='//TRIM(ADJUSTL(TimeDir))
	CALL Read_File()
	Uall(:,:,:,0)=Utmp
	
	TimeDir=''
	WRITE(TimeDir,'(F10.3)')TIME_LATEST1
	TimeDir=TRIM(ADJUSTL(PathOutput))//'/'//'Time='//TRIM(ADJUSTL(TimeDir))
	CALL Read_File
	Uall(:,:,:,-1)=Utmp
	
	TimeDir=''
	WRITE(TimeDir,'(F10.3)')TIME_LATEST2
	TimeDir=TRIM(ADJUSTL(PathOutput))//'/'//'Time='//TRIM(ADJUSTL(TimeDir))
	CALL Read_File
	Uall(:,:,:,-2)=Utmp
ENDIF

ENDSUBROUTINE InitFields





!=================================================================================================================================
!> Directly set solutions to zero except for boundary
!=================================================================================================================================
SUBROUTINE TIME0
USE MOD_Global
USE MOD_MPI
USE MOD_FileIO, ONLY: Print_Result
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER::i,K
!=================================================================================================================================

! Initialize solutions 
Utmp=0.
Utmp(0 ,0:Ny,4)=-0.5
Utmp(Nx,0:Ny,4)= 0.5

! Write solutions to file 
Time=0.D0
IF(MyID==0)THEN
	WRITE(TimeDir,'(F10.3)')Time
	TimeDir=TRIM(ADJUSTL(PathOutput))//'/'//'Time='//TRIM(ADJUSTL(TimeDir))
	CALL system('mkdir "'//TRIM(ADJUSTL(TimeDir))//'"')
	U0=Utmp
	CALL Print_Result()
ENDIF

! Initialize solutions at the require previous three time steps
DO i=-2,0
	Uall(:,:,:,i)=Utmp
ENDDO

END SUBROUTINE





