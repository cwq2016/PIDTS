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
!> Routines to update the solutions at boundies.
!=================================================================================================================================
MODULE MOD_BC
USE MOD_Global
IMPLICIT NONE
PRIVATE
SAVE
!---------------------------------------------------------------------------------------------------------------------------------
! PUBLIC VARIABLES
PUBLIC::bc_adj
PUBLIC::Init_BC
!---------------------------------------------------------------------------------------------------------------------------------
! PRIVATE VARIABLES
REAL(RP) :: Ainvx(0:1,0:1),Ainvy(0:1,0:1)   !< the inverse of coefficient matrix 
                                            !< for boundary equations. 
											!< See (W.Q. Chen et al. 2018) 
REAL(RP),ALLOCATABLE,DIMENSION(:,:,:  ):: ResiBoundx,ResiBoundy
									        !< The residual of boundary equations. 
!=================================================================================================================================

CONTAINS
SUBROUTINE bc_adj
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER :: xid(2),yid(2) ! The index of boundary in x and y direction
INTEGER :: i,j,n
REAL(RP):: source(0:1),xx(0:1)
!=================================================================================================================================


xid = (/0,Nx/)
yid = (/0,Ny/)
CALL Residual_Interface

!-------------- update boundary pressure ----------------------
n=1 ! pressure
! bounday in x direction
! dp_dx+d2v_dxdy=0
DO j=1,Ny-1
	source=-ResiBoundx(0:1,j,n)
	CALL SolveBC(Ainvx,source,xx)
	Mid(xid,j,n)=Mid(xid,j,n)+xx;
ENDDO

! bounday in y direction
! dp_dy+d2u_dxdy+T=0
DO i=0,Nx         
	source=-ResiBoundy(i,0:1,n)
	CALL SolveBC(Ainvy,source,xx)
	Mid(i,yid,n)=Mid(i,yid,n)+xx;
	!******************************
ENDDO

!-------------- update boundary temperature -------------------
n=4 ! temperature
! bounday in y direction 
! dT_dy=0
DO i=0,Nx
	source=-ResiBoundy(i,0:1,n)
	CALL SolveBC(Ainvy,source,xx)
	Mid(i,yid,n)=Mid(i,yid,n)+xx;
	!******************************
ENDDO
          
! bounday in x direction
! Dirichlet boundary condition
Mid(0 ,0:Ny,4)= -0.5
Mid(Nx,0:Ny,4)=  0.5+Amplitude*SIN(2.d0*PI*Time/Period)

END SUBROUTINE 


!=================================================================================================================================
!> Solve boundary equations
!=================================================================================================================================
SUBROUTINE SolveBC(Ainv,b,x)
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(RP),INTENT(IN )::Ainv(0:1,0:1),b(0:1)
REAL(RP),INTENT(OUT)::x(0:1)
!=================================================================================================================================
x(0)= Ainv(0,0)*b(0) + Ainv(0,1)*b(1)
x(1)= Ainv(1,0)*b(0) + Ainv(1,1)*b(1)
END SUBROUTINE



!=================================================================================================================================
!> Calculate the residual of boundary equations
!=================================================================================================================================
SUBROUTINE Residual_Interface
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER ::loc
INTEGER i,j,n
REAL(RP)::ux(0:Nx),vy(0:Ny)
!=================================================================================================================================

ResiBoundx=0.
ResiBoundy=0.

!-------------- pressure boundary residual ----------------------
n=1
! Add velocity derivatives to the pressure bounday residual
DO loc=1,4
	CALL pbound(loc)
ENDDO

! Add temperature term to the pressure bounday residual
DO i=0,Nx
	ResiBoundy(i,0,1)=ResiBoundy(i,0,1)+Mid(i,0 ,4)
	ResiBoundy(i,1,1)=ResiBoundy(i,1,1)+Mid(i,Ny,4)
ENDDO

!-------------- temperature boundary residual ----------------------
n=4
! dth_dy=0 at top and lower sidewall
DO i=0,Nx
	ResiBoundy(i,0,4)=sum( DY(0 ,0:Ny) * Mid(i,0:Ny,4) )
	ResiBoundy(i,1,4)=sum( DY(Ny,0:Ny) * Mid(i,0:Ny,4) )
ENDDO
END SUBROUTINE

!=================================================================================================================================
!> Calculate part of  pressure boundary residual,
!> namely, d2v_dxdy and d2u_dxdy
!=================================================================================================================================
SUBROUTINE pbound(loc)
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER i,j,n
INTEGER,intent(in)::loc
REAL(RP)::ux(0:Nx),vy(0:Ny)
REAL(RP):: alpha=1.D0,beta=0.D0
!=================================================================================================================================

!loc is the index of sidewall: 2--left£¬1--right£¬ 3--top£¬ 4--bottom

n=1
IF (loc==1)then
	! Right boundary--dp_dx+d2v_dxdy
	DO i=0,Nx
		!vy
		CALL dgemv('n',Ny+1,Ny+1,alpha,DY(0:Ny,0:Ny),Ny+1,Mid(i,0:Ny,3),1,beta,vy(0:Ny),1)
		
		ResiBoundx(0,0:Ny,n)=ResiBoundx(0,0:Ny,n)+DX(0,i)*(Kappa(2)*vy(0:Ny)+Mid(i,0:Ny,n))
	ENDDO
elseif(loc==2)then
	! Right boundary--dp_dx+d2v_dxdy
	DO i=0,Nx
		!vy
		CALL dgemv('n',Ny+1,Ny+1,alpha,DY(0:Ny,0:Ny),Ny+1,Mid(i,0:Ny,3),1,beta,vy(0:Ny),1)
		
		ResiBoundx(1,0:Ny,n)=ResiBoundx(1,0:Ny,n)+DX(Nx,i)*(Kappa(2)*vy(0:Ny)+Mid(i,0:Ny,n))
	ENDDO

elseif(loc==3)then
	!Top side wall--dp_dy+d2u_dxdy
	DO j=0,Ny
		!ux
		CALL dgemv('n',Nx+1,Nx+1,alpha,DX(0:Nx,0:Nx),Nx+1,Mid(0:Nx,j,2),1,beta,ux(0:Nx),1)
		
		ResiBoundy(0:Nx,0,n)=ResiBoundy(0:Nx,0,n)+DY(0,j)*(Kappa(3)*ux(0:Nx)+Mid(0:Nx,j,n))
	ENDDO
elseif(loc==4)then
	!Bottom side wall--dp_dy+d2u_dxdy
	DO j=0,Ny
		!ux
		CALL dgemv('n',Nx+1,Nx+1,alpha,DX(0:Nx,0:Nx),Nx+1,Mid(0:Nx,j,2),1,beta,ux(0:Nx),1)
		
		ResiBoundy(0:Nx,1,n)=ResiBoundy(0:Nx,1,n)+DY(Ny,j)*(Kappa(3)*ux(0:Nx)+Mid(0:Nx,j,n))
	ENDDO
endif
END SUBROUTINE

!=================================================================================================================================
!> Precompute the inverse of coefficient matrix for boundary equations
!=================================================================================================================================
SUBROUTINE Init_BC
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
REAL(RP)::Ax(   0:1,0:1),Ay(   0:1,0:1)
REAL(RP)::detAx,detAy
!=================================================================================================================================

! dphi_dx -->  Ax
Ax(0,0)=DX(0 ,0 )
Ax(0,1)=DX(0 ,Nx)
Ax(1,0)=DX(Nx,0 )
Ax(1,1)=DX(Nx,Nx)

! Derive the inverse of Ax 
detAx = Ax(0,0)*Ax(1,1) - Ax(1,0)*Ax(0,1)
Ainvx(0,0)= Ax(1,1)
Ainvx(0,1)=-Ax(0,1)
Ainvx(1,0)=-Ax(1,0)
Ainvx(1,1)= Ax(0,0)
Ainvx = Ainvx/detAx

! dphi_dy -->  Ay
Ay(0,0)=DY(0 ,0 )
Ay(0,1)=DY(0 ,Ny)
Ay(1,0)=DY(Ny,0 )
Ay(1,1)=DY(Ny,Ny)

! Derive the inverse of Ay 
detAy = Ay(0,0)*Ay(1,1) - Ay(1,0)*Ay(0,1)
Ainvy(0,0)= Ay(1,1)
Ainvy(0,1)=-Ay(0,1)
Ainvy(1,0)=-Ay(1,0)
Ainvy(1,1)= Ay(0,0)
Ainvy = Ainvy/detAy

! Allocate space for boundary residual
ALLOCATE(ResiBoundx(0:1,0:Ny,nVar),ResiBoundy(0:Nx,0:1,nVar))	
END SUBROUTINE

END MODULE
