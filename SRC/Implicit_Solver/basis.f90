! This file is part of FLEXI.
! The original copyright notice follows:
!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. IF not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================
!==================================================================================================================================
!> Routines to provide and evaluate basis FUNCTION coefficients, or provide fast interpolation coefficients
!==================================================================================================================================
MODULE MOD_Basis
! MODULES
USE MOD_Global, ONLY:RP
IMPLICIT NONE
PRIVATE
SAVE
!----------------------------------------------------------------------------------------------------------------------------------
! MOD_Global VARIABLES
!----------------------------------------------------------------------------------------------------------------------------------

PUBLIC::ChebyLob2Legendre
PUBLIC::ChebyGaussLobNodesAndWeights
PUBLIC::InitializeVandermonde
PUBLIC::LegendreGaussNodesAndWeights
PUBLIC::BarycentricWeights
!==================================================================================================================================

CONTAINS

!===================================================================================================================================
!> Build a 1D Vandermonde matrix using the Lagrange basis functions of degree
!> N_In, evaluated at the interpolation points xi_Out
!===================================================================================================================================
PURE SUBROUTINE InitializeVandermonde(N_In,N_Out,wBary_In,xi_In,xi_Out,Vdm)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_In                           !< (IN)  input polynomial degree
INTEGER,INTENT(IN) :: N_Out                          !< (IN)  output polynomial degree
REAL(RP),INTENT(IN)    :: xi_In(0:N_In)              !< (IN)  input nodal positions [-1,1]
REAL(RP),INTENT(IN)    :: xi_Out(0:N_Out)            !< (IN)  outout nodal positions [-1,1]
REAL(RP),INTENT(IN)    :: wBary_In(0:N_In)           !< (IN)  input interpolation weights
REAL(RP),INTENT(OUT)   :: Vdm(0:N_Out,0:N_In)        !< (OUT) nodal Vandermonde from N_In to N_out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iXi
!==================================================================================================================================
DO iXi=0,N_Out
  CALL LagrangeInterpolationPolys(xi_Out(iXi),N_In,xi_In,wBary_In,Vdm(iXi,:)) !l(0:N_In)
END DO
END SUBROUTINE InitializeVandermonde


!==================================================================================================================================
!> Compute Chebychev-Gauss-Lobatto nodes and integration weights (algorithm 27, Kopriva book)
!==================================================================================================================================
PURE SUBROUTINE ChebyGaussLobNodesAndWeights(N_in,xGP,wGP)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in         !< polynomial degree, (N_in+1) CLpoints
REAL(RP),INTENT(OUT)          :: xGP(0:N_in)  !< Gauss point positions for the reference interval [-1,1]
REAL(RP),INTENT(OUT),OPTIONAL :: wGP(0:N_in)  !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP
!==================================================================================================================================
DO iGP=0,N_in
  xGP(iGP)=-COS(iGP/DBLE(N_in)*ACOS(-1.))
END DO
IF(PRESENT(wGP))THEN
  DO iGP=0,N_in
    wGP(iGP)=ACOS(-1.)/DBLE(N_in)
  END DO
  wGP(0)=wGP(0)*0.5
  wGP(N_in)=wGP(N_in)*0.5
END IF
END SUBROUTINE ChebyGaussLobNodesAndWeights


!==================================================================================================================================
!> @brief Compute Legendre-Gauss nodes and integration weights (algorithm 23, Kopriva book)
!>
!> Starting with Chebychev point positions, a Newton method is used to find the roots
!> of the Legendre Polynomial L_(N_in+1), which are the positions of Gausspoints
!> uses LegendrePolynomialAndDerivative SUBROUTINE
!==================================================================================================================================
SUBROUTINE LegendreGaussNodesAndWeights(N_in,xGP,wGP)
!MODULES
!USE MOD_Globals
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)        :: N_in              !< polynomial degree, (N_in+1) Gausspoints
REAL(RP),INTENT(OUT)          :: xGP(0:N_in)       !< Gauss point positions for the reference interval [-1,1]
REAL(RP),INTENT(OUT),OPTIONAL :: wGP(0:N_in)       !< Gauss point weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER                   :: nIter = 10        ! max. number of newton iterations
REAL(RP)                      :: Tol   = 1.E-15    ! tolerance for Newton iteration: TODO: use variable tolerance here!
INTEGER                   :: iGP,iter
REAL(RP)                      :: L_Np1,Lder_Np1    ! L_{N_in+1},Lder_{N_in+1}
REAL(RP)                      :: DX                ! Newton step
REAL(RP)                      :: cheb_tmp          ! temporary variable for evaluation of chebychev node positions
!==================================================================================================================================
IF(N_in .EQ. 0) THEN
  xGP=0.
  IF(PRESENT(wGP))wGP=2.
  RETURN
ELSEIF(N_in.EQ.1)THEN
  xGP(0)=-sqrt(1./3.)
  xGP(N_in)=-xGP(0)
  IF(PRESENT(wGP))wGP=1.
  RETURN
ELSE ! N_in>1
  cheb_tmp=2.*atan(1.)/DBLE(N_in+1) ! PI/(2N+2)
  DO iGP=0,(N_in+1)/2-1 !since points are symmetric, only left side is computed
    xGP(iGP)=-cos(cheb_tmp*DBLE(2*iGP+1)) !initial guess
    ! Newton iteration
    DO iter=0,nIter
      CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
      DX=-L_Np1/Lder_Np1
      xGP(iGP)=xGP(iGP)+DX
      IF(abs(DX).LT.Tol*abs(xGP(iGP))) EXIT
    END DO ! iter
    IF(iter.GT.nIter) THEN
      WRITE(*,*) 'maximum iteration steps >10 in Newton iteration for Legendre Gausspoint'
      xGP(iGP)=-cos(cheb_tmp*DBLE(2*iGP+1)) !initial guess
      ! Newton iteration
      DO iter=0,nIter
        WRITE(*,*)iter,xGP(iGP)    !DEBUG
        CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
        DX=-L_Np1/Lder_Np1
        xGP(iGP)=xGP(iGP)+DX
        IF(abs(DX).LT.Tol*abs(xGP(iGP))) EXIT
      END DO !iter
      WRITE(*,*) 'ERROR: Legendre Gauss nodes could not be computed up to desired precision. Code stopped!'
      STOP
    END IF ! (iter.GT.nIter)
    CALL LegendrePolynomialAndDerivative(N_in+1,xGP(iGP),L_Np1,Lder_Np1)
    xGP(N_in-iGP)=-xGP(iGP)
    IF(PRESENT(wGP))THEN
      !wGP(iGP)=2./((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1) !IF Legendre not normalized
      wGP(iGP)=(2.*N_in+3)/((1.-xGP(iGP)*xGP(iGP))*Lder_Np1*Lder_Np1)
      wGP(N_in-iGP)=wGP(iGP)
    END IF
  END DO !iGP
END IF ! N_in
IF(mod(N_in,2) .EQ. 0) THEN
  xGP(N_in/2)=0.
  CALL LegendrePolynomialAndDerivative(N_in+1,xGP(N_in/2),L_Np1,Lder_Np1)
  !IF(PRESENT(wGP))wGP(N_in/2)=2./(Lder_Np1*Lder_Np1) !IF Legendre not normalized
  IF(PRESENT(wGP))wGP(N_in/2)=(2.*N_in+3)/(Lder_Np1*Lder_Np1)
END IF ! (mod(N_in,2) .EQ. 0)
END SUBROUTINE LegendreGaussNodesAndWeights

!===================================================================================================================================
!> Evaluate the Legendre polynomial L_N and its derivative at position x[-1,1]
!> recursive algorithm using the N_in-1 N_in-2 Legendre polynomials
!> algorithm 22, Kopriva book
!===================================================================================================================================
ELEMENTAL SUBROUTINE LegendrePolynomialAndDerivative(N_in,x,L,Lder)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in   !< (IN)  polynomial degree, (n+1) CLpoints
REAL(RP),INTENT(IN)    :: x      !< (IN)  coordinate value in the interval [-1,1]
REAL(RP),INTENT(OUT)   :: L      !< (OUT) Legedre polynomial evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
REAL(RP),INTENT(OUT)   :: Lder   !< (OUT) Legedre polynomial deriv. evaluated at \f$ \xi: L_N(\xi), \partial/\partial\xi L_N(\xi) \f$
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER     :: iLegendre
REAL(RP)    :: L_Nm1,L_Nm2 ! L_{N_in-2},L_{N_in-1}
REAL(RP)    :: Lder_Nm1,Lder_Nm2 ! Lder_{N_in-2},Lder_{N_in-1}
!==================================================================================================================================
IF(N_in .EQ. 0)THEN
  L=1.
  Lder=0.
ELSEIF(N_in .EQ. 1) THEN
  L=x
  Lder=1.
ELSE ! N_in > 1
  L_Nm2=1.
  L_Nm1=x
  Lder_Nm2=0.
  Lder_Nm1=1.
  DO iLegendre=2,N_in
    L=(REAL(2*iLegendre-1)*x*L_Nm1 - REAL(iLegendre-1)*L_Nm2)/REAL(iLegendre)
    Lder=Lder_Nm2 + REAL(2*iLegendre-1)*L_Nm1
    L_Nm2=L_Nm1
    L_Nm1=L
    Lder_Nm2=Lder_Nm1
    Lder_Nm1=Lder
  END DO !iLegendre=2,N_in
END IF ! N_in
!normalize
L=L*SQRT(REAL(N_in)+0.5)
Lder=Lder*SQRT(REAL(N_in)+0.5)
END SUBROUTINE LegendrePolynomialAndDerivative




!==================================================================================================================================
!> Computes barycentric (interpolation) weights for interpolation polynomial given by set of nodes. (Algorithm 30, Kopriva book)
!==================================================================================================================================
PURE SUBROUTINE BarycentricWeights(N_in,xGP,wBary)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN) :: N_in               !< polynomial degree
REAL(RP),INTENT(IN)    :: xGP(0:N_in)        !< Gauss point positions for the reference interval [-1,1]
REAL(RP),INTENT(OUT)   :: wBary(0:N_in)      !< barycentric weights
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP,jGP
!==================================================================================================================================
wBary(:)=1.
DO iGP=1,N_in
  DO jGP=0,iGP-1
    wBary(jGP)=wBary(jGP)*(xGP(jGP)-xGP(iGP))
    wBary(iGP)=wBary(iGP)*(xGP(iGP)-xGP(jGP))
  END DO ! jGP
END DO ! iGP
wBary(:)=1./wBary(:)
END SUBROUTINE BarycentricWeights

!============================================================================================================================
!> Computes all Lagrange functions evaluated at position x in [-1,1]
!> For details see paper Barycentric Lagrange Interpolation by Berrut and Trefethen (SIAM 2004)
!> Uses FUNCTION ALMOSTEQUAL
!> Algorithm 34, Kopriva book
!============================================================================================================================
PURE SUBROUTINE LagrangeInterpolationPolys(x,N_in,xGP,wBary,L)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(RP), INTENT(IN)   :: x                !< Coordinate
INTEGER,INTENT(IN) :: N_in             !< polynomial degree
REAL(RP),INTENT(IN)    :: xGP(0:N_in)      !< Gauss point positions for the reference interval [-1,1]
REAL(RP),INTENT(IN)    :: wBary(0:N_in)    !< Barycentric weights
REAL(RP),INTENT(OUT)   :: L(0:N_in)        !< Lagrange basis functions evaluated at x
!----------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER            :: iGP
LOGICAL            :: xEqualGP         ! is x equal to a Gauss Point
REAL(RP)               :: DummySum
!============================================================================================================================
xEqualGP=.FALSE.
DO iGP=0,N_in
  L(iGP)=0.
  IF(ALMOSTEQUAL(x,xGP(iGP))) THEN
    L(iGP)=1.
    xEqualGP=.TRUE.
  END IF
END DO

! IF x is equal to a Gauss point, L=(0,....,1,....0)
IF(xEqualGP) RETURN
DummySum=0.
DO iGP=0, N_in
  L(iGP)=wBary(iGP)/(x-xGP(iGP))
  DummySum=DummySum+L(iGP)
END DO

DO iGP=0,N_in
  L(iGP)=L(iGP)/DummySum
END DO
END SUBROUTINE LagrangeInterpolationPolys



!==================================================================================================================================
!> Interpolate an array of order NIn at 2D Chebyshev-Gauss-Lobatto points to that
!> of order NOut at Legendre points using the 1D Vandermonde matrix Vdm.
!> IF X_Out is present, the output will be written to X_Out, otherwise X_In will be overwritten.
!==================================================================================================================================
SUBROUTINE ChebyLob2Legendre(NIn,NOut,VDM,X_In,X_Out)
! MODULES
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
INTEGER,INTENT(IN)          :: NIn                   !< Input polynomial degree, no. of points = NIn+1
INTEGER,INTENT(IN)          :: NOut                  !< Output polynomial degree, no. of points = NOut+1
REAL(RP),INTENT(INOUT)        :: X_In( 0:NIn ,0:NIn )    !< Input field, dimensions must match Dim1,NIn
REAL(RP),INTENT(OUT),OPTIONAL :: X_Out(0:NOut,0:NOut)  !< Output field, dimensions must match Dim1,NOut, OPTIONAL
REAL(RP),INTENT(IN)           :: Vdm(  0:NOut,0:NIn )     !< 1D Vandermonde In -> Out
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
INTEGER             :: iI,jI,iO,jO
REAL(RP)              :: X_Buf1(0:NOut,0:NIn)       ! first intermediate results from 1D interpolations
!==================================================================================================================================
X_buf1=0.
! first direction iI
DO jI=0,NIn; DO iI=0,NIn
  DO iO=0,NOut
    X_Buf1(iO,jI)=X_Buf1(iO,jI)+Vdm(iO,iI)*X_In(iI,jI)
  END DO
END DO; END DO
IF (PRESENT(X_OUT)) THEN
  X_OUT=0.
  ! SECOND DIRECTION JI
  DO JI=0,NIN; DO JO=0,NOUT
    DO IO=0,NOUT
      X_OUT(IO,JO)=X_OUT(IO,JO)+VDM(JO,JI)*X_BUF1(IO,JI)
    END DO
  END DO; END DO
ELSE
  X_IN=0.
  ! SECOND DIRECTION JI
  DO JI=0,NIN; DO JO=0,NOUT
    DO IO=0,NOUT
      X_IN(IO,JO)=X_IN(IO,JO)+VDM(JO,JI)*X_BUF1(IO,JI)
    END DO
  END DO; END DO
END IF
END SUBROUTINE ChebyLob2Legendre










!==================================================================================================================================
!> Determines IF two REAL(RP) numbers are equal up to a specified tolerance (=PP_RealTolerance, normaly set to machine precision)
!> Takes into account that x,y are located in-between [-1;1] for additional accuracy
!> Based on Algorithm 139, Kopriva
!==================================================================================================================================
ELEMENTAL FUNCTION ALMOSTEQUAL(x,y)
IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
REAL(RP),INTENT(IN) :: x                !< (IN)  first scalar to be compared
REAL(RP),INTENT(IN) :: y                !< (IN)  second scalar to be compared
LOGICAL         :: AlmostEqual      !< (OUT) TRUE IF |x-y| < 2*PP_RealTolerance
!----------------------------------------------------------------------------------------------------------------------------------
! LOCAL VARIABLES
!==================================================================================================================================
AlmostEqual=.FALSE.
IF((x.EQ.0.).OR.(y.EQ.0.)) THEN
  IF(ABS(x-y).LE.2.*1.D-15) AlmostEqual=.TRUE.
ELSE ! x, y not zero
  IF((ABS(x-y).LE.2.*1.D-15).AND.((ABS(x-y).LE. 1.D-15*ABS(y)))) AlmostEqual=.TRUE.
END IF ! x,y zero
END FUNCTION ALMOSTEQUAL






END MODULE MOD_Basis