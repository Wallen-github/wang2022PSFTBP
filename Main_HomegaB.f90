    PROGRAM MAIN
    USE GlobalDef
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER(4) :: time_begin, time_end
    REAL(8) :: OmegaB,OD,TT,X(8),T0,X0(8),Xp(8),flag
    REAL(8) :: MOMENTUM, ENERGY
    INTEGER(4) :: Num,I
    REAL(8) :: OmegaBSta,OmegaBEnd,OmegaBStp

104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    CALL system_clock(time_begin)
    CALL Instantiate_Global_Values()
    
    OPEN(13,file='POT.DAT', status='REPLACE')
    
    Num = 500
    OmegaBSta = 1.55D0!1.0123D0!2.1998D0
    OmegaBEnd = 1.47D0!2.1998D0
    OmegaBStp = (OmegaBEnd-OmegaBSta)/Num
    
    X0(1)=1.0D0
	X0(2)=0D0*PI/2D0
	X0(3)=0D0*PI/2D0
	X0(4)=0D0*PI/2D0
    X0(5)=0D0
	X0(6)=0D0
	X0(7)=OmegaBSta-1D0
	X0(8)=1D0
    
    DO I = 0,Num
        
        OmegaB = OmegaBSta + I*OmegaBStp
        OD=1D0-OmegaB
	    T0=PI2/DABS(OD)
        WRITE(*,*) OmegaB,T0
    
        CALL PeriodT(T0,X0,Tp,Xp)
        CALL EigMon(Tp,Xp,flag)
        !MOMENTUM=(IzA+IzB+MASS*Xp(1)**2)*Xp(8)+IzB*Xp(7)+MASS*Xp(1)**2*Xp(6)
        MOMENTUM = Mass + IzA*Xp(8) + IzB*(Xp(7)+Xp(8))
        ENERGY=-MASS*(1D0/Xp(1)+1D0/Xp(1)**3*(A1+A2*COS(2D0*Xp(2))+A3*COS(2D0*(Xp(2)-Xp(3)))))+ &
    	    & MASS/2D0*((Xp(1)*(Xp(6)+Xp(8)))**2+Xp(5)**2)+IzA/2D0*Xp(8)**2+IzB/2D0*(Xp(7)+Xp(8))**2
        
        WRITE(13,104) OmegaB,MOMENTUM,ENERGY,flag
        !CALL PeriodS(T0,X0,TT,X)
        
    
    ENDDO
    
    CALL system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    END
    
    SUBROUTINE PeriodT(T0,X0,TT,X)
    USE GlobalDef
    IMPLICIT NONE 
	REAL*8 X0(8),X(8),Y(8),A(8,8),C1(4,4),C2(4,4),BB1(4),BB2(4)
    REAL*8 OD,TT,T0,DD,HH,EE
    REAL*8 EPS,FACTOR
    INTEGER I,FLAG,J

    X = X0
    TT = T0
	EPS=1D0
    BB2=0D0
	DO WHILE(EPS.GT.1D-012)
		FACTOR=1D0!+1D0*EPS
		X(1)=X(1)+BB2(1)/FACTOR
	  	X(6)=X(6)+BB2(2)/FACTOR
	  	X(7)=X(7)+BB2(3)/FACTOR
	  	X(8)=X(8)+BB2(4)/FACTOR
	  	
          
    	CALL MONODROMY3(X,TT/2D0,Y,A)     !  STM METHOD
! 		CALL MONODROMY4(X,TT,Y,A)     !  DIFFERENCE METHOD    
		C1(1,1)=A(2,1)
	  	C1(1,2)=A(2,6)
	 	C1(1,3)=A(2,7)
	    C1(1,4)=A(2,8)

	  	C1(2,1)=A(3,1)	  	
	  	C1(2,2)=A(3,6)
	  	C1(2,3)=A(3,7)
	  	C1(2,4)=A(3,8)

	  	C1(4,1)=A(5,1)
	  	C1(4,2)=A(5,6)
	  	C1(4,3)=A(5,7)
	  	C1(4,4)=A(5,8)

	  	BB1(1)=X(2)-Y(2)
        BB1(2)=X(3)-Y(3)+PI     ! omegaB>1的周期轨道族
	  	!BB1(2)=X(3)-Y(3)-PI    ! omegaB<1的周期轨道族
	  	BB1(4)=-Y(5)

	  	CALL MONODROMY3(X,PI2,Y,A)

	  	C1(3,1)=A(4,1)	  	
	  	C1(3,2)=A(4,6)
	  	C1(3,3)=A(4,7)
	  	C1(3,4)=A(4,8)  	

	  	BB1(3)=-Y(4)+PI2

	  	CALL ANTIMATRIX(C1,C2,4,FLAG)
        BB2=MATMUL(C2,BB1)
	  	EPS=DSQRT(BB1(1)**2+BB1(2)**2+BB1(3)**2+BB1(4)**2)

	  	!PRINT*, EPS
! 	  	PAUSE
    ENDDO
    WRITE(*,*) 'Success !'
	!DO I=1,8
	!WRITE(*,'(E25.14)'), X(I)
	!ENDDO

    END
    
    SUBROUTINE PeriodS(T0,X0,TT,X)
    USE GlobalDef
    IMPLICIT NONE 
	REAL*8 X0(8),X(8),Y(8),A(8,8),C1(4,4),C2(4,4),B1(4),B2(4)
    REAL*8 OD,TT,T0,DD,HH,EE,OmegaB
    REAL*8 EPS,FACTOR
    INTEGER I,FLAG,J

    X = X0
    TT = T0
	EPS=1D0
    B2=0D0
	DO WHILE(EPS.GT.1D-012)
		FACTOR=1D0+1D0*EPS
		TT=TT+B2(1)/FACTOR
	  	X(6)=X(6)+B2(2)/FACTOR
	  	X(7)=X(7)+B2(3)/FACTOR
	  	X(8)=X(8)+B2(4)/FACTOR
	  	
          
    	CALL MONODROMY3(X,TT/2D0,Y,A)     !  STM METHOD
! 		CALL MONODROMY2(X,TT,Y,A)     !  DIFFERENCE METHOD    
		C1(1,1)=Y(6)
	  	C1(1,2)=A(2,6)
	 	C1(1,3)=A(2,7)
	    C1(1,4)=A(2,8)

	  	C1(2,1)=Y(7)	  	
	  	C1(2,2)=A(3,6)
	  	C1(2,3)=A(3,7)
	  	C1(2,4)=A(3,8)

	  	C1(4,1)=+Y(1)*(Y(6)+Y(8))**2-(1D0/Y(1)**2+3D0/Y(1)**4*&
  			(A1+A2*COS(Y(2))+A3*COS(Y(2)-Y(3))))
	  	C1(4,2)=A(5,6)
	  	C1(4,3)=A(5,7)
	  	C1(4,4)=A(5,8)

	  	B1(1)=X(2)-Y(2)
	  	B1(2)=X(3)-Y(3)+PI
	  	B1(4)=-Y(5)

	  	CALL MONODROMY3(X,PI2,Y,A)

	  	C1(3,1)=0D0	  	
	  	C1(3,2)=A(4,6)
	  	C1(3,3)=A(4,7)
	  	C1(3,4)=A(4,8)  	

	  	B1(3)=-Y(4)+PI2

	  	CALL ANTIMATRIX(C1,C2,4,FLAG)
        !B2=MATMUL(C2,B1)
        CALL MATRIXVECTOR(C2,B1,B2,4,4)
	  	EPS=DSQRT(B1(1)**2+B1(2)**2+B1(3)**2+B1(4)**2)

	  	PRINT*, EPS
! 	  	PAUSE
	ENDDO

	OmegaB=PI2/TT+1D0

	WRITE(*,'(A10,F25.15,A10)') '----------',OmegaB,'----------'
	DO I=1,8
		WRITE(*,'(D25.14)') X(I)
	ENDDO

    END
    
    SUBROUTINE EigMon(TT,X,flag)
    
    IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION WR(8),WI(8)
    REAL(8) :: X(8),TT,A(8,8),XT(8),Y(8),flag
    Integer(4) :: I
    REAL(8) :: DIS(8)
    
    CALL MONODROMY3(X,TT,Y,A)
    CALL EIGENVALUE(8,A,WR,WI)
    
    DIS = DSQRT(WR**2D0+WI**2D0)
    
    flag = 1D0
    DO I = 1,8
        IF (ABS(DIS(I)-1D0)>1D-10) THEN
            flag = -1D0
        ENDIF
        !WRITE(*,*) WR(I),WI(I),DIS(I)
    ENDDO
    
    END
    
    
    SUBROUTINE MONODROMY3(X,TT,Y,A)
    IMPLICIT NONE
    EXTERNAL :: YHCM
    REAL*8 PI,PI2
    REAL*8 X(8),Y(8),A(8,8)
    REAL*8 Y1(72)
	PARAMETER(PI=3.14159265358979323846D0,PI2=PI*2D0)
    REAL*8 TT,DD,HH,T0,EE
    INTEGER I
    Y1=0D0
    Y1(1:8)=X
    Y1(9)=1D0
    Y1(18)=1D0
    Y1(27)=1D0
    Y1(36)=1D0
    Y1(45)=1D0
    Y1(54)=1D0
    Y1(63)=1D0
    Y1(72)=1D0
	DD=1D-012
	T0=0D0
	HH=1D-006
	DO WHILE(T0.LT.TT)
	  	CALL RKF78(YHCM,HH,T0,Y1,EE,DD,72)
	ENDDO
	CALL RKF78(YHCM,TT-T0,T0,Y1,EE,1D0,72)   
    Y=Y1(1:8)
    DO I=1,8
        A(I,:)=Y1(I*8+1:(I+1)*8)
    ENDDO
	RETURN
	END