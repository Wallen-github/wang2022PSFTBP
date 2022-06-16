    PROGRAM MAIN
    USE GlobalDef
    IMPLICIT REAL(8) (A-H,O-Z)
    INTEGER(4) :: time_begin, time_end
    REAL(8) :: OmegaB,OD,T0,X0(8),Xp(8),flag, Tmid
    REAL(8) :: MOMENTUM, ENERGY,Index, ss


104 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
105 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
110 format(1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X, E25.18,1X,E25.18,1X,E25.18,1X,E25.18,1X,E25.18)
    
    CALL system_clock(time_begin)
    CALL Instantiate_Global_Values()
    
    OPEN(13,file='PO.DAT', status='REPLACE')
    
    OmegaB =  1.55D0
    !OmegaB = 1.44123!1999KW4 1.3~1.5共振段初始值
    !OmegaB = 1.4523 !Didymos 1.3~1.5共振段初始值
    OD=1D0-OmegaB
	T0=PI2/DABS(OD)
    
    X0(1)=1D0
	X0(2)=0D0*PI/2D0
	X0(3)=0D0*PI/2D0
	X0(4)=0D0*PI/2D0
    X0(5)=0D0
	X0(6)=0D0
	X0(7)=OmegaB-1D0
	X0(8)=1D0
    
    Index = 1D0
    
    write(*,*) ' '
    WRITE(*,*) '======================== Peroid Orbit Family =========================='
    write(*,*) ' '
    
    DO I = 0,5000
        
        IF (Index==1D0) THEN
            OD=1D0-OmegaB
	        T0=PI2/DABS(OD)
    
            CALL PeriodT(T0,X0,Tp,Xp,Index)
            OmegaB = SIGN(1D0,OmegaB-1D0)*PI2/Tp+1D0        ! 智能判断
            ss = ABS(X0(1)-Xp(1))/X0(1) + 1D0
            !OmegaB = OmegaB + 3.123D-5/ss
            OmegaB = OmegaB - 3.123D-5/ss
            IF (Index<0D0) THEN
                Index = 2D0
                T0 = Tmid
                !CYCLE
                write(*,*) ' '
                WRITE(*,*) '======================== Stop Point =========================='
                write(*,*) ' '
                WRITE(*,'(E25.18)') T0
                WRITE(*,104) X0
                PAUSE
                CYCLE
                !EXIT
            ENDIF
            CALL EigMon(Tp,Xp,flag)
            MOMENTUM=(IzA+IzB+MASS*Xp(1)**2)*Xp(8)+IzB*Xp(7)+MASS*Xp(1)**2*Xp(6)
            ENERGY=-MASS*(1D0/Xp(1)+1D0/Xp(1)**3*(A1+A2*COS(2D0*Xp(2))+A3*COS(2D0*(Xp(2)-Xp(3)))))+ &
    	        & MASS/2D0*((Xp(1)*(Xp(6)+Xp(8)))**2+Xp(5)**2)+IzA/2D0*Xp(8)**2+IzB/2D0*(Xp(7)+Xp(8))**2
            WRITE(13,105) OmegaB,MOMENTUM,ENERGY,flag,Xp(1)
            WRITE(*,*) OmegaB,Tp,flag
            WRITE(*,104) Xp
            
            X0 = Xp
            Tmid = Tp
            !T0 = Tp + 1D-1
            
        ELSEIF (Index==2D0) THEN
            
            CALL PeriodS(T0,X0,Tp,Xp,Index)
            OmegaB = SIGN(1D0,OmegaB-1D0)*PI2/Tp+1D0        ! 智能判断
            !OmegaB = PI2/Tp+1D0        ! omegaB>1的周期轨道族
            !OmegaB = -PI2/Tp+1D0         ! omegaB>1的周期轨道族
            IF (Index<0D0) THEN
                !Index = 1D0
                write(*,*) ' '
                WRITE(*,*) '======================== Stop Point =========================='
                write(*,*) ' '
                WRITE(*,'(E25.18)') T0
                WRITE(*,104) X0
                X0(1) = X0(1) + 1.22D-4/ss
                ss = ss+1D1
                X0(1) = X0(1) - 1.22D-4/ss
                PAUSE
                CYCLE
                IF (ss>1D2) THEN
                    EXIT
                ENDIF
            ENDIF
            CALL EigMon(Tp,Xp,flag)
            MOMENTUM=(IzA+IzB+MASS*Xp(1)**2)*Xp(8)+IzB*Xp(7)+MASS*Xp(1)**2*Xp(6)
            ENERGY=-MASS*(1D0/Xp(1)+1D0/Xp(1)**3*(A1+A2*COS(2D0*Xp(2))+A3*COS(2D0*(Xp(2)-Xp(3)))))+ &
    	        & MASS/2D0*((Xp(1)*(Xp(6)+Xp(8)))**2+Xp(5)**2)+IzA/2D0*Xp(8)**2+IzB/2D0*(Xp(7)+Xp(8))**2
            WRITE(13,105) OmegaB,MOMENTUM,ENERGY,flag,Xp(1)
            WRITE(*,*) OmegaB,T0,flag
            WRITE(*,104) Xp
            ss = ABS(Tp-T0)/T0+1D0
            X0 = Xp
            T0 = Tp
            X0(1) = X0(1) - 1.22D-4/ss
        ENDIF
        
    ENDDO
    
    write(*,*) ' '
    WRITE(*,*) '======================== End Point =========================='
    write(*,*) ' '
    WRITE(*,'(E25.18)') Tp
    WRITE(*,104) Xp
    
    CALL system_clock(time_end)
    PRINT *, 'Time of operation was ', (time_end - time_begin)/1D3, ' seconds'
    END
    