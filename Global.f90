
    MODULE GlobalDef
    
    REAL(8),PARAMETER :: GG = 6.67D-11, AU=1.4959787D011,SecOfYear=31536D3
    REAL(8),PARAMETER :: PI=3.14159265358979323846D0,PI2=PI*2D0
    REAL(8) :: IzA,IzB,A1,A2,A3,AlphaA,AlphaB,Mass,MassRatio
    REAL(8) :: B1,B2,B3,B4,B5,B6,B7
    REAL(8) :: Tunit,Runit,Munit
    REAL(8) :: Axis,ECC,OmegaA_Con
    
    END MODULE
    
    SUBROUTINE Instantiate_Global_Values()
    
    USE GlobalDef
     
    REAL(8) :: Axis_aA,Axis_bA,Axis_cA,Axis_aB,Axis_bB,Axis_cB,rhoA,rhoB
    REAL(8) :: RadiiA,RadiiB,SemiAxis,MassA,MassB
    REAL(8) :: J2A,J2B,J22A,J22B,J4A,J42A,J44A,J4B,J42B,J44B
    
    !Axis_aA =   5D2
    !Axis_bA =   48D1
    !Axis_cA =   48D1
    !Axis_aB =   1D3
    !Axis_bB =   96D1
    !Axis_cB =   96D1
    !SemiAxis=   3D3
    !RhoA     =   2D3
    !RhoB     =   2D3
    
    ! Case 2
    Axis_aA =   5D2
    Axis_bA =   5D2
    Axis_cA =   5D2
    Axis_aB =   15D2
    Axis_bB =   1450D0
    Axis_cB =   1450D0
    SemiAxis=   3D3
    RhoA     =   2D3
    RhoB     =   2D3
    
    ECC     =   0D-1
    RadiiA  =   (Axis_aA*Axis_bA*Axis_cA)**(1D0/3D0)
    RadiiB  =   (Axis_aB*Axis_bB*Axis_cB)**(1D0/3D0)
    MassA   =   4D0/3D0*PI*RadiiA**3D0*RhoA
    MassB   =   4D0/3D0*PI*RadiiB**3D0*RhoB
    MassRatio = MassB/(MassA+MassB)
    Mass = MassRatio*(1D0-MassRatio)
    
    Runit   =   SemiAxis
    Munit   =   MassA+MassB
    Tunit   =   DSQRT(Runit**3/(GG * Munit))
    
    OmegaA_Con = 1D0
    Axis = SemiAxis/Runit
    AlphaA  =   RadiiA/Runit
    AlphaB  =   RadiiB/Runit
    J2A = (Axis_aA**2D0+Axis_bA**2D0-2D0*Axis_cA**2D0)/(1D1*RadiiA**2D0)
    J2B = (Axis_aB**2D0+Axis_bB**2D0-2D0*Axis_cB**2D0)/(1D1*RadiiB**2D0)
    J22A = (Axis_aA**2D0-Axis_bA**2D0)/(2D1*RadiiA**2D0)
    J22B = (Axis_aB**2D0-Axis_bB**2D0)/(2D1*RadiiB**2D0)
    J4A=-15D0/7D0*(J2A**2+2D0*J22A**2)
	J42A=-5D0/7D0*J2A*J22A
	J44A=5D0/28D0*J22A**2
	J4B=-15D0/7D0*(J2B**2+2D0*J22B**2)
	J42B=-5D0/7D0*J2B*J22B
	J44B=5D0/28D0*J22B**2
    A1 = 1D0/2D0*(J2A*AlphaA**2D0+J2B*AlphaB**2D0)
    A2 = 3D0*J22A*AlphaA**2D0
    A3 = 3D0*J22B*AlphaB**2D0
    B1=-3D0/8D0*(J4A*AlphaA**4+J4B*AlphaB**4)+9D0/4D0*J2A*J2B*AlphaA**2*AlphaB**2
	B2=-15D0/2D0*(J42A*AlphaA**4-J22A*J2B*AlphaA**2*AlphaB**2)
	B3=105D0*J44A*AlphaA**4
	B4=-15D0/2D0*(J42B*AlphaB**4-J2A*J22B*AlphaA**2*AlphaB**2)
	B5=105D0*J44B*AlphaB**4
	B6=9D0/2D0*J22A*J22B*AlphaA**2*AlphaB**2
	B7=105D0/2D0*J22A*J22B*AlphaA**2*AlphaB**2
    IzA = (1D0-MassRatio)*(Axis_aA**2D0+Axis_bA**2D0)/(5D0*Runit**2D0)
    IzB = MassRatio*(Axis_aB**2D0+Axis_bB**2D0)/(5D0*Runit**2D0)
    
    
    write(*,*) ' '
    WRITE(*,*) '======================== System Parameters =========================='
    write(*,*) ' '
    WRITE(*,*) 'aA  (m)                 =   ', Axis_aA
    WRITE(*,*) 'bA  (m)                 =   ', Axis_bA
    WRITE(*,*) 'cA  (m)                 =   ', Axis_cA
    WRITE(*,*) 'aB  (m)                 =   ', Axis_aB
    WRITE(*,*) 'bB  (m)                 =   ', Axis_bB
    WRITE(*,*) 'cB  (m)                 =   ', Axis_cB
    WRITE(*,*) 'MassRatio               =   ', MassRatio
    WRITE(*,*) 'MassRatio(q)            =   ', MassA/MassB
    WRITE(*,*) 'OrbitPeroid             =   ', PI2*DSQRT(Axis**3D0)*Tunit/36D2
    WRITE(*,*) ' '
    WRITE(*,*) 'Runit                   =   ', Runit
    WRITE(*,*) 'Munit                   =   ', Munit
    WRITE(*,*) 'Tunit                   =   ', Tunit/36D2
        
    END SUBROUTINE
    
    ! Test
    !Axis_aB =   658.5D0
    !Axis_bB =   Axis_aB/1.01D0
    !Axis_cB =   Axis_bB/1D0
    !Axis_aA =   295D0
    !Axis_bA =   Axis_aA/1D0
    !Axis_cA =   Axis_bA/1D0
    !SemiAxis=   2548D0
    !RhoA    =   1300D0
    !RhoB    =   1970D0
        
    ! Case 1
    !Axis_aA =   15D2
    !Axis_bA =   15D2
    !Axis_cA =   15D2
    !Axis_aB =   5D2
    !Axis_bB =   48D1
    !Axis_cB =   48D1
    !SemiAxis=   3D3
    !RhoA     =   2D3
    !RhoB     =   2D3
    !
    ! Case 2
    !Axis_aA =   5D2
    !Axis_bA =   5D2
    !Axis_cA =   5D2
    !Axis_aB =   15D2
    !Axis_bB =   1493D0
    !Axis_cB =   1493D0
    !SemiAxis=   3D3
    !RhoA     =   2D3
    !RhoB     =   2D3
    
    !Axis_aA =   5D2
    !Axis_bA =   48D1
    !Axis_cA =   48D1
    !Axis_aB =   1D3
    !Axis_bB =   96D1
    !Axis_cB =   96D1
    !SemiAxis=   3D3
    !RhoA     =   2D3
    !RhoB     =   2D3
    
    ! Unperturbed
    !Axis_aA =   5D2
    !Axis_bA =   5D2
    !Axis_cA =   5D2
    !Axis_aB =   1D3
    !Axis_bB =   1D3
    !Axis_cB =   1D3
    !SemiAxis=   3D3
    !RhoA     =   2D3
    !RhoB     =   2D3
    
    ! 1999 KW4 / (66391) Moshup and Squannit
    !Axis_aB =   658.5D0
    !Axis_bB =   Axis_aB/1.02D0
    !Axis_cB =   Axis_bB/1.11D0
    !Axis_aA =   295D0
    !Axis_bA =   Axis_aA/1.23D0
    !Axis_cA =   Axis_bA/1.33D0
    !SemiAxis=   2548D0
    !RhoA    =   1300D0
    !RhoB    =   1970D0
    
    ! Didymos
    !Axis_aA =   103D0
    !Axis_bA =   79D0
    !Axis_cA =   66D0
    !Axis_aB =   399D0
    !Axis_bB =   392D0
    !Axis_cB =   380D0
    !SemiAxis=   1180D0
    !RhoA     =   2D3
    !RhoB     =   2D3
    
    ! 1999 VO123 + 2002 ES90
    !Axis_aB =   450D0
    !Axis_bB =   Axis_aB/1.2D0
    !Axis_cB =   Axis_aB/1.2D0
    !Axis_aA =   Axis_aB*0.32D0
    !Axis_bA =   Axis_aA/1.5D0
    !Axis_cA =   Axis_aA/1.5D0
    !SemiAxis=   Axis_aB*3.1D0
    !RhoA     =   2D3
    !RhoB     =   2D3