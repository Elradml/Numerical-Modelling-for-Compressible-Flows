
! *************************************************************************************
! Lab for Msc course 2016
! Programmer V.A. Titarev/P. Tsoutsanis
! *************************************************************************************

! Finite-volume scheme with TVD Runge-Kutta time stepping
! for one-dimensional compressible Euler equations
! Use third order TVD Runge-Kutta in time

! Reference solution for test problem one is given in ref1.out


 ! Declaration of variables
 
 IMPLICIT NONE

 ! Spatial order
 Integer :: SpatialOrder ! 1 or 2
 
 ! Flux type
 Integer FluxType

 ! Courant number
 Real  CFL

 ! Type of initial condition, number of spatial cells
 Integer  IvType, N

 ! Vector of conservative variables CSV(1:4,-6+n+5), first index is the Runge-Kutta stage
 Real, ALLOCATABLE::  CSV(:,:,:) 
 ! Primitive variables, no RK stage, W = (rho,u,P)
 Real, ALLOCATABLE::  PV(:,:,:)
 ! Intercell fluxes
 Real, ALLOCATABLE::  IntercellFlux(:,:,:)
 

 ! other variables
 Integer  ::  OutFreq = 25, Rkstage=0
 Real     LB, RB, h,Time
 Real, ALLOCATABLE::  CELL_CENTERS(:)
 Real, parameter :: GM=1.4
 Real DT,T_
 Integer IT

 ! ---------------------------------- START THE PROGRAM----------------------

 
 ! initialise of variables, grid etc
  CALL INITALL
  
 ! START TIME CYCLE
 IT=0
 do  !
   ! Compute a stable time step
   CALL ComputeTimeStep(dt)
   ! Use third-order TVD Runge-Kutta
   Call ThirdOrderTVD
   ! advance time and time counter
   T_ = T_  + DT
   IT = IT+1

   If ( MOD(IT,OutFreq) ==0) PRINT*,' it= ',it,'  t= ',T_
   If (mod(it,2) .eq. 0)  Call OutputTecplot 
   
   ! check whether we reached the output time
   If ( ABS(T_ - TIME)/TIME .LE. 1D-8) GOTO 101
 enddo

 101 CONTINUE
 Call Output
 Call OutputTecplot 
 close(111) ! close the movie file
 print*,' Job finished.'
 Print*,' Number of time steps : ',it
 PRINT*,' The end' 


 ! //////// code's subroutines ///////////////
 
 Contains 

 
 !%%%%%%%%%% initialization of the run %%%%%%%%%%
 Subroutine InitAll
  Integer I
  Real X,U1,U2,U3

 ! read the input file
  Open(1,file='euler.ini')
   Read(1,*) SpatialOrder
   Read(1,*) IvType
   Read(1,*) N
   Read(1,*) CFL
   Read(1,*) FluxType
  Close(1)


 SELECT Case(IVTYPE)
 Case(1) 
  LB = 0.
  RB = +1. 
  TIME = 0.2
 Case(2)
  LB = -5.
  RB = +5.
 Time = 5.d0   
 Case Default
  print*,' Wrong test problem number. Stop the code.'
  stop
 END SELECT

 ! Spatial Cell size
 h = (RB-LB)/N

 ! Vector of conservative variables QC(1:4,-6+n+5), first index is the Runge-Kutta stage
 ALLOCATE(CSV(1:4,3,-7:n+6),PV(4,1:4,-7:n+6),InterCellFlux(1:4,3,-1:n+1))
 ALLOCATE(CELL_CENTERS(1:N))

 ! calculate cell centers
 do I=1,n
   CELL_CENTERS(I) =  LB+I*H - H/2
 enddo

 ! initialise the vector of conserved quantities
  do I=1,N
    X = CELL_CENTERS(I)
    CALL  U0(x,CSV(1,:,i))
  enddo
 
  ! calculate primitive variables from conservative
  Rkstage = 1
  do I=1,N
   PV(Rkstage,1,i) = CSV(Rkstage,1,I)
   PV(Rkstage,2,I) = CSV(Rkstage,2,I)/CSV(Rkstage,1,I)
   PV(Rkstage,3,I) = (GM-1)*( CSV(Rkstage,3,I) - 0.5*CSV(Rkstage,2,I)*PV(Rkstage,2,I))
   PV(Rkstage,4,I) = sqrt(GM*PV(Rkstage,3,I)/PV(Rkstage,1,I))
  enddo
  
  ! set flow time to zero
   T_=0.D0
   
  OPEN(UNIT = 111, FILE = 'movie.dat', STATUS = 'UNKNOWN')
  WRITE(111,*)'TITLE="Solution" '
  WRITE(111,*)'VARIABLES="X" "rho" "u" "p"'
  Call OutputTecplot 
   
 End subroutine


 !%%%%%%%%%%%%%% Set up boundary conditions for given stage of the Runge Kutta marching %%%%%%%%%%
 
 Subroutine SetBC(Rkstage)
   Integer k,i,Rkstage
   
   ! set up ghost cells 
   
   Do i=-5,0
    do k=1,3
    CSV(Rkstage,k,i)  = CSV(Rkstage,k,abs(i)+1)
	enddo
    do k=1,4
     PV(Rkstage,k,i)  = PV(Rkstage,k,abs(i)+1)
	enddo
   Enddo

   Do i=1,6
    do k=1,3
     CSV(Rkstage,k,N+i)  = CSV(Rkstage,k,N-1-I)    
	enddo
    do k=1,4
     PV(Rkstage,k,N+i)  = PV(Rkstage,k,N-1-I)    
	enddo
   Enddo
 End subroutine   


 !%%%%%%%%%%% Compute initial data at t=0 for given spatial position 'x' %%%%%%%%%%%%%
 Subroutine U0(X,Q)
  ! U1 = RHO, U2 = RHOU, U3 = E
   Real X,U1,U2,U3,Q(3)
   Real DL,DR,UL,UR,PL,PR,x0
   Real :: pi= 3.141592653589793
   ! IvType = 1 : Sod' Shock Tube  Problem
   ! IvType = 2    Shock - turbulence interaction

  SELECT Case(IVTYPE)
   Case(1) 
    DL=1.0 ;   UL=0.0 ;  PL=1.0 
    DR=0.125 ;  UR=0.0 ;  PR=0.1 

    If (X .LE. 0.5)  THEN
     U1 = DL
     U2 = DL*UL
     U3 = PL/(GM-1) + 0.5*DL*UL**2
    Else
     U1 = DR
     U2 = DR*UR
     U3 = PR/(GM-1) + 0.5*DR*UR**2
    Endif
    
 Case(2)    
  ! Long time shock/turbulence interaction
  ! Mach number 1.1, S=1.5
  DL=   1.51569506726457     
  UL=   0.523345519274197     
  PL=   1.80500000000000     
  
  IF (X .LE. -4.5)  THEN
   U1 = DL
   U2 = DL*UL
   U3 = PL/(GM-1) + 0.5*DL*UL**2
  ELSE
   U1 = 1 + 0.1d0*sin(20*pi*x)
   U2 = 0.
   U3 = 1./(GM-1)  ! U = 0
  ENDIF
   
 End Select

 Q(1) = U1
 Q(2) = U2
 Q(3) = U3
end subroutine


!%%%%%%%%%%%%%%%% Time marching algorithm, which uses third order TVD RK method  %%%%%%%
!%%%%%%  Jiang G.S. and Shu C.W. Efficient Implementation of  weighted ENO schemes //J. Comput. Phys. 1996.  V. 126.  pp.202-212.

  Subroutine  ThirdOrderTVD
   Integer i,k

  ! loop stages from 1 to 3
   do Rkstage=1,3
     ! set up boundary conditions
     Call SetBc(Rkstage)
     ! calculate intercell fluxes
     Call ComputeFlux(Rkstage)
     ! perform the update
     CALL Update(Rkstage)
     Do i=1,n
       PV(Rkstage+1,1,i) = CSV(Rkstage+1,1,i)
       PV(Rkstage+1,2,i) = CSV(Rkstage+1,2,i)/CSV(Rkstage+1,1,i)
       PV(Rkstage+1,3,i) = (gm-1)*( CSV(Rkstage+1,3,i) - 0.5*PV(Rkstage+1,1,i)*PV(Rkstage+1,2,i)**2)
       PV(Rkstage+1,4,I) = sqrt(GM*PV(Rkstage+1,3,I)/PV(Rkstage+1,1,I))
     Enddo
   enddo
  
   ! re-assign the flow variables to stage 1 of RK method
   Do i=1,n
   do k=1,3
     CSV(1,k,i) = CSV(4,k,i)
   enddo
   do k=1,4
     PV(1,k,i) = PV(4,k,i)
   enddo
   Enddo

  End subroutine

 !%%%%%%%%%%% Solution update for each stage of TVD RK method %%%%%%%%%%
  Subroutine UPDATE(Rkstage)
   Integer I,K,Rkstage

   SELECT Case(Rkstage)
    Case(1)
     do i=1,n
      do K=1,3
       CSV(2,K,i)  =  CSV(1,K,i)  - (Dt/H)*(InterCellFlux(1,K,i) - InterCellFlux(1,K,i-1))
      enddo
     enddo

    Case(2)
     do i=1,n
      do K=1,3
       CSV(3,k,i) =  0.75*Csv(1,k,i) + 0.25*CSV(2,k,i)     - (0.25*Dt/H)*(InterCellFlux(2,k,i) - InterCellFlux(2,k,i-1))
      enddo
     enddo

    Case(3)
     do i=1,n
      do k=1,3
       CSV(4,K,i) =  (1./3)*CSV(1,K,i) + (2./3)*CSV(3,K,i)     - (2./3)*(Dt/H)*(InterCellFlux(3,K,i) - InterCellFlux(3,K,i-1))
	  enddo
     enddo
    END SELECT
  end subroutine


 !%%%%%%%%% write the output file %%%%%%%%%
  Subroutine    Output
  Integer i
   202 format(6(2x,e11.4))
   open(1,file='results.dat')
   WRITE(1,*)'TITLE="Solution" '
  WRITE(1,*)'VARIABLES="X" "rho""u""P" "T"' ! "u" "p"
   WRITE(1,*)'ZONE ',',I=',n, ',F="POINT"'
   do i=1,n   
    write(1,202) cell_centers(i),CSV(1,1,i),PV(1,2,I),PV(1,3,I),PV(1,3,I)/PV(1,1,I)
   enddo
   close(1)
  End subroutine   


  !%%%%%%%%%% Evaluation of the physical flux function from the conserved vector CDS =(rho,rho*u,E)
  Subroutine FluEval(CDS,Flux)
    Real cds(3),p,u,flux(3)
    
    u = cds(2)/cds(1)
    p = (gM-1)*(Cds(3) - 0.5*cds(1)*u**2)
    Flux(1) = cds(2)
    Flux(2) = cds(2)*u + p
    Flux(3) = (cds(3)+p)*u 
    
   End subroutine


  !%%%%%%%%%% Calculation of a stable time step %%%%%
  Subroutine ComputeTimeStep(dt)
   Integer i
   Real Umax,dt, a
  
   umax  = 0.0 
   Do i=1,n
    ! compute the sound speed
    a = ComputeSoundSpeed(CSV(1,:,i))
    umax = max(umax, a + abs(PV(1,2,i)))
   Enddo

   ! reduce the time step for first 10 time steps
   If ( IT<10) THEN
    dt = MIN(0.1*H/UMAX, TIME-T_)
   Else
    dt = MIN(CFL*H/UMAX, TIME-T_)
   Endif
  End subroutine

 
 !%%%%%%%%%%%%%%%% Calculation of the sound speed  on the conserved vector CDS
  Real function ComputeSoundSpeed(cds)
    Real cds(3),p,u  
    u = cds(2)/cds(1)
    p = (gm-1)*(cds(3) - 0.5*cds(2)*u)
    ComputeSoundSpeed=sqrt(gm*p/cds(1))  
  End function


  !%%%%%%%%%%%%%%%%% minmod slope limiter %%%%%%%%%%%%%%%%%
  Real function minmod(x,y)
   Real x,y
	if (x * y.ge.0.0d0)then
	minmod = sign(1.0, x) * min(abs(x), abs(y))
	else
	minmod = 0.0
	end if
  End function

 ! Compute the numerical flux
 Subroutine ComputeFlux(Rkstage)
  Integer i,k,Rkstage
  Real CDL(3),CDR(3),LocalFlux(3) 

  ! Loop over the spatial index i
  do I=0,N
   
   ! call reconstruction procedure at RK stage Rkstage to compute left CDL and right CDR
   ! values of the conserved vector between cells i and i+1
   CALL Reconstruction(CSV(Rkstage,:,i-2:i+3),CDL,CDR)
  
   ! calculate the numerical flux using reconstructed conserved vectors CDL, CDR
   ! Left   initial data  for the local Riemann problem is given by CDL
   ! Right  initial data  for the local Riemann problem is given by CDR

   Select Case(FluxType)
   Case(1)
      CALL LxF(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(2)
      CALL Rusanov(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(3)
      CALL HLL(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux

   Case(4)
      CALL HLLC(CDL,CDR,LocalFlux)
      InterCellFlux(Rkstage,:,i) = LocalFlux
      
   Case default
    print*,' the flux is not defined. stop the program'
	read*
	stop
  End select	 

  Enddo
 end subroutine


  !%%%%%%%%% Lax Friedrich flux %%%%%%%%
   Subroutine LxF(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3)
      
        CALL FLUEVAL(CDL,FL)
        CALL FLUEVAL(CDR,FR)
        
        Flux = 0.5*(FL+FR) - 0.5*(h/dt)*(CDR-CDL) 
  End subroutine
  
  
 !%%%%%%%% Rusanov flux %%%%%%%%%%%%%%%%%%
   Subroutine Rusanov(CDL,CDR,Flux)
      Real  FL(3), FR(3),CDL(3),CDR(3),Flux(3),Speed, speed_l, speed_r

	speed_l = abs(cdl(2)/cdl(1)) + ComputeSoundSpeed(cdl)
	speed_r = abs(cdr(2)/cdr(1)) + ComputeSoundSpeed(cdr)
	Speed = max(speed_l, speed_r)
	
	Flux = 0.5 * (FL+FR) - 0.5 * (Speed)*(CDR-CDL)
      
  End subroutine 

 !%%%%%%%%%%%%%% HLL flux %%%%%%%%%%%%%%%%%%%%%
   Subroutine HLL(CDL,CDR,Flux)
      Real  CDL(3),CDR(3),Flux(3), fl(3), fr(3), speed_l, speed_r

	CALL FLUEVAL(CDL, FL)
	CALL FLUEVAL(CDR, FR)

	speed_l = min((cdl(2)/cdl(1)) - ComputeSoundSpeed(cdl), (cdr(2)/cdr(1)) - ComputeSoundSpeed(cdr))
	speed_r = max((cdr(2)/cdr(1)) + ComputeSoundSpeed(cdr), (cdr(2)/cdr(1)) + ComputeSoundSpeed(cdr))

	if (speed_l.ge.0.0d0)then
	Flux = fl
	end if
	if (speed_r.le.0.0d0)then
	Flux = fr
	end if
	if ((speed_l.le.0.0d0).and.(speed_r.ge.0.0d0))then
	Flux = ((speed_r*fl) - (speed_l*fr) + (speed_l*speed_r*(cdr-cdl))) / (speed_r-speed_l)
	end if

  End subroutine 


 !%%%%%%%%%%%% HLLC flux %%%%%%%%%%%%%%%%%%%%%
   Subroutine HLLC(CDL,CDR,Flux)
      Real  CDL(3),CDR(3),Flux(3), fl(3), fr(3), speed_l, speed_r, speed_m, uln, urn
      Real  rho_l, rho_r, PL, PR, EL, ER, ul, ur, al, ar
      Real  Ustarl(3), Ustarr(3)

	! setting left variables
	rho_l = cdl(1)
	ul = cdl(2) / rho_l
	EL =  cdl(3)
	PL = (GM - 1.0) * (EL - 0.5 * rho_l * ul**2)
	al = sqrt(GM * PL / rho_l)

	! setting right variables
	rho_r = cdr(1)
	ur = cdr(2) / rho_r
	ER = cdr(3)
	PR = (GM - 1.0) * (ER - 0.5 * rho_r * ur**2)
	ar = sqrt(GM * PR / rho_r)
	
	CALL FLUEVAL(CDL, FL)
	CALL FLUEVAL(CDR, FR)
	
	speed_l = min(ul - al, ur - ar)
	speed_r = max(ul + al, ur + ar)

	!middle wave speed
	speed_m = ((PR - PL) + rho_l * ul * (speed_l * - ul) - rho_r * ur * (speed_r - ur)) / (rho_l * (speed_l - ul) - rho_r *(speed_r - ur))	
	
	! left star states
	Ustarl(1) = rho_l * ((speed_l - ul) / (speed_l - speed_m))
	Ustarl(2) = speed_m * (rho_l * ((speed_l - ul) / (speed_l - speed_m)))
	Ustarl(3) = (rho_l * ((speed_l - ul) / (speed_l - speed_m))) * ((EL / rho_l) + (speed_m - ul) *(speed_m + (PL / rho_l *(speed_l - ul))))

	! right star states
	Ustarr(1) = rho_r * ((speed_r - ur) / (speed_r - speed_m))
	Ustarr(2) = speed_m * (rho_r * ((speed_r - ur) / (speed_r - speed_m)))
	Ustarr(3) = (rho_r * ((speed_r - ur) / (speed_r - speed_m))) * ((ER / rho_r) + (speed_m - ur) *(speed_m + (PR / rho_r *(speed_r - ur))))

	
	if (speed_l.ge.0.0d0)then
	Flux = fl
	end if
	if ((speed_l.le.0.0d0).and.(speed_m.ge.0.0d0))then
	Flux = fl + speed_l * (Ustarl - cdl)
	end if
	if ((speed_m.le.0.0d0).and.(speed_r.ge.0.0d0))then
	Flux = fr + speed_r * (Ustarr - cdr)
	end if
	if (speed_r.le.0.0d0)then	
	Flux = fr
	end if

  End subroutine 

 !%%%%%%%%%%%%%%%% Reconstruction procedure%%%%%%%%%%%%%%
 ! Input: one-dimensional array U1D of flow quantities near cell interface i+1/2
 ! Output: left CDL and right CDR values at interface
 Subroutine Reconstruction(U1D,CDL,CDR)
   Integer I, F
   Real U1d(3,-2:3),CDL(3),CDR(3)
   Real r, lim_slope(3), slope_l, slope_r

   select Case(SpatialOrder)
  
	 Case(1)
	 
	 ! First order 
	  Do f=1,3
	   CDL(f) = U1D(f,0)  
	   CDR(f) = U1D(f,1)  	
	  Enddo 
	 
	 ! second order TVD 
	 Case(2)
	 Do F = 1, 3
	 slope_l = U1D(F,0) - U1D(F, -1)
	 slope_r = U1d(F,1) - U1D(F, 0)

	 lim_slope(F) = minmod(slope_l, slope_r)
	
	 cdl(F) = U1D(F, 0) + 0.5 * lim_slope(F)
	 cdr(F) = U1D(F, 1) - 0.5 * lim_slope(F)
	 End Do

	 Case default
	 print*,' Wrong spatial accuracy. Stop the code'
	 stop  
  end select

 end subroutine

 Subroutine OutPutTecplot
    Integer i,j
	real(8) x
    
   55 Format (4(2x,e11.4))
    WRITE(111,*)'ZONE ',',I=',n, ',F="POINT"'
    WRITE(111,*) ', SOLUTIONTIME=',T_
    DO  I = 1,n
	  x = lb + i*h-h/2
      WRITE(111,55)X,pv(1,1,i),pv(1,2,i),pv(1,3,i)
    Enddo

 End subroutine
  

 END 
