program rh_special_function_probe
  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: i128 = selected_int_kind(38)
  integer :: chain, child_length, ios, i, ifact, jbot, mmax
  real(dp) :: p, sp, tpp, tpm, c, gam, zet(3:13)
  real(dp) :: a1, xr, fracL, sx, con1, con2, con3, sav1, sav2, z
  real(dp) :: psi_args(6), psi_values(6)
  complex(kind=16) :: coeff(7,0:6), zlarcoeff(-1:6), zsmalcoeff(-1:6)
  complex(kind=16) :: c1, c2, c3, c4, c5, c6, c7, epi4
  complex(kind=16) :: endpoint, erf_values(2), erf_weights(2)
  complex(kind=16) :: fn, cr2, cr3, source_t1, source_t2, source_t4
  character(len=32) :: h_a1, h_xr, h_frac, h_endpoint_re, h_endpoint_im

  common /PARMS/ coeff,zlarcoeff,zsmalcoeff
  common /PARS2/ gam,zet
  common /PARS4/ p,sp,tpp,tpm,epi4

! Set p=pi, tpp=2*pi and tpm=-2*pi  

  C=1.0
  p=4*ATAN(C)
  mmax=10
  tpp=2*p
  tpm=-tpp

  !       Standard constants sqrt(pi) exp(pi*i/4).

  C=1.0
  sp=sqrt(p)
  epi4=(1.0,1.0)/sqrt(2*C)

  !       More standard constants relating to the zeta function and gamma function for real values.
  !       Double precision is fine for the values of zeta(3-9).

  GAM=0.57721566490153286060D0

  zet(3)=1.20205690315959428540d0
  zet(5)=1.03692775514336992633d0
  zet(7)=1.00834927738192282684d0
  zet(9)=1.00200839282608221442d0
  zet(11)=1.0004941886041194646d0
  zet(13)=1.0001227133475784891d0


  !     ** Coefficients used to calculate the complex error function. If you can make the
  !     intrinsic error function erf work for complex arguments all the below is unnecessary

  C1=((1.0,0.0)*cos(1.0d0)+(0.0,1.0)*sin(1.0d0))*sqrt(2.0/p) 
  C2=((1.0,0.0)*cos(4.0d0)+(0.0,1.0)*sin(4.0d0))*sqrt(2.0/p)        
  C3=((1.0,0.0)*cos(2.25d0)+(0.0,1.0)*sin(2.25d0))*sqrt(2.0/p)
  C4=((1.0,0.0)*COS(6.25d0)+(0.0,1.0)*SIN(6.25d0))*sqrt(2.0/p)
  C5=((1.0,0.0)*cos(1.5625d0)+(0.0,1.0)*sin(1.5625d0))*sqrt(2.0/p)        
  C6=((1.0,0.0)*cos(3.0625d0)+(0.0,1.0)*sin(3.0625d0))*sqrt(2.0/p)
  C7=((1.0,0.0)*COS(5.0625d0)+(0.0,1.0)*SIN(5.0625d0))*sqrt(2.0/p)
  !     error function coefficients start here.

  !     COEFF(1,..) WILL APPLY FOR Z=1

  COEFF(1,0)=(0.969264211944215930381490d0,-0.474147636640994245161680d0)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1)

  COEFF(1,1)=C1*(1.0,-1.0)
  COEFF(1,2)=C1*(1.0,1.0)
  COEFF(1,3)=-C1*(1.0,-3.0)/3.0d0
  COEFF(1,4)=-C1*(5.0,-1.0)/6.0d0
  COEFF(1,5)=-C1*(11.0,13.0)/30.0d0
  COEFF(1,6)=C1*(9.0d0,-31.0d0)/90.0d0

  !     NOW DO SAME FOR ERF(EXP(-i*PI/4)*2)

  COEFF(2,0)=(1.010311712025489491642626d0,0.273925759463539899021137d0)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*2)

  COEFF(2,1)=C2*(1.0,-1.0)
  COEFF(2,2)=C2*(2.0,2.0)
  COEFF(2,3)=C2*(-7.0,9.0)/3.0d0
  COEFF(2,4)=-C2*(11.0,5.0)/3.0d0
  COEFF(2,5)=C2*(13.0,-109.0)/30.0d0
  COEFF(2,6)=C2*(129.0,-31.0)/45.0d0

  COEFF(3,0)=(1.338389640116239225341021d0,-0.096501782737190712909169d0)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1.5)

  COEFF(3,1)=C3*(1.0,-1.0)
  COEFF(3,2)=C3*(1.5,1.5)
  COEFF(3,3)=C3*(-7.0,11.0)/6.0d0
  COEFF(3,4)=C3*(-15.0,-3.0)/8.0d0
  COEFF(3,5)=C3*(-13.0,-59.0)/40.0d0
  COEFF(3,6)=C3*(67.0,-53.0)/80.0d0

  COEFF(4,0)=(8.264692402664926098554228d-1,-1.394623184589031908854228d-1)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*2.5)

  COEFF(4,1)=C4*(1.0,-1.0)
  COEFF(4,2)=C4*(2.5,2.5)
  COEFF(4,3)=C4*(-23.0,27.0)/6.0d0
  COEFF(4,4)=C4*(-155.0,-95.0)/24.0d0
  COEFF(4,5)=C4*(313.0,-913.0)/120.0d0
  COEFF(4,6)=C4*(1065.0,65.0)/144.0d0

  COEFF(5,0)=(1.215497305362876593269559d0,-0.344267547866857083642667d0)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1.25)

  COEFF(5,1)=C5*(1.0,-1.0)
  COEFF(5,2)=C5*(1.25,1.25)
  COEFF(5,3)=C5*(-17.0,33.0)/24.0d0
  COEFF(5,4)=C5*(-245.0,-5.0)/192.0d0
  COEFF(5,5)=C5*(-767.0,-1633.0)/1920.0d0
  COEFF(5,6)=C5*(1665.0,-2335.0)/4608.0d0

  COEFF(6,0)=(1.260051027257371974124166d0,1.664741404414525369657245d-1)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*1.75)

  COEFF(6,1)=C6*(1.0,-1.0)
  COEFF(6,2)=C6*(1.75,1.75)
  COEFF(6,3)=C6*(-41.0,57.0)/24.0d0
  COEFF(6,4)=C6*(-511.0,-175.0)/192.0d0
  COEFF(6,5)=C6*(-143.0,-4561.0)/1920.0d0
  COEFF(6,6)=C6*(37527.0,-17353.0)/23040.0d0

  COEFF(7,0)=(7.873277503318070969692577d-1,1.234468235979886882660918d-1)

  !     EXACT VALUE FOR ERF(EXP(-i*PI/4)*2.25)

  COEFF(7,1)=C7*(1.0,-1.0)
  COEFF(7,2)=C7*(2.25,2.25)
  COEFF(7,3)=C7*(-73.0,89.0)/24.0d0
  COEFF(7,4)=C7*(-315.0,-171.0)/64.0d0
  COEFF(7,5)=C7*(827.0,-3419.0)/640.0d0
  COEFF(7,6)=C7*(12081.0,-879.0)/2560.0d0


  !     NEXT COMPUTE THE COEFFICIENTS RELEVANT TO THE ASYMPTOTIC APPROXIMATION
  !     OF ERF(Z*EXP(-i*PI/4)) WHEN ABS(Z)>>1.

  ZLARCOEFF(-1)=(1.0,1.0)/SQRT(2*p)
  ZLARCOEFF(0)=(1.0,0.0)
  DO I=1,6
     ZLARCOEFF(I)=(2*I-1)*ZLARCOEFF(I-1)/(0.0,2.0)
  ENDDO

  !     NEXT COMPUTE THE COEFFICIENTS RELEVANT TO THE APPROXIMATION
  !     OF ERF(Z*EXP(-i*PI/4)) FOR SMALL ABS(Z)<<1

  ZSMALCOEFF(-1)=(1.0,-1.0)/SQRT(0.5*p)
  ZSMALCOEFF(0)=(1.0,0.0)
  C1=ZSMALCOEFF(0)
  IFACT=1
  DO I=1,6
     IFACT=IFACT*I
     C1=C1*(0.0,1.0)
     ZSMALCOEFF(I)=C1/((2*I+1.0)*IFACT*1.0)
  ENDDO

  
  write(*,'(A)',advance='no') '#constants'
  write(*,'(1X,A)',advance='no') real_hex(p)
  write(*,'(1X,A)',advance='no') real_hex(sp)
  write(*,'(1X,A)',advance='no') real_hex(tpp)
  write(*,'(1X,A)',advance='no') real_hex(tpm)
  write(*,'(1X,A)',advance='no') real_hex(real(epi4,kind=dp))
  write(*,'(1X,A)') real_hex(aimag(epi4))

  do
     read(*,*,iostat=ios) chain, child_length, h_a1, h_xr, h_frac, h_endpoint_re, h_endpoint_im
     if (ios.lt.0) exit
     if (ios.ne.0) error stop 'special-function probe input failure'
     call real_from_hex(h_a1,a1)
     call real_from_hex(h_xr,xr)
     call real_from_hex(h_frac,fracL)
     call complex_from_hex(h_endpoint_re,h_endpoint_im,endpoint)

     jbot=ceiling(a1)
     psi_args(1)=child_length+fracL+1.0
     psi_args(2)=1.0-fracL
     psi_args(3)=a1+1.0
     psi_args(4)=child_length+1.0-a1
     psi_args(5)=jbot-a1
     psi_args(6)=jbot-(child_length+fracL)
     do i=1,6
        call psi(psi_args(i),psi_values(i))
     enddo

     sav1=psi_values(2)
     con2=psi_values(1)-psi_values(2)
     sav2=psi_values(4)
     con3=psi_values(3)-psi_values(4)
     source_t1=(0.0,1.0)*(conjg(endpoint)*con2-con3)/tpp

     sx=sqrt(xr)
     con1=1.0-fracL
     z=sp*con1/sx
     call erf(z,6,erf_values(1))
     cr2=1.0-erf_values(1)
     fn=-(0.0,1.0)*p*(con1**2)/xr
     cr3=p*exp(fn)*cr2/(sx*epi4)-1/con1-1/(child_length*1.0+1.0)
     source_t2=(0.0,1.0)*conjg(endpoint)*cr3/tpp
     erf_weights(1)=conjg(-(0.0,1.0)*conjg(endpoint)*p*exp(fn)/(sx*epi4*tpp))

     if (a1.gt.0.0) then
        con1=a1
     else
        con1=a1+1.0
     endif
     z=sp*con1/sx
     call erf(z,6,erf_values(2))
     cr2=1.0-erf_values(2)
     fn=-(0.0,1.0)*p*(con1**2)/xr
     cr3=(0.0,1.0)*exp(fn)*cr2/(2*sx*epi4)
     if (a1.gt.0.0) then
        source_t4=-(0.0,1.0)*conjg(endpoint)/(tpp*(child_length+fracL))+cr3
     else
        source_t4=(0.0,1.0)*(a1/con1-1.0)/tpp+cr3
     endif
     erf_weights(2)=conjg(-(0.0,1.0)*exp(fn)/(2*sx*epi4))

     write(*,'(I0)',advance='no') chain
     do i=1,6
        write(*,'(1X,A,1X,A)',advance='no') real_hex(psi_args(i)),real_hex(psi_values(i))
     enddo
     write(*,'(1X,A,1X,A,1X,A)',advance='no') real_hex(sp*(1.0-fracL)/sx), &
          real_hex(real(erf_values(1),kind=dp)),real_hex(aimag(erf_values(1)))
     write(*,'(1X,A,1X,A,1X,A)',advance='no') real_hex(z), &
          real_hex(real(erf_values(2),kind=dp)),real_hex(aimag(erf_values(2)))
     do i=1,2
        write(*,'(1X,A,1X,A)',advance='no') real_hex(real(erf_weights(i),kind=dp)), &
             real_hex(aimag(erf_weights(i)))
     enddo
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(source_t1,kind=dp)),real_hex(aimag(source_t1))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(source_t2,kind=dp)),real_hex(aimag(source_t2))
     write(*,'(1X,A,1X,A)') real_hex(real(source_t4,kind=dp)),real_hex(aimag(source_t4))
  enddo

contains

  function real_hex(value) result(text)
    real(dp), intent(in) :: value
    character(len=32) :: text
    integer(i128) :: bits
    bits=transfer(value,bits)
    write(text,'(Z32.32)') bits
  end function real_hex

  subroutine real_from_hex(text,value)
    character(len=32), intent(in) :: text
    real(dp), intent(out) :: value
    integer(i128) :: bits
    read(text,'(Z32)') bits
    value=transfer(bits,value)
  end subroutine real_from_hex

  subroutine complex_from_hex(real_text,imag_text,value)
    character(len=32), intent(in) :: real_text,imag_text
    complex(kind=16), intent(out) :: value
    real(dp) :: real_part,imag_part
    call real_from_hex(real_text,real_part)
    call real_from_hex(imag_text,imag_part)
    value=cmplx(real_part,imag_part,kind=16)
  end subroutine complex_from_hex

end program rh_special_function_probe

SUBROUTINE ERF(Z,N,ER)

  !     COMPUTES ERF(EXP(-i*PI/4)*Z) FOR Z REAL.

  !     IT USES JUST N TERMS ONLY.
  !     USES TAYLOR SERIES IF ABS(Z) NEAR 1, POWER SERIES IF ABS(Z)<<1 AND
  !     ASYMPTOTIC SERIES IF ABS(Z)>>1. THE COEFFICIENTS OF THESE SERIES ARE
  !     STORED IN CEF,ZSMAL AND ZLAR RESPECTIVELY. THREE TAYLOR SERIES ARE USED
  !     CENTRED ON Z=1, 1.5 AND 2. BELOW ABS(Z)<0.8 THE POWER SERIES ABOUT Z=0
  !     IS EMPLOYED AND ABOVE ABS(Z)>2.25 THE ASYMPTOTIC SERIES IS USED.

  implicit none

  integer, parameter :: dp = selected_real_kind(33) 

  integer       :: N,I
  real(dp)      :: Z,P,E,ZA,SP,TPP,TPM
  complex(kind=16)   :: ER,CEF(7,0:6),ZLAR(-1:6),ZSMAL(-1:6),V,EPI4

  COMMON/PARMS/ CEF,ZLAR,ZSMAL
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4

  IF (N.GT.6) THEN
     WRITE(6,*) 'INCREASE NUMBER OF COEFFICIENTS IN SUBROUTINE ERF.'
     STOP
  ENDIF

  !     DO ALL MAIN CALCULATIONS ASSUMING Z>0.

  ZA=ABS(Z)
  IF (ZA.LT.1E-14) THEN
     ER=ZSMAL(-1)*ZA
     GOTO 1000
  ENDIF
  ER=(0.0,0.0)
  IF (ZA.GE.0.8.AND.ZA.LT.1.125) THEN

     !     ZA CLOSE TO 1, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1)
     !     TO GET ACCURATE SOLUTION

     E=ZA-1.0
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(1,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(1,0)
  ELSE IF (ZA.GE.1.125.AND.ZA.LT.1.375) THEN

     !     ZA CLOSE TO 1.25, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1.25)
     !     TO GET ACCURATE SOLUTION

     E=ZA-1.25
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(5,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(5,0)
  ELSE IF (ZA.GE.1.375.AND.ZA.LT.1.625) THEN

     !     ZA CLOSE TO 1.5, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1.5)
     !     TO GET ACCURATE SOLUTION


     E=ZA-1.5
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(3,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(3,0)
  ELSE IF (ZA.GE.1.625.AND.ZA.LT.1.875) THEN

     !     ZA CLOSE TO 1.75, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*1.75)
     !     TO GET ACCURATE SOLUTION


     E=ZA-1.75
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(6,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(6,0)
  ELSE IF (ZA.GE.1.875.AND.ZA.LT.2.125) THEN

     !     ZA CLOSE TO 2.0, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*2)
     !     TO GET ACCURATE SOLUTION


     E=ZA-2.0
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(2,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(2,0) 
  ELSE IF (ZA.GE.2.125.AND.ZA.LT.2.375) THEN

     !     ZA CLOSE TO 2.25, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*2.25)
     !     TO GET ACCURATE SOLUTION


     E=ZA-2.25
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(7,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(7,0)

  ELSE IF (ZA.GE.2.375.AND.ZA.LT.2.625) THEN

     !     ZA CLOSE TO 2.5, USE TAYLOR APPROX OF ERF(EXP(-I*PI/4)*2.5)
     !     TO GET ACCURATE SOLUTION


     E=ZA-2.5
     IF (ABS(E).GT.1E-14) THEN
        DO I=1,N
           ER=ER+CEF(4,I)*(E**I)
        ENDDO
     ENDIF
     ER=ER+CEF(4,0)     
  ELSE IF (ZA.GE.2.625) THEN

     !    ZA LARGE, USE ASMPTOTIC APPROXIMATION. ACTUALLY WE COMPUTE
     !    ERFC(EXP(-i*PI/4)*Z) FIRST AND THEN USE ERF(Z)=1-ERFC(Z).
     !    THE APPROXIMATION IS VALID WHEN ABS(ARG(EXP(-i*PI/4)*Z))=PI/4,
     !    WHICH APPLIES WHEN Z>0, IS LESS THAN 3*PI/4. THE CASE Z<0 IS DEALT
     !    WITH AT THE END.

     E=ZA*ZA
     V=((1.0,0.0)*COS(E)+(0.0,1.0)*SIN(E))/ZA
     DO I=0,N
        ER=ER+ZLAR(I)/(E**I)
     ENDDO

     !    CONVERT FROM ERFC TO ERF

     ER=(1.0D0,0.0D0)-ZLAR(-1)*V*ER
  ELSE

     !    ZA SMALL<1, USE THE SERIES SOLUTION WHICH CONVERGES RAPIDLY.

     E=ZA*ZA
     DO I=0,N
        ER=ER+ZSMAL(I)*(E**I)
     ENDDO
     ER=ZA*ZSMAL(-1)*ER
  ENDIF

  !     NOW IF Z<0 ORIGINALLY USE ERF(-Z)=-ERF(Z)

1000 IF (Z.LT.0.0D0) THEN
     ER=-ER
  ENDIF
  RETURN
END SUBROUTINE ERF

SUBROUTINE PSI(X,PS)

  IMPLICIT NONE

  integer, parameter :: dp = selected_real_kind(33) 

  real(dp) X,PS,P,GAM,ZET(3:13),XX,Z,SUM,CA,Q,SP,TPP,TPM
  complex(kind=16)   :: EPI4
  INTEGER  :: I

  COMMON/PARS2/ GAM,ZET
  COMMON/PARS4/ P,SP,TPP,TPM,EPI4

  XX=ABS(X)

  IF (XX.LT.1e-12) THEN
     WRITE(6,*) 'PSI FUNCTION UNDEFINED AT X=0'
     STOP
  ENDIF

  IF (X.LT.0.0.AND.ABS(X-NINT(X)).LT.1e-12) THEN
     WRITE(6,*) 'PSI FUNCTION UNDEFINED AT NEGATIVE INTEGERS'
     STOP
  ENDIF

  !      SOME AWKWARD SPECIAL VALUES PSI(1)=-GAM, PSI(2), PSI(3), PSI(0.5).

  IF (ABS(XX-1.0).LT.1e-12) THEN
     PS=-GAM
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  IF (ABS(XX-2.0).LT.1e-12) THEN
     PS=0.422784335098467d0
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  IF (ABS(XX-3.0).LT.1e-12) THEN
     PS=0.922784335098467d0
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  !      PSI(1/2)=-GAM-2*LOG(2)

  IF (ABS(XX-0.5).LT.1e-12) THEN
     PS=-GAM-2*LOG(2.0)
     IF (X.LT.0.0) THEN
        GOTO 4020
     ENDIF
     RETURN
  ENDIF

  !      WORK WITH POSITIVE X VALUES VIZ. XX=ABS(X) FIRST. WORK OUT PSI(XX)

  IF (XX.LT.4.0) THEN

     !   magnitude of xx<4 use series solution valid for small xx. For positive xx=abs(x)

     IF (XX.LT.0.5) THEN
        Z=XX
        Q=-1/Z
     ELSE IF (XX.GT.0.5.AND.XX.LE.1.5) THEN
        Z=XX-1.0
        Q=0.0
     ELSE IF (XX.GT.1.5.AND.XX.LE.2.5) THEN
        Z=XX-2.0
        Q=1/(Z+1.0)
     ELSE IF (XX.GT.2.5.AND.XX.LE.3.5) THEN
        Z=XX-3.0   
        Q=1/(Z+1.0)+1/(Z+2.0)
     ELSE
        Z=XX-4.0
        Q=1/(Z+1.0)+1/(Z+2.0)+1/(Z+3.0)
     ENDIF

     !   series for psi(1+z) equation 5.7.5. psi(2+z)=psi(1+z)+1/(1+z); psi(3+z)=psi(1+z)+1/(1+z)+1/(2+z); with z an element [-0.5,0.5].

     SUM=0.0
     DO I=1,6
        SUM=SUM+(ZET(2*I+1)-1.0)*(Z**(2*I))
     ENDDO
     IF (ABS((Z-0.5)-NINT(Z-0.5)).LT.1e-12) THEN
        CA=0.0
     ELSE
        CA=0.5*P/TAN(P*Z)
     ENDIF
     PS=1/(2*Z)-CA+1/(Z*Z-1.0)+1-GAM-SUM+Q
  ELSE

     ! magnitude of xx>4 asymptotic series converges fast

     PS=LOG(XX)-0.5/XX-1/(12.0*XX*XX)+1/(120.0*(XX**4))-1/(252.0*(XX**6))+1/(240.0*(XX**8))
  ENDIF

  !  if x<0 then use reflection formula

4020 IF (X.LT.0.0)THEN
     IF (ABS((XX-0.5)-NINT(XX-0.5)).LT.1e-12) THEN
        CA=0.0
     ELSE
        CA=P/TAN(P*XX)
     ENDIF
     PS=PS+1/XX+CA
  ENDIF
  RETURN
END SUBROUTINE PSI

