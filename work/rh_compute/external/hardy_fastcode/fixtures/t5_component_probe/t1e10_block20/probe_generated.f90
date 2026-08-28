program rh_t5_component_probe
  implicit none

  integer, parameter :: dp = selected_real_kind(33)
  integer, parameter :: i128 = selected_int_kind(38)
  integer :: chain, child_length, parent_length, ip, degree, ios
  integer :: i, j, ifact, jbot, mmax, ip1, slot_count
  integer :: slot_active(3), slot_region(3), slot_index(3)
  real(dp) :: p, sp, tpp, tpm, c, gam, zet(3:13)
  real(dp) :: a1, a2, a3, xr, fracL, sx, con1, con2, con3, wm
  real(dp) :: beta, xbeta, c1q, ecor, gc, sav1, sav2, ps1, ps2, tol
  real(dp) :: erfc_args(3), erfc_values(3), phi(3)
  complex(kind=16) :: coeff(7,0:6), zlarcoeff(-1:6), zsmalcoeff(-1:6)
  complex(kind=16) :: c1, c2, c3, c4, c5, c6, c7, epi4
  complex(kind=16) :: endpoint, cr1, cr2, cr3, sum0, sum1, t5, t5psi
  complex(kind=16) :: erfc_weights(3), erfc_terms(3)
  character(len=32) :: h_a1, h_a2, h_a3, h_xr, h_frac
  character(len=32) :: h_endpoint_re, h_endpoint_im

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
     read(*,*,iostat=ios) chain, child_length, parent_length, ip, degree, &
          h_a1, h_a2, h_a3, h_xr, h_frac, h_endpoint_re, h_endpoint_im
     if (ios.lt.0) exit
     if (ios.ne.0) error stop 't5 component probe input failure'
     if (degree.ne.3) error stop 't5 component probe degree drift'
     call real_from_hex(h_a1,a1)
     call real_from_hex(h_a2,a2)
     call real_from_hex(h_a3,a3)
     call real_from_hex(h_xr,xr)
     call real_from_hex(h_frac,fracL)
     call complex_from_hex(h_endpoint_re,h_endpoint_im,endpoint)
     phi(1)=a1
     phi(2)=a2
     phi(3)=a3

     ! Replay the source state inherited by the t5 block.
     sx=sqrt(xr)
     con1=1.0-fracL
     con2=child_length+fracL+1.0
     call psi(con2,ps2)
     call psi(con1,ps1)
     sav1=ps1
     con2=ps2-ps1
     call psi(a1+1.0,ps2)
     call psi(child_length+1.0-a1,ps1)
     sav2=ps1
     con3=ps2-ps1

     slot_count=0
     slot_active=0
     slot_region=0
     slot_index=0
     erfc_args=0.0
     erfc_values=0.0
     erfc_weights=(0.0,0.0)
     erfc_terms=(0.0,0.0)

     sum0=(0.0,0.0)
     jbot=ceiling(a1)

     con1=jbot-a1
     call psi(con1,ps1)
     sav2=sav2-ps1
     con2=jbot-(child_length+fracL)
     call psi(con2,ps2)
     sav1=ps2-sav1
     t5=(0.0,1.0)*(sav2+endpoint*sav1)/tpp
     t5psi=t5

     ip1=ip
     if (ip1.gt.(0.5*child_length)) then
        ip1=max1(0.2*child_length,1.0)
     endif
     do i=jbot,ip1
        con1=i*1.0-a1
        wm=3*a3*(con1**2)/(xr**3)
        c=con1/xr-wm
        sav1=sp*sx*c*sqrt(2.0)
        slot_count=slot_count+1
        if (slot_count.gt.3) error stop 't5 component probe slot overflow'
        slot_active(slot_count)=1
        slot_region(slot_count)=0
        slot_index(slot_count)=i
        erfc_args(slot_count)=sav1
        sav1=erfc(sav1)
        erfc_values(slot_count)=sav1
        gc=i*c
        do j=1,degree
           gc=gc-phi(j)*(c**j)
        enddo
        cr1=(0.0,1.0)*tpp*gc
        cr1=exp(cr1)
        sav2=exp(tpm*con1*c)
        cr2=(0.0,1.0)*a2/p
        cr3=(1/con1+cr2*(1+tpp*c*con1*(1.0+p*con1*c))/(con1**3))
        erfc_weights(slot_count)=-cr1/(2*sx*epi4)
        erfc_terms(slot_count)=-cr1*sav1/(2*sx*epi4)
        sum0=sum0+(-cr1*sav1/(2*sx*EPI4)-(0.0,1.0)*(sav2*cr3-cr2/(con1**3))/tpp)
     enddo

     sum1=(0.0,0.0)
     c1q=child_length+fracL
     ecor=0.0
     do i=child_length,child_length-ip1+1,-1
        con1=i*1.0-a1
        wm=3*a3*(con1**2)/(xr**3)
        c=con1/xr-wm
        beta=c/parent_length
        xbeta=parent_length*(1-beta)
        sav1=sp*sx*xbeta*sqrt(2.0)
        slot_count=slot_count+1
        if (slot_count.gt.3) error stop 't5 component probe slot overflow'
        slot_active(slot_count)=1
        slot_region(slot_count)=1
        slot_index(slot_count)=i
        erfc_args(slot_count)=sav1
        sav1=erfc(sav1)
        erfc_values(slot_count)=sav1
        gc=i*c
        do j=1,degree
           gc=gc-phi(j)*(c**j)
        enddo
        cr1=(0.0,1.0)*tpp*gc
        cr1=exp(cr1)
        cr2=(0.0,1.0)*a2/p
        sav2=exp(tpm*xbeta*(c1q-i))
        cr3=(1/(c1q-i)+cr2*(1+tpp*xbeta*(c1q-i)*(1.0+p*(c1q-i)*xbeta))/((c1q-i)**3))
        erfc_weights(slot_count)=-cr1/(2*sx*epi4)
        erfc_terms(slot_count)=(-sav1)*cr1/(2*sx*epi4)
        sum1=sum1+(-sav1)*cr1/(2*sx*EPI4)
        sum1=sum1-((0.0,1.0)*endpoint*(sav2*cr3-cr2/((c1q-i)**3))/tpp)
     enddo

     t5=t5+(sum0+sum1)

     write(*,'(I0)',advance='no') chain
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(t5psi,kind=dp)),real_hex(aimag(t5psi))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(sum0,kind=dp)),real_hex(aimag(sum0))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(sum1,kind=dp)),real_hex(aimag(sum1))
     write(*,'(1X,A,1X,A)',advance='no') real_hex(real(t5,kind=dp)),real_hex(aimag(t5))
     do i=1,3
        if (i.lt.3) then
           write(*,'(1X,I0,1X,I0,1X,I0,6(1X,A))',advance='no') &
                slot_active(i),slot_region(i),slot_index(i),real_hex(erfc_args(i)), &
                real_hex(erfc_values(i)),real_hex(real(erfc_weights(i),kind=dp)), &
                real_hex(aimag(erfc_weights(i))),real_hex(real(erfc_terms(i),kind=dp)), &
                real_hex(aimag(erfc_terms(i)))
        else
           write(*,'(1X,I0,1X,I0,1X,I0,6(1X,A))') &
                slot_active(i),slot_region(i),slot_index(i),real_hex(erfc_args(i)), &
                real_hex(erfc_values(i)),real_hex(real(erfc_weights(i),kind=dp)), &
                real_hex(aimag(erfc_weights(i))),real_hex(real(erfc_terms(i),kind=dp)), &
                real_hex(aimag(erfc_terms(i)))
        endif
     enddo
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

end program rh_t5_component_probe

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

