      SUBROUTINE FOUR2 (DATA,N,NDIM,ISIGN,IFORM)   
C               
C     COOLEY-TUKEY FAST FOURIER TRANSFORM IN USASI BASIC FORTRAN.       
C     MULTI-DIMENSIONAL TRANSFORM, EACH DIMENSION A POWER OF TWO,       
C     COMPLEX OR REAL DATA.          
C     TRANSFORM(K1,K2,...) = SUM(DATA(J1,J2,...)*EXP(ISIGN*2*PI*SQRT(-1)
C     *((J1-1)*(K1-1)/N(1)+(J2-1)*(K2-1)/N(2)+...))), SUMMED FOR ALL    
C     J1 AND K1 FROM 1 TO N(1), J2 AND K2 FROM 1 TO N(2), ETC. FOR ALL  
C     NDIM SUBSCRIPTS.  NDIM MUST BE POSITIVE AND EACH N(IDIM) MUST BE  
C     A POWER OF TWO.  ISIGN IS +1 OR -1.  LET NTOT = N(1)*N(2)...      
C     ...*N(NDIM).  THEN A -1 TRANSFORM FOLLOWED BY A +1 ONE     
C     (OR VICE VERSA) RETURNS NTOT (NTOT/2 IF IFORM = 0 OR       
C     -1) TIMES THE ORIGINAL DATA.  IFORM = 1, 0 OR -1, AS DATA IS      
C     COMPLEX, REAL OR THE FIRST HALF OF A COMPLEX ARRAY.  TRANSFORM    
C     VALUES ARE RETURNED TO ARRAY DATA.  THEY ARE COMPLEX, REAL OR     
C     THE FIRST HALF OF A COMPLEX ARRAY, AS IFORM = 1, -1 OR 0.  
C     THE TRANSFORM OF A REAL ARRAY (IFORM = 0) DIMENSIONED N(1) BY N(2)
C     BY ... WILL BE RETURNED IN THE SAME ARRAY, NOW CONSIDERED TO      
C     BE COMPLEX OF DIMENSIONS N(1)/2+1 BY N(2) BY ....  NOTE THAT IF   
C     IFORM = 0 OR -1, N(1) MUST BE EVEN, AND ENOUGH ROOM MUST BE       
C     RESERVED.  THE MISSING VALUES MAY BE OBTAINED BY COMPLEX CONJU-   
C     GATION.  THE REVERSE TRANSFORMATION, OF A HALF COMPLEX ARRAY      
C     DIMENSIONED N(1)/2+1 BY N(2) BY ..., IS ACCOMPLISHED SETTING IFORM
C     TO -1.  IN THE N ARRAY, N(1) MUST BE THE TRUE N(1), NOT N(1)/2+1. 
C     THE TRANSFORM WILL BE REAL AND RETURNED TO THE INPUT ARRAY.       
C     RUNNING TIME IS PROPORTIONAL TO NTOT*LOG2(NTOT), RATHER THAN      
C     THE NAIVE NTOT**2.      
C     ORIGINAL BY NORMAN BRENNER OF MIT LINCOLN LABORATORY, JUNE 1968.  
C     SEE-- IEEE AUDIO TRANSACTIONS (JUNE 1967), SPECIAL ISSUE ON FFT.  
      DIMENSION DATA(*), N(*)        
      NTOT=1           
      DO 10 IDIM=1,NDIM       
 10   NTOT=NTOT*N(IDIM)       
      IF (IFORM) 70,20,20            
 20   NREM=NTOT        
      DO 60 IDIM=1,NDIM       
      NREM=NREM/N(IDIM)       
      NPREV=NTOT/(N(IDIM)*NREM)      
      NCURR=N(IDIM)           
      IF (IDIM-1+IFORM) 30,30,40     
 30   NCURR=NCURR/2           
 40   CALL BITRV (DATA,NPREV,NCURR,NREM)    
      CALL COOL2 (DATA,NPREV,NCURR,NREM,ISIGN)     
      IF (IDIM-1+IFORM) 50,50,60     
 50   CALL FIXRL (DATA,N(1),NREM,ISIGN,IFORM)      
      NTOT=(NTOT/N(1))*(N(1)/2+1)           

 60   CONTINUE         
      RETURN           
 70   NTOT=(NTOT/N(1))*(N(1)/2+1)           
      NREM=1           
      DO 100 JDIM=1,NDIM      
      IDIM=NDIM+1-JDIM        
      NCURR=N(IDIM)           
      IF (IDIM-1) 80,80,90           
 80   NCURR=NCURR/2           
      CALL FIXRL (DATA,N(1),NREM,ISIGN,IFORM)      
      NTOT=NTOT/(N(1)/2+1)*N(1)      
 90   NPREV=NTOT/(N(IDIM)*NREM)      
      CALL BITRV (DATA,NPREV,NCURR,NREM)    
      CALL COOL2 (DATA,NPREV,NCURR,NREM,ISIGN)     
 100  NREM=NREM*N(IDIM)       
      RETURN           
      END              
C
C
      SUBROUTINE BITRV (DATA,NPREV,N,NREM)         
C     SHUFFLE THE DATA BY @BIT REVERSAL@.          
C     DIMENSION DATA(NPREV,N,NREM)          
C     DATA(I1,I2REV,I3) = DATA(I1,I2,I3), FOR ALL I1 FROM 1 TO NPREV,   
C     ALL I2 FROM 1 TO N (WHICH MUST BE A POWER OF TWO), AND ALL I3     
C     FROM 1 TO NREM, WHERE I2REV-1 IS THE BITWISE REVERSAL OF I2-1.    
C     FOR EXAMPLE, N = 32, I2-1 = 10011 AND I2REV-1 = 11001.     
      DIMENSION DATA(*)       
      IP0=2            
      IP1=IP0*NPREV           
      IP4=IP1*N        
      IP5=IP4*NREM            
      I4REV=1          
      DO 60 I4=1,IP4,IP1      
      IF (I4-I4REV) 10,30,30         
 10   I1MAX=I4+IP1-IP0        
      DO 20 I1=I4,I1MAX,IP0          
      DO 20 I5=I1,IP5,IP4            
      I5REV=I4REV+I5-I4       
      TEMPR=DATA(I5)          
      TEMPI=DATA(I5+1)        
      DATA(I5)=DATA(I5REV)           
      DATA(I5+1)=DATA(I5REV+1)       
      DATA(I5REV)=TEMPR       
 20   DATA(I5REV+1)=TEMPI            
 30   IP2=IP4/2        
 40   IF (I4REV-IP2) 60,60,50        
 50   I4REV=I4REV-IP2         
      IP2=IP2/2        
      IF (IP2-IP1) 60,40,40          
 60   I4REV=I4REV+IP2         
      RETURN           
      END              
C
C
      SUBROUTINE COOL2 (DATA,NPREV,N,NREM,ISIGN)   
C     FOURIER TRANSFORM OF LENGTH N BY THE COOLEY-TUKEY ALGORITHM.      
C     BIT-REVERSED TO NORMAL ORDER.         
C     DIMENSION DATA(NPREV,N,NREM)          
C     COMPLEX DATA            
C     DATA(I1,J2,I3) = SUM(DATA(I1,I2,I3)*EXP(ISIGN*2*PI*I*((I2-1)*     
C     (J2-1)/N))), SUMMED OVER I2 = 1 TO N FOR ALL I1 FROM 1 TO NPREV,  
C     J2 FROM 1 TO N AND I3 FROM 1 TO NREM.  N MUST BE A POWER OF TWO.  
C     FACTORING N BY FOURS SAVES ABOUT TWENTY FIVE PERCENT OVER FACTOR- 
C     ING BY TWOS.            
C     NOTE--IT IS NOT NECESSARY TO REWRITE THIS SUBROUTINE INTO COMPLEX 
C     NOTATION SO LONG AS THE FORTRAN COMPILER USED STORES REAL AND     
C     IMAGINARY PARTS IN ADJACENT STORAGE LOCATIONS.  IT MUST ALSO      
C     STORE ARRAYS WITH THE FIRST SUBSCRIPT INCREASING FASTEST.  
      DIMENSION DATA(*)       
      TWOPI=6.283185307179586*FLOAT(ISIGN)
      IP0=2
      IP1=IP0*NPREV           
      IP4=IP1*N        
      IP5=IP4*NREM            
      IP2=IP1          
      NPART=N          
 10   IF (NPART-2) 50,30,20          
 20   NPART=NPART/4           
      GO TO 10         
C     DO A FOURIER TRANSFORM OF LENGTH TWO         
 30   IP3=IP2*2        
      DO 40 I1=1,IP1,IP0      
      DO 40 I5=I1,IP5,IP3            
      J0=I5            
      J1=J0+IP2        
      TEMPR=DATA(J1)          
      TEMPI=DATA(J1+1)        
      DATA(J1)=DATA(J0)-TEMPR        
      DATA(J1+1)=DATA(J0+1)-TEMPI           
      DATA(J0)=DATA(J0)+TEMPR        
 40   DATA(J0+1)=DATA(J0+1)+TEMPI           
      GO TO 140        
C     DO A FOURIER TRANSFORM OF LENGTH FOUR (FROM BIT REVERSED ORDER)   
 50   IP3=IP2*4        
      THETA=TWOPI/FLOAT(IP3/IP1)     
      SINTH=SIN(THETA/2.)            
      WSTPR=-2.*SINTH*SINTH          
C     COS(THETA)-1, FOR ACCURACY.           
      WSTPI=SIN(THETA)        
      WR=1.            
      WI=0.            
      DO 130 I2=1,IP2,IP1            
      IF (I2-1) 70,70,60      
 60   W2R=WR*WR-WI*WI         
      W2I=2.*WR*WI            
      W3R=W2R*WR-W2I*WI       
      W3I=W2R*WI+W2I*WR       
 70   I1MAX=I2+IP1-IP0        
      DO 120 I1=I2,I1MAX,IP0         
      DO 120 I5=I1,IP5,IP3           
      J0=I5            
      J1=J0+IP2        
      J2=J1+IP2        
      J3=J2+IP2        
      IF (I2-1) 90,90,80      
C     APPLY THE PHASE SHIFT FACTORS         
 80   TEMPR=DATA(J1)          
      DATA(J1)=W2R*TEMPR-W2I*DATA(J1+1)     
      DATA(J1+1)=W2R*DATA(J1+1)+W2I*TEMPR          
      TEMPR=DATA(J2)          
      DATA(J2)=WR*TEMPR-WI*DATA(J2+1)       
      DATA(J2+1)=WR*DATA(J2+1)+WI*TEMPR     
      TEMPR=DATA(J3)          
      DATA(J3)=W3R*TEMPR-W3I*DATA(J3+1)     
      DATA(J3+1)=W3R*DATA(J3+1)+W3I*TEMPR          
 90   T0R=DATA(J0)+DATA(J1)          
      T0I=DATA(J0+1)+DATA(J1+1)      
      T1R=DATA(J0)-DATA(J1)          
      T1I=DATA(J0+1)-DATA(J1+1)      
      T2R=DATA(J2)+DATA(J3)          
      T2I=DATA(J2+1)+DATA(J3+1)      
      T3R=DATA(J2)-DATA(J3)          
      T3I=DATA(J2+1)-DATA(J3+1)      
      DATA(J0)=T0R+T2R        
      DATA(J0+1)=T0I+T2I      
      DATA(J2)=T0R-T2R        
      DATA(J2+1)=T0I-T2I      
      IF (ISIGN) 100,100,110         
 100  T3R=-T3R         
      T3I=-T3I         
 110  DATA(J1)=T1R-T3I        
      DATA(J1+1)=T1I+T3R      
      DATA(J3)=T1R+T3I        
 120  DATA(J3+1)=T1I-T3R      
      TEMPR=WR         
      WR=WSTPR*TEMPR-WSTPI*WI+TEMPR         
 130  WI=WSTPR*WI+WSTPI*TEMPR+WI     
 140  IP2=IP3          
      IF (IP3-IP4) 50,150,150        
 150  RETURN           
      END              
C
C
      SUBROUTINE FIXRL (DATA,N,NREM,ISIGN,IFORM)   
C     FOR IFORM = 0, CONVERT THE TRANSFORM OF A DOUBLED-UP REAL ARRAY,  
C     CONSIDERED COMPLEX, INTO ITS TRUE TRANSFORM.  SUPPLY ONLY THE     
C     FIRST HALF OF THE COMPLEX TRANSFORM, AS THE SECOND HALF HAS       
C     CONJUGATE SYMMETRY.  FOR IFORM = -1, CONVERT THE FIRST HALF       
C     OF THE TRUE TRANSFORM INTO THE TRANSFORM OF A DOUBLED-UP REAL     
C     ARRAY.  N MUST BE EVEN.        
C     USING COMPLEX NOTATION AND SUBSCRIPTS STARTING AT ZERO, THE       
C     TRANSFORMATION IS--            
C     DIMENSION DATA(N,NREM)         
C     ZSTP = EXP(ISIGN*2*PI*I/N)     
C     DO 10 I2=0,NREM-1       
C     DATA(0,I2) = CONJ(DATA(0,I2))*(1+I)/(1-IFORM)       
C     DO 10 I1=1,N/4          
C     Z = (1+(2*IFORM+1)*I*ZSTP**I1)/2      
C     I1CNJ = N/2-I1          
C     DIF = DATA(I1,I2)-CONJ(DATA(I1CNJ,I2))       
C     TEMP = Z*DIF            
C     DATA(I1,I2) = DATA(I1,I2)-TEMP        
C     10  DATA(I1CNJ,I2) = DATA(I1CNJ,I2)+CONJ(TEMP)      
C     IF I1=I1CNJ, THE CALCULATION FOR THAT VALUE COLLAPSES INTO 
C     A SIMPLE CONJUGATION OF DATA(I1,I2).         
      DIMENSION DATA(*)       
      TWOPI=6.283185307179586*FLOAT(ISIGN)
      IP0=2            
      IP1=IP0*(N/2)           
      IP2=IP1*NREM            
      IF (IFORM) 10,60,60            
 10   J1=IP1+1         
      DATA(2)=DATA(J1)        
      IF(NREM-1)60,60,15      
 15   J1=J1+IP0        
      I2MIN=IP1+1             
      DO 50 I2=I2MIN,IP2,IP1         
      DATA(I2)=DATA(J1)       
      J1=J1+IP0        
      IF (N-2) 40,40,20       
 20   I1MIN=I2+IP0            
      I1MAX=I2+IP1-IP0        
      DO 30 I1=I1MIN,I1MAX,IP0       
      DATA(I1)=DATA(J1)       
      DATA(I1+1)=DATA(J1+1)          
 30   J1=J1+IP0        
 40   DATA(I2+1)=DATA(J1)            
 50   J1=J1+IP0        
 60   DO 80 I2=1,IP2,IP1      
      TEMPR=DATA(I2)          
      DATA(I2)=DATA(I2)+DATA(I2+1)          
      DATA(I2+1)=TEMPR-DATA(I2+1)           
      IF (IFORM) 70,80,80            
 70   DATA(I2)=DATA(I2)/2.           
      DATA(I2+1)=DATA(I2+1)/2.       
 80   CONTINUE         
      IF (N-2) 170,170,90            
 90   THETA=TWOPI/FLOAT(N)           
      SINTH=SIN(THETA/2.)            
      ZSTPR=-2.*SINTH*SINTH          
      ZSTPI=SIN(THETA)        
      ZR=(1.-ZSTPI)/2.        
      ZI=(1.+ZSTPR)/2.        
      IF (IFORM) 100,110,110         
 100  ZR=1.-ZR         
      ZI=-ZI           
 110  I1MIN=IP0+1             
      I1MAX=IP0*(N/4)+1       
      DO 160 I1=I1MIN,I1MAX,IP0      
      DO 150 I2=I1,IP2,IP1           
      I2CNJ=N+IP0-2*I1+I2            
      IF (I2-I2CNJ) 140,120,120      
 120  IF (ISIGN*(2*IFORM+1)) 130,150,150    
 130  DATA(I2+1)=-DATA(I2+1)         
      GO TO 150        
 140  DIFR=DATA(I2)-DATA(I2CNJ)      
      DIFI=DATA(I2+1)+DATA(I2CNJ+1)         
      TEMPR=DIFR*ZR-DIFI*ZI          
      TEMPI=DIFR*ZI+DIFI*ZR          
      DATA(I2)=DATA(I2)-TEMPR        
      DATA(I2+1)=DATA(I2+1)-TEMPI           
      DATA(I2CNJ)=DATA(I2CNJ)+TEMPR         
      DATA(I2CNJ+1)=DATA(I2CNJ+1)-TEMPI     
 150  CONTINUE         
      TEMPR=ZR-.5             
      ZR=ZSTPR*TEMPR-ZSTPI*ZI+ZR     
 160  ZI=ZSTPR*ZI+ZSTPI*TEMPR+ZI     
 170  IF (IFORM) 240,180,180         
 180  I2=IP2+1         
      I1=I2            
      J1=IP0*(N/2+1)*NREM+1          
      GO TO 220        
 190  DATA(J1)=DATA(I1)       
      DATA(J1+1)=DATA(I1+1)          
      I1=I1-IP0        
      J1=J1-IP0        
 200  IF (I2-I1) 190,210,210         
 210  DATA(J1)=DATA(I1)       
      DATA(J1+1)=0.           
 220  I2=I2-IP1        
      J1=J1-IP0        
      DATA(J1)=DATA(I2+1)            
      DATA(J1+1)=0.           
      I1=I1-IP0        
      J1=J1-IP0        
      IF (I2-1) 230,230,200          
 230  DATA(2)=0.       
 240  RETURN           
      END              
