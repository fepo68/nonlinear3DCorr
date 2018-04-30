	SUBROUTINE DXPOSE(A, N1, N2, N12, MOVED, NWORK)
C  TRANSPOSITION OF A RECTANGULAR MATRIX IN SITU.
C  BY NORMAN BRENNER, MIT, 1/72.  CF. ALG. 380, CACM, 5/70.
C  TRANSPOSITION OF THE N1 BY N2 MATRIX A AMOUNTS TO
C  REPLACING THE ELEMENT AT VECTOR POSITION I (0-ORIGIN)
C  WITH THE ELEMENT AT POSITION N1*I (MOD N1*N2-I).
C  EACH SUBCYCLE OF THIS PERMUTATION IS COMPLETED IN ORDER.
C  MOVED IS A LOGICAL WORK ARRAY OF LENGTH NWORK.
	LOGICAL MOVED
	DOUBLE PRECISION A
	DOUBLE PRECISION ATEMP
	DOUBLE PRECISION BTEMP
	DIMENSION A(N12), MOVED(NWORK)
C  REALLY A(N1,N2), BUT N12 = N1*N2
	DIMENSION IFACT(8), IPOWER(8), NEXP(8), IEXP(8)
	IF (N1.LT.2 .OR. N2.LT.2) RETURN
	N = N1
	M = N1*N2 - 1
	IF (N1.NE.N2) GO TO 30
C  SQUARE MATRICES ARE DONE SEPARATELY FOR SPEED
	I1MIN = 2
	DO 20 I1MAX=N,M,N
	I2 = I1MIN + N - 1
	DO 10 I1=I1MIN,I1MAX
	ATEMP = A(I1)
	A(I1) = A(I2)
	A(I2) = ATEMP
	I2 = I2 + N
10	CONTINUE
	I1MIN = I1MIN + N + 1
20	CONTINUE
	RETURN
C  MODULUS M IS FACTORED INTO PRIME POWERS.  EIGHT FACTORS
C  SUFFICE UP TO M = 2*3*5*7*11*13*17*19 = 9,767,520.
30	CALL FACTOR(M, IFACT, IPOWER, NEXP, NPOWER)
	DO 40 IP=1,NPOWER
	IEXP(IP) = 0
40	CONTINUE
C  GENERATE EVERY DIVISOR OF M LESS THAN M/2
	IDIV = 1
50	IF (IDIV.GE.M/2) GO TO 190
C  THE NUMBER OF ELEMENTS WHOSE INDEX IS DIVISIBLE BY IDIV
C  AND BY NO OTHER DIVISOR OF M IS THE EULER TOTIENT
C  FUNCTION, PHI(M/IDIV).
	NCOUNT = M/IDIV
	DO 60 IP=1,NPOWER
	IF (IEXP(IP).EQ.NEXP(IP)) GO TO 60
	NCOUNT = (NCOUNT/IFACT(IP))*(IFACT(IP)-1)
60	CONTINUE
	DO 70 I=1,NWORK
	MOVED(I) = .FALSE.
70	CONTINUE
C  THE STARTING POINT OF A SUBCYCLE IS DIVISIBLE ONLY BY IDIV
C  AND MUST NOT APPEAR IN ANY OTHER SUBCYCLE.
	ISTART = IDIV
80	MMIST = M - ISTART
	IF (ISTART.EQ.IDIV) GO TO 120
	IF (ISTART.GT.NWORK) GO TO 90
	IF (MOVED(ISTART)) GO TO 160
90	ISOID = ISTART/IDIV
	DO 100 IP=1,NPOWER
	IF (IEXP(IP).EQ.NEXP(IP)) GO TO 100
	IF (MOD(ISOID,IFACT(IP)).EQ.0) GO TO 160
100	CONTINUE
	IF (ISTART.LE.NWORK) GO TO 120
	ITEST = ISTART
110	ITEST = MOD(N*ITEST,M)
	IF (ITEST.LT.ISTART .OR. ITEST.GT.MMIST) GO TO 160
	IF (ITEST.GT.ISTART .AND. ITEST.LT.MMIST) GO TO 110
120	ATEMP = A(ISTART+1)
	BTEMP = A(MMIST+I)
	IA1 = ISTART
130	IA2 = MOD(N*IA1,M)
	MMIA1 = M - IA1
	MMIA2 = M - IA2
	IF (IA1.LE.NWORK) MOVED(IA1) = .TRUE.
	IF (MMIA1.LE.NWORK) MOVED(MMIA1) = .TRUE.
	NCOUNT = NCOUNT - 2
C  MOVE TWO ELEMENTS, THE SECOND FROM THE NEGATIVE
C  SUBCYCLE.  CHECK FIRST FOR SUBCYCLE CLOSURE.
	IF (IA2.EQ.ISTART) GO TO 140
	IF (MMIA2.EQ.ISTART) GO TO 150
	A(IA1+1) = A(IA2+1)
	A(MMIA1+1) = A(MMIA2+1)
	IA1 = IA2
	GO TO 130
140	A(IA1+1) = ATEMP
	A(MMIA1+1) = BTEMP
	GO TO 160
150	A(IA1+1) = BTEMP
	A(MMIA1+1) = ATEMP
160	ISTART = ISTART + IDIV
	IF (NCOUNT.GT.0) GO TO 80
	DO 180 IP=1,NPOWER
	IF (IEXP(IP).EQ.NEXP(IP)) GO TO 170
	IEXP(IP) = IEXP(IP) + 1
	IDIV = IDIV*IFACT(IP)
	GO TO 50
170	IEXP(IP) = 0
	IDIV = IDIV/IPOWER(IP)
180	CONTINUE
190	RETURN
	END

	SUBROUTINE FACTOR(N, IFACT, IPOWER, NEXP, NPOWER)
C  FACTOR N INTO ITS PRIME POWERS, NPOWER IN NUMBER.
C  E.G., FOR N=1970=2**3 *5 *7**2, NPOWER=3, IFACT=3,5,7,
C  IPOWER=8,5,49, AND NEXP=3,1,2.
	DIMENSION IFACT(8), IPOWER(8), NEXP(8)
	IP = 0
	IFCUR = 0
	NPART = N
	IDIV = 2
10	IQUOT = NPART/IDIV
	IF (NPART-IDIV*IQUOT) 60, 20, 60
20	IF (IDIV-IFCUR) 40, 40, 30
30	IP = IP + 1
	IFACT(IP) = IDIV
	IPOWER(IP) = IDIV
	IFCUR = IDIV
	NEXP(IP) = 1
	GO TO 50
40	IPOWER(IP) = IDIV*IPOWER(IP)
	NEXP(IP) = NEXP(IP) + 1
50	NPART = IQUOT
	GO TO 10
60	IF (IQUOT-IDIV) 100, 100, 70
70	IF (IDIV-2) 80, 80, 90
80	IDIV = 3
	GO TO 10
 90	IDIV = IDIV + 2
	GO TO 10
100	IF (NPART-1) 140, 140, 110
110	IF (NPART-IFCUR) 130, 130, 120
120	IP = IP + 1
	IFACT(IP) = NPART
	IPOWER(IP) = NPART
	NEXP(IP) = 1
	GO TO 140
130	IPOWER(IP) = NPART*IPOWER(IP)
	NEXP(IP) = NEXP(IP) + 1
140	NPOWER = IP
	RETURN
	END
