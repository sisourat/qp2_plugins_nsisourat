      FUNCTION PYTHAG(A,B)
      REAL*16 PYTHAG,A,B
C
C     FINDS QSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      REAL*16 P,R,S,T,U
      P = MAX1(ABS(A),ABS(B))
      IF (P .EQ. 0.0Q0) GO TO 20
      R = (MIN1(ABS(A),ABS(B))/P)**2
   10 CONTINUE
         T = 4.0Q0 + R
         IF (T .EQ. 4.0Q0) GO TO 20
         S = R/T
         U = 1.0Q0 + 2.0Q0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END

