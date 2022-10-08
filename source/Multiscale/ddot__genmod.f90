        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec  8 17:06:28 2016
        MODULE DDOT__genmod
          INTERFACE 
            FUNCTION DDOT(N,DX,INCX,DY,INCY)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: DX(*)
              INTEGER(KIND=4), INTENT(IN) :: INCX
              REAL(KIND=8), INTENT(IN) :: DY(*)
              INTEGER(KIND=4), INTENT(IN) :: INCY
              REAL(KIND=8) :: DDOT
            END FUNCTION DDOT
          END INTERFACE 
        END MODULE DDOT__genmod
