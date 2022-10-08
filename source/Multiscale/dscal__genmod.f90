        !COMPILER-GENERATED INTERFACE MODULE: Thu Dec  8 17:06:28 2016
        MODULE DSCAL__genmod
          INTERFACE 
            SUBROUTINE DSCAL(N,SA,X,INCX)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: SA
              REAL(KIND=8), INTENT(INOUT) :: X(*)
              INTEGER(KIND=4), INTENT(IN) :: INCX
            END SUBROUTINE DSCAL
          END INTERFACE 
        END MODULE DSCAL__genmod
