!> Routines needed to perform compuations with complex matrices
  module MatFunc
    !
    implicit none
    public
    !

  contains

  subroutine cholesky_dec ( UPLO, N, A, LDA, INFO )
   !======================================================================
   ! DESCRIPTION OF THE VARIABLES
   !
   !       UPLO is CHARACTER*1
   !          = 'U':  Upper triangle of A is stored;
   !          = 'L':  Lower triangle of A is stored.
   !
   !       N is INTEGER
   !          The order of the matrix A.  N >= 0.
   !
   !       A is COMPLEX*16 array, dimension (LDA,N)
   !          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
   !          N-by-N upper triangular part of A contains the upper
   !          triangular part of the matrix A, and the strictly lower
   !          triangular part of A is not referenced.  If UPLO = 'L', the
   !          leading N-by-N lower triangular part of A contains the lower
   !          triangular part of the matrix A, and the strictly upper
   !          triangular part of A is not referenced.
   !
   !          On exit, if INFO = 0, the factor U or L from the Cholesky
   !          factorization A = U**H*U or A = L*L**H.
   !
   !       LDA is INTEGER
   !          The leading dimension of the array A.  LDA >= max(1,N).
   !      INFO is INTEGER
   !          = 0:  successful exit
   !          < 0:  if INFO = -i, the i-th argument had an illegal value
   !          > 0:  if INFO = i, the leading minor of order i is not
   !                positive definite, and the factorization could not be
   !                completed.
   !======================================================================
   ! .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
   !    ..
   !    .. Array Arguments ..
       COMPLEX*16            A( lda, * )
   !   ..
   ! =====================================================================
   !
   !     .. Parameters ..
         DOUBLE PRECISION   ONE
         COMPLEX*16         CONE
         parameter( one = 1.0d0, cone = ( 1.0d0, 0.0d0 ) )
   !     ..
   !     .. Local Scalars ..
         LOGICAL            UPPER
         INTEGER            J, JB, NB
   !    ..
   !    .. External Functions ..
    !     LOGICAL            LSAME
    !     INTEGER            ILAENV
    !     EXTERNAL           lsame, ilaenv
   !     ..
   !     .. External Subroutines ..
     !    EXTERNAL           zgemm, zpotf2, zherk, ztrsm, xerbla
   !     ..
   !     .. Intrinsic Functions ..
         INTRINSIC          max, min
   !    ..
   !    .. Executable Statements ..
   !
   !     Test the input parameters.
   !
        info = 0
        upper = lsame( uplo, 'U' )
        IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
            info = -1
        ELSE IF( n.LT.0 ) THEN
            info = -2
        ELSE IF( lda.LT.max( 1, n ) ) THEN
            info = -4
        END IF
        IF( info.NE.0 ) THEN
           CALL xerbla( 'ZPOTRF', -info )
            RETURN
        END IF
   !
   !    Quick return if possible
   !
        IF( n.EQ.0 )   RETURN
   !
   !    Determine the block size for this environment.
   !
         nb = ilaenv( 1, 'ZPOTRF', uplo, n, -1, -1, -1 )
         IF( nb.LE.1 .OR. nb.GE.n ) THEN
   !
   !        Use unblocked code.
   !
            CALL zpotf2( uplo, n, a, lda, info )
         ELSE
   !
   !        Use blocked code.
   !
         IF( upper ) THEN
   !
   !          Compute the Cholesky factorization A = U'*U.
   !
           DO 10 j = 1, n, nb
   !
   !              Update and factorize the current diagonal block and test
   !              for non-positive-definiteness.
   !
                jb = min( nb, n-j+1 )
   !
                CALL zpotf2( 'Upper', jb, a( j, j ), lda, info )
   !
                IF( info.NE.0 )  GO TO 30
   !
                IF( j+jb.LE.n ) THEN
   !
   !                 Updating the trailing submatrix.
   !
                CALL ztrsm( 'Left', 'Upper', 'Conjugate Transpose', &
                               'Non-unit', jb, n-j-jb+1, cone, a( j, j ), &
                                lda, a( j, j+jb ), lda )
                CALL zherk( 'Upper', 'Conjugate transpose', n-j-jb+1, &
                               jb, -one, a( j, j+jb ), lda, &
                               one, a( j+jb, j+jb ), lda )
                END IF
     10       CONTINUE
   !
                ELSE
   !
   !           Compute the Cholesky factorization A = L*L'.
   !
                DO 20 j = 1, n, nb
   !
   !              Update and factorize the current diagonal block and test
   !              for non-positive-definiteness.
   !
                  jb = min( nb, n-j+1 )
   !
                  CALL zpotf2( 'Lower', jb, a( j, j ), lda, info )
   !
                  IF( info.NE.0 )    GO TO 30
   !
                  IF( j+jb.LE.n ) THEN
   !
   !                Updating the trailing submatrix.
   !
                  CALL ztrsm( 'Right', 'Lower', 'Conjugate Transpose', &
                             'Non-unit', n-j-jb+1, jb, cone, a( j, j ), &
                             lda, a( j+jb, j ), lda )
   !
                  CALL zherk( 'Lower', 'No Transpose', n-j-jb+1, jb, &
                             -one, a( j+jb, j ), lda, &
                             one, a( j+jb, j+jb ), lda )
                 END IF
     20       CONTINUE
            END IF
         END IF
         GO TO 40
  !
     30 CONTINUE
         info = info + j - 1
  !
     40 CONTINUE
        RETURN
  !
  !     End of cholesky_dec
  !
  end subroutine cholesky_dec

! ===================================================================
  subroutine zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!  ZGEMM  performs one of the matrix-matrix operations
!
!    C := alpha*op( A )*op( B ) + beta*C,
!
! where  op( X ) is one of
!
!    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
!
! alpha and beta are scalars, and A, B and C are matrices, with op( A )
! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!   TRANSA is CHARACTER*1
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A**T.
!
!              TRANSA = 'C' or 'c',  op( A ) = A**H.
!
!   TRANSB is CHARACTER*1
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B**T.
!
!              TRANSB = 'C' or 'c',  op( B ) = B**H.
!   M is INTEGER
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!   N is INTEGER
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!   K is INTEGER
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!   ALPHA is COMPLEX*16
!           On entry, ALPHA specifies the scalar alpha.
!                    A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!   A is COMPLEX*16 array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!   LDA is INTEGER
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!   B is COMPLEX*16 array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!   LDB is INTEGER
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!   BETA is COMPLEX*16
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!   C is COMPLEX*16 array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!   LDC is INTEGER
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!
 !
 !  -- Reference BLAS level3 routine (version 3.6.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2015
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA,BETA
       INTEGER K,LDA,LDB,LDC,M,N
       CHARACTER TRANSA,TRANSB
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),B(ldb,*),C(ldc,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
  !     LOGICAL LSAME
  !     EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
  !     EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
       LOGICAL CONJA,CONJB,NOTA,NOTB
 !     ..
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d0,0.0d0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !
 !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
 !     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
 !     B  respectively are to be  transposed but  not conjugated  and set
 !     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
 !     and the number of rows of  B  respectively.
 !
       nota = lsame(transa,'N')
       notb = lsame(transb,'N')
       conja = lsame(transa,'C')
       conjb = lsame(transb,'C')
       IF (nota) THEN
           nrowa = m
           ncola = k
       ELSE
           nrowa = k
           ncola = m
       END IF
       IF (notb) THEN
           nrowb = k
       ELSE
           nrowb = n
       END IF
 !
 !     Test the input parameters.
 !
       info = 0
       IF ((.NOT.nota) .AND. (.NOT.conja) .AND. &
           (.NOT.lsame(transa,'T'))) THEN
           info = 1
       ELSE IF ((.NOT.notb) .AND. (.NOT.conjb) .AND. &
           (.NOT.lsame(transb,'T'))) THEN
           info = 2
       ELSE IF (m.LT.0) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (k.LT.0) THEN
           info = 5
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 8
       ELSE IF (ldb.LT.max(1,nrowb)) THEN
           info = 10
       ELSE IF (ldc.LT.max(1,m)) THEN
           info = 13
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZGEMM ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (((alpha.EQ.zero).OR. &
          (k.EQ.0)).AND. (beta.EQ.one))) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           IF (beta.EQ.zero) THEN
               DO 20 j = 1,n
                   DO 10 i = 1,m
                       c(i,j) = zero
    10             CONTINUE
    20         CONTINUE
           ELSE
               DO 40 j = 1,n
                   DO 30 i = 1,m
                       c(i,j) = beta*c(i,j)
    30             CONTINUE
    40         CONTINUE
           END IF
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (notb) THEN
           IF (nota) THEN
 !
 !           Form  C := alpha*A*B + beta*C.
 !
               DO 90 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 50 i = 1,m
                           c(i,j) = zero
    50                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 60 i = 1,m
                           c(i,j) = beta*c(i,j)
    60                 CONTINUE
                   END IF
                   DO 80 l = 1,k
                       temp = alpha*b(l,j)
                       DO 70 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
    70                 CONTINUE
    80             CONTINUE
    90         CONTINUE
           ELSE IF (conja) THEN
 !
 !           Form  C := alpha*A**H*B + beta*C.
 !
               DO 120 j = 1,n
                   DO 110 i = 1,m
                       temp = zero
                       DO 100 l = 1,k
                           temp = temp + dconjg(a(l,i))*b(l,j)
   100                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   110             CONTINUE
   120         CONTINUE
           ELSE
 !
 !           Form  C := alpha*A**T*B + beta*C
 !
               DO 150 j = 1,n
                   DO 140 i = 1,m
                       temp = zero
                       DO 130 l = 1,k
                           temp = temp + a(l,i)*b(l,j)
   130                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   140             CONTINUE
   150         CONTINUE
           END IF
       ELSE IF (nota) THEN
           IF (conjb) THEN
 !
 !           Form  C := alpha*A*B**H + beta*C.
 !
               DO 200 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 160 i = 1,m
                           c(i,j) = zero
   160                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 170 i = 1,m
                           c(i,j) = beta*c(i,j)
   170                 CONTINUE
                   END IF
                   DO 190 l = 1,k
                       temp = alpha*dconjg(b(j,l))
                       DO 180 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
   180                 CONTINUE
   190             CONTINUE
   200         CONTINUE
           ELSE
 !
 !           Form  C := alpha*A*B**T + beta*C
 !
               DO 250 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 210 i = 1,m
                           c(i,j) = zero
   210                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 220 i = 1,m
                           c(i,j) = beta*c(i,j)
   220                 CONTINUE
                   END IF
                   DO 240 l = 1,k
                       temp = alpha*b(j,l)
                       DO 230 i = 1,m
                           c(i,j) = c(i,j) + temp*a(i,l)
   230                 CONTINUE
   240             CONTINUE
   250         CONTINUE
           END IF
       ELSE IF (conja) THEN
           IF (conjb) THEN
 !
 !           Form  C := alpha*A**H*B**H + beta*C.
 !
               DO 280 j = 1,n
                   DO 270 i = 1,m
                       temp = zero
                       DO 260 l = 1,k
                           temp = temp + dconjg(a(l,i))*dconjg(b(j,l))
   260                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   270             CONTINUE
   280         CONTINUE
           ELSE
 !
 !           Form  C := alpha*A**H*B**T + beta*C
 !
               DO 310 j = 1,n
                   DO 300 i = 1,m
                       temp = zero
                       DO 290 l = 1,k
                           temp = temp + dconjg(a(l,i))*b(j,l)
   290                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   300             CONTINUE
   310         CONTINUE
           END IF
       ELSE
           IF (conjb) THEN
 !
 !           Form  C := alpha*A**T*B**H + beta*C
 !
               DO 340 j = 1,n
                   DO 330 i = 1,m
                       temp = zero
                       DO 320 l = 1,k
                           temp = temp + a(l,i)*dconjg(b(j,l))
   320                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   330             CONTINUE
   340         CONTINUE
           ELSE
 !
 !           Form  C := alpha*A**T*B**T + beta*C
 !
               DO 370 j = 1,n
                   DO 360 i = 1,m
                       temp = zero
                       DO 350 l = 1,k
                           temp = temp + a(l,i)*b(j,l)
   350                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   360             CONTINUE
   370         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZGEMM .
 !
  end subroutine zgemm

 !===================================================================
  subroutine xerbla( SRNAME, INFO )
 !
 !  -- LAPACK test routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER*(*)      SRNAME
       INTEGER            INFO
 !     ..
 !
 !  =====================================================================
 !
 !     .. Scalars in Common ..
       LOGICAL            LERR, OK
       CHARACTER*32       SRNAMT
       INTEGER            INFOT, NOUT
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          len_trim
 !     ..
 !     .. Common blocks ..
       COMMON             / infoc / infot, nout, ok, lerr
       COMMON             / srnamc / srnamt
 !     ..
 !     .. Executable Statements ..
 !
       lerr = .true.
       IF( info.NE.infot ) THEN
          IF( infot.NE.0 ) THEN
             WRITE( nout, fmt = 9999 ) &
           srnamt( 1:len_trim( srnamt ) ), info, infot
          ELSE
             WRITE( nout, fmt = 9997 ) &
           srname( 1:len_trim( srname ) ), info
          END IF
          ok = .false.
       END IF
       IF( srname.NE.srnamt ) THEN
          WRITE( nout, fmt = 9998 ) &
           srname( 1:len_trim( srname ) ), &
           srnamt( 1:len_trim( srnamt ) )
          ok = .false.
       END IF
       RETURN
 !
  9999 FORMAT( ' *** XERBLA was called from ', a, ' with INFO = ', i6, &
            ' instead of ', i2, ' ***' )
  9998 FORMAT( ' *** XERBLA was called with SRNAME = ', a, &
            ' instead of ', a6, ' ***' )
  9997 FORMAT( ' *** On entry to ', a, ' parameter number ', i6, &
            ' had an illegal value ***' )
 !
 !     End of XERBLA
 !
 end subroutine xerbla


 !===================================================================
  logical function lsame(CA,CB)
 !
 !  -- Reference BLAS level1 routine (version 3.1) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER CA,CB
 !     ..
 !
 ! =====================================================================
 !
 !     .. Intrinsic Functions ..
       INTRINSIC ichar
 !     ..
 !     .. Local Scalars ..
       INTEGER INTA,INTB,ZCODE
 !     ..
 !
 !     Test if the characters are equal
 !
       lsame = ca .EQ. cb
       IF (lsame) RETURN
 !
 !     Now test for equivalence if both characters are alphabetic.
 !
       zcode = ichar('Z')
 !
 !     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
 !     machines, on which ICHAR returns a value with bit 8 set.
 !     ICHAR('A') on Prime machines returns 193 which is the same as
 !     ICHAR('A') on an EBCDIC machine.
 !
       inta = ichar(ca)
       intb = ichar(cb)
 !
       IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
 !
 !        ASCII is assumed - ZCODE is the ASCII code of either lower or
 !        upper case 'Z'.
 !
           IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
           IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
 !
       ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
 !
 !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
 !        upper case 'Z'.
 !
           IF (inta.GE.129 .AND. inta.LE.137 .OR. &
               inta.GE.145 .AND. inta.LE.153 .OR. &
               inta.GE.162 .AND. inta.LE.169) inta = inta + 64
           IF (intb.GE.129 .AND. intb.LE.137 .OR. &
               intb.GE.145 .AND. intb.LE.153 .OR. &
               intb.GE.162 .AND. intb.LE.169) intb = intb + 64
 !
       ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
 !
 !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
 !        plus 128 of either lower or upper case 'Z'.
 !
           IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
           IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
       END IF
       lsame = inta .EQ. intb
 !
 !     RETURN
 !
 !     End of LSAME
 !
  end function lsame

! ====================================================================
  integer function ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
  !
  !  -- LAPACK auxiliary routine (version 3.6.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     November 2015
  !
  !     .. Scalar Arguments ..
        CHARACTER*( * )    NAME, OPTS
        INTEGER            ISPEC, N1, N2, N3, N4
  !     ..
  !
  !  =====================================================================
  !
  !     .. Local Scalars ..
        INTEGER            I, IC, IZ, NB, NBMIN, NX
        LOGICAL            CNAME, SNAME
        CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*6
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          char, ichar, int, min, real
  !     ..
  !     .. External Functions ..
   !     INTEGER            IEEECK, IPARMQ
   !     EXTERNAL           ieeeck, iparmq
  !     ..
  !     .. Executable Statements ..
  !
        GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, &
                130, 140, 150, 160, 160, 160, 160, 160 )ispec
  !
  !     Invalid value for ISPEC
  !
        ilaenv = -1
        RETURN
  !
     10 CONTINUE
  !
  !     Convert NAME to upper case if the first character is lower case.
  !
        ilaenv = 1
        subnam = name
        ic = ichar( subnam( 1: 1 ) )
        iz = ichar( 'Z' )
        IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
  !
  !        ASCII character set
  !
           IF( ic.GE.97 .AND. ic.LE.122 ) THEN
              subnam( 1: 1 ) = char( ic-32 )
              DO 20 i = 2, 6
                 ic = ichar( subnam( i: i ) )
                 IF( ic.GE.97 .AND. ic.LE.122 ) &
                    subnam( i: i ) = char( ic-32 )
     20       CONTINUE
           END IF
  !
        ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
  !
  !        EBCDIC character set
  !
           IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
               ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
               ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
              subnam( 1: 1 ) = char( ic+64 )
              DO 30 i = 2, 6
                 ic = ichar( subnam( i: i ) )
                 IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                     ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                     ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i: &
                     i ) = char( ic+64 )
     30       CONTINUE
           END IF
  !
        ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
  !
  !        Prime machines:  ASCII+128
  !
           IF( ic.GE.225 .AND. ic.LE.250 ) THEN
              subnam( 1: 1 ) = char( ic-32 )
              DO 40 i = 2, 6
                 ic = ichar( subnam( i: i ) )
                 IF( ic.GE.225 .AND. ic.LE.250 ) &
                    subnam( i: i ) = char( ic-32 )
     40       CONTINUE
           END IF
        END IF
  !
        c1 = subnam( 1: 1 )
        sname = c1.EQ.'S' .OR. c1.EQ.'D'
        cname = c1.EQ.'C' .OR. c1.EQ.'Z'
        IF( .NOT.( cname .OR. sname ) ) &
           RETURN
        c2 = subnam( 2: 3 )
        c3 = subnam( 4: 6 )
        c4 = c3( 2: 3 )
  !
        GO TO ( 50, 60, 70 )ispec
  !
     50 CONTINUE
  !
  !     ISPEC = 1:  block size
  !
  !     In these examples, separate code is provided for setting NB for
  !     real and complex.  We assume that NB will take the same value in
  !     single or double precision.
  !
        nb = 1
  !
        IF( c2.EQ.'GE' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              IF( sname ) THEN
                 nb = 64
              ELSE
                 nb = 64
              END IF
           ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. &
                    c3.EQ.'QLF' ) THEN
              IF( sname ) THEN
                 nb = 32
              ELSE
                 nb = 32
              END IF
           ELSE IF( c3.EQ.'HRD' ) THEN
              IF( sname ) THEN
                 nb = 32
              ELSE
                 nb = 32
              END IF
           ELSE IF( c3.EQ.'BRD' ) THEN
              IF( sname ) THEN
                 nb = 32
              ELSE
                 nb = 32
              END IF
           ELSE IF( c3.EQ.'TRI' ) THEN
              IF( sname ) THEN
                 nb = 64
              ELSE
                 nb = 64
              END IF
           END IF
        ELSE IF( c2.EQ.'PO' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              IF( sname ) THEN
                 nb = 64
              ELSE
                 nb = 64
              END IF
           END IF
        ELSE IF( c2.EQ.'SY' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              IF( sname ) THEN
                 nb = 64
              ELSE
                 nb = 64
              END IF
           ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
              nb = 32
           ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
              nb = 64
           END IF
        ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              nb = 64
           ELSE IF( c3.EQ.'TRD' ) THEN
              nb = 32
           ELSE IF( c3.EQ.'GST' ) THEN
              nb = 64
           END IF
        ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
           IF( c3( 1: 1 ).EQ.'G' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nb = 32
              END IF
           ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nb = 32
              END IF
           END IF
        ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
           IF( c3( 1: 1 ).EQ.'G' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nb = 32
              END IF
           ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nb = 32
              END IF
           END IF
        ELSE IF( c2.EQ.'GB' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              IF( sname ) THEN
                 IF( n4.LE.64 ) THEN
                    nb = 1
                 ELSE
                    nb = 32
                 END IF
              ELSE
                 IF( n4.LE.64 ) THEN
                    nb = 1
                 ELSE
                    nb = 32
                 END IF
              END IF
           END IF
        ELSE IF( c2.EQ.'PB' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              IF( sname ) THEN
                 IF( n2.LE.64 ) THEN
                    nb = 1
                 ELSE
                    nb = 32
                 END IF
              ELSE
                 IF( n2.LE.64 ) THEN
                    nb = 1
                 ELSE
                    nb = 32
                 END IF
              END IF
           END IF
        ELSE IF( c2.EQ.'TR' ) THEN
           IF( c3.EQ.'TRI' ) THEN
              IF( sname ) THEN
                 nb = 64
              ELSE
                 nb = 64
              END IF
           END IF
        ELSE IF( c2.EQ.'LA' ) THEN
           IF( c3.EQ.'UUM' ) THEN
              IF( sname ) THEN
                 nb = 64
              ELSE
                 nb = 64
              END IF
           END IF
        ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
           IF( c3.EQ.'EBZ' ) THEN
              nb = 1
           END IF
        ELSE IF( c2.EQ.'GG' ) THEN
           nb = 32
           IF( c3.EQ.'HD3' ) THEN
              IF( sname ) THEN
                 nb = 32
              ELSE
                 nb = 32
              END IF
           END IF
        END IF
        ilaenv = nb
        RETURN
  !
     60 CONTINUE
  !
  !     ISPEC = 2:  minimum block size
  !
        nbmin = 2
        IF( c2.EQ.'GE' ) THEN
           IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. &
              'QLF' ) THEN
              IF( sname ) THEN
                 nbmin = 2
              ELSE
                 nbmin = 2
              END IF
           ELSE IF( c3.EQ.'HRD' ) THEN
              IF( sname ) THEN
                 nbmin = 2
              ELSE
                 nbmin = 2
              END IF
           ELSE IF( c3.EQ.'BRD' ) THEN
              IF( sname ) THEN
                 nbmin = 2
              ELSE
                 nbmin = 2
              END IF
           ELSE IF( c3.EQ.'TRI' ) THEN
              IF( sname ) THEN
                 nbmin = 2
              ELSE
                 nbmin = 2
              END IF
           END IF
        ELSE IF( c2.EQ.'SY' ) THEN
           IF( c3.EQ.'TRF' ) THEN
              IF( sname ) THEN
                 nbmin = 8
              ELSE
                 nbmin = 8
              END IF
           ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
              nbmin = 2
           END IF
        ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
           IF( c3.EQ.'TRD' ) THEN
              nbmin = 2
           END IF
        ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
           IF( c3( 1: 1 ).EQ.'G' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nbmin = 2
              END IF
           ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nbmin = 2
              END IF
           END IF
        ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
           IF( c3( 1: 1 ).EQ.'G' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nbmin = 2
              END IF
           ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nbmin = 2
              END IF
           END IF
        ELSE IF( c2.EQ.'GG' ) THEN
           nbmin = 2
           IF( c3.EQ.'HD3' ) THEN
              nbmin = 2
           END IF
        END IF
        ilaenv = nbmin
        RETURN
  !
     70 CONTINUE
  !
  !     ISPEC = 3:  crossover point
  !
        nx = 0
        IF( c2.EQ.'GE' ) THEN
           IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ. &
              'QLF' ) THEN
              IF( sname ) THEN
                 nx = 128
              ELSE
                 nx = 128
              END IF
           ELSE IF( c3.EQ.'HRD' ) THEN
              IF( sname ) THEN
                 nx = 128
              ELSE
                 nx = 128
              END IF
           ELSE IF( c3.EQ.'BRD' ) THEN
              IF( sname ) THEN
                 nx = 128
              ELSE
                 nx = 128
              END IF
           END IF
        ELSE IF( c2.EQ.'SY' ) THEN
           IF( sname .AND. c3.EQ.'TRD' ) THEN
              nx = 32
           END IF
        ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
           IF( c3.EQ.'TRD' ) THEN
              nx = 32
           END IF
        ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
           IF( c3( 1: 1 ).EQ.'G' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nx = 128
              END IF
           END IF
        ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
           IF( c3( 1: 1 ).EQ.'G' ) THEN
              IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ. &
                 'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' ) &
                  THEN
                 nx = 128
              END IF
           END IF
        ELSE IF( c2.EQ.'GG' ) THEN
           nx = 128
           IF( c3.EQ.'HD3' ) THEN
              nx = 128
           END IF
        END IF
        ilaenv = nx
        RETURN
  !
     80 CONTINUE
  !
  !     ISPEC = 4:  number of shifts (used by xHSEQR)
  !
        ilaenv = 6
        RETURN
  !
     90 CONTINUE
  !
  !     ISPEC = 5:  minimum column dimension (not used)
  !
        ilaenv = 2
        RETURN
  !
    100 CONTINUE
  !
  !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
  !
        ilaenv = int( REAL( MIN( N1, N2 ) )*1.6e0 )
        RETURN
  !
    110 CONTINUE
  !
  !     ISPEC = 7:  number of processors (not used)
  !
        ilaenv = 1
        RETURN
  !
    120 CONTINUE
  !
  !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
  !
        ilaenv = 50
        RETURN
  !
    130 CONTINUE
  !
  !     ISPEC = 9:  maximum size of the subproblems at the bottom of the
  !                 computation tree in the divide-and-conquer algorithm
  !                 (used by xGELSD and xGESDD)
  !
        ilaenv = 25
        RETURN
  !
    140 CONTINUE
  !
  !     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
  !
  !     ILAENV = 0
        ilaenv = 1
        IF( ilaenv.EQ.1 ) THEN
           ilaenv = ieeeck( 1, 0.0, 1.0 )
        END IF
        RETURN
  !
    150 CONTINUE
  !
  !     ISPEC = 11: infinity arithmetic can be trusted not to trap
  !
  !     ILAENV = 0
        ilaenv = 1
        IF( ilaenv.EQ.1 ) THEN
           ilaenv = ieeeck( 0, 0.0, 1.0 )
        END IF
        RETURN
  !
    160 CONTINUE
  !
  !     12 <= ISPEC <= 16: xHSEQR or related subroutines.
  !
        ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
        RETURN
  !
  !     End of ILAENV
  !
  end function ilaenv

  !===================================================================
  integer function  ieeeck( ISPEC, ZERO, ONE )
  !
  !  -- LAPACK auxiliary routine (version 3.4.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     November 2011
  !
  !     .. Scalar Arguments ..
        INTEGER            ISPEC
        REAL               ONE, ZERO
  !     ..
  !
  !  =====================================================================
  !
  !     .. Local Scalars ..
        REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF, &
                           negzro, newzro, posinf
  !     ..
  !     .. Executable Statements ..
        ieeeck = 1
  !
        posinf = one / zero
        IF( posinf.LE.one ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        neginf = -one / zero
        IF( neginf.GE.zero ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        negzro = one / ( neginf+one )
        IF( negzro.NE.zero ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        neginf = one / negzro
        IF( neginf.GE.zero ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        newzro = negzro + zero
        IF( newzro.NE.zero ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        posinf = one / newzro
        IF( posinf.LE.one ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        neginf = neginf*posinf
        IF( neginf.GE.zero ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        posinf = posinf*posinf
        IF( posinf.LE.one ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
  !
  !
  !
  !     Return if we were only asked to check infinity arithmetic
  !
        IF( ispec.EQ.0 ) RETURN
  !
        nan1 = posinf + neginf
  !
        nan2 = posinf / neginf
  !
        nan3 = posinf / posinf
  !
        nan4 = posinf*zero
  !
        nan5 = neginf*negzro
  !
        nan6 = nan5*zero
  !
        IF( nan1.EQ.nan1 ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        IF( nan2.EQ.nan2 ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        IF( nan3.EQ.nan3 ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        IF( nan4.EQ.nan4 ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        IF( nan5.EQ.nan5 ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        IF( nan6.EQ.nan6 ) THEN
           ieeeck = 0
           RETURN
        END IF
  !
        RETURN
  end function ieeeck

 !  =====================================================================
  integer function iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
 !
 !  -- LAPACK auxiliary routine (version 3.6.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2015
 !
 !     .. Scalar Arguments ..
       INTEGER            IHI, ILO, ISPEC, LWORK, N
       CHARACTER          NAME*( * ), OPTS*( * )
 !
 !  ================================================================
 !     .. Parameters ..
       INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
       parameter( inmin = 12, inwin = 13, inibl = 14, ishfts = 15, iacc22 = 16 )
       INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
       parameter( nmin = 75, k22min = 14, kacmin = 14, nibble = 14, knwswp = 500 )
       REAL               TWO
       parameter( two = 2.0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            NH, NS
       INTEGER            I, IC, IZ
       CHARACTER          SUBNAM*6
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          log, max, mod, nint, real
 !     ..
 !     .. Executable Statements ..
       IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR. &
           ( ispec.EQ.iacc22 ) ) THEN
 !
 !        ==== Set the number simultaneous shifts ====
 !
          nh = ihi - ilo + 1
          ns = 2
          IF( nh.GE.30 ) &
             ns = 4
          IF( nh.GE.60 ) &
             ns = 10
          IF( nh.GE.150 ) &
             ns = max( 10, nh / nint( log( REAL( NH ) ) / log( TWO ) ) )
          IF( nh.GE.590 ) &
             ns = 64
          IF( nh.GE.3000 ) &
             ns = 128
          IF( nh.GE.6000 ) &
             ns = 256
          ns = max( 2, ns-mod( ns, 2 ) )
       END IF
 !
       IF( ispec.EQ.inmin ) THEN
 !
 !
 !        ===== Matrices of order smaller than NMIN get sent
 !        .     to xLAHQR, the classic double shift algorithm.
 !        .     This must be at least 11. ====
 !
          iparmq = nmin
 !
       ELSE IF( ispec.EQ.inibl ) THEN
 !
 !        ==== INIBL: skip a multi-shift qr iteration and
 !        .    whenever aggressive early deflation finds
 !        .    at least (NIBBLE*(window size)/100) deflations. ====
 !
          iparmq = nibble
 !
       ELSE IF( ispec.EQ.ishfts ) THEN
 !
 !        ==== NSHFTS: The number of simultaneous shifts =====
 !
          iparmq = ns
 !
       ELSE IF( ispec.EQ.inwin ) THEN
 !
 !        ==== NW: deflation window size.  ====
 !
          IF( nh.LE.knwswp ) THEN
             iparmq = ns
          ELSE
             iparmq = 3*ns / 2
          END IF
 !
       ELSE IF( ispec.EQ.iacc22 ) THEN
 !
 !        ==== IACC22: Whether to accumulate reflections
 !        .     before updating the far-from-diagonal elements
 !        .     and whether to use 2-by-2 block structure while
 !        .     doing it.  A small amount of work could be saved
 !        .     by making this choice dependent also upon the
 !        .     NH=IHI-ILO+1.
 !
 !
 !        Convert NAME to upper case if the first character is lower case.
 !
          iparmq = 0
          subnam = name
          ic = ichar( subnam( 1: 1 ) )
          iz = ichar( 'Z' )
          IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
 !
 !           ASCII character set
 !
             IF( ic.GE.97 .AND. ic.LE.122 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.97 .AND. ic.LE.122 ) &
                      subnam( i: i ) = char( ic-32 )
                END DO
             END IF
 !
          ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
 !
 !           EBCDIC character set
 !
             IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                 ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                 ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
                subnam( 1: 1 ) = char( ic+64 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR. &
                       ( ic.GE.145 .AND. ic.LE.153 ) .OR. &
                       ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i: &
                       i ) = char( ic+64 )
                END DO
             END IF
 !
          ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
 !
 !           Prime machines:  ASCII+128
 !
             IF( ic.GE.225 .AND. ic.LE.250 ) THEN
                subnam( 1: 1 ) = char( ic-32 )
                DO i = 2, 6
                   ic = ichar( subnam( i: i ) )
                   IF( ic.GE.225 .AND. ic.LE.250 ) &
                      subnam( i: i ) = char( ic-32 )
                END DO
             END IF
          END IF
 !
          IF( subnam( 2:6 ).EQ.'GGHRD' .OR. &
              subnam( 2:6 ).EQ.'GGHD3' ) THEN
             iparmq = 1
             IF( nh.GE.k22min ) &
                iparmq = 2
          ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
             IF( nh.GE.kacmin ) &
                iparmq = 1
             IF( nh.GE.k22min ) &
                iparmq = 2
          ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR. &
                    subnam( 2:5 ).EQ.'LAQR' ) THEN
             IF( ns.GE.kacmin ) &
                iparmq = 1
             IF( ns.GE.k22min ) &
                iparmq = 2
          END IF
 !
       ELSE
 !        ===== invalid value of ispec =====
          iparmq = -1
 !
       END IF
 !
 !     ==== End of IPARMQ ====
 !
  end function iparmq

  !====================================================================
 subroutine zpotf2( UPLO, N, A, LDA, INFO )
  !
  !  -- LAPACK computational routine (version 3.4.2) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     September 2012
  !
  !     .. Scalar Arguments ..
        CHARACTER          UPLO
        INTEGER            INFO, LDA, N
  !     ..
  !     .. Array Arguments ..
        COMPLEX*16         A( lda, * )
  !     ..
  !
  !  =====================================================================
  !
  !     .. Parameters ..
        DOUBLE PRECISION   ONE, ZERO
        parameter( one = 1.0d0, zero = 0.0d0 )
        COMPLEX*16         CONE
        parameter( cone = ( 1.0d0, 0.0d0 ) )
  !     ..
  !     .. Local Scalars ..
        LOGICAL            UPPER
        INTEGER            J
        DOUBLE PRECISION   AJJ
  !     ..
  !     .. External Functions ..
   !     LOGICAL            LSAME, DISNAN
   !     COMPLEX*16         ZDOTC
   !     EXTERNAL           lsame, zdotc, disnan
  !     ..
  !     .. External Subroutines ..
   !     EXTERNAL           xerbla, zdscal, zgemv, zlacgv
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          dble, max, sqrt
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
        info = 0
        upper = lsame( uplo, 'U' )
        IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
           info = -1
        ELSE IF( n.LT.0 ) THEN
           info = -2
        ELSE IF( lda.LT.max( 1, n ) ) THEN
           info = -4
        END IF
        IF( info.NE.0 ) THEN
           CALL xerbla( 'ZPOTF2', -info )
           RETURN
        END IF
  !
  !     Quick return if possible
  !
        IF( n.EQ.0 )   RETURN
  !
        IF( upper ) THEN
  !
  !        Compute the Cholesky factorization A = U**H *U.
  !
           DO 10 j = 1, n
  !
  !           Compute U(J,J) and test for non-positive-definiteness.
  !
              ajj = dble( a( j, j ) ) - zdotc( j-1, a( 1, j ), 1, &
                    a( 1, j ), 1 )
              IF( ajj.LE.zero.OR.disnan( ajj ) ) THEN
                 a( j, j ) = ajj
                 GO TO 30
              END IF
              ajj = sqrt( ajj )
              a( j, j ) = ajj
  !
  !           Compute elements J+1:N of row J.
  !
              IF( j.LT.n ) THEN
                 CALL zlacgv( j-1, a( 1, j ), 1 )
                 CALL zgemv( 'Transpose', j-1, n-j, -cone, a( 1, j+1 ), &
                             lda, a( 1, j ), 1, cone, a( j, j+1 ), lda )
                 CALL zlacgv( j-1, a( 1, j ), 1 )
                 CALL zdscal( n-j, one / ajj, a( j, j+1 ), lda )
              END IF
     10    CONTINUE
        ELSE
  !
  !        Compute the Cholesky factorization A = L*L**H.
  !
           DO 20 j = 1, n
  !
  !           Compute L(J,J) and test for non-positive-definiteness.
  !
              ajj = dble( a( j, j ) ) - zdotc( j-1, a( j, 1 ), lda, &
                    a( j, 1 ), lda )
              IF( ajj.LE.zero.OR.disnan( ajj ) ) THEN
                 a( j, j ) = ajj
                 GO TO 30
              END IF
              ajj = sqrt( ajj )
              a( j, j ) = ajj
  !
  !           Compute elements J+1:N of column J.
  !
              IF( j.LT.n ) THEN
                 CALL zlacgv( j-1, a( j, 1 ), lda )
                 CALL zgemv( 'No transpose', n-j, j-1, -cone, a( j+1, 1 ), &
                             lda, a( j, 1 ), lda, cone, a( j+1, j ), 1 )
                 CALL zlacgv( j-1, a( j, 1 ), lda )
                 CALL zdscal( n-j, one / ajj, a( j+1, j ), 1 )
              END IF
     20    CONTINUE
        END IF
        GO TO 40
  !
     30 CONTINUE
        info = j
  !
     40 CONTINUE
        RETURN
  !
  !     End of ZPOTF2
  !
  end subroutine zpotf2

  !====================================================================
  complex*16 function zdotc(N,ZX,INCX,ZY,INCY)
  !
  !  -- Reference BLAS level1 routine (version 3.6.0) --
  !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     November 2015
  !
  !     .. Scalar Arguments ..
        INTEGER INCX,INCY,N
  !     ..
  !     .. Array Arguments ..
        COMPLEX*16 ZX(*),ZY(*)
  !     ..
  !
  !  =====================================================================
  !
  !     .. Local Scalars ..
        COMPLEX*16 ZTEMP
        INTEGER I,IX,IY
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC dconjg
  !     ..
        ztemp = (0.0d0,0.0d0)
        zdotc = (0.0d0,0.0d0)
        IF (n.LE.0) RETURN
        IF (incx.EQ.1 .AND. incy.EQ.1) THEN
  !
  !        code for both increments equal to 1
  !
           DO i = 1,n
              ztemp = ztemp + dconjg(zx(i))*zy(i)
           END DO
        ELSE
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
           ix = 1
           iy = 1
           IF (incx.LT.0) ix = (-n+1)*incx + 1
           IF (incy.LT.0) iy = (-n+1)*incy + 1
           DO i = 1,n
              ztemp = ztemp + dconjg(zx(ix))*zy(iy)
              ix = ix + incx
              iy = iy + incy
           END DO
        END IF
        zdotc = ztemp
        RETURN
  end function zdotc

 !===================================================================
  logical function disnan( DIN )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   DIN
 !     ..
 !
 !  =====================================================================
 !
 !  .. External Functions ..
  !     LOGICAL DLAISNAN
  !     EXTERNAL dlaisnan
 !  ..
 !  .. Executable Statements ..
       disnan = dlaisnan(din,din)
       RETURN
  end function disnan

 !=====================================================================
  logical function dlaisnan( DIN1, DIN2 )

 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   DIN1, DIN2
 !     ..
 !
 !  =====================================================================
 !
 !  .. Executable Statements ..
       dlaisnan = (din1.NE.din2)
       RETURN
  end function dlaisnan

 !  =====================================================================
  subroutine zdscal(N,DA,ZX,INCX)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION DA
       INTEGER INCX,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,NINCX
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dcmplx
 !     ..
       IF (n.LE.0 .OR. incx.LE.0) RETURN
       IF (incx.EQ.1) THEN
 !
 !        code for increment equal to 1
 !
          DO i = 1,n
             zx(i) = dcmplx(da,0.0d0)*zx(i)
          END DO
       ELSE
 !
 !        code for increment not equal to 1
 !
          nincx = n*incx
          DO i = 1,nincx,incx
             zx(i) = dcmplx(da,0.0d0)*zx(i)
          END DO
       END IF
       RETURN
  end subroutine zdscal

 !  =====================================================================
  subroutine zgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
 !
 !  -- Reference BLAS level2 routine (version 3.6.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2015
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA,BETA
       INTEGER INCX,INCY,LDA,M,N
       CHARACTER TRANS
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),X(*),Y(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d0,0.0d0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
       LOGICAL NOCONJ
 !     ..
 !     .. External Functions ..
  !     LOGICAL LSAME
  !     EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
  !     EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. &
           .NOT.lsame(trans,'C')) THEN
           info = 1
       ELSE IF (m.LT.0) THEN
           info = 2
       ELSE IF (n.LT.0) THEN
           info = 3
       ELSE IF (lda.LT.max(1,m)) THEN
           info = 6
       ELSE IF (incx.EQ.0) THEN
           info = 8
       ELSE IF (incy.EQ.0) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZGEMV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. &
           ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
 !
       noconj = lsame(trans,'T')
 !
 !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
 !     up the start points in  X  and  Y.
 !
       IF (lsame(trans,'N')) THEN
           lenx = n
           leny = m
       ELSE
           lenx = m
           leny = n
       END IF
       IF (incx.GT.0) THEN
           kx = 1
       ELSE
           kx = 1 - (lenx-1)*incx
       END IF
       IF (incy.GT.0) THEN
           ky = 1
       ELSE
           ky = 1 - (leny-1)*incy
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
 !     First form  y := beta*y.
 !
       IF (beta.NE.one) THEN
           IF (incy.EQ.1) THEN
               IF (beta.EQ.zero) THEN
                   DO 10 i = 1,leny
                       y(i) = zero
    10             CONTINUE
               ELSE
                   DO 20 i = 1,leny
                       y(i) = beta*y(i)
    20             CONTINUE
               END IF
           ELSE
               iy = ky
               IF (beta.EQ.zero) THEN
                   DO 30 i = 1,leny
                       y(iy) = zero
                       iy = iy + incy
    30             CONTINUE
               ELSE
                   DO 40 i = 1,leny
                       y(iy) = beta*y(iy)
                       iy = iy + incy
    40             CONTINUE
               END IF
           END IF
       END IF
       IF (alpha.EQ.zero) RETURN
       IF (lsame(trans,'N')) THEN
 !
 !        Form  y := alpha*A*x + y.
 !
           jx = kx
           IF (incy.EQ.1) THEN
               DO 60 j = 1,n
                   temp = alpha*x(jx)
                   DO 50 i = 1,m
                       y(i) = y(i) + temp*a(i,j)
    50             CONTINUE
                   jx = jx + incx
    60         CONTINUE
           ELSE
               DO 80 j = 1,n
                   temp = alpha*x(jx)
                   iy = ky
                   DO 70 i = 1,m
                       y(iy) = y(iy) + temp*a(i,j)
                       iy = iy + incy
    70             CONTINUE
                   jx = jx + incx
    80         CONTINUE
           END IF
       ELSE
 !
 !        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
 !
           jy = ky
           IF (incx.EQ.1) THEN
               DO 110 j = 1,n
                   temp = zero
                   IF (noconj) THEN
                       DO 90 i = 1,m
                           temp = temp + a(i,j)*x(i)
    90                 CONTINUE
                   ELSE
                       DO 100 i = 1,m
                           temp = temp + dconjg(a(i,j))*x(i)
   100                 CONTINUE
                   END IF
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
   110         CONTINUE
           ELSE
               DO 140 j = 1,n
                   temp = zero
                   ix = kx
                   IF (noconj) THEN
                       DO 120 i = 1,m
                           temp = temp + a(i,j)*x(ix)
                           ix = ix + incx
   120                 CONTINUE
                   ELSE
                       DO 130 i = 1,m
                           temp = temp + dconjg(a(i,j))*x(ix)
                           ix = ix + incx
   130                 CONTINUE
                   END IF
                   y(jy) = y(jy) + alpha*temp
                   jy = jy + incy
   140         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZGEMV .
 !
  end subroutine zgemv

  !====================================================================
  subroutine zlacgv( N, X, INCX )
  !
  !  -- LAPACK auxiliary routine (version 3.4.2) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     September 2012
  !
  !     .. Scalar Arguments ..
        INTEGER            INCX, N
  !     ..
  !     .. Array Arguments ..
        COMPLEX*16         X( * )
  !     ..
  !
  ! =====================================================================
  !
  !     .. Local Scalars ..
        INTEGER            I, IOFF
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          dconjg
  !     ..
  !     .. Executable Statements ..
  !
        IF( incx.EQ.1 ) THEN
           DO 10 i = 1, n
              x( i ) = dconjg( x( i ) )
     10    CONTINUE
        ELSE
           ioff = 1
           IF( incx.LT.0 ) &
              ioff = 1 - ( n-1 )*incx
           DO 20 i = 1, n
              x( ioff ) = dconjg( x( ioff ) )
              ioff = ioff + incx
     20    CONTINUE
        END IF
        RETURN
  !
  !     End of ZLACGV
  !
  end subroutine zlacgv

 !  =====================================================================
  subroutine zherk(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
 !
 !  -- Reference BLAS level3 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION ALPHA,BETA
       INTEGER K,LDA,LDC,N
       CHARACTER TRANS,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),C(ldc,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
  !     LOGICAL LSAME
  !     EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
  !     EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dble,dcmplx,dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       DOUBLE PRECISION RTEMP
       INTEGER I,INFO,J,L,NROWA
       LOGICAL UPPER
 !     ..
 !     .. Parameters ..
       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d0,zero=0.0d0)
 !     ..
 !
 !     Test the input parameters.
 !
       IF (lsame(trans,'N')) THEN
           nrowa = n
       ELSE
           nrowa = k
       END IF
       upper = lsame(uplo,'U')
 !
       info = 0
       IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 1
       ELSE IF ((.NOT.lsame(trans,'N')) .AND. &
                (.NOT.lsame(trans,'C'))) THEN
           info = 2
       ELSE IF (n.LT.0) THEN
           info = 3
       ELSE IF (k.LT.0) THEN
           info = 4
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 7
       ELSE IF (ldc.LT.max(1,n)) THEN
           info = 10
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZHERK ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR. &
           (k.EQ.0)).AND. (beta.EQ.one))) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           IF (upper) THEN
               IF (beta.EQ.zero) THEN
                   DO 20 j = 1,n
                       DO 10 i = 1,j
                           c(i,j) = zero
    10                 CONTINUE
    20             CONTINUE
               ELSE
                   DO 40 j = 1,n
                       DO 30 i = 1,j - 1
                           c(i,j) = beta*c(i,j)
    30                 CONTINUE
                       c(j,j) = beta*dble(c(j,j))
    40             CONTINUE
               END IF
           ELSE
               IF (beta.EQ.zero) THEN
                   DO 60 j = 1,n
                       DO 50 i = j,n
                           c(i,j) = zero
    50                 CONTINUE
    60             CONTINUE
               ELSE
                   DO 80 j = 1,n
                       c(j,j) = beta*dble(c(j,j))
                       DO 70 i = j + 1,n
                           c(i,j) = beta*c(i,j)
    70                 CONTINUE
    80             CONTINUE
               END IF
           END IF
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (lsame(trans,'N')) THEN
 !
 !        Form  C := alpha*A*A**H + beta*C.
 !
           IF (upper) THEN
               DO 130 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 90 i = 1,j
                           c(i,j) = zero
    90                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 100 i = 1,j - 1
                           c(i,j) = beta*c(i,j)
   100                 CONTINUE
                       c(j,j) = beta*dble(c(j,j))
                   ELSE
                       c(j,j) = dble(c(j,j))
                   END IF
                   DO 120 l = 1,k
                       IF (a(j,l).NE.dcmplx(zero)) THEN
                           temp = alpha*dconjg(a(j,l))
                           DO 110 i = 1,j - 1
                               c(i,j) = c(i,j) + temp*a(i,l)
   110                     CONTINUE
                           c(j,j) = dble(c(j,j)) + dble(temp*a(i,l))
                       END IF
   120             CONTINUE
   130         CONTINUE
           ELSE
               DO 180 j = 1,n
                   IF (beta.EQ.zero) THEN
                       DO 140 i = j,n
                           c(i,j) = zero
   140                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       c(j,j) = beta*dble(c(j,j))
                       DO 150 i = j + 1,n
                           c(i,j) = beta*c(i,j)
   150                 CONTINUE
                   ELSE
                       c(j,j) = dble(c(j,j))
                   END IF
                   DO 170 l = 1,k
                       IF (a(j,l).NE.dcmplx(zero)) THEN
                           temp = alpha*dconjg(a(j,l))
                           c(j,j) = dble(c(j,j)) + dble(temp*a(j,l))
                           DO 160 i = j + 1,n
                               c(i,j) = c(i,j) + temp*a(i,l)
   160                     CONTINUE
                       END IF
   170             CONTINUE
   180         CONTINUE
           END IF
       ELSE
 !
 !        Form  C := alpha*A**H*A + beta*C.
 !
           IF (upper) THEN
               DO 220 j = 1,n
                   DO 200 i = 1,j - 1
                       temp = zero
                       DO 190 l = 1,k
                           temp = temp + dconjg(a(l,i))*a(l,j)
   190                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   200             CONTINUE
                   rtemp = zero
                   DO 210 l = 1,k
                       rtemp = rtemp + dconjg(a(l,j))*a(l,j)
   210             CONTINUE
                   IF (beta.EQ.zero) THEN
                       c(j,j) = alpha*rtemp
                   ELSE
                       c(j,j) = alpha*rtemp + beta*dble(c(j,j))
                   END IF
   220         CONTINUE
           ELSE
               DO 260 j = 1,n
                   rtemp = zero
                   DO 230 l = 1,k
                       rtemp = rtemp + dconjg(a(l,j))*a(l,j)
   230             CONTINUE
                   IF (beta.EQ.zero) THEN
                       c(j,j) = alpha*rtemp
                   ELSE
                       c(j,j) = alpha*rtemp + beta*dble(c(j,j))
                   END IF
                   DO 250 i = j + 1,n
                       temp = zero
                       DO 240 l = 1,k
                           temp = temp + dconjg(a(l,i))*a(l,j)
   240                 CONTINUE
                       IF (beta.EQ.zero) THEN
                           c(i,j) = alpha*temp
                       ELSE
                           c(i,j) = alpha*temp + beta*c(i,j)
                       END IF
   250             CONTINUE
   260         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZHERK .
 !
  end subroutine zherk

 !  =====================================================================
  subroutine ztrsm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
 !
 !  -- Reference BLAS level3 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       INTEGER LDA,LDB,M,N
       CHARACTER DIAG,SIDE,TRANSA,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),B(ldb,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
  !     LOGICAL LSAME
  !     EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
  !     EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,J,K,NROWA
       LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
 !     ..
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d0,0.0d0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !
 !     Test the input parameters.
 !
       lside = lsame(side,'L')
       IF (lside) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF
       noconj = lsame(transa,'T')
       nounit = lsame(diag,'N')
       upper = lsame(uplo,'U')
 !
       info = 0
       IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
           info = 1
       ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 2
       ELSE IF ((.NOT.lsame(transa,'N')) .AND. &
                (.NOT.lsame(transa,'T')) .AND. &
                (.NOT.lsame(transa,'C'))) THEN
           info = 3
       ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
           info = 4
       ELSE IF (m.LT.0) THEN
           info = 5
       ELSE IF (n.LT.0) THEN
           info = 6
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 9
       ELSE IF (ldb.LT.max(1,m)) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZTRSM ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (m.EQ.0 .OR. n.EQ.0) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           DO 20 j = 1,n
               DO 10 i = 1,m
                   b(i,j) = zero
    10         CONTINUE
    20     CONTINUE
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (lside) THEN
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*inv( A )*B.
 !
               IF (upper) THEN
                   DO 60 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 30 i = 1,m
                               b(i,j) = alpha*b(i,j)
    30                     CONTINUE
                       END IF
                       DO 50 k = m,1,-1
                           IF (b(k,j).NE.zero) THEN
                               IF (nounit) b(k,j) = b(k,j)/a(k,k)
                               DO 40 i = 1,k - 1
                                   b(i,j) = b(i,j) - b(k,j)*a(i,k)
    40                         CONTINUE
                           END IF
    50                 CONTINUE
    60             CONTINUE
               ELSE
                   DO 100 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 70 i = 1,m
                               b(i,j) = alpha*b(i,j)
    70                     CONTINUE
                       END IF
                       DO 90 k = 1,m
                           IF (b(k,j).NE.zero) THEN
                               IF (nounit) b(k,j) = b(k,j)/a(k,k)
                               DO 80 i = k + 1,m
                                   b(i,j) = b(i,j) - b(k,j)*a(i,k)
    80                         CONTINUE
                           END IF
    90                 CONTINUE
   100             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*inv( A**T )*B
 !           or    B := alpha*inv( A**H )*B.
 !
               IF (upper) THEN
                   DO 140 j = 1,n
                       DO 130 i = 1,m
                           temp = alpha*b(i,j)
                           IF (noconj) THEN
                               DO 110 k = 1,i - 1
                                   temp = temp - a(k,i)*b(k,j)
   110                         CONTINUE
                               IF (nounit) temp = temp/a(i,i)
                           ELSE
                               DO 120 k = 1,i - 1
                                   temp = temp - dconjg(a(k,i))*b(k,j)
   120                         CONTINUE
                               IF (nounit) temp = temp/dconjg(a(i,i))
                           END IF
                           b(i,j) = temp
   130                 CONTINUE
   140             CONTINUE
               ELSE
                   DO 180 j = 1,n
                       DO 170 i = m,1,-1
                           temp = alpha*b(i,j)
                           IF (noconj) THEN
                               DO 150 k = i + 1,m
                                   temp = temp - a(k,i)*b(k,j)
   150                         CONTINUE
                               IF (nounit) temp = temp/a(i,i)
                           ELSE
                               DO 160 k = i + 1,m
                                   temp = temp - dconjg(a(k,i))*b(k,j)
   160                         CONTINUE
                               IF (nounit) temp = temp/dconjg(a(i,i))
                           END IF
                           b(i,j) = temp
   170                 CONTINUE
   180             CONTINUE
               END IF
           END IF
       ELSE
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*B*inv( A ).
 !
               IF (upper) THEN
                   DO 230 j = 1,n
                       IF (alpha.NE.one) THEN
                           DO 190 i = 1,m
                               b(i,j) = alpha*b(i,j)
   190                     CONTINUE
                       END IF
                       DO 210 k = 1,j - 1
                           IF (a(k,j).NE.zero) THEN
                               DO 200 i = 1,m
                                   b(i,j) = b(i,j) - a(k,j)*b(i,k)
   200                         CONTINUE
                           END IF
   210                 CONTINUE
                       IF (nounit) THEN
                           temp = one/a(j,j)
                           DO 220 i = 1,m
                               b(i,j) = temp*b(i,j)
   220                     CONTINUE
                       END IF
   230             CONTINUE
               ELSE
                   DO 280 j = n,1,-1
                       IF (alpha.NE.one) THEN
                           DO 240 i = 1,m
                               b(i,j) = alpha*b(i,j)
   240                     CONTINUE
                       END IF
                       DO 260 k = j + 1,n
                           IF (a(k,j).NE.zero) THEN
                               DO 250 i = 1,m
                                   b(i,j) = b(i,j) - a(k,j)*b(i,k)
   250                         CONTINUE
                           END IF
   260                 CONTINUE
                       IF (nounit) THEN
                           temp = one/a(j,j)
                           DO 270 i = 1,m
                               b(i,j) = temp*b(i,j)
   270                     CONTINUE
                       END IF
   280             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*B*inv( A**T )
 !           or    B := alpha*B*inv( A**H ).
 !
               IF (upper) THEN
                   DO 330 k = n,1,-1
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = one/a(k,k)
                           ELSE
                               temp = one/dconjg(a(k,k))
                           END IF
                           DO 290 i = 1,m
                               b(i,k) = temp*b(i,k)
   290                     CONTINUE
                       END IF
                       DO 310 j = 1,k - 1
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = a(j,k)
                               ELSE
                                   temp = dconjg(a(j,k))
                               END IF
                               DO 300 i = 1,m
                                   b(i,j) = b(i,j) - temp*b(i,k)
   300                         CONTINUE
                           END IF
   310                 CONTINUE
                       IF (alpha.NE.one) THEN
                           DO 320 i = 1,m
                               b(i,k) = alpha*b(i,k)
   320                     CONTINUE
                       END IF
   330             CONTINUE
               ELSE
                   DO 380 k = 1,n
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = one/a(k,k)
                           ELSE
                               temp = one/dconjg(a(k,k))
                           END IF
                           DO 340 i = 1,m
                               b(i,k) = temp*b(i,k)
   340                     CONTINUE
                       END IF
                       DO 360 j = k + 1,n
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = a(j,k)
                               ELSE
                                   temp = dconjg(a(j,k))
                               END IF
                               DO 350 i = 1,m
                                   b(i,j) = b(i,j) - temp*b(i,k)
   350                         CONTINUE
                           END IF
   360                 CONTINUE
                       IF (alpha.NE.one) THEN
                           DO 370 i = 1,m
                               b(i,k) = alpha*b(i,k)
   370                     CONTINUE
                       END IF
   380             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRSM .
 !
  end subroutine ztrsm

 !  =====================================================================
  subroutine zheev( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, info )
 !
 ! ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
 ! complex Hermitian matrix A.
 !
 ![in]    JOBZ
 !          JOBZ is CHARACTER*1
 !          = 'N':  Compute eigenvalues only;
 !          = 'V':  Compute eigenvalues and eigenvectors.
 ![in]    UPLO
 !          UPLO is CHARACTER*1
 !          = 'U':  Upper triangle of A is stored;
 !          = 'L':  Lower triangle of A is stored.
 ![in]    N
 !          N is INTEGER
 !          The order of the matrix A.  N >= 0.
 ! [in,out]    A
 !         A is COMPLEX*16 array, dimension (LDA, N)
 !         On entry, the Hermitian matrix A.  If UPLO = 'U', the
 !         leading N-by-N upper triangular part of A contains the
 !         upper triangular part of the matrix A.  If UPLO = 'L',
 !         the leading N-by-N lower triangular part of A contains
 !         the lower triangular part of the matrix A.
 !         On exit, if JOBZ = 'V', then if INFO = 0, A contains the
 !         orthonormal eigenvectors of the matrix A.
 !         If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
 !         or the upper triangle (if UPLO='U') of A, including the
 !         diagonal, is destroyed.
 ! [in]    LDA
 !         LDA is INTEGER
 !         The leading dimension of the array A.  LDA >= max(1,N).
 ! [out]   W
 !         W is DOUBLE PRECISION array, dimension (N)
 !         If INFO = 0, the eigenvalues in ascending order.
 ! [out]   WORK
 !         WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
 !         On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 ! [in]    LWORK
 !         LWORK is INTEGER
 !         The length of the array WORK.  LWORK >= max(1,2*N-1).
 !         For optimal efficiency, LWORK >= (NB+1)*N,
 !         where NB is the blocksize for ZHETRD returned by ILAENV.
 !
 !         If LWORK = -1, then a workspace query is assumed; the routine
 !         only calculates the optimal size of the WORK array, returns
 !         this value as the first entry of the WORK array, and no error
 !         message related to LWORK is issued by XERBLA.
 ! [out]   RWORK
 !         RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))
 ! [out]   INFO
 !        INFO is INTEGER
 !         = 0:  successful exit
 !         < 0:  if INFO = -i, the i-th argument had an illegal value
 !         > 0:  if INFO = i, the algorithm failed to converge; i
 !               off-diagonal elements of an intermediate tridiagonal
 !                form did not converge to zero.
 !
 !
 !  -- LAPACK driver routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER          JOBZ, UPLO
       INTEGER            INFO, LDA, LWORK, N
 !     ..
 !     .. Array Arguments ..

       DOUBLE PRECISION   RWORK( * ), W( * )
       COMPLEX*16         A( lda, * ), WORK( *)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO, ONE
       parameter( zero = 0.0d0, one = 1.0d0 )
       COMPLEX*16         CONE
       parameter( cone = ( 1.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            LOWER, LQUERY, WANTZ
       INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE, &
                          llwork, lwkopt, nb
       DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, &
                          smlnum
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      INTEGER            ILAENV
 !      DOUBLE PRECISION   DLAMCH, ZLANHE
 !      EXTERNAL           lsame, ilaenv, dlamch, zlanhe
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           dscal, dsterf, xerbla, zhetrd, zlascl, zsteqr, &
 !                         zungtr
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, sqrt
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       wantz = lsame( jobz, 'V' )
       lower = lsame( uplo, 'L' )
       lquery = ( lwork.EQ.-1 )
 !
       info = 0
       IF( .NOT.( wantz .OR. lsame( jobz, 'N' ) ) ) THEN
          info = -1
       ELSE IF( .NOT.( lower .OR. lsame( uplo, 'U' ) ) ) THEN
          info = -2
       ELSE IF( n.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       END IF
 !
       IF( info.EQ.0 ) THEN
          nb = ilaenv( 1, 'ZHETRD', uplo, n, -1, -1, -1 )
          lwkopt = max( 1, ( nb+1 )*n )
          work( 1 ) = lwkopt
 !
          IF( lwork.LT.max( 1, 2*n-1 ) .AND. .NOT.lquery ) &
             info = -8
       END IF
 !
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZHEEV ', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 ) THEN
          RETURN
       END IF
 !
       IF( n.EQ.1 ) THEN
          w( 1 ) = a( 1, 1 )
          work( 1 ) = 1
          IF( wantz ) &
             a( 1, 1 ) = cone
          RETURN
       END IF
 !
 !     Get machine constants.
 !
       safmin = dlamch( 'Safe minimum' )
       eps = dlamch( 'Precision' )
       smlnum = safmin / eps
       bignum = one / smlnum
       rmin = sqrt( smlnum )
       rmax = sqrt( bignum )
 !
 !     Scale matrix to allowable range, if necessary.
 !
       anrm = zlanhe( 'M', uplo, n, a, lda, rwork )
       iscale = 0
       IF( anrm.GT.zero .AND. anrm.LT.rmin ) THEN
          iscale = 1
          sigma = rmin / anrm
       ELSE IF( anrm.GT.rmax ) THEN
          iscale = 1
          sigma = rmax / anrm
       END IF
       IF( iscale.EQ.1 ) &
          CALL zlascl( uplo, 0, 0, one, sigma, n, n, a, lda, info )
 !
 !     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
 !
       inde = 1
       indtau = 1
       indwrk = indtau + n
       llwork = lwork - indwrk + 1
       CALL zhetrd( uplo, n, a, lda, w, rwork( inde ), work( indtau ), &
                    work( indwrk ), llwork, iinfo )
 !
 !     For eigenvalues only, call DSTERF.  For eigenvectors, first call
 !     ZUNGTR to generate the unitary matrix, then call ZSTEQR.
 !
       IF( .NOT.wantz ) THEN
          CALL dsterf( n, w, rwork( inde ), info )
       ELSE
          CALL zungtr( uplo, n, a, lda, work( indtau ), work( indwrk ), &
                       llwork, iinfo )
          indwrk = inde + n
          CALL zsteqr( jobz, n, w, rwork( inde ), a, lda, &
                       rwork( indwrk ), info )
       END IF

 !
 !     If matrix was scaled, then rescale eigenvalues appropriately.
 !
       IF( iscale.EQ.1 ) THEN
          IF( info.EQ.0 ) THEN
             imax = n
          ELSE
             imax = info - 1
          END IF
          CALL dscal( imax, one / sigma, w, 1 )
       END IF
 !
 !     Set WORK(1) to optimal complex workspace size.
 !
       work( 1 ) = lwkopt
 !
       RETURN
 !
 !     End of ZHEEV
 !
  end subroutine zheev

 !  =====================================================================
  double precision function dlamch( CMACH )
 !
 !  -- LAPACK auxiliary routine (version 3.6.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2015
 !
 !     .. Scalar Arguments ..
       CHARACTER          CMACH
 !     ..
 !
 ! =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE, ZERO
       parameter( one = 1.0d0, zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION   RND, EPS, SFMIN, SMALL, RMACH
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      EXTERNAL           lsame
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          digits, epsilon, huge, maxexponent, &
                          minexponent, radix, tiny
 !     ..
 !     .. Executable Statements ..
 !
 !
 !     Assume rounding, not chopping. Always.
 !
       rnd = one
 !
       IF( one.EQ.rnd ) THEN
          eps = epsilon(zero) * 0.5
       ELSE
          eps = epsilon(zero)
       END IF
 !
       IF( lsame( cmach, 'E' ) ) THEN
          rmach = eps
       ELSE IF( lsame( cmach, 'S' ) ) THEN
          sfmin = tiny(zero)
          small = one / huge(zero)
          IF( small.GE.sfmin ) THEN
 !
 !           Use SMALL plus a bit, to avoid the possibility of rounding
 !           causing overflow when computing  1/sfmin.
 !
             sfmin = small*( one+eps )
          END IF
          rmach = sfmin
       ELSE IF( lsame( cmach, 'B' ) ) THEN
          rmach = radix(zero)
       ELSE IF( lsame( cmach, 'P' ) ) THEN
          rmach = eps * radix(zero)
       ELSE IF( lsame( cmach, 'N' ) ) THEN
          rmach = digits(zero)
       ELSE IF( lsame( cmach, 'R' ) ) THEN
          rmach = rnd
       ELSE IF( lsame( cmach, 'M' ) ) THEN
          rmach = minexponent(zero)
       ELSE IF( lsame( cmach, 'U' ) ) THEN
          rmach = tiny(zero)
       ELSE IF( lsame( cmach, 'L' ) ) THEN
          rmach = maxexponent(zero)
       ELSE IF( lsame( cmach, 'O' ) ) THEN
          rmach = huge(zero)
       ELSE
          rmach = zero
       END IF
 !
       dlamch = rmach
       RETURN
 !
 !     End of DLAMCH
 !
  end function dlamch

 !  =====================================================================
  double precision function zlanhe( NORM, UPLO, N, A, LDA, WORK )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          NORM, UPLO
       INTEGER            LDA, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   WORK( * )
       COMPLEX*16         A( lda, * )
 !     ..
 !
 ! =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE, ZERO
       parameter( one = 1.0d0, zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, J
       DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME, DISNAN
 !      EXTERNAL           lsame, disnan
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           zlassq
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, dble, sqrt
 !     ..
 !     .. Executable Statements ..
 !
       IF( n.EQ.0 ) THEN
          VALUE = zero
       ELSE IF( lsame( norm, 'M' ) ) THEN
 !
 !        Find max(abs(A(i,j))).
 !
          VALUE = zero
          IF( lsame( uplo, 'U' ) ) THEN
             DO 20 j = 1, n
                DO 10 i = 1, j - 1
                   sum = abs( a( i, j ) )
                   IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
    10          CONTINUE
                sum = abs( dble( a( j, j ) ) )
                IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
    20       CONTINUE
          ELSE
             DO 40 j = 1, n
                sum = abs( dble( a( j, j ) ) )
                IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
                DO 30 i = j + 1, n
                   sum = abs( a( i, j ) )
                   IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
    30          CONTINUE
    40       CONTINUE
          END IF
       ELSE IF( ( lsame( norm, 'I' ) ) .OR. ( lsame( norm, 'O' ) ) .OR. &
               ( norm.EQ.'1' ) ) THEN
 !
 !        Find normI(A) ( = norm1(A), since A is hermitian).
 !
          VALUE = zero
          IF( lsame( uplo, 'U' ) ) THEN
             DO 60 j = 1, n
                sum = zero
                DO 50 i = 1, j - 1
                   absa = abs( a( i, j ) )
                   sum = sum + absa
                   work( i ) = work( i ) + absa
    50          CONTINUE
                work( j ) = sum + abs( dble( a( j, j ) ) )
    60       CONTINUE
             DO 70 i = 1, n
                sum = work( i )
                IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
    70       CONTINUE
          ELSE
             DO 80 i = 1, n
                work( i ) = zero
    80       CONTINUE
             DO 100 j = 1, n
                sum = work( j ) + abs( dble( a( j, j ) ) )
                DO 90 i = j + 1, n
                   absa = abs( a( i, j ) )
                   sum = sum + absa
                   work( i ) = work( i ) + absa
    90          CONTINUE
                IF( VALUE .LT. sum .OR. disnan( sum ) ) VALUE = sum
   100       CONTINUE
          END IF
       ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
 !
 !        Find normF(A).
 !
          scale = zero
          sum = one
          IF( lsame( uplo, 'U' ) ) THEN
             DO 110 j = 2, n
                CALL zlassq( j-1, a( 1, j ), 1, scale, sum )
   110       CONTINUE
          ELSE
             DO 120 j = 1, n - 1
                CALL zlassq( n-j, a( j+1, j ), 1, scale, sum )
   120       CONTINUE
          END IF
          sum = 2*sum
          DO 130 i = 1, n
             IF( dble( a( i, i ) ).NE.zero ) THEN
                absa = abs( dble( a( i, i ) ) )
                IF( scale.LT.absa ) THEN
                   sum = one + sum*( scale / absa )**2
                   scale = absa
                ELSE
                   sum = sum + ( absa / scale )**2
                END IF
             END IF
   130    CONTINUE
          VALUE = scale*sqrt( sum )
       END IF
 !
       zlanhe = VALUE
       RETURN
 !
 !     End of ZLANHE
 !
  end function zlanhe

 !  =====================================================================
  subroutine zlassq( N, X, INCX, SCALE, SUMSQ )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       INTEGER            INCX, N
       DOUBLE PRECISION   SCALE, SUMSQ
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         X( * )
 !     ..
 !
 ! =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            IX
       DOUBLE PRECISION   TEMP1
 !     ..
 !     .. External Functions ..
 !      LOGICAL            DISNAN
 !      EXTERNAL           disnan
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, dble, dimag
 !     ..
 !     .. Executable Statements ..
 !
       IF( n.GT.0 ) THEN
          DO 10 ix = 1, 1 + ( n-1 )*incx, incx
             temp1 = abs( dble( x( ix ) ) )
             IF( temp1.GT.zero.OR.disnan( temp1 ) ) THEN
                IF( scale.LT.temp1 ) THEN
                   sumsq = 1 + sumsq*( scale / temp1 )**2
                   scale = temp1
                ELSE
                   sumsq = sumsq + ( temp1 / scale )**2
                END IF
             END IF
             temp1 = abs( dimag( x( ix ) ) )
             IF( temp1.GT.zero.OR.disnan( temp1 ) ) THEN
                IF( scale.LT.temp1 ) THEN
                   sumsq = 1 + sumsq*( scale / temp1 )**2
                   scale = temp1
                ELSE
                   sumsq = sumsq + ( temp1 / scale )**2
                END IF
             END IF
    10    CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZLASSQ
 !
  end subroutine zlassq

 !  =====================================================================
  subroutine dscal(N,DA,DX,INCX)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION DA
       INTEGER INCX,N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION DX(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,M,MP1,NINCX
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC mod
 !     ..
       IF (n.LE.0 .OR. incx.LE.0) RETURN
       IF (incx.EQ.1) THEN
 !
 !        code for increment equal to 1
 !
 !
 !        clean-up loop
 !
          m = mod(n,5)
          IF (m.NE.0) THEN
             DO i = 1,m
                dx(i) = da*dx(i)
             END DO
             IF (n.LT.5) RETURN
          END IF
          mp1 = m + 1
          DO i = mp1,n,5
             dx(i) = da*dx(i)
             dx(i+1) = da*dx(i+1)
             dx(i+2) = da*dx(i+2)
             dx(i+3) = da*dx(i+3)
             dx(i+4) = da*dx(i+4)
          END DO
       ELSE
 !
 !        code for increment not equal to 1
 !
          nincx = n*incx
          DO i = 1,nincx,incx
             dx(i) = da*dx(i)
          END DO
       END IF
       RETURN
  end subroutine dscal

 !  =====================================================================
  subroutine dsterf( N, D, E, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   D( * ), E( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO, ONE, TWO, THREE
       parameter( zero = 0.0d0, one = 1.0d0, two = 2.0d0, &
                          three = 3.0d0 )
       INTEGER            MAXIT
       parameter( maxit = 30 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, &
                          nmaxit
       DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC, &
                          oldgam, p, r, rt1, rt2, rte, s, safmax, safmin, &
                          sigma, ssfmax, ssfmin, rmax
 !     ..
 !     .. External Functions ..
 !      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
 !      EXTERNAL           dlamch, dlanst, dlapy2
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           dlae2, dlascl, dlasrt, xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, sign, sqrt
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       info = 0
 !
 !     Quick return if possible
 !
       IF( n.LT.0 ) THEN
          info = -1
          CALL xerbla( 'DSTERF', -info )
          RETURN
       END IF
       IF( n.LE.1 ) &
          RETURN
 !
 !     Determine the unit roundoff for this environment.
 !
       eps = dlamch( 'E' )
       eps2 = eps**2
       safmin = dlamch( 'S' )
       safmax = one / safmin
       ssfmax = sqrt( safmax ) / three
       ssfmin = sqrt( safmin ) / eps2
       rmax = dlamch( 'O' )
 !
 !     Compute the eigenvalues of the tridiagonal matrix.
 !
       nmaxit = n*maxit
       sigma = zero
       jtot = 0
 !
 !     Determine where the matrix splits and choose QL or QR iteration
 !     for each block, according to whether top or bottom diagonal
 !     element is smaller.
 !
       l1 = 1
 !
    10 CONTINUE
       IF( l1.GT.n ) &
          GO TO 170
       IF( l1.GT.1 ) &
          e( l1-1 ) = zero
       DO 20 m = l1, n - 1
          IF( abs( e( m ) ).LE.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+ &
              1 ) ) ) )*eps ) THEN
             e( m ) = zero
             GO TO 30
          END IF
    20 CONTINUE
       m = n
 !
    30 CONTINUE
       l = l1
       lsv = l
       lend = m
       lendsv = lend
       l1 = m + 1
       IF( lend.EQ.l ) &
          GO TO 10
 !
 !     Scale submatrix in rows and columns L to LEND
 !
       anorm = dlanst( 'M', lend-l+1, d( l ), e( l ) )
       iscale = 0
       IF( anorm.EQ.zero ) &
          GO TO 10
       IF( (anorm.GT.ssfmax) ) THEN
          iscale = 1
          CALL dlascl( 'G', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n, &
                       info )
          CALL dlascl( 'G', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n, &
                       info )
       ELSE IF( anorm.LT.ssfmin ) THEN
          iscale = 2
          CALL dlascl( 'G', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n, &
                       info )
          CALL dlascl( 'G', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n, &
                       info )
       END IF
 !
       DO 40 i = l, lend - 1
          e( i ) = e( i )**2
    40 CONTINUE
 !
 !     Choose between QL and QR iteration
 !
       IF( abs( d( lend ) ).LT.abs( d( l ) ) ) THEN
          lend = lsv
          l = lendsv
       END IF
 !
       IF( lend.GE.l ) THEN
 !
 !        QL Iteration
 !
 !        Look for small subdiagonal element.
 !
    50    CONTINUE
          IF( l.NE.lend ) THEN
             DO 60 m = l, lend - 1
                IF( abs( e( m ) ).LE.eps2*abs( d( m )*d( m+1 ) ) ) &
                   GO TO 70
    60       CONTINUE
          END IF
          m = lend
 !
    70    CONTINUE
          IF( m.LT.lend ) &
             e( m ) = zero
          p = d( l )
          IF( m.EQ.l ) &
             GO TO 90
 !
 !        If remaining matrix is 2 by 2, use DLAE2 to compute its
 !        eigenvalues.
 !
          IF( m.EQ.l+1 ) THEN
             rte = sqrt( e( l ) )
             CALL dlae2( d( l ), rte, d( l+1 ), rt1, rt2 )
             d( l ) = rt1
             d( l+1 ) = rt2
             e( l ) = zero
             l = l + 2
             IF( l.LE.lend ) &
                GO TO 50
             GO TO 150
          END IF
 !
          IF( jtot.EQ.nmaxit ) &
             GO TO 150
          jtot = jtot + 1
 !
 !        Form shift.
 !
          rte = sqrt( e( l ) )
          sigma = ( d( l+1 )-p ) / ( two*rte )
          r = dlapy2( sigma, one )
          sigma = p - ( rte / ( sigma+sign( r, sigma ) ) )
 !
          c = one
          s = zero
          gamma = d( m ) - sigma
          p = gamma*gamma
 !
 !        Inner loop
 !
          DO 80 i = m - 1, l, -1
             bb = e( i )
             r = p + bb
             IF( i.NE.m-1 ) &
                e( i+1 ) = s*r
             oldc = c
             c = p / r
             s = bb / r
             oldgam = gamma
             alpha = d( i )
             gamma = c*( alpha-sigma ) - s*oldgam
             d( i+1 ) = oldgam + ( alpha-gamma )
             IF( c.NE.zero ) THEN
                p = ( gamma*gamma ) / c
             ELSE
                p = oldc*bb
             END IF
    80    CONTINUE
 !
          e( l ) = s*p
          d( l ) = sigma + gamma
          GO TO 50
 !
 !        Eigenvalue found.
 !
    90    CONTINUE
          d( l ) = p
 !
          l = l + 1
          IF( l.LE.lend ) &
             GO TO 50
          GO TO 150
 !
       ELSE
 !
 !        QR Iteration
 !
 !        Look for small superdiagonal element.
 !
   100    CONTINUE
          DO 110 m = l, lend + 1, -1
             IF( abs( e( m-1 ) ).LE.eps2*abs( d( m )*d( m-1 ) ) ) &
                GO TO 120
   110    CONTINUE
          m = lend
 !
   120    CONTINUE
          IF( m.GT.lend ) &
             e( m-1 ) = zero
          p = d( l )
          IF( m.EQ.l ) &
             GO TO 140
 !
 !        If remaining matrix is 2 by 2, use DLAE2 to compute its
 !        eigenvalues.
 !
          IF( m.EQ.l-1 ) THEN
             rte = sqrt( e( l-1 ) )
             CALL dlae2( d( l ), rte, d( l-1 ), rt1, rt2 )
             d( l ) = rt1
             d( l-1 ) = rt2
             e( l-1 ) = zero
             l = l - 2
             IF( l.GE.lend ) &
                GO TO 100
             GO TO 150
          END IF
 !
          IF( jtot.EQ.nmaxit ) &
             GO TO 150
          jtot = jtot + 1
 !
 !        Form shift.
 !
          rte = sqrt( e( l-1 ) )
          sigma = ( d( l-1 )-p ) / ( two*rte )
          r = dlapy2( sigma, one )
          sigma = p - ( rte / ( sigma+sign( r, sigma ) ) )
 !
          c = one
          s = zero
          gamma = d( m ) - sigma
          p = gamma*gamma
 !
 !        Inner loop
 !
          DO 130 i = m, l - 1
             bb = e( i )
             r = p + bb
             IF( i.NE.m ) &
                e( i-1 ) = s*r
             oldc = c
             c = p / r
             s = bb / r
             oldgam = gamma
             alpha = d( i+1 )
             gamma = c*( alpha-sigma ) - s*oldgam
             d( i ) = oldgam + ( alpha-gamma )
             IF( c.NE.zero ) THEN
                p = ( gamma*gamma ) / c
             ELSE
                p = oldc*bb
             END IF
   130    CONTINUE
 !
          e( l-1 ) = s*p
          d( l ) = sigma + gamma
          GO TO 100
 !
 !        Eigenvalue found.
 !
   140    CONTINUE
          d( l ) = p
 !
          l = l - 1
          IF( l.GE.lend ) &
             GO TO 100
          GO TO 150
 !
       END IF
 !
 !     Undo scaling if necessary
 !
   150 CONTINUE
       IF( iscale.EQ.1 ) &
          CALL dlascl( 'G', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1, &
                       d( lsv ), n, info )
       IF( iscale.EQ.2 ) &
          CALL dlascl( 'G', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1, &
                       d( lsv ), n, info )
 !
 !     Check for no convergence to an eigenvalue after a total
 !     of N*MAXIT iterations.
 !
       IF( jtot.LT.nmaxit ) &
          GO TO 10
       DO 160 i = 1, n - 1
          IF( e( i ).NE.zero ) &
             info = info + 1
   160 CONTINUE
       GO TO 180
 !
 !     Sort eigenvalues in increasing order.
 !
   170 CONTINUE
       CALL dlasrt( 'I', n, d, info )
 !
   180 CONTINUE
       RETURN
 !
 !     End of DSTERF
 !
  end subroutine dsterf

  !  =====================================================================
  double precision function dlanst( NORM, N, D, E )
  !
  !  -- LAPACK auxiliary routine (version 3.4.2) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     September 2012
  !
  !     .. Scalar Arguments ..
        CHARACTER          NORM
        INTEGER            N
  !     ..
  !     .. Array Arguments ..
        DOUBLE PRECISION   D( * ), E( * )
  !     ..
  !
  !  =====================================================================
  !
  !     .. Parameters ..
        DOUBLE PRECISION   ONE, ZERO
        parameter( one = 1.0d0, zero = 0.0d0 )
  !     ..
  !     .. Local Scalars ..
        INTEGER            I
        DOUBLE PRECISION   ANORM, SCALE, SUM
  !     ..
  !     .. External Functions ..
  !      LOGICAL            LSAME, DISNAN
  !      EXTERNAL           lsame, disnan
  !     ..
  !     .. External Subroutines ..
  !      EXTERNAL           dlassq
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC          abs, sqrt
  !     ..
  !     .. Executable Statements ..
  !
        IF( n.LE.0 ) THEN
           anorm = zero
        ELSE IF( lsame( norm, 'M' ) ) THEN
  !
  !        Find max(abs(A(i,j))).
  !
           anorm = abs( d( n ) )
           DO 10 i = 1, n - 1
              sum = abs( d( i ) )
              IF( anorm .LT. sum .OR. disnan( sum ) ) anorm = sum
              sum = abs( e( i ) )
              IF( anorm .LT. sum .OR. disnan( sum ) ) anorm = sum
     10    CONTINUE
        ELSE IF( lsame( norm, 'O' ) .OR. norm.EQ.'1' .OR. &
                 lsame( norm, 'I' ) ) THEN
  !
  !        Find norm1(A).
  !
           IF( n.EQ.1 ) THEN
              anorm = abs( d( 1 ) )
           ELSE
              anorm = abs( d( 1 ) )+abs( e( 1 ) )
              sum = abs( e( n-1 ) )+abs( d( n ) )
              IF( anorm .LT. sum .OR. disnan( sum ) ) anorm = sum
              DO 20 i = 2, n - 1
                 sum = abs( d( i ) )+abs( e( i ) )+abs( e( i-1 ) )
                 IF( anorm .LT. sum .OR. disnan( sum ) ) anorm = sum
     20       CONTINUE
           END IF
        ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
  !
  !        Find normF(A).
  !
           scale = zero
           sum = one
           IF( n.GT.1 ) THEN
              CALL dlassq( n-1, e, 1, scale, sum )
              sum = 2*sum
           END IF
           CALL dlassq( n, d, 1, scale, sum )
           anorm = scale*sqrt( sum )
        END IF
  !
        dlanst = anorm
        RETURN
  !
  !     End of DLANST
  !
  end function dlanst

 !  =====================================================================
  subroutine dlassq( N, X, INCX, SCALE, SUMSQ )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       INTEGER            INCX, N
       DOUBLE PRECISION   SCALE, SUMSQ
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   X( * )
 !     ..
 !
 ! =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            IX
       DOUBLE PRECISION   ABSXI
 !     ..
 !     .. External Functions ..
 !      LOGICAL            DISNAN
 !      EXTERNAL           disnan
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs
 !     ..
 !     .. Executable Statements ..
 !
       IF( n.GT.0 ) THEN
          DO 10 ix = 1, 1 + ( n-1 )*incx, incx
             absxi = abs( x( ix ) )
             IF( absxi.GT.zero.OR.disnan( absxi ) ) THEN
                IF( scale.LT.absxi ) THEN
                   sumsq = 1 + sumsq*( scale / absxi )**2
                   scale = absxi
                ELSE
                   sumsq = sumsq + ( absxi / scale )**2
                END IF
             END IF
    10    CONTINUE
       END IF
       RETURN
 !
 !     End of DLASSQ
 !
  end subroutine dlassq

 !  =====================================================================
  double precision function dlapy2( X, Y )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   X, Y
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d0 )
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION   W, XABS, YABS, Z
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, max, min, sqrt
 !     ..
 !     .. Executable Statements ..
 !
       xabs = abs( x )
       yabs = abs( y )
       w = max( xabs, yabs )
       z = min( xabs, yabs )
       IF( z.EQ.zero ) THEN
          dlapy2 = w
       ELSE
          dlapy2 = w*sqrt( one+( z / w )**2 )
       END IF
       RETURN
 !
 !     End of DLAPY2
 !
  end function dlapy2

 !  =====================================================================
  subroutine dlae2( A, B, C, RT1, RT2 )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   A, B, C, RT1, RT2
 !     ..
 !
 ! =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d0 )
       DOUBLE PRECISION   TWO
       parameter( two = 2.0d0 )
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
       DOUBLE PRECISION   HALF
       parameter( half = 0.5d0 )
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION   AB, ACMN, ACMX, ADF, DF, RT, SM, TB
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, sqrt
 !     ..
 !     .. Executable Statements ..
 !
 !     Compute the eigenvalues
 !
       sm = a + c
       df = a - c
       adf = abs( df )
       tb = b + b
       ab = abs( tb )
       IF( abs( a ).GT.abs( c ) ) THEN
          acmx = a
          acmn = c
       ELSE
          acmx = c
          acmn = a
       END IF
       IF( adf.GT.ab ) THEN
          rt = adf*sqrt( one+( ab / adf )**2 )
       ELSE IF( adf.LT.ab ) THEN
          rt = ab*sqrt( one+( adf / ab )**2 )
       ELSE
 !
 !        Includes case AB=ADF=0
 !
          rt = ab*sqrt( two )
       END IF
       IF( sm.LT.zero ) THEN
          rt1 = half*( sm-rt )
 !
 !        Order of execution important.
 !        To get fully accurate smaller eigenvalue,
 !        next line needs to be executed in higher precision.
 !
          rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
       ELSE IF( sm.GT.zero ) THEN
          rt1 = half*( sm+rt )
 !
 !        Order of execution important.
 !        To get fully accurate smaller eigenvalue,
 !        next line needs to be executed in higher precision.
 !
          rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
       ELSE
 !
 !        Includes case RT1 = RT2 = 0
 !
          rt1 = half*rt
          rt2 = -half*rt
       END IF
       RETURN
 !
 !     End of DLAE2
 !
  end subroutine dlae2

 !  =====================================================================
  subroutine dlascl( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          TYPE
       INTEGER            INFO, KL, KU, LDA, M, N
       DOUBLE PRECISION   CFROM, CTO
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO, ONE
       parameter( zero = 0.0d0, one = 1.0d0 )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            DONE
       INTEGER            I, ITYPE, J, K1, K2, K3, K4
       DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME, DISNAN
 !      DOUBLE PRECISION   DLAMCH
 !      EXTERNAL           lsame, dlamch, disnan
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, max, min
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
 !
       IF( lsame( TYPE, 'G' ) ) then
          itype = 0
       ELSE IF( lsame( TYPE, 'L' ) ) then
          itype = 1
       ELSE IF( lsame( TYPE, 'U' ) ) then
          itype = 2
       ELSE IF( lsame( TYPE, 'H' ) ) then
          itype = 3
       ELSE IF( lsame( TYPE, 'B' ) ) then
          itype = 4
       ELSE IF( lsame( TYPE, 'Q' ) ) then
          itype = 5
       ELSE IF( lsame( TYPE, 'Z' ) ) then
          itype = 6
       ELSE
          itype = -1
       END IF
 !
       IF( itype.EQ.-1 ) THEN
          info = -1
       ELSE IF( cfrom.EQ.zero .OR. disnan(cfrom) ) THEN
          info = -4
       ELSE IF( disnan(cto) ) THEN
          info = -5
       ELSE IF( m.LT.0 ) THEN
          info = -6
       ELSE IF( n.LT.0 .OR. ( itype.EQ.4 .AND. n.NE.m ) .OR. &
                ( itype.EQ.5 .AND. n.NE.m ) ) THEN
          info = -7
       ELSE IF( itype.LE.3 .AND. lda.LT.max( 1, m ) ) THEN
          info = -9
       ELSE IF( itype.GE.4 ) THEN
          IF( kl.LT.0 .OR. kl.GT.max( m-1, 0 ) ) THEN
             info = -2
          ELSE IF( ku.LT.0 .OR. ku.GT.max( n-1, 0 ) .OR. &
                   ( ( itype.EQ.4 .OR. itype.EQ.5 ) .AND. kl.NE.ku ) ) &
                    THEN
             info = -3
          ELSE IF( ( itype.EQ.4 .AND. lda.LT.kl+1 ) .OR. &
                   ( itype.EQ.5 .AND. lda.LT.ku+1 ) .OR. &
                   ( itype.EQ.6 .AND. lda.LT.2*kl+ku+1 ) ) THEN
             info = -9
          END IF
       END IF
 !
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DLASCL', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 .OR. m.EQ.0 ) &
          RETURN
 !
 !     Get machine parameters
 !
       smlnum = dlamch( 'S' )
       bignum = one / smlnum
 !
       cfromc = cfrom
       ctoc = cto
 !
    10 CONTINUE
       cfrom1 = cfromc*smlnum
       IF( cfrom1.EQ.cfromc ) THEN
 !        CFROMC is an inf.  Multiply by a correctly signed zero for
 !        finite CTOC, or a NaN if CTOC is infinite.
          mul = ctoc / cfromc
          done = .true.
          cto1 = ctoc
       ELSE
          cto1 = ctoc / bignum
          IF( cto1.EQ.ctoc ) THEN
 !           CTOC is either 0 or an inf.  In both cases, CTOC itself
 !           serves as the correct multiplication factor.
             mul = ctoc
             done = .true.
             cfromc = one
          ELSE IF( abs( cfrom1 ).GT.abs( ctoc ) .AND. ctoc.NE.zero ) THEN
             mul = smlnum
             done = .false.
             cfromc = cfrom1
          ELSE IF( abs( cto1 ).GT.abs( cfromc ) ) THEN
             mul = bignum
             done = .false.
             ctoc = cto1
          ELSE
             mul = ctoc / cfromc
             done = .true.
          END IF
       END IF
 !
       IF( itype.EQ.0 ) THEN
 !
 !        Full matrix
 !
          DO 30 j = 1, n
             DO 20 i = 1, m
                a( i, j ) = a( i, j )*mul
    20       CONTINUE
    30    CONTINUE
 !
       ELSE IF( itype.EQ.1 ) THEN
 !
 !        Lower triangular matrix
 !
          DO 50 j = 1, n
             DO 40 i = j, m
                a( i, j ) = a( i, j )*mul
    40       CONTINUE
    50    CONTINUE
 !
       ELSE IF( itype.EQ.2 ) THEN
 !
 !        Upper triangular matrix
 !
          DO 70 j = 1, n
             DO 60 i = 1, min( j, m )
                a( i, j ) = a( i, j )*mul
    60       CONTINUE
    70    CONTINUE
 !
       ELSE IF( itype.EQ.3 ) THEN
 !
 !        Upper Hessenberg matrix
 !
          DO 90 j = 1, n
             DO 80 i = 1, min( j+1, m )
                a( i, j ) = a( i, j )*mul
    80       CONTINUE
    90    CONTINUE
 !
       ELSE IF( itype.EQ.4 ) THEN
 !
 !        Lower half of a symmetric band matrix
 !
          k3 = kl + 1
          k4 = n + 1
          DO 110 j = 1, n
             DO 100 i = 1, min( k3, k4-j )
                a( i, j ) = a( i, j )*mul
   100       CONTINUE
   110    CONTINUE
 !
       ELSE IF( itype.EQ.5 ) THEN
 !
 !        Upper half of a symmetric band matrix
 !
          k1 = ku + 2
          k3 = ku + 1
          DO 130 j = 1, n
             DO 120 i = max( k1-j, 1 ), k3
                a( i, j ) = a( i, j )*mul
   120       CONTINUE
   130    CONTINUE
 !
       ELSE IF( itype.EQ.6 ) THEN
 !
 !        Band matrix
 !
          k1 = kl + ku + 2
          k2 = kl + 1
          k3 = 2*kl + ku + 1
          k4 = kl + ku + 1 + m
          DO 150 j = 1, n
             DO 140 i = max( k1-j, k2 ), min( k3, k4-j )
                a( i, j ) = a( i, j )*mul
   140       CONTINUE
   150    CONTINUE
 !
       END IF
 !
       IF( .NOT.done ) &
          GO TO 10
 !
       RETURN
 !
 !     End of DLASCL
 !
  end subroutine dlascl

 !  =====================================================================
  subroutine dlasrt( ID, N, D, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          ID
       INTEGER            INFO, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   D( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       INTEGER            SELECT
       parameter( SELECT = 20 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            DIR, ENDD, I, J, START, STKPNT
       DOUBLE PRECISION   D1, D2, D3, DMNMX, TMP
 !     ..
 !     .. Local Arrays ..
       INTEGER            STACK( 2, 32 )
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      EXTERNAL           lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input paramters.
 !
       info = 0
       dir = -1
       IF( lsame( id, 'D' ) ) THEN
          dir = 0
       ELSE IF( lsame( id, 'I' ) ) THEN
          dir = 1
       END IF
       IF( dir.EQ.-1 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DLASRT', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.LE.1 ) &
          RETURN
 !
       stkpnt = 1
       stack( 1, 1 ) = 1
       stack( 2, 1 ) = n
    10 CONTINUE
       start = stack( 1, stkpnt )
       ENDD = STACK( 2, STKPNT )
       stkpnt = stkpnt - 1
       IF( endd-start.LE.SELECT .AND. endd-start.GT.0 ) THEN
 !
 !        Do Insertion sort on D( START:ENDD )
 !
          IF( dir.EQ.0 ) THEN
 !
 !           Sort into decreasing order
 !
             DO 30 i = start + 1, endd
                DO 20 j = i, start + 1, -1
                   IF( d( j ).GT.d( j-1 ) ) THEN
                      dmnmx = d( j )
                      d( j ) = d( j-1 )
                      d( j-1 ) = dmnmx
                   ELSE
                      GO TO 30
                   END IF
    20          CONTINUE
    30       CONTINUE
 !
          ELSE
 !
 !           Sort into increasing order
 !
             DO 50 i = start + 1, endd
                DO 40 j = i, start + 1, -1
                   IF( d( j ).LT.d( j-1 ) ) THEN
                      dmnmx = d( j )
                      d( j ) = d( j-1 )
                      d( j-1 ) = dmnmx
                   ELSE
                      GO TO 50
                   END IF
    40          CONTINUE
    50       CONTINUE
 !
          END IF
 !
       ELSE IF( endd-start.GT.SELECT ) THEN
 !
 !        Partition D( START:ENDD ) and stack parts, largest one first
 !
 !        Choose partition entry as median of 3
 !
          d1 = d( start )
          d2 = d( endd )
          i = ( start+endd ) / 2
          d3 = d( i )
          IF( d1.LT.d2 ) THEN
             IF( d3.LT.d1 ) THEN
                dmnmx = d1
             ELSE IF( d3.LT.d2 ) THEN
                dmnmx = d3
             ELSE
                dmnmx = d2
             END IF
          ELSE
             IF( d3.LT.d2 ) THEN
                dmnmx = d2
             ELSE IF( d3.LT.d1 ) THEN
                dmnmx = d3
             ELSE
                dmnmx = d1
             END IF
          END IF
 !
          IF( dir.EQ.0 ) THEN
 !
 !           Sort into decreasing order
 !
             i = start - 1
             j = endd + 1
    60       CONTINUE
    70       CONTINUE
             j = j - 1
             IF( d( j ).LT.dmnmx ) &
                GO TO 70
    80       CONTINUE
             i = i + 1
             IF( d( i ).GT.dmnmx ) &
                GO TO 80
             IF( i.LT.j ) THEN
                tmp = d( i )
                d( i ) = d( j )
                d( j ) = tmp
                GO TO 60
             END IF
             IF( j-start.GT.endd-j-1 ) THEN
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
             ELSE
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
             END IF
          ELSE
 !
 !           Sort into increasing order
 !
             i = start - 1
             j = endd + 1
    90       CONTINUE
   100       CONTINUE
             j = j - 1
             IF( d( j ).GT.dmnmx ) &
                GO TO 100
   110       CONTINUE
             i = i + 1
             IF( d( i ).LT.dmnmx ) &
                GO TO 110
             IF( i.LT.j ) THEN
                tmp = d( i )
                d( i ) = d( j )
                d( j ) = tmp
                GO TO 90
             END IF
             IF( j-start.GT.endd-j-1 ) THEN
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
             ELSE
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = j + 1
                stack( 2, stkpnt ) = endd
                stkpnt = stkpnt + 1
                stack( 1, stkpnt ) = start
                stack( 2, stkpnt ) = j
             END IF
          END IF
       END IF
       IF( stkpnt.GT.0 ) &
          GO TO 10
       RETURN
 !
 !     End of DLASRT
 !
  end subroutine dlasrt

 !  =====================================================================
  subroutine zhetrd( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LWORK, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   D( * ), E( * )
       COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d0 )
       COMPLEX*16         CONE
       parameter( cone = ( 1.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            LQUERY, UPPER
       INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, &
                          nbmin, nx
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zher2k, zhetd2, zlatrd
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      INTEGER            ILAENV
 !      EXTERNAL           lsame, ilaenv
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters
 !
       info = 0
       upper = lsame( uplo, 'U' )
       lquery = ( lwork.EQ.-1 )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -4
       ELSE IF( lwork.LT.1 .AND. .NOT.lquery ) THEN
          info = -9
       END IF
 !
       IF( info.EQ.0 ) THEN
 !
 !        Determine the block size.
 !
          nb = ilaenv( 1, 'ZHETRD', uplo, n, -1, -1, -1 )
          lwkopt = n*nb
          work( 1 ) = lwkopt
       END IF
 !
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZHETRD', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 ) THEN
          work( 1 ) = 1
          RETURN
       END IF
 !
       nx = n
       iws = 1
       IF( nb.GT.1 .AND. nb.LT.n ) THEN
 !
 !        Determine when to cross over from blocked to unblocked code
 !        (last block is always handled by unblocked code).
 !
          nx = max( nb, ilaenv( 3, 'ZHETRD', uplo, n, -1, -1, -1 ) )
          IF( nx.LT.n ) THEN
 !
 !           Determine if workspace is large enough for blocked code.
 !
             ldwork = n
             iws = ldwork*nb
             IF( lwork.LT.iws ) THEN
 !
 !              Not enough workspace to use optimal NB:  determine the
 !              minimum value of NB, and reduce NB or force use of
 !              unblocked code by setting NX = N.
 !
                nb = max( lwork / ldwork, 1 )
                nbmin = ilaenv( 2, 'ZHETRD', uplo, n, -1, -1, -1 )
                IF( nb.LT.nbmin ) &
                   nx = n
             END IF
          ELSE
             nx = n
          END IF
       ELSE
          nb = 1
       END IF
 !
       IF( upper ) THEN
 !
 !        Reduce the upper triangle of A.
 !        Columns 1:kk are handled by the unblocked method.
 !
          kk = n - ( ( n-nx+nb-1 ) / nb )*nb
          DO 20 i = n - nb + 1, kk + 1, -nb
 !
 !           Reduce columns i:i+nb-1 to tridiagonal form and form the
 !           matrix W which is needed to update the unreduced part of
 !           the matrix
 !
             CALL zlatrd( uplo, i+nb-1, nb, a, lda, e, tau, work, &
                          ldwork )
 !
 !           Update the unreduced submatrix A(1:i-1,1:i-1), using an
 !           update of the form:  A := A - V*W**H - W*V**H
 !
             CALL zher2k( uplo, 'No transpose', i-1, nb, -cone, &
                          a( 1, i ), lda, work, ldwork, one, a, lda )
 !
 !           Copy superdiagonal elements back into A, and diagonal
 !           elements into D
 !
             DO 10 j = i, i + nb - 1
                a( j-1, j ) = e( j-1 )
                d( j ) = a( j, j )
    10       CONTINUE
    20    CONTINUE
 !
 !        Use unblocked code to reduce the last or only block
 !
          CALL zhetd2( uplo, kk, a, lda, d, e, tau, iinfo )
       ELSE
 !
 !        Reduce the lower triangle of A
 !
          DO 40 i = 1, n - nx, nb
 !
 !           Reduce columns i:i+nb-1 to tridiagonal form and form the
 !           matrix W which is needed to update the unreduced part of
 !           the matrix
 !
             CALL zlatrd( uplo, n-i+1, nb, a( i, i ), lda, e( i ), &
                          tau( i ), work, ldwork )
 !
 !           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
 !           an update of the form:  A := A - V*W**H - W*V**H
 !
             CALL zher2k( uplo, 'No transpose', n-i-nb+1, nb, -cone, &
                          a( i+nb, i ), lda, work( nb+1 ), ldwork, one, &
                          a( i+nb, i+nb ), lda )
 !
 !           Copy subdiagonal elements back into A, and diagonal
 !           elements into D
 !
             DO 30 j = i, i + nb - 1
                a( j+1, j ) = e( j )
                d( j ) = a( j, j )
    30       CONTINUE
    40    CONTINUE
 !
 !        Use unblocked code to reduce the last or only block
 !
          CALL zhetd2( uplo, n-i+1, a( i, i ), lda, d( i ), e( i ), &
                       tau( i ), iinfo )
       END IF
 !
       work( 1 ) = lwkopt
       RETURN
 !
 !     End of ZHETRD
 !
  end subroutine zhetrd

 !  =====================================================================
  subroutine zher2k(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
 !
 !  -- Reference BLAS level3 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       DOUBLE PRECISION BETA
       INTEGER K,LDA,LDB,LDC,N
       CHARACTER TRANS,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),B(ldb,*),C(ldc,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
 !      LOGICAL LSAME
 !      EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dble,dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP1,TEMP2
       INTEGER I,INFO,J,L,NROWA
       LOGICAL UPPER
 !     ..
 !     .. Parameters ..
       DOUBLE PRECISION ONE
       parameter(one=1.0d0)
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !
 !     Test the input parameters.
 !
       IF (lsame(trans,'N')) THEN
           nrowa = n
       ELSE
           nrowa = k
       END IF
       upper = lsame(uplo,'U')
 !
       info = 0
       IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 1
       ELSE IF ((.NOT.lsame(trans,'N')) .AND. &
                (.NOT.lsame(trans,'C'))) THEN
           info = 2
       ELSE IF (n.LT.0) THEN
           info = 3
       ELSE IF (k.LT.0) THEN
           info = 4
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 7
       ELSE IF (ldb.LT.max(1,nrowa)) THEN
           info = 9
       ELSE IF (ldc.LT.max(1,n)) THEN
           info = 12
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZHER2K',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((n.EQ.0) .OR. (((alpha.EQ.zero).OR. &
           (k.EQ.0)).AND. (beta.EQ.one))) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           IF (upper) THEN
               IF (beta.EQ.dble(zero)) THEN
                   DO 20 j = 1,n
                       DO 10 i = 1,j
                           c(i,j) = zero
    10                 CONTINUE
    20             CONTINUE
               ELSE
                   DO 40 j = 1,n
                       DO 30 i = 1,j - 1
                           c(i,j) = beta*c(i,j)
    30                 CONTINUE
                       c(j,j) = beta*dble(c(j,j))
    40             CONTINUE
               END IF
           ELSE
               IF (beta.EQ.dble(zero)) THEN
                   DO 60 j = 1,n
                       DO 50 i = j,n
                           c(i,j) = zero
    50                 CONTINUE
    60             CONTINUE
               ELSE
                   DO 80 j = 1,n
                       c(j,j) = beta*dble(c(j,j))
                       DO 70 i = j + 1,n
                           c(i,j) = beta*c(i,j)
    70                 CONTINUE
    80             CONTINUE
               END IF
           END IF
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (lsame(trans,'N')) THEN
 !
 !        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
 !                   C.
 !
           IF (upper) THEN
               DO 130 j = 1,n
                   IF (beta.EQ.dble(zero)) THEN
                       DO 90 i = 1,j
                           c(i,j) = zero
    90                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 100 i = 1,j - 1
                           c(i,j) = beta*c(i,j)
   100                 CONTINUE
                       c(j,j) = beta*dble(c(j,j))
                   ELSE
                       c(j,j) = dble(c(j,j))
                   END IF
                   DO 120 l = 1,k
                       IF ((a(j,l).NE.zero) .OR. (b(j,l).NE.zero)) THEN
                           temp1 = alpha*dconjg(b(j,l))
                           temp2 = dconjg(alpha*a(j,l))
                           DO 110 i = 1,j - 1
                               c(i,j) = c(i,j) + a(i,l)*temp1 + &
                                        b(i,l)*temp2
   110                     CONTINUE
                           c(j,j) = dble(c(j,j)) + &
                                    dble(a(j,l)*temp1+b(j,l)*temp2)
                       END IF
   120             CONTINUE
   130         CONTINUE
           ELSE
               DO 180 j = 1,n
                   IF (beta.EQ.dble(zero)) THEN
                       DO 140 i = j,n
                           c(i,j) = zero
   140                 CONTINUE
                   ELSE IF (beta.NE.one) THEN
                       DO 150 i = j + 1,n
                           c(i,j) = beta*c(i,j)
   150                 CONTINUE
                       c(j,j) = beta*dble(c(j,j))
                   ELSE
                       c(j,j) = dble(c(j,j))
                   END IF
                   DO 170 l = 1,k
                       IF ((a(j,l).NE.zero) .OR. (b(j,l).NE.zero)) THEN
                           temp1 = alpha*dconjg(b(j,l))
                           temp2 = dconjg(alpha*a(j,l))
                           DO 160 i = j + 1,n
                               c(i,j) = c(i,j) + a(i,l)*temp1 + &
                                        b(i,l)*temp2
   160                     CONTINUE
                           c(j,j) = dble(c(j,j)) + &
                                    dble(a(j,l)*temp1+b(j,l)*temp2)
                       END IF
   170             CONTINUE
   180         CONTINUE
           END IF
       ELSE
 !
 !        Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
 !                   C.
 !
           IF (upper) THEN
               DO 210 j = 1,n
                   DO 200 i = 1,j
                       temp1 = zero
                       temp2 = zero
                       DO 190 l = 1,k
                           temp1 = temp1 + dconjg(a(l,i))*b(l,j)
                           temp2 = temp2 + dconjg(b(l,i))*a(l,j)
   190                 CONTINUE
                       IF (i.EQ.j) THEN
                           IF (beta.EQ.dble(zero)) THEN
                               c(j,j) = dble(alpha*temp1+ &
                                       dconjg(alpha)*temp2)
                           ELSE
                               c(j,j) = beta*dble(c(j,j)) + &
                                        dble(alpha*temp1+ &
                                        dconjg(alpha)*temp2)
                           END IF
                       ELSE
                           IF (beta.EQ.dble(zero)) THEN
                               c(i,j) = alpha*temp1 + dconjg(alpha)*temp2
                           ELSE
                               c(i,j) = beta*c(i,j) + alpha*temp1 + &
                                        dconjg(alpha)*temp2
                           END IF
                       END IF
   200             CONTINUE
   210         CONTINUE
           ELSE
               DO 240 j = 1,n
                   DO 230 i = j,n
                       temp1 = zero
                       temp2 = zero
                       DO 220 l = 1,k
                           temp1 = temp1 + dconjg(a(l,i))*b(l,j)
                           temp2 = temp2 + dconjg(b(l,i))*a(l,j)
   220                 CONTINUE
                       IF (i.EQ.j) THEN
                           IF (beta.EQ.dble(zero)) THEN
                               c(j,j) = dble(alpha*temp1+ &
                                        dconjg(alpha)*temp2)
                           ELSE
                               c(j,j) = beta*dble(c(j,j)) + &
                                        dble(alpha*temp1+ &
                                        dconjg(alpha)*temp2)
                           END IF
                       ELSE
                           IF (beta.EQ.dble(zero)) THEN
                               c(i,j) = alpha*temp1 + dconjg(alpha)*temp2
                           ELSE
                               c(i,j) = beta*c(i,j) + alpha*temp1 + &
                                        dconjg(alpha)*temp2
                           END IF
                       END IF
   230             CONTINUE
   240         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZHER2K.
 !
 end subroutine zher2k

 !  =====================================================================
  subroutine zhetd2( UPLO, N, A, LDA, D, E, TAU, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   D( * ), E( * )
       COMPLEX*16         A( lda, * ), TAU( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO, HALF
       parameter( one = ( 1.0d0, 0.0d0 ), &
                          zero = ( 0.0d0, 0.0d0 ), &
                          half = ( 0.5d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            UPPER
       INTEGER            I
       COMPLEX*16         ALPHA, TAUI
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zaxpy, zhemv, zher2, zlarfg
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      COMPLEX*16         ZDOTC
 !      EXTERNAL           lsame, zdotc
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          dble, max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters
 !
       info = 0
       upper = lsame( uplo, 'U')
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -4
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZHETD2', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.LE.0 ) &
          RETURN
 !
       IF( upper ) THEN
 !
 !        Reduce the upper triangle of A
 !
          a( n, n ) = dble( a( n, n ) )
          DO 10 i = n - 1, 1, -1
 !
 !           Generate elementary reflector H(i) = I - tau * v * v**H
 !           to annihilate A(1:i-1,i+1)
 !
             alpha = a( i, i+1 )
             CALL zlarfg( i, alpha, a( 1, i+1 ), 1, taui )
             e( i ) = alpha
 !
             IF( taui.NE.zero ) THEN
 !
 !              Apply H(i) from both sides to A(1:i,1:i)
 !
                a( i, i+1 ) = one
 !
 !              Compute  x := tau * A * v  storing x in TAU(1:i)
 !
                CALL zhemv( uplo, i, taui, a, lda, a( 1, i+1 ), 1, zero, &
                            tau, 1 )
 !
 !              Compute  w := x - 1/2 * tau * (x**H * v) * v
 !
                alpha = -half*taui*zdotc( i, tau, 1, a( 1, i+1 ), 1 )
                CALL zaxpy( i, alpha, a( 1, i+1 ), 1, tau, 1 )
 !
 !              Apply the transformation as a rank-2 update:
 !                 A := A - v * w**H - w * v**H
 !
                CALL zher2( uplo, i, -one, a( 1, i+1 ), 1, tau, 1, a, &
                            lda )
 !
             ELSE
                a( i, i ) = dble( a( i, i ) )
             END IF
             a( i, i+1 ) = e( i )
             d( i+1 ) = a( i+1, i+1 )
             tau( i ) = taui
    10    CONTINUE
          d( 1 ) = a( 1, 1 )
       ELSE
 !
 !        Reduce the lower triangle of A
 !
          a( 1, 1 ) = dble( a( 1, 1 ) )
          DO 20 i = 1, n - 1
 !
 !           Generate elementary reflector H(i) = I - tau * v * v**H
 !           to annihilate A(i+2:n,i)
 !
             alpha = a( i+1, i )
             CALL zlarfg( n-i, alpha, a( min( i+2, n ), i ), 1, taui )
             e( i ) = alpha
 !
             IF( taui.NE.zero ) THEN
 !
 !              Apply H(i) from both sides to A(i+1:n,i+1:n)
 !
                a( i+1, i ) = one
 !
 !              Compute  x := tau * A * v  storing y in TAU(i:n-1)
 !
                CALL zhemv( uplo, n-i, taui, a( i+1, i+1 ), lda, &
                            a( i+1, i ), 1, zero, tau( i ), 1 )
 !
 !              Compute  w := x - 1/2 * tau * (x**H * v) * v
 !
                alpha = -half*taui*zdotc( n-i, tau( i ), 1, a( i+1, i ), 1 )
                CALL zaxpy( n-i, alpha, a( i+1, i ), 1, tau( i ), 1 )
 !
 !              Apply the transformation as a rank-2 update:
 !                 A := A - v * w**H - w * v**H
 !
                CALL zher2( uplo, n-i, -one, a( i+1, i ), 1, tau( i ), 1, &
                            a( i+1, i+1 ), lda )
 !
             ELSE
                a( i+1, i+1 ) = dble( a( i+1, i+1 ) )
             END IF
             a( i+1, i ) = e( i )
             d( i ) = a( i, i )
             tau( i ) = taui
    20    CONTINUE
          d( n ) = a( n, n )
       END IF
 !
       RETURN
 !
 !     End of ZHETD2
 !
  end subroutine zhetd2

 !  =====================================================================
  subroutine zaxpy(N,ZA,ZX,INCX,ZY,INCY)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ZA
       INTEGER INCX,INCY,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*),ZY(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,IX,IY
 !     ..
 !     .. External Functions ..
 !      DOUBLE PRECISION DCABS1
 !      EXTERNAL dcabs1
 !     ..
       IF (n.LE.0) RETURN
       IF (dcabs1(za).EQ.0.0d0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 !
 !        code for both increments equal to 1
 !
          DO i = 1,n
             zy(i) = zy(i) + za*zx(i)
          END DO
       ELSE
 !
 !        code for unequal increments or equal increments
 !          not equal to 1
 !
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             zy(iy) = zy(iy) + za*zx(ix)
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
 !
       RETURN
  end subroutine zaxpy

 !  =====================================================================
  double precision FUNCTION dcabs1(Z)
 !
 !  -- Reference BLAS level1 routine (version 3.6.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2015
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 Z
 !     ..
 !     ..
 !  =====================================================================
 !
 !     .. Intrinsic Functions ..
       INTRINSIC abs,dble,dimag
 !
       dcabs1 = abs(dble(z)) + abs(dimag(z))
       RETURN
  end function dcabs1

 !  =====================================================================
  subroutine zhemv(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
 !
 !  -- Reference BLAS level2 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA,BETA
       INTEGER INCX,INCY,LDA,N
       CHARACTER UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),X(*),Y(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d0,0.0d0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP1,TEMP2
       INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
 !     ..
 !     .. External Functions ..
 !      LOGICAL LSAME
 !      EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dble,dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
           info = 1
       ELSE IF (n.LT.0) THEN
           info = 2
       ELSE IF (lda.LT.max(1,n)) THEN
           info = 5
       ELSE IF (incx.EQ.0) THEN
           info = 7
       ELSE IF (incy.EQ.0) THEN
           info = 10
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZHEMV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((n.EQ.0) .OR. ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
 !
 !     Set up the start points in  X  and  Y.
 !
       IF (incx.GT.0) THEN
           kx = 1
       ELSE
           kx = 1 - (n-1)*incx
       END IF
       IF (incy.GT.0) THEN
           ky = 1
       ELSE
           ky = 1 - (n-1)*incy
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through the triangular part
 !     of A.
 !
 !     First form  y := beta*y.
 !
       IF (beta.NE.one) THEN
           IF (incy.EQ.1) THEN
               IF (beta.EQ.zero) THEN
                   DO 10 i = 1,n
                       y(i) = zero
    10             CONTINUE
               ELSE
                   DO 20 i = 1,n
                       y(i) = beta*y(i)
    20             CONTINUE
               END IF
           ELSE
               iy = ky
               IF (beta.EQ.zero) THEN
                   DO 30 i = 1,n
                       y(iy) = zero
                       iy = iy + incy
    30             CONTINUE
               ELSE
                   DO 40 i = 1,n
                       y(iy) = beta*y(iy)
                       iy = iy + incy
    40             CONTINUE
               END IF
           END IF
       END IF
       IF (alpha.EQ.zero) RETURN
       IF (lsame(uplo,'U')) THEN
 !
 !        Form  y  when A is stored in upper triangle.
 !
           IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
               DO 60 j = 1,n
                   temp1 = alpha*x(j)
                   temp2 = zero
                   DO 50 i = 1,j - 1
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + dconjg(a(i,j))*x(i)
    50             CONTINUE
                   y(j) = y(j) + temp1*dble(a(j,j)) + alpha*temp2
    60         CONTINUE
           ELSE
               jx = kx
               jy = ky
               DO 80 j = 1,n
                   temp1 = alpha*x(jx)
                   temp2 = zero
                   ix = kx
                   iy = ky
                   DO 70 i = 1,j - 1
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + dconjg(a(i,j))*x(ix)
                       ix = ix + incx
                       iy = iy + incy
    70             CONTINUE
                   y(jy) = y(jy) + temp1*dble(a(j,j)) + alpha*temp2
                   jx = jx + incx
                   jy = jy + incy
    80         CONTINUE
           END IF
       ELSE
 !
 !        Form  y  when A is stored in lower triangle.
 !
           IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
               DO 100 j = 1,n
                   temp1 = alpha*x(j)
                   temp2 = zero
                   y(j) = y(j) + temp1*dble(a(j,j))
                   DO 90 i = j + 1,n
                       y(i) = y(i) + temp1*a(i,j)
                       temp2 = temp2 + dconjg(a(i,j))*x(i)
    90             CONTINUE
                   y(j) = y(j) + alpha*temp2
   100         CONTINUE
           ELSE
               jx = kx
               jy = ky
               DO 120 j = 1,n
                   temp1 = alpha*x(jx)
                   temp2 = zero
                   y(jy) = y(jy) + temp1*dble(a(j,j))
                   ix = jx
                   iy = jy
                   DO 110 i = j + 1,n
                       ix = ix + incx
                       iy = iy + incy
                       y(iy) = y(iy) + temp1*a(i,j)
                       temp2 = temp2 + dconjg(a(i,j))*x(ix)
   110             CONTINUE
                   y(jy) = y(jy) + alpha*temp2
                   jx = jx + incx
                   jy = jy + incy
   120         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZHEMV .
 !
  end subroutine zhemv

 !  =====================================================================
  subroutine zher2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
 !
 !  -- Reference BLAS level2 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       INTEGER INCX,INCY,LDA,N
       CHARACTER UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),X(*),Y(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP1,TEMP2
       INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY
 !     ..
 !     .. External Functions ..
 !      LOGICAL LSAME
 !      EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dble,dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
           info = 1
       ELSE IF (n.LT.0) THEN
           info = 2
       ELSE IF (incx.EQ.0) THEN
           info = 5
       ELSE IF (incy.EQ.0) THEN
           info = 7
       ELSE IF (lda.LT.max(1,n)) THEN
           info = 9
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZHER2 ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
 !
 !     Set up the start points in X and Y if the increments are not both
 !     unity.
 !
       IF ((incx.NE.1) .OR. (incy.NE.1)) THEN
           IF (incx.GT.0) THEN
               kx = 1
           ELSE
               kx = 1 - (n-1)*incx
           END IF
           IF (incy.GT.0) THEN
               ky = 1
           ELSE
               ky = 1 - (n-1)*incy
           END IF
           jx = kx
           jy = ky
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through the triangular part
 !     of A.
 !
       IF (lsame(uplo,'U')) THEN
 !
 !        Form  A  when A is stored in the upper triangle.
 !
           IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
               DO 20 j = 1,n
                   IF ((x(j).NE.zero) .OR. (y(j).NE.zero)) THEN
                       temp1 = alpha*dconjg(y(j))
                       temp2 = dconjg(alpha*x(j))
                       DO 10 i = 1,j - 1
                           a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
    10                 CONTINUE
                       a(j,j) = dble(a(j,j)) + &
                                dble(x(j)*temp1+y(j)*temp2)
                   ELSE
                       a(j,j) = dble(a(j,j))
                   END IF
    20         CONTINUE
           ELSE
               DO 40 j = 1,n
                   IF ((x(jx).NE.zero) .OR. (y(jy).NE.zero)) THEN
                       temp1 = alpha*dconjg(y(jy))
                       temp2 = dconjg(alpha*x(jx))
                       ix = kx
                       iy = ky
                       DO 30 i = 1,j - 1
                           a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
                           ix = ix + incx
                           iy = iy + incy
    30                 CONTINUE
                       a(j,j) = dble(a(j,j)) + &
                                dble(x(jx)*temp1+y(jy)*temp2)
                   ELSE
                       a(j,j) = dble(a(j,j))
                   END IF
                   jx = jx + incx
                   jy = jy + incy
    40         CONTINUE
           END IF
       ELSE
 !
 !        Form  A  when A is stored in the lower triangle.
 !
           IF ((incx.EQ.1) .AND. (incy.EQ.1)) THEN
               DO 60 j = 1,n
                   IF ((x(j).NE.zero) .OR. (y(j).NE.zero)) THEN
                       temp1 = alpha*dconjg(y(j))
                       temp2 = dconjg(alpha*x(j))
                       a(j,j) = dble(a(j,j)) + &
                                dble(x(j)*temp1+y(j)*temp2)
                       DO 50 i = j + 1,n
                           a(i,j) = a(i,j) + x(i)*temp1 + y(i)*temp2
    50                 CONTINUE
                   ELSE
                       a(j,j) = dble(a(j,j))
                   END IF
    60         CONTINUE
           ELSE
               DO 80 j = 1,n
                   IF ((x(jx).NE.zero) .OR. (y(jy).NE.zero)) THEN
                       temp1 = alpha*dconjg(y(jy))
                       temp2 = dconjg(alpha*x(jx))
                       a(j,j) = dble(a(j,j)) + &
                                dble(x(jx)*temp1+y(jy)*temp2)
                       ix = jx
                       iy = jy
                       DO 70 i = j + 1,n
                           ix = ix + incx
                           iy = iy + incy
                           a(i,j) = a(i,j) + x(ix)*temp1 + y(iy)*temp2
    70                 CONTINUE
                   ELSE
                       a(j,j) = dble(a(j,j))
                   END IF
                   jx = jx + incx
                   jy = jy + incy
    80         CONTINUE
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZHER2 .
 !
  end subroutine zher2

 !  =====================================================================
  subroutine zlarfg( N, ALPHA, X, INCX, TAU )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       INTEGER            INCX, N
       COMPLEX*16         ALPHA, TAU
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         X( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE, ZERO
       parameter( one = 1.0d0, zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            J, KNT
       DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
 !     ..
 !     .. External Functions ..
 !      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
 !      COMPLEX*16         ZLADIV
 !      EXTERNAL           dlamch, dlapy3, dznrm2, zladiv
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, dble, dcmplx, dimag, sign
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           zdscal, zscal
 !     ..
 !     .. Executable Statements ..
 !
       IF( n.LE.0 ) THEN
          tau = zero
          RETURN
       END IF
 !
       xnorm = dznrm2( n-1, x, incx )
       alphr = dble( alpha )
       alphi = dimag( alpha )
 !
       IF( xnorm.EQ.zero .AND. alphi.EQ.zero ) THEN
 !
 !        H  =  I
 !
          tau = zero
       ELSE
 !
 !        general case
 !
          beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
          safmin = dlamch( 'S' ) / dlamch( 'E' )
          rsafmn = one / safmin
 !
          knt = 0
          IF( abs( beta ).LT.safmin ) THEN
 !
 !           XNORM, BETA may be inaccurate; scale X and recompute them
 !
    10       CONTINUE
             knt = knt + 1
             CALL zdscal( n-1, rsafmn, x, incx )
             beta = beta*rsafmn
             alphi = alphi*rsafmn
             alphr = alphr*rsafmn
             IF( abs( beta ).LT.safmin ) &
               GO TO 10
 !
 !           New BETA is at most 1, at least SAFMIN
 !
             xnorm = dznrm2( n-1, x, incx )
             alpha = dcmplx( alphr, alphi )
             beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
          END IF
          tau = dcmplx( ( beta-alphr ) / beta, -alphi / beta )
          alpha = zladiv( dcmplx( one ), alpha-beta )
          CALL zscal( n-1, alpha, x, incx )
 !
 !        If ALPHA is subnormal, it may lose relative accuracy
 !
          DO 20 j = 1, knt
             beta = beta*safmin
  20      CONTINUE
          alpha = beta
       END IF
 !
       RETURN
 !
 !     End of ZLARFG
 !
  end subroutine zlarfg

 !  =====================================================================
  double precision function dlapy3( X, Y, Z )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   X, Y, Z
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION   W, XABS, YABS, ZABS
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, max, sqrt
 !     ..
 !     .. Executable Statements ..
 !
       xabs = abs( x )
       yabs = abs( y )
       zabs = abs( z )
       w = max( xabs, yabs, zabs )
       IF( w.EQ.zero ) THEN
 !     W can be zero for max(0,nan,0)
 !     adding all three entries together will make sure
 !     NaN will not disappear.
          dlapy3 =  xabs + yabs + zabs
       ELSE
          dlapy3 = w*sqrt( ( xabs / w )**2+( yabs / w )**2+ &
                   ( zabs / w )**2 )
       END IF
       RETURN
 !
 !     End of DLAPY3
 !
  end function dlapy3

 !  =====================================================================
  double precision function dznrm2(N,X,INCX)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 X(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION ONE,ZERO
       parameter(one=1.0d0,zero=0.0d0)
 !     ..
 !     .. Local Scalars ..
       DOUBLE PRECISION NORM,SCALE,SSQ,TEMP
       INTEGER IX
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC abs,dble,dimag,sqrt
 !     ..
       IF (n.LT.1 .OR. incx.LT.1) THEN
           norm = zero
       ELSE
           scale = zero
           ssq = one
 !        The following loop is equivalent to this call to the LAPACK
 !        auxiliary routine:
 !        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
 !
           DO 10 ix = 1,1 + (n-1)*incx,incx
               IF (dble(x(ix)).NE.zero) THEN
                   temp = abs(dble(x(ix)))
                   IF (scale.LT.temp) THEN
                       ssq = one + ssq* (scale/temp)**2
                       scale = temp
                   ELSE
                       ssq = ssq + (temp/scale)**2
                   END IF
               END IF
               IF (dimag(x(ix)).NE.zero) THEN
                   temp = abs(dimag(x(ix)))
                   IF (scale.LT.temp) THEN
                       ssq = one + ssq* (scale/temp)**2
                       scale = temp
                   ELSE
                       ssq = ssq + (temp/scale)**2
                   END IF
               END IF
    10     CONTINUE
           norm = scale*sqrt(ssq)
       END IF
 !
       dznrm2 = norm
       RETURN
 !
 !     End of DZNRM2.
 !
  end function dznrm2

 !  =====================================================================
  complex*16 function zladiv( X, Y )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       COMPLEX*16         X, Y
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       DOUBLE PRECISION   ZI, ZR
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           dladiv
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          dble, dcmplx, dimag
 !     ..
 !     .. Executable Statements ..
 !
       CALL dladiv( dble( x ), dimag( x ), dble( y ), dimag( y ), zr, zi )
       zladiv = dcmplx( zr, zi )
 !
       RETURN
 !
 !     End of ZLADIV
 !
  end function zladiv

 !  =====================================================================
  subroutine zscal(N,ZA,ZX,INCX)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ZA
       INTEGER INCX,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,NINCX
 !     ..
       IF (n.LE.0 .OR. incx.LE.0) RETURN
       IF (incx.EQ.1) THEN
 !
 !        code for increment equal to 1
 !
          DO i = 1,n
             zx(i) = za*zx(i)
          END DO
       ELSE
 !
 !        code for increment not equal to 1
 !
          nincx = n*incx
          DO i = 1,nincx,incx
             zx(i) = za*zx(i)
          END DO
       END IF
       RETURN
  end subroutine zscal

 !  =====================================================================
  subroutine zlatrd( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            LDA, LDW, N, NB
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   E( * )
       COMPLEX*16         A( lda, * ), TAU( * ), W( ldw, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ZERO, ONE, HALF
       parameter( zero = ( 0.0d0, 0.0d0 ), &
                          one = ( 1.0d0, 0.0d0 ), &
                          half = ( 0.5d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, IW
       COMPLEX*16         ALPHA
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           zaxpy, zgemv, zhemv, zlacgv, zlarfg, zscal
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      COMPLEX*16         ZDOTC
 !      EXTERNAL           lsame, zdotc
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          dble, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Quick return if possible
 !
       IF( n.LE.0 ) &
          RETURN
 !
       IF( lsame( uplo, 'U' ) ) THEN
 !
 !        Reduce last NB columns of upper triangle
 !
          DO 10 i = n, n - nb + 1, -1
             iw = i - n + nb
             IF( i.LT.n ) THEN
 !
 !              Update A(1:i,i)
 !
                a( i, i ) = dble( a( i, i ) )
                CALL zlacgv( n-i, w( i, iw+1 ), ldw )
                CALL zgemv( 'No transpose', i, n-i, -one, a( 1, i+1 ), &
                            lda, w( i, iw+1 ), ldw, one, a( 1, i ), 1 )
                CALL zlacgv( n-i, w( i, iw+1 ), ldw )
                CALL zlacgv( n-i, a( i, i+1 ), lda )
                CALL zgemv( 'No transpose', i, n-i, -one, w( 1, iw+1 ), &
                            ldw, a( i, i+1 ), lda, one, a( 1, i ), 1 )
                CALL zlacgv( n-i, a( i, i+1 ), lda )
                a( i, i ) = dble( a( i, i ) )
             END IF
             IF( i.GT.1 ) THEN
 !
 !              Generate elementary reflector H(i) to annihilate
 !              A(1:i-2,i)
 !
                alpha = a( i-1, i )
                CALL zlarfg( i-1, alpha, a( 1, i ), 1, tau( i-1 ) )
                e( i-1 ) = alpha
                a( i-1, i ) = one
 !
 !              Compute W(1:i-1,i)
 !
                CALL zhemv( 'Upper', i-1, one, a, lda, a( 1, i ), 1, &
                            zero, w( 1, iw ), 1 )
                IF( i.LT.n ) THEN
                   CALL zgemv( 'Conjugate transpose', i-1, n-i, one, &
                               w( 1, iw+1 ), ldw, a( 1, i ), 1, zero, &
                               w( i+1, iw ), 1 )
                   CALL zgemv( 'No transpose', i-1, n-i, -one, &
                               a( 1, i+1 ), lda, w( i+1, iw ), 1, one, &
                               w( 1, iw ), 1 )
                   CALL zgemv( 'Conjugate transpose', i-1, n-i, one, &
                               a( 1, i+1 ), lda, a( 1, i ), 1, zero, &
                               w( i+1, iw ), 1 )
                   CALL zgemv( 'No transpose', i-1, n-i, -one, &
                               w( 1, iw+1 ), ldw, w( i+1, iw ), 1, one, &
                               w( 1, iw ), 1 )
                END IF
                CALL zscal( i-1, tau( i-1 ), w( 1, iw ), 1 )
                alpha = -half*tau( i-1 )*zdotc( i-1, w( 1, iw ), 1, &
                        a( 1, i ), 1 )
                CALL zaxpy( i-1, alpha, a( 1, i ), 1, w( 1, iw ), 1 )
             END IF
 !
    10    CONTINUE
       ELSE
 !
 !        Reduce first NB columns of lower triangle
 !
          DO 20 i = 1, nb
 !
 !           Update A(i:n,i)
 !
             a( i, i ) = dble( a( i, i ) )
             CALL zlacgv( i-1, w( i, 1 ), ldw )
             CALL zgemv( 'No transpose', n-i+1, i-1, -one, a( i, 1 ), &
                         lda, w( i, 1 ), ldw, one, a( i, i ), 1 )
             CALL zlacgv( i-1, w( i, 1 ), ldw )
             CALL zlacgv( i-1, a( i, 1 ), lda )
             CALL zgemv( 'No transpose', n-i+1, i-1, -one, w( i, 1 ), &
                         ldw, a( i, 1 ), lda, one, a( i, i ), 1 )
             CALL zlacgv( i-1, a( i, 1 ), lda )
             a( i, i ) = dble( a( i, i ) )
             IF( i.LT.n ) THEN
 !
 !              Generate elementary reflector H(i) to annihilate
 !              A(i+2:n,i)
 !
                alpha = a( i+1, i )
                CALL zlarfg( n-i, alpha, a( min( i+2, n ), i ), 1, &
                             tau( i ) )
                e( i ) = alpha
                a( i+1, i ) = one
 !
 !              Compute W(i+1:n,i)
 !
                CALL zhemv( 'Lower', n-i, one, a( i+1, i+1 ), lda, &
                            a( i+1, i ), 1, zero, w( i+1, i ), 1 )
                CALL zgemv( 'Conjugate transpose', n-i, i-1, one, &
                            w( i+1, 1 ), ldw, a( i+1, i ), 1, zero, &
                            w( 1, i ), 1 )
                CALL zgemv( 'No transpose', n-i, i-1, -one, a( i+1, 1 ), &
                            lda, w( 1, i ), 1, one, w( i+1, i ), 1 )
                CALL zgemv( 'Conjugate transpose', n-i, i-1, one, &
                            a( i+1, 1 ), lda, a( i+1, i ), 1, zero, &
                            w( 1, i ), 1 )
                CALL zgemv( 'No transpose', n-i, i-1, -one, w( i+1, 1 ), &
                            ldw, w( 1, i ), 1, one, w( i+1, i ), 1 )
                CALL zscal( n-i, tau( i ), w( i+1, i ), 1 )
                alpha = -half*tau( i )*zdotc( n-i, w( i+1, i ), 1, &
                        a( i+1, i ), 1 )
                CALL zaxpy( n-i, alpha, a( i+1, i ), 1, w( i+1, i ), 1 )
             END IF
 !
    20    CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZLATRD
 !
  end subroutine zlatrd

 !  =====================================================================
  subroutine zlascl( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          TYPE
       INTEGER            INFO, KL, KU, LDA, M, N
       DOUBLE PRECISION   CFROM, CTO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO, ONE
       parameter( zero = 0.0d0, one = 1.0d0 )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            DONE
       INTEGER            I, ITYPE, J, K1, K2, K3, K4
       DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME, DISNAN
 !      DOUBLE PRECISION   DLAMCH
 !      EXTERNAL           lsame, dlamch, disnan
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, max, min
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
 !
       IF( lsame( TYPE, 'G' ) ) then
          itype = 0
       ELSE IF( lsame( TYPE, 'L' ) ) then
          itype = 1
       ELSE IF( lsame( TYPE, 'U' ) ) then
          itype = 2
       ELSE IF( lsame( TYPE, 'H' ) ) then
          itype = 3
       ELSE IF( lsame( TYPE, 'B' ) ) then
          itype = 4
       ELSE IF( lsame( TYPE, 'Q' ) ) then
          itype = 5
       ELSE IF( lsame( TYPE, 'Z' ) ) then
          itype = 6
       ELSE
          itype = -1
       END IF
 !
       IF( itype.EQ.-1 ) THEN
          info = -1
       ELSE IF( cfrom.EQ.zero .OR. disnan(cfrom) ) THEN
          info = -4
       ELSE IF( disnan(cto) ) THEN
          info = -5
       ELSE IF( m.LT.0 ) THEN
          info = -6
       ELSE IF( n.LT.0 .OR. ( itype.EQ.4 .AND. n.NE.m ) .OR. &
                ( itype.EQ.5 .AND. n.NE.m ) ) THEN
          info = -7
       ELSE IF( itype.LE.3 .AND. lda.LT.max( 1, m ) ) THEN
          info = -9
       ELSE IF( itype.GE.4 ) THEN
          IF( kl.LT.0 .OR. kl.GT.max( m-1, 0 ) ) THEN
             info = -2
          ELSE IF( ku.LT.0 .OR. ku.GT.max( n-1, 0 ) .OR. &
                   ( ( itype.EQ.4 .OR. itype.EQ.5 ) .AND. kl.NE.ku ) ) &
                    THEN
             info = -3
          ELSE IF( ( itype.EQ.4 .AND. lda.LT.kl+1 ) .OR. &
                   ( itype.EQ.5 .AND. lda.LT.ku+1 ) .OR. &
                   ( itype.EQ.6 .AND. lda.LT.2*kl+ku+1 ) ) THEN
             info = -9
          END IF
       END IF
 !
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZLASCL', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 .OR. m.EQ.0 ) &
          RETURN
 !
 !     Get machine parameters
 !
       smlnum = dlamch( 'S' )
       bignum = one / smlnum
 !
       cfromc = cfrom
       ctoc = cto
 !
    10 CONTINUE
       cfrom1 = cfromc*smlnum
       IF( cfrom1.EQ.cfromc ) THEN
 !        CFROMC is an inf.  Multiply by a correctly signed zero for
 !        finite CTOC, or a NaN if CTOC is infinite.
          mul = ctoc / cfromc
          done = .true.
          cto1 = ctoc
       ELSE
          cto1 = ctoc / bignum
          IF( cto1.EQ.ctoc ) THEN
 !           CTOC is either 0 or an inf.  In both cases, CTOC itself
 !           serves as the correct multiplication factor.
             mul = ctoc
             done = .true.
             cfromc = one
          ELSE IF( abs( cfrom1 ).GT.abs( ctoc ) .AND. ctoc.NE.zero ) THEN
             mul = smlnum
             done = .false.
             cfromc = cfrom1
          ELSE IF( abs( cto1 ).GT.abs( cfromc ) ) THEN
             mul = bignum
             done = .false.
             ctoc = cto1
          ELSE
             mul = ctoc / cfromc
             done = .true.
          END IF
       END IF
 !
       IF( itype.EQ.0 ) THEN
 !
 !        Full matrix
 !
          DO 30 j = 1, n
             DO 20 i = 1, m
                a( i, j ) = a( i, j )*mul
    20       CONTINUE
    30    CONTINUE
 !
       ELSE IF( itype.EQ.1 ) THEN
 !
 !        Lower triangular matrix
 !
          DO 50 j = 1, n
             DO 40 i = j, m
                a( i, j ) = a( i, j )*mul
    40       CONTINUE
    50    CONTINUE
 !
       ELSE IF( itype.EQ.2 ) THEN
 !
 !        Upper triangular matrix
 !
          DO 70 j = 1, n
             DO 60 i = 1, min( j, m )
                a( i, j ) = a( i, j )*mul
    60       CONTINUE
    70    CONTINUE
 !
       ELSE IF( itype.EQ.3 ) THEN
 !
 !        Upper Hessenberg matrix
 !
          DO 90 j = 1, n
             DO 80 i = 1, min( j+1, m )
                a( i, j ) = a( i, j )*mul
    80       CONTINUE
    90    CONTINUE
 !
       ELSE IF( itype.EQ.4 ) THEN
 !
 !        Lower half of a symmetric band matrix
 !
          k3 = kl + 1
          k4 = n + 1
          DO 110 j = 1, n
             DO 100 i = 1, min( k3, k4-j )
                a( i, j ) = a( i, j )*mul
   100       CONTINUE
   110    CONTINUE
 !
       ELSE IF( itype.EQ.5 ) THEN
 !
 !        Upper half of a symmetric band matrix
 !
          k1 = ku + 2
          k3 = ku + 1
          DO 130 j = 1, n
             DO 120 i = max( k1-j, 1 ), k3
                a( i, j ) = a( i, j )*mul
   120       CONTINUE
   130    CONTINUE
 !
       ELSE IF( itype.EQ.6 ) THEN
 !
 !        Band matrix
 !
          k1 = kl + ku + 2
          k2 = kl + 1
          k3 = 2*kl + ku + 1
          k4 = kl + ku + 1 + m
          DO 150 j = 1, n
             DO 140 i = max( k1-j, k2 ), min( k3, k4-j )
                a( i, j ) = a( i, j )*mul
   140       CONTINUE
   150    CONTINUE
 !
       END IF
 !
       IF( .NOT.done ) &
          GO TO 10
 !
       RETURN
 !
 !     End of ZLASCL
 !
  end subroutine zlascl

 !  =====================================================================
  subroutine zsteqr( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER          COMPZ
       INTEGER            INFO, LDZ, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   D( * ), E( * ), WORK( * )
       COMPLEX*16         Z( ldz, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO, ONE, TWO, THREE
       parameter( zero = 0.0d0, one = 1.0d0, two = 2.0d0, three = 3.0d0 )
       COMPLEX*16         CZERO, CONE
       parameter( czero = ( 0.0d0, 0.0d0 ), cone = ( 1.0d0, 0.0d0 ) )
       INTEGER            MAXIT
       parameter( maxit = 30 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                          lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1, &
                          nm1, nmaxit
       DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                          s, safmax, safmin, ssfmax, ssfmin, tst
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
 !      EXTERNAL           lsame, dlamch, dlanst, dlapy2
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           dlae2, dlaev2, dlartg, dlascl, dlasrt, xerbla, &
 !                         zlaset, zlasr, zswap
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, max, sign, sqrt
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       info = 0
 !
       IF( lsame( compz, 'N' ) ) THEN
          icompz = 0
       ELSE IF( lsame( compz, 'V' ) ) THEN
          icompz = 1
       ELSE IF( lsame( compz, 'I' ) ) THEN
          icompz = 2
       ELSE
          icompz = -1
       END IF
       IF( icompz.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( ( ldz.LT.1 ) .OR. ( icompz.GT.0 .AND. ldz.LT.max( 1, &
                n ) ) ) THEN
          info = -6
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZSTEQR', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 ) &
          RETURN
 !
       IF( n.EQ.1 ) THEN
          IF( icompz.EQ.2 ) &
             z( 1, 1 ) = cone
          RETURN
       END IF
 !
 !     Determine the unit roundoff and over/underflow thresholds.
 !
       eps = dlamch( 'E' )
       eps2 = eps**2
       safmin = dlamch( 'S' )
       safmax = one / safmin
       ssfmax = sqrt( safmax ) / three
       ssfmin = sqrt( safmin ) / eps2
 !
 !     Compute the eigenvalues and eigenvectors of the tridiagonal
 !     matrix.
 !
       IF( icompz.EQ.2 ) &
          CALL zlaset( 'Full', n, n, czero, cone, z, ldz )
 !
       nmaxit = n*maxit
       jtot = 0
 !
 !     Determine where the matrix splits and choose QL or QR iteration
 !     for each block, according to whether top or bottom diagonal
 !     element is smaller.
 !
       l1 = 1
       nm1 = n - 1
 !
    10 CONTINUE
       IF( l1.GT.n ) &
          GO TO 160
       IF( l1.GT.1 ) &
          e( l1-1 ) = zero
       IF( l1.LE.nm1 ) THEN
          DO 20 m = l1, nm1
             tst = abs( e( m ) )
             IF( tst.EQ.zero ) &
                GO TO 30
             IF( tst.LE.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+ &
                 1 ) ) ) )*eps ) THEN
                e( m ) = zero
                GO TO 30
             END IF
    20    CONTINUE
       END IF
       m = n
 !
    30 CONTINUE
       l = l1
       lsv = l
       lend = m
       lendsv = lend
       l1 = m + 1
       IF( lend.EQ.l ) &
          GO TO 10
 !
 !     Scale submatrix in rows and columns L to LEND
 !
       anorm = dlanst( 'I', lend-l+1, d( l ), e( l ) )
       iscale = 0
       IF( anorm.EQ.zero ) &
          GO TO 10
       IF( anorm.GT.ssfmax ) THEN
          iscale = 1
          CALL dlascl( 'G', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n, &
                       info )
          CALL dlascl( 'G', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n, &
                       info )
       ELSE IF( anorm.LT.ssfmin ) THEN
          iscale = 2
          CALL dlascl( 'G', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n, &
                       info )
          CALL dlascl( 'G', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n, &
                       info )
       END IF
 !
 !     Choose between QL and QR iteration
 !
       IF( abs( d( lend ) ).LT.abs( d( l ) ) ) THEN
          lend = lsv
          l = lendsv
       END IF
 !
       IF( lend.GT.l ) THEN
 !
 !        QL Iteration
 !
 !        Look for small subdiagonal element.
 !
    40    CONTINUE
          IF( l.NE.lend ) THEN
             lendm1 = lend - 1
             DO 50 m = l, lendm1
                tst = abs( e( m ) )**2
                IF( tst.LE.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+ &
                    safmin )GO TO 60
    50       CONTINUE
          END IF
 !
          m = lend
 !
    60    CONTINUE
          IF( m.LT.lend ) &
             e( m ) = zero
          p = d( l )
          IF( m.EQ.l ) &
             GO TO 80
 !
 !        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
 !        to compute its eigensystem.
 !
          IF( m.EQ.l+1 ) THEN
             IF( icompz.GT.0 ) THEN
                CALL dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
                work( l ) = c
                work( n-1+l ) = s
                CALL zlasr( 'R', 'V', 'B', n, 2, work( l ), &
                            work( n-1+l ), z( 1, l ), ldz )
             ELSE
                CALL dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
             END IF
             d( l ) = rt1
             d( l+1 ) = rt2
             e( l ) = zero
             l = l + 2
             IF( l.LE.lend ) &
                GO TO 40
             GO TO 140
          END IF
 !
          IF( jtot.EQ.nmaxit ) &
             GO TO 140
          jtot = jtot + 1
 !
 !        Form shift.
 !
          g = ( d( l+1 )-p ) / ( two*e( l ) )
          r = dlapy2( g, one )
          g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
 !
          s = one
          c = one
          p = zero
 !
 !        Inner loop
 !
          mm1 = m - 1
          DO 70 i = mm1, l, -1
             f = s*e( i )
             b = c*e( i )
             CALL dlartg( g, f, c, s, r )
             IF( i.NE.m-1 ) &
                e( i+1 ) = r
             g = d( i+1 ) - p
             r = ( d( i )-g )*s + two*c*b
             p = s*r
             d( i+1 ) = g + p
             g = c*r - b
 !
 !           If eigenvectors are desired, then save rotations.
 !
             IF( icompz.GT.0 ) THEN
                work( i ) = c
                work( n-1+i ) = -s
             END IF
 !
    70    CONTINUE
 !
 !        If eigenvectors are desired, then apply saved rotations.
 !
          IF( icompz.GT.0 ) THEN
             mm = m - l + 1
             CALL zlasr( 'R', 'V', 'B', n, mm, work( l ), work( n-1+l ), &
                         z( 1, l ), ldz )
          END IF
 !
          d( l ) = d( l ) - p
          e( l ) = g
          GO TO 40
 !
 !        Eigenvalue found.
 !
    80    CONTINUE
          d( l ) = p
 !
          l = l + 1
          IF( l.LE.lend ) &
             GO TO 40
          GO TO 140
 !
       ELSE
 !
 !        QR Iteration
 !
 !        Look for small superdiagonal element.
 !
    90    CONTINUE
          IF( l.NE.lend ) THEN
             lendp1 = lend + 1
             DO 100 m = l, lendp1, -1
                tst = abs( e( m-1 ) )**2
                IF( tst.LE.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+ &
                    safmin )GO TO 110
   100       CONTINUE
          END IF
 !
          m = lend
 !
   110    CONTINUE
          IF( m.GT.lend ) &
             e( m-1 ) = zero
          p = d( l )
          IF( m.EQ.l ) &
             GO TO 130
 !
 !        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
 !        to compute its eigensystem.
 !
          IF( m.EQ.l-1 ) THEN
             IF( icompz.GT.0 ) THEN
                CALL dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
                work( m ) = c
                work( n-1+m ) = s
                CALL zlasr( 'R', 'V', 'F', n, 2, work( m ), &
                            work( n-1+m ), z( 1, l-1 ), ldz )
             ELSE
                CALL dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
             END IF
             d( l-1 ) = rt1
             d( l ) = rt2
             e( l-1 ) = zero
             l = l - 2
             IF( l.GE.lend ) &
                GO TO 90
             GO TO 140
          END IF
 !
          IF( jtot.EQ.nmaxit ) &
             GO TO 140
          jtot = jtot + 1
 !
 !        Form shift.
 !
          g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
          r = dlapy2( g, one )
          g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
 !
          s = one
          c = one
          p = zero
 !
 !        Inner loop
 !
          lm1 = l - 1
          DO 120 i = m, lm1
             f = s*e( i )
             b = c*e( i )
             CALL dlartg( g, f, c, s, r )
             IF( i.NE.m ) &
                e( i-1 ) = r
             g = d( i ) - p
             r = ( d( i+1 )-g )*s + two*c*b
             p = s*r
             d( i ) = g + p
             g = c*r - b
 !
 !           If eigenvectors are desired, then save rotations.
 !
             IF( icompz.GT.0 ) THEN
                work( i ) = c
                work( n-1+i ) = s
             END IF
 !
   120    CONTINUE
 !
 !        If eigenvectors are desired, then apply saved rotations.
 !
          IF( icompz.GT.0 ) THEN
             mm = l - m + 1
             CALL zlasr( 'R', 'V', 'F', n, mm, work( m ), work( n-1+m ), &
                         z( 1, m ), ldz )
          END IF
 !
          d( l ) = d( l ) - p
          e( lm1 ) = g
          GO TO 90
 !
 !        Eigenvalue found.
 !
   130    CONTINUE
          d( l ) = p
 !
          l = l - 1
          IF( l.GE.lend ) &
             GO TO 90
          GO TO 140
 !
       END IF
 !
 !     Undo scaling if necessary
 !
   140 CONTINUE
       IF( iscale.EQ.1 ) THEN
          CALL dlascl( 'G', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1,  &
                       d( lsv ), n, info )
          CALL dlascl( 'G', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ), &
                      n, info )
       ELSE IF( iscale.EQ.2 ) THEN
          CALL dlascl( 'G', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1, &
                       d( lsv ), n, info )
          CALL dlascl( 'G', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ), &
                       n, info )
       END IF
 !
 !     Check for no convergence to an eigenvalue after a total
 !     of N*MAXIT iterations.
 !
       IF( jtot.EQ.nmaxit ) THEN
          DO 150 i = 1, n - 1
             IF( e( i ).NE.zero ) &
                info = info + 1
   150    CONTINUE
          RETURN
       END IF
       GO TO 10
 !
 !     Order eigenvalues and eigenvectors.
 !
   160 CONTINUE
       IF( icompz.EQ.0 ) THEN
 !
 !        Use Quick Sort
 !
          CALL dlasrt( 'I', n, d, info )
 !
       ELSE
 !
 !        Use Selection Sort to minimize swaps of eigenvectors
 !
          DO 180 ii = 2, n
             i = ii - 1
             k = i
             p = d( i )
             DO 170 j = ii, n
                IF( d( j ).LT.p ) THEN
                   k = j
                   p = d( j )
                END IF
   170       CONTINUE
             IF( k.NE.i ) THEN
                d( k ) = d( i )
                d( i ) = p
                CALL zswap( n, z( 1, i ), 1, z( 1, k ), 1 )
             END IF
   180    CONTINUE
       END IF
       RETURN
 !
 !     End of ZSTEQR
 !
  end subroutine zsteqr


 !  =====================================================================
  subroutine dlaev2( A, B, C, RT1, RT2, CS1, SN1 )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   A, B, C, CS1, RT1, RT2, SN1
 !     ..
 !
 ! =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d0 )
       DOUBLE PRECISION   TWO
       parameter( two = 2.0d0 )
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
       DOUBLE PRECISION   HALF
       parameter( half = 0.5d0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            SGN1, SGN2
       DOUBLE PRECISION   AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, tb, tn
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, sqrt
 !     ..
 !     .. Executable Statements ..
 !
 !     Compute the eigenvalues
 !
       sm = a + c
       df = a - c
       adf = abs( df )
       tb = b + b
       ab = abs( tb )
       IF( abs( a ).GT.abs( c ) ) THEN
          acmx = a
          acmn = c
       ELSE
          acmx = c
          acmn = a
       END IF
       IF( adf.GT.ab ) THEN
          rt = adf*sqrt( one+( ab / adf )**2 )
       ELSE IF( adf.LT.ab ) THEN
          rt = ab*sqrt( one+( adf / ab )**2 )
       ELSE
 !
 !        Includes case AB=ADF=0
 !
          rt = ab*sqrt( two )
       END IF
       IF( sm.LT.zero ) THEN
          rt1 = half*( sm-rt )
          sgn1 = -1
 !
 !        Order of execution important.
 !        To get fully accurate smaller eigenvalue,
 !        next line needs to be executed in higher precision.
 !
          rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
       ELSE IF( sm.GT.zero ) THEN
          rt1 = half*( sm+rt )
          sgn1 = 1
 !
 !        Order of execution important.
 !        To get fully accurate smaller eigenvalue,
 !        next line needs to be executed in higher precision.
 !
          rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
       ELSE
 !
 !        Includes case RT1 = RT2 = 0
 !
          rt1 = half*rt
          rt2 = -half*rt
          sgn1 = 1
       END IF
 !
 !     Compute the eigenvector
 !
       IF( df.GE.zero ) THEN
          cs = df + rt
          sgn2 = 1
       ELSE
          cs = df - rt
          sgn2 = -1
       END IF
       acs = abs( cs )
       IF( acs.GT.ab ) THEN
          ct = -tb / cs
          sn1 = one / sqrt( one+ct*ct )
          cs1 = ct*sn1
       ELSE
          IF( ab.EQ.zero ) THEN
             cs1 = one
             sn1 = zero
          ELSE
             tn = -cs / tb
             cs1 = one / sqrt( one+tn*tn )
             sn1 = tn*cs1
          END IF
       END IF
       IF( sgn1.EQ.sgn2 ) THEN
          tn = cs1
          cs1 = -sn1
          sn1 = tn
       END IF
       RETURN
 !
 !     End of DLAEV2
 !
  end subroutine dlaev2

 !  =====================================================================
  subroutine dlartg( F, G, CS, SN, R )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       DOUBLE PRECISION   CS, F, G, R, SN
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ZERO
       parameter( zero = 0.0d0 )
       DOUBLE PRECISION   ONE
       parameter( one = 1.0d0 )
       DOUBLE PRECISION   TWO
       parameter( two = 2.0d0 )
 !     ..
 !     .. Local Scalars ..
 !     LOGICAL            FIRST
       INTEGER            COUNT, I
       DOUBLE PRECISION   EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
 !     ..
 !     .. External Functions ..
 !      DOUBLE PRECISION   DLAMCH
 !      EXTERNAL           dlamch
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          abs, int, log, max, sqrt
 !     ..
 !     .. Save statement ..
 !     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
 !     ..
 !     .. Data statements ..
 !     DATA               FIRST / .TRUE. /
 !     ..
 !     .. Executable Statements ..
 !
 !     IF( FIRST ) THEN
          safmin = dlamch( 'S' )
          eps = dlamch( 'E' )
          safmn2 = dlamch( 'B' )**int( log( safmin / eps ) / &
                   log( dlamch( 'B' ) ) / two )
          safmx2 = one / safmn2
 !        FIRST = .FALSE.
 !     END IF
       IF( g.EQ.zero ) THEN
          cs = one
          sn = zero
          r = f
       ELSE IF( f.EQ.zero ) THEN
          cs = zero
          sn = one
          r = g
       ELSE
          f1 = f
          g1 = g
          scale = max( abs( f1 ), abs( g1 ) )
          IF( scale.GE.safmx2 ) THEN
             count = 0
    10       CONTINUE
             count = count + 1
             f1 = f1*safmn2
             g1 = g1*safmn2
             scale = max( abs( f1 ), abs( g1 ) )
             IF( scale.GE.safmx2 ) &
                GO TO 10
             r = sqrt( f1**2+g1**2 )
             cs = f1 / r
             sn = g1 / r
             DO 20 i = 1, count
                r = r*safmx2
    20       CONTINUE
          ELSE IF( scale.LE.safmn2 ) THEN
             count = 0
    30       CONTINUE
             count = count + 1
             f1 = f1*safmx2
             g1 = g1*safmx2
             scale = max( abs( f1 ), abs( g1 ) )
             IF( scale.LE.safmn2 ) &
                GO TO 30
             r = sqrt( f1**2+g1**2 )
             cs = f1 / r
             sn = g1 / r
             DO 40 i = 1, count
                r = r*safmn2
    40       CONTINUE
          ELSE
             r = sqrt( f1**2+g1**2 )
             cs = f1 / r
             sn = g1 / r
          END IF
          IF( abs( f ).GT.abs( g ) .AND. cs.LT.zero ) THEN
             cs = -cs
             sn = -sn
             r = -r
          END IF
       END IF
       RETURN
 !
 !     End of DLARTG
 !
  end subroutine dlartg

 !  =====================================================================
  subroutine zlaset( UPLO, M, N, ALPHA, BETA, A, LDA )
 !
 !  -- LAPACK auxiliary routine (version 3.6.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2015
 !
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            LDA, M, N
       COMPLEX*16         ALPHA, BETA
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER            I, J
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      EXTERNAL           lsame
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          min
 !     ..
 !     .. Executable Statements ..
 !
       IF( lsame( uplo, 'U' ) ) THEN
 !
 !        Set the diagonal to BETA and the strictly upper triangular
 !        part of the array to ALPHA.
 !
          DO 20 j = 2, n
             DO 10 i = 1, min( j-1, m )
                a( i, j ) = alpha
    10       CONTINUE
    20    CONTINUE
          DO 30 i = 1, min( n, m )
             a( i, i ) = beta
    30    CONTINUE
 !
       ELSE IF( lsame( uplo, 'L' ) ) THEN
 !
 !        Set the diagonal to BETA and the strictly lower triangular
 !        part of the array to ALPHA.
 !
          DO 50 j = 1, min( m, n )
             DO 40 i = j + 1, m
                a( i, j ) = alpha
    40       CONTINUE
    50    CONTINUE
          DO 60 i = 1, min( n, m )
             a( i, i ) = beta
    60    CONTINUE
 !
       ELSE
 !
 !        Set the array to BETA on the diagonal and ALPHA on the
 !        offdiagonal.
 !
          DO 80 j = 1, n
             DO 70 i = 1, m
                a( i, j ) = alpha
    70       CONTINUE
    80    CONTINUE
          DO 90 i = 1, min( m, n )
             a( i, i ) = beta
    90    CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZLASET
 !
  end subroutine zlaset

 !  =====================================================================
  subroutine zlasr( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          DIRECT, PIVOT, SIDE
       INTEGER            LDA, M, N
 !     ..
 !     .. Array Arguments ..
       DOUBLE PRECISION   C( * ), S( * )
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       DOUBLE PRECISION   ONE, ZERO
       parameter( one = 1.0d0, zero = 0.0d0 )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, INFO, J
       DOUBLE PRECISION   CTEMP, STEMP
       COMPLEX*16         TEMP
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      EXTERNAL           lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters
 !
       info = 0
       IF( .NOT.( lsame( side, 'L' ) .OR. lsame( side, 'R' ) ) ) THEN
          info = 1
       ELSE IF( .NOT.( lsame( pivot, 'V' ) .OR. lsame( pivot, &
                'T' ) .OR. lsame( pivot, 'B' ) ) ) THEN
          info = 2
       ELSE IF( .NOT.( lsame( direct, 'F' ) .OR. lsame( direct, 'B' ) ) ) &
                 THEN
          info = 3
       ELSE IF( m.LT.0 ) THEN
          info = 4
       ELSE IF( n.LT.0 ) THEN
          info = 5
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = 9
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZLASR ', info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( ( m.EQ.0 ) .OR. ( n.EQ.0 ) ) &
          RETURN
       IF( lsame( side, 'L' ) ) THEN
 !
 !        Form  P * A
 !
          IF( lsame( pivot, 'V' ) ) THEN
             IF( lsame( direct, 'F' ) ) THEN
                DO 20 j = 1, m - 1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 10 i = 1, n
                         temp = a( j+1, i )
                         a( j+1, i ) = ctemp*temp - stemp*a( j, i )
                         a( j, i ) = stemp*temp + ctemp*a( j, i )
    10                CONTINUE
                   END IF
    20          CONTINUE
             ELSE IF( lsame( direct, 'B' ) ) THEN
                DO 40 j = m - 1, 1, -1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 30 i = 1, n
                         temp = a( j+1, i )
                         a( j+1, i ) = ctemp*temp - stemp*a( j, i )
                         a( j, i ) = stemp*temp + ctemp*a( j, i )
    30                CONTINUE
                   END IF
    40          CONTINUE
             END IF
          ELSE IF( lsame( pivot, 'T' ) ) THEN
             IF( lsame( direct, 'F' ) ) THEN
                DO 60 j = 2, m
                   ctemp = c( j-1 )
                   stemp = s( j-1 )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 50 i = 1, n
                         temp = a( j, i )
                         a( j, i ) = ctemp*temp - stemp*a( 1, i )
                         a( 1, i ) = stemp*temp + ctemp*a( 1, i )
    50                CONTINUE
                   END IF
    60          CONTINUE
             ELSE IF( lsame( direct, 'B' ) ) THEN
                DO 80 j = m, 2, -1
                   ctemp = c( j-1 )
                   stemp = s( j-1 )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 70 i = 1, n
                         temp = a( j, i )
                         a( j, i ) = ctemp*temp - stemp*a( 1, i )
                         a( 1, i ) = stemp*temp + ctemp*a( 1, i )
    70                CONTINUE
                   END IF
    80          CONTINUE
             END IF
          ELSE IF( lsame( pivot, 'B' ) ) THEN
             IF( lsame( direct, 'F' ) ) THEN
                DO 100 j = 1, m - 1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 90 i = 1, n
                         temp = a( j, i )
                         a( j, i ) = stemp*a( m, i ) + ctemp*temp
                         a( m, i ) = ctemp*a( m, i ) - stemp*temp
    90                CONTINUE
                   END IF
   100          CONTINUE
             ELSE IF( lsame( direct, 'B' ) ) THEN
                DO 120 j = m - 1, 1, -1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 110 i = 1, n
                         temp = a( j, i )
                         a( j, i ) = stemp*a( m, i ) + ctemp*temp
                         a( m, i ) = ctemp*a( m, i ) - stemp*temp
   110                CONTINUE
                   END IF
   120          CONTINUE
             END IF
          END IF
       ELSE IF( lsame( side, 'R' ) ) THEN
 !
 !        Form A * P**T
 !
          IF( lsame( pivot, 'V' ) ) THEN
             IF( lsame( direct, 'F' ) ) THEN
                DO 140 j = 1, n - 1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 130 i = 1, m
                         temp = a( i, j+1 )
                         a( i, j+1 ) = ctemp*temp - stemp*a( i, j )
                         a( i, j ) = stemp*temp + ctemp*a( i, j )
   130                CONTINUE
                   END IF
   140          CONTINUE
             ELSE IF( lsame( direct, 'B' ) ) THEN
                DO 160 j = n - 1, 1, -1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 150 i = 1, m
                         temp = a( i, j+1 )
                         a( i, j+1 ) = ctemp*temp - stemp*a( i, j )
                         a( i, j ) = stemp*temp + ctemp*a( i, j )
   150                CONTINUE
                   END IF
   160          CONTINUE
             END IF
          ELSE IF( lsame( pivot, 'T' ) ) THEN
             IF( lsame( direct, 'F' ) ) THEN
                DO 180 j = 2, n
                   ctemp = c( j-1 )
                   stemp = s( j-1 )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 170 i = 1, m
                         temp = a( i, j )
                         a( i, j ) = ctemp*temp - stemp*a( i, 1 )
                         a( i, 1 ) = stemp*temp + ctemp*a( i, 1 )
   170                CONTINUE
                   END IF
   180          CONTINUE
             ELSE IF( lsame( direct, 'B' ) ) THEN
                DO 200 j = n, 2, -1
                   ctemp = c( j-1 )
                   stemp = s( j-1 )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 190 i = 1, m
                         temp = a( i, j )
                         a( i, j ) = ctemp*temp - stemp*a( i, 1 )
                         a( i, 1 ) = stemp*temp + ctemp*a( i, 1 )
   190                CONTINUE
                   END IF
   200          CONTINUE
             END IF
          ELSE IF( lsame( pivot, 'B' ) ) THEN
             IF( lsame( direct, 'F' ) ) THEN
                DO 220 j = 1, n - 1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 210 i = 1, m
                         temp = a( i, j )
                         a( i, j ) = stemp*a( i, n ) + ctemp*temp
                         a( i, n ) = ctemp*a( i, n ) - stemp*temp
   210                CONTINUE
                   END IF
   220          CONTINUE
             ELSE IF( lsame( direct, 'B' ) ) THEN
                DO 240 j = n - 1, 1, -1
                   ctemp = c( j )
                   stemp = s( j )
                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
                      DO 230 i = 1, m
                         temp = a( i, j )
                         a( i, j ) = stemp*a( i, n ) + ctemp*temp
                         a( i, n ) = ctemp*a( i, n ) - stemp*temp
   230                CONTINUE
                   END IF
   240          CONTINUE
             END IF
          END IF
       END IF
 !
       RETURN
 !
 !     End of ZLASR
 !
  end subroutine zlasr

 !  =====================================================================
  subroutine zswap(N,ZX,INCX,ZY,INCY)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,INCY,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*),ZY(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       COMPLEX*16 ZTEMP
       INTEGER I,IX,IY
 !     ..
       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 !
 !       code for both increments equal to 1
          DO i = 1,n
             ztemp = zx(i)
             zx(i) = zy(i)
             zy(i) = ztemp
          END DO
       ELSE
 !
 !       code for unequal increments or equal increments not equal
 !         to 1
 !
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             ztemp = zx(ix)
             zx(ix) = zy(iy)
             zy(iy) = ztemp
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       RETURN
  end subroutine zswap

 !  =====================================================================
  subroutine zungtr( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER          UPLO
       INTEGER            INFO, LDA, LWORK, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ZERO, ONE
       parameter( zero = ( 0.0d0, 0.0d0 ), &
                          one = ( 1.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            LQUERY, UPPER
       INTEGER            I, IINFO, J, LWKOPT, NB
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      INTEGER            ILAENV
 !      EXTERNAL           lsame, ilaenv
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zungql, zungqr
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
       lquery = ( lwork.EQ.-1 )
       upper = lsame( uplo, 'U' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -4
       ELSE IF( lwork.LT.max( 1, n-1 ) .AND. .NOT.lquery ) THEN
          info = -7
       END IF
 !
       IF( info.EQ.0 ) THEN
          IF( upper ) THEN
             nb = ilaenv( 1, 'ZUNGQL', ' ', n-1, n-1, n-1, -1 )
          ELSE
             nb = ilaenv( 1, 'ZUNGQR', ' ', n-1, n-1, n-1, -1 )
          END IF
          lwkopt = max( 1, n-1 )*nb
          work( 1 ) = lwkopt
       END IF
 !
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZUNGTR', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 ) THEN
          work( 1 ) = 1
          RETURN
       END IF
 !
       IF( upper ) THEN
 !
 !        Q was determined by a call to ZHETRD with UPLO = 'U'
 !
 !        Shift the vectors which define the elementary reflectors one
 !        column to the left, and set the last row and column of Q to
 !        those of the unit matrix
 !
          DO 20 j = 1, n - 1
             DO 10 i = 1, j - 1
                a( i, j ) = a( i, j+1 )
    10       CONTINUE
             a( n, j ) = zero
    20    CONTINUE
          DO 30 i = 1, n - 1
             a( i, n ) = zero
    30    CONTINUE
          a( n, n ) = one
 !
 !        Generate Q(1:n-1,1:n-1)
 !
          CALL zungql( n-1, n-1, n-1, a, lda, tau, work, lwork, iinfo )
 !
       ELSE
 !
 !        Q was determined by a call to ZHETRD with UPLO = 'L'.
 !
 !        Shift the vectors which define the elementary reflectors one
 !        column to the right, and set the first row and column of Q to
 !        those of the unit matrix
 !
          DO 50 j = n, 2, -1
             a( 1, j ) = zero
             DO 40 i = j + 1, n
                a( i, j ) = a( i, j-1 )
    40       CONTINUE
    50    CONTINUE
          a( 1, 1 ) = one
          DO 60 i = 2, n
             a( i, 1 ) = zero
    60    CONTINUE
          IF( n.GT.1 ) THEN
 !
 !           Generate Q(2:n,2:n)
 !
             CALL zungqr( n-1, n-1, n-1, a( 2, 2 ), lda, tau, work, &
                          lwork, iinfo )
          END IF
       END IF
       work( 1 ) = lwkopt
       RETURN
 !
 !     End of ZUNGTR
 !
  end subroutine zungtr


 !  =====================================================================
  subroutine zungql( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, K, LDA, LWORK, M, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ZERO
       parameter( zero = ( 0.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            LQUERY
       INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT, &
                          nb, nbmin, nx
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zlarfb, zlarft, zung2l
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. External Functions ..
 !      INTEGER            ILAENV
 !      EXTERNAL           ilaenv
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
       lquery = ( lwork.EQ.-1 )
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
          info = -2
       ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -5
       END IF
 !
       IF( info.EQ.0 ) THEN
          IF( n.EQ.0 ) THEN
             lwkopt = 1
          ELSE
             nb = ilaenv( 1, 'ZUNGQL', ' ', m, n, k, -1 )
             lwkopt = n*nb
          END IF
          work( 1 ) = lwkopt
 !
          IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
             info = -8
          END IF
       END IF
 !
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZUNGQL', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.LE.0 ) THEN
          RETURN
       END IF
 !
       nbmin = 2
       nx = 0
       iws = n
       IF( nb.GT.1 .AND. nb.LT.k ) THEN
 !
 !        Determine when to cross over from blocked to unblocked code.
 !
          nx = max( 0, ilaenv( 3, 'ZUNGQL', ' ', m, n, k, -1 ) )
          IF( nx.LT.k ) THEN
 !
 !           Determine if workspace is large enough for blocked code.
 !
             ldwork = n
             iws = ldwork*nb
             IF( lwork.LT.iws ) THEN
 !
 !              Not enough workspace to use optimal NB:  reduce NB and
 !              determine the minimum value of NB.
 !
                nb = lwork / ldwork
                nbmin = max( 2, ilaenv( 2, 'ZUNGQL', ' ', m, n, k, -1 ) )
             END IF
          END IF
       END IF
 !
       IF( nb.GE.nbmin .AND. nb.LT.k .AND. nx.LT.k ) THEN
 !
 !        Use blocked code after the first block.
 !        The last kk columns are handled by the block method.
 !
          kk = min( k, ( ( k-nx+nb-1 ) / nb )*nb )
 !
 !        Set A(m-kk+1:m,1:n-kk) to zero.
 !
          DO 20 j = 1, n - kk
             DO 10 i = m - kk + 1, m
                a( i, j ) = zero
    10       CONTINUE
    20    CONTINUE
       ELSE
          kk = 0
       END IF
 !
 !     Use unblocked code for the first or only block.
 !
       CALL zung2l( m-kk, n-kk, k-kk, a, lda, tau, work, iinfo )
 !
       IF( kk.GT.0 ) THEN
 !
 !        Use blocked code
 !
          DO 50 i = k - kk + 1, k, nb
             ib = min( nb, k-i+1 )
             IF( n-k+i.GT.1 ) THEN
 !
 !              Form the triangular factor of the block reflector
 !              H = H(i+ib-1) . . . H(i+1) H(i)
 !
                CALL zlarft( 'Backward', 'Columnwise', m-k+i+ib-1, ib, &
                             a( 1, n-k+i ), lda, tau( i ), work, ldwork )
 !
 !              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
 !
                CALL zlarfb( 'Left', 'No transpose', 'Backward', &
                             'Columnwise', m-k+i+ib-1, n-k+i-1, ib, &
                             a( 1, n-k+i ), lda, work, ldwork, a, lda, &
                             work( ib+1 ), ldwork )
             END IF
 !
 !           Apply H to rows 1:m-k+i+ib-1 of current block
 !
             CALL zung2l( m-k+i+ib-1, ib, ib, a( 1, n-k+i ), lda, &
                          tau( i ), work, iinfo )
 !
 !           Set rows m-k+i+ib:m of current block to zero
 !
             DO 40 j = n - k + i, n - k + i + ib - 1
                DO 30 l = m - k + i + ib, m
                   a( l, j ) = zero
    30          CONTINUE
    40       CONTINUE
    50    CONTINUE
       END IF
 !
       work( 1 ) = iws
       RETURN
 !
 !     End of ZUNGQL
 !
  end subroutine zungql

 !  =====================================================================
  subroutine zlarfb( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, &
                          t, ldt, c, ldc, work, ldwork )
 !
 !  -- LAPACK auxiliary routine (version 3.5.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     June 2013
 !
 !     .. Scalar Arguments ..
       CHARACTER          DIRECT, SIDE, STOREV, TRANS
       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         C( ldc, * ), T( ldt, * ), V( ldv, * ), &
                          work( ldwork, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE
       parameter( one = ( 1.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       CHARACTER          TRANST
       INTEGER            I, J
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      EXTERNAL           lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           zcopy, zgemm, zlacgv, ztrmm
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          dconjg
 !     ..
 !     .. Executable Statements ..
 !
 !     Quick return if possible
 !
       IF( m.LE.0 .OR. n.LE.0 ) &
          RETURN
 !
       IF( lsame( trans, 'N' ) ) THEN
          transt = 'C'
       ELSE
          transt = 'N'
       END IF
 !
       IF( lsame( storev, 'C' ) ) THEN
 !
          IF( lsame( direct, 'F' ) ) THEN
 !
 !           Let  V =  ( V1 )    (first K rows)
 !                     ( V2 )
 !           where  V1  is unit lower triangular.
 !
             IF( lsame( side, 'L' ) ) THEN
 !
 !              Form  H * C  or  H**H * C  where  C = ( C1 )
 !                                                    ( C2 )
 !
 !              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
 !
 !              W := C1**H
 !
                DO 10 j = 1, k
                   CALL zcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
                   CALL zlacgv( n, work( 1, j ), 1 )
    10          CONTINUE
 !
 !              W := W * V1
 !
                CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', n, &
                            k, one, v, ldv, work, ldwork )
                IF( m.GT.k ) THEN
 !
 !                 W := W + C2**H * V2
 !
                   CALL zgemm( 'Conjugate transpose', 'No transpose', n, &
                               k, m-k, one, c( k+1, 1 ), ldc, &
                               v( k+1, 1 ), ldv, one, work, ldwork )
                END IF
 !
 !              W := W * T**H  or  W * T
 !
                CALL ztrmm( 'Right', 'Upper', transt, 'Non-unit', n, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - V * W**H
 !
                IF( m.GT.k ) THEN
 !
 !                 C2 := C2 - V2 * W**H
 !
                   CALL zgemm( 'No transpose', 'Conjugate transpose', &
                               m-k, n, k, -one, v( k+1, 1 ), ldv, work, &
                               ldwork, one, c( k+1, 1 ), ldc )
                END IF
 !
 !              W := W * V1**H
 !
                CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose', &
                            'Unit', n, k, one, v, ldv, work, ldwork )
 !
 !              C1 := C1 - W**H
 !
                DO 30 j = 1, k
                   DO 20 i = 1, n
                      c( j, i ) = c( j, i ) - dconjg( work( i, j ) )
    20             CONTINUE
    30          CONTINUE
 !
             ELSE IF( lsame( side, 'R' ) ) THEN
 !
 !              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
 !
 !              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
 !
 !              W := C1
 !
                DO 40 j = 1, k
                   CALL zcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
    40          CONTINUE
 !
 !              W := W * V1
 !
                CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', m, &
                            k, one, v, ldv, work, ldwork )
                IF( n.GT.k ) THEN
 !
 !                 W := W + C2 * V2
 !
                   CALL zgemm( 'No transpose', 'No transpose', m, k, n-k, &
                               one, c( 1, k+1 ), ldc, v( k+1, 1 ), ldv, &
                               one, work, ldwork )
                END IF
 !
 !              W := W * T  or  W * T**H
 !
                CALL ztrmm( 'Right', 'Upper', trans, 'Non-unit', m, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - W * V**H
 !
                IF( n.GT.k ) THEN
 !
 !                 C2 := C2 - W * V2**H
 !
                   CALL zgemm( 'No transpose', 'Conjugate transpose', m, &
                               n-k, k, -one, work, ldwork, v( k+1, 1 ), &
                               ldv, one, c( 1, k+1 ), ldc )
                END IF
 !
 !              W := W * V1**H
 !
                CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose', &
                            'Unit', m, k, one, v, ldv, work, ldwork )
 !
 !              C1 := C1 - W
 !
                DO 60 j = 1, k
                   DO 50 i = 1, m
                      c( i, j ) = c( i, j ) - work( i, j )
    50             CONTINUE
    60          CONTINUE
             END IF
 !
          ELSE
 !
 !           Let  V =  ( V1 )
 !                     ( V2 )    (last K rows)
 !           where  V2  is unit upper triangular.
 !
             IF( lsame( side, 'L' ) ) THEN
 !
 !              Form  H * C  or  H**H * C  where  C = ( C1 )
 !                                                    ( C2 )
 !
 !              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
 !
 !              W := C2**H
 !
                DO 70 j = 1, k
                   CALL zcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1 )
                   CALL zlacgv( n, work( 1, j ), 1 )
    70          CONTINUE
 !
 !              W := W * V2
 !
                CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', n, &
                            k, one, v( m-k+1, 1 ), ldv, work, ldwork )
                IF( m.GT.k ) THEN
 !
 !                 W := W + C1**H * V1
 !
                   CALL zgemm( 'Conjugate transpose', 'No transpose', n, &
                               k, m-k, one, c, ldc, v, ldv, one, work, &
                               ldwork )
                END IF
 !
 !              W := W * T**H  or  W * T
 !
                CALL ztrmm( 'Right', 'Lower', transt, 'Non-unit', n, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - V * W**H
 !
                IF( m.GT.k ) THEN
 !
 !                 C1 := C1 - V1 * W**H
 !
                   CALL zgemm( 'No transpose', 'Conjugate transpose', &
                               m-k, n, k, -one, v, ldv, work, ldwork, &
                               one, c, ldc )
                END IF
 !
 !              W := W * V2**H
 !
                CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose', &
                            'Unit', n, k, one, v( m-k+1, 1 ), ldv, work, &
                            ldwork )
 !
 !              C2 := C2 - W**H
 !
                DO 90 j = 1, k
                   DO 80 i = 1, n
                      c( m-k+j, i ) = c( m-k+j, i ) - &
                                      dconjg( work( i, j ) )
    80             CONTINUE
    90          CONTINUE
 !
             ELSE IF( lsame( side, 'R' ) ) THEN
 !
 !              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
 !
 !              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
 !
 !              W := C2
 !
                DO 100 j = 1, k
                   CALL zcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
   100          CONTINUE
 !
 !              W := W * V2
 !
                CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', m, &
                            k, one, v( n-k+1, 1 ), ldv, work, ldwork )
                IF( n.GT.k ) THEN
 !
 !                 W := W + C1 * V1
 !
                   CALL zgemm( 'No transpose', 'No transpose', m, k, n-k, &
                               one, c, ldc, v, ldv, one, work, ldwork )
                END IF
 !
 !              W := W * T  or  W * T**H
 !
                CALL ztrmm( 'Right', 'Lower', trans, 'Non-unit', m, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - W * V**H
 !
                IF( n.GT.k ) THEN
 !
 !                 C1 := C1 - W * V1**H
 !
                   CALL zgemm( 'No transpose', 'Conjugate transpose', m, &
                               n-k, k, -one, work, ldwork, v, ldv, one, &
                               c, ldc )
                END IF
 !
 !              W := W * V2**H
 !
                CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose', &
                            'Unit', m, k, one, v( n-k+1, 1 ), ldv, work, &
                            ldwork )
 !
 !              C2 := C2 - W
 !
                DO 120 j = 1, k
                   DO 110 i = 1, m
                      c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
   110             CONTINUE
   120          CONTINUE
             END IF
          END IF
 !
       ELSE IF( lsame( storev, 'R' ) ) THEN
 !
          IF( lsame( direct, 'F' ) ) THEN
 !
 !           Let  V =  ( V1  V2 )    (V1: first K columns)
 !           where  V1  is unit upper triangular.
 !
             IF( lsame( side, 'L' ) ) THEN
 !
 !              Form  H * C  or  H**H * C  where  C = ( C1 )
 !                                                    ( C2 )
 !
 !              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
 !
 !              W := C1**H
 !
                DO 130 j = 1, k
                   CALL zcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
                   CALL zlacgv( n, work( 1, j ), 1 )
   130          CONTINUE
 !
 !              W := W * V1**H
 !
                CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose', &
                            'Unit', n, k, one, v, ldv, work, ldwork )
                IF( m.GT.k ) THEN
 !
 !                 W := W + C2**H * V2**H
 !
                   CALL zgemm( 'Conjugate transpose', &
                               'Conjugate transpose', n, k, m-k, one, &
                               c( k+1, 1 ), ldc, v( 1, k+1 ), ldv, one, &
                               work, ldwork )
                END IF
 !
 !              W := W * T**H  or  W * T
 !
                CALL ztrmm( 'Right', 'Upper', transt, 'Non-unit', n, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - V**H * W**H
 !
                IF( m.GT.k ) THEN
 !
 !                 C2 := C2 - V2**H * W**H
 !
                   CALL zgemm( 'Conjugate transpose', &
                               'Conjugate transpose', m-k, n, k, -one, &
                               v( 1, k+1 ), ldv, work, ldwork, one, &
                               c( k+1, 1 ), ldc )
                END IF
 !
 !              W := W * V1
 !
                CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', n, &
                            k, one, v, ldv, work, ldwork )
 !
 !              C1 := C1 - W**H
 !
                DO 150 j = 1, k
                   DO 140 i = 1, n
                      c( j, i ) = c( j, i ) - dconjg( work( i, j ) )
   140             CONTINUE
   150          CONTINUE
 !
             ELSE IF( lsame( side, 'R' ) ) THEN
 !
 !              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
 !
 !              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
 !
 !              W := C1
 !
                DO 160 j = 1, k
                   CALL zcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
   160          CONTINUE
 !
 !              W := W * V1**H
 !
                CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose', &
                            'Unit', m, k, one, v, ldv, work, ldwork )
                IF( n.GT.k ) THEN
 !
 !                 W := W + C2 * V2**H
 !
                   CALL zgemm( 'No transpose', 'Conjugate transpose', m, &
                               k, n-k, one, c( 1, k+1 ), ldc, &
                               v( 1, k+1 ), ldv, one, work, ldwork )
                END IF
 !
 !              W := W * T  or  W * T**H
 !
                CALL ztrmm( 'Right', 'Upper', trans, 'Non-unit', m, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - W * V
 !
                IF( n.GT.k ) THEN
 !
 !                 C2 := C2 - W * V2
 !
                   CALL zgemm( 'No transpose', 'No transpose', m, n-k, k, &
                               -one, work, ldwork, v( 1, k+1 ), ldv, one, &
                               c( 1, k+1 ), ldc )
                END IF
 !
 !              W := W * V1
 !
                CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', m, &
                            k, one, v, ldv, work, ldwork )
 !
 !              C1 := C1 - W
 !
                DO 180 j = 1, k
                   DO 170 i = 1, m
                      c( i, j ) = c( i, j ) - work( i, j )
   170             CONTINUE
   180          CONTINUE
 !
             END IF
 !
          ELSE
 !
 !           Let  V =  ( V1  V2 )    (V2: last K columns)
 !           where  V2  is unit lower triangular.
 !
             IF( lsame( side, 'L' ) ) THEN
 !
 !              Form  H * C  or  H**H * C  where  C = ( C1 )
 !                                                    ( C2 )
 !
 !              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
 !
 !              W := C2**H
 !
                DO 190 j = 1, k
                   CALL zcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1 )
                   CALL zlacgv( n, work( 1, j ), 1 )
   190          CONTINUE
 !
 !              W := W * V2**H
 !
                CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose', &
                            'Unit', n, k, one, v( 1, m-k+1 ), ldv, work, &
                            ldwork )
                IF( m.GT.k ) THEN
 !
 !                 W := W + C1**H * V1**H
 !
                   CALL zgemm( 'Conjugate transpose', &
                               'Conjugate transpose', n, k, m-k, one, c, &
                               ldc, v, ldv, one, work, ldwork )
                END IF
 !
 !              W := W * T**H  or  W * T
 !
                CALL ztrmm( 'Right', 'Lower', transt, 'Non-unit', n, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - V**H * W**H
 !
                IF( m.GT.k ) THEN
 !
 !                 C1 := C1 - V1**H * W**H
 !
                   CALL zgemm( 'Conjugate transpose', &
                               'Conjugate transpose', m-k, n, k, -one, v, &
                               ldv, work, ldwork, one, c, ldc )
                END IF
 !
 !              W := W * V2
 !
                CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', n, &
                            k, one, v( 1, m-k+1 ), ldv, work, ldwork )
 !
 !              C2 := C2 - W**H
 !
                DO 210 j = 1, k
                   DO 200 i = 1, n
                      c( m-k+j, i ) = c( m-k+j, i ) - &
                                      dconjg( work( i, j ) )
   200             CONTINUE
   210          CONTINUE
 !
             ELSE IF( lsame( side, 'R' ) ) THEN
 !
 !              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
 !
 !              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
 !
 !              W := C2
 !
                DO 220 j = 1, k
                   CALL zcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
   220          CONTINUE
 !
 !              W := W * V2**H
 !
                CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose', &
                            'Unit', m, k, one, v( 1, n-k+1 ), ldv, work, &
                            ldwork )
                IF( n.GT.k ) THEN
 !
 !                 W := W + C1 * V1**H
 !
                   CALL zgemm( 'No transpose', 'Conjugate transpose', m, &
                               k, n-k, one, c, ldc, v, ldv, one, work, &
                               ldwork )
                END IF
 !
 !              W := W * T  or  W * T**H
 !
                CALL ztrmm( 'Right', 'Lower', trans, 'Non-unit', m, k, &
                            one, t, ldt, work, ldwork )
 !
 !              C := C - W * V
 !
                IF( n.GT.k ) THEN
 !
 !                 C1 := C1 - W * V1
 !
                   CALL zgemm( 'No transpose', 'No transpose', m, n-k, k, &
                               -one, work, ldwork, v, ldv, one, c, ldc )
                END IF
 !
 !              W := W * V2
 !
                CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', m, &
                            k, one, v( 1, n-k+1 ), ldv, work, ldwork )
 !
 !              C1 := C1 - W
 !
                DO 240 j = 1, k
                   DO 230 i = 1, m
                      c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
   230             CONTINUE
   240          CONTINUE
 !
             END IF
 !
          END IF
       END IF
 !
       RETURN
 !
 !     End of ZLARFB
 !
  end subroutine zlarfb

 !  =====================================================================
  subroutine zcopy(N,ZX,INCX,ZY,INCY)
 !
 !  -- Reference BLAS level1 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,INCY,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 ZX(*),ZY(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Local Scalars ..
       INTEGER I,IX,IY
 !     ..
       IF (n.LE.0) RETURN
       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
 !
 !        code for both increments equal to 1
 !
          DO i = 1,n
           zy(i) = zx(i)
          END DO
       ELSE
 !
 !        code for unequal increments or equal increments
 !          not equal to 1
 !
          ix = 1
          iy = 1
          IF (incx.LT.0) ix = (-n+1)*incx + 1
          IF (incy.LT.0) iy = (-n+1)*incy + 1
          DO i = 1,n
             zy(iy) = zx(ix)
             ix = ix + incx
             iy = iy + incy
          END DO
       END IF
       RETURN
  end subroutine zcopy


 !  =====================================================================
  subroutine ztrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
 !
 !  -- Reference BLAS level3 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       INTEGER LDA,LDB,M,N
       CHARACTER DIAG,SIDE,TRANSA,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),B(ldb,*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. External Functions ..
 !      LOGICAL LSAME
 !      EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,J,K,NROWA
       LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
 !     ..
 !     .. Parameters ..
       COMPLEX*16 ONE
       parameter(one= (1.0d0,0.0d0))
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !
 !     Test the input parameters.
 !
       lside = lsame(side,'L')
       IF (lside) THEN
           nrowa = m
       ELSE
           nrowa = n
       END IF
       noconj = lsame(transa,'T')
       nounit = lsame(diag,'N')
       upper = lsame(uplo,'U')
 !
       info = 0
       IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
           info = 1
       ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
           info = 2
       ELSE IF ((.NOT.lsame(transa,'N')) .AND. &
                (.NOT.lsame(transa,'T')) .AND. &
                (.NOT.lsame(transa,'C'))) THEN
           info = 3
       ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
           info = 4
       ELSE IF (m.LT.0) THEN
           info = 5
       ELSE IF (n.LT.0) THEN
           info = 6
       ELSE IF (lda.LT.max(1,nrowa)) THEN
           info = 9
       ELSE IF (ldb.LT.max(1,m)) THEN
           info = 11
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZTRMM ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (m.EQ.0 .OR. n.EQ.0) RETURN
 !
 !     And when  alpha.eq.zero.
 !
       IF (alpha.EQ.zero) THEN
           DO 20 j = 1,n
               DO 10 i = 1,m
                   b(i,j) = zero
    10         CONTINUE
    20     CONTINUE
           RETURN
       END IF
 !
 !     Start the operations.
 !
       IF (lside) THEN
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*A*B.
 !
               IF (upper) THEN
                   DO 50 j = 1,n
                       DO 40 k = 1,m
                           IF (b(k,j).NE.zero) THEN
                               temp = alpha*b(k,j)
                               DO 30 i = 1,k - 1
                                   b(i,j) = b(i,j) + temp*a(i,k)
    30                         CONTINUE
                               IF (nounit) temp = temp*a(k,k)
                               b(k,j) = temp
                           END IF
    40                 CONTINUE
    50             CONTINUE
               ELSE
                   DO 80 j = 1,n
                       DO 70 k = m,1,-1
                           IF (b(k,j).NE.zero) THEN
                               temp = alpha*b(k,j)
                               b(k,j) = temp
                               IF (nounit) b(k,j) = b(k,j)*a(k,k)
                               DO 60 i = k + 1,m
                                   b(i,j) = b(i,j) + temp*a(i,k)
    60                         CONTINUE
                           END IF
    70                 CONTINUE
    80             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
 !
               IF (upper) THEN
                   DO 120 j = 1,n
                       DO 110 i = m,1,-1
                           temp = b(i,j)
                           IF (noconj) THEN
                               IF (nounit) temp = temp*a(i,i)
                               DO 90 k = 1,i - 1
                                   temp = temp + a(k,i)*b(k,j)
    90                         CONTINUE
                           ELSE
                               IF (nounit) temp = temp*dconjg(a(i,i))
                               DO 100 k = 1,i - 1
                                   temp = temp + dconjg(a(k,i))*b(k,j)
   100                         CONTINUE
                           END IF
                           b(i,j) = alpha*temp
   110                 CONTINUE
   120             CONTINUE
               ELSE
                   DO 160 j = 1,n
                       DO 150 i = 1,m
                           temp = b(i,j)
                           IF (noconj) THEN
                               IF (nounit) temp = temp*a(i,i)
                               DO 130 k = i + 1,m
                                   temp = temp + a(k,i)*b(k,j)
   130                         CONTINUE
                           ELSE
                               IF (nounit) temp = temp*dconjg(a(i,i))
                               DO 140 k = i + 1,m
                                   temp = temp + dconjg(a(k,i))*b(k,j)
   140                         CONTINUE
                           END IF
                           b(i,j) = alpha*temp
   150                 CONTINUE
   160             CONTINUE
               END IF
           END IF
       ELSE
           IF (lsame(transa,'N')) THEN
 !
 !           Form  B := alpha*B*A.
 !
               IF (upper) THEN
                   DO 200 j = n,1,-1
                       temp = alpha
                       IF (nounit) temp = temp*a(j,j)
                       DO 170 i = 1,m
                           b(i,j) = temp*b(i,j)
   170                 CONTINUE
                       DO 190 k = 1,j - 1
                           IF (a(k,j).NE.zero) THEN
                               temp = alpha*a(k,j)
                               DO 180 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   180                         CONTINUE
                           END IF
   190                 CONTINUE
   200             CONTINUE
               ELSE
                   DO 240 j = 1,n
                       temp = alpha
                       IF (nounit) temp = temp*a(j,j)
                       DO 210 i = 1,m
                           b(i,j) = temp*b(i,j)
   210                 CONTINUE
                       DO 230 k = j + 1,n
                           IF (a(k,j).NE.zero) THEN
                               temp = alpha*a(k,j)
                               DO 220 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   220                         CONTINUE
                           END IF
   230                 CONTINUE
   240             CONTINUE
               END IF
           ELSE
 !
 !           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
 !
               IF (upper) THEN
                   DO 280 k = 1,n
                       DO 260 j = 1,k - 1
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = alpha*a(j,k)
                               ELSE
                                   temp = alpha*dconjg(a(j,k))
                               END IF
                               DO 250 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   250                         CONTINUE
                           END IF
   260                 CONTINUE
                       temp = alpha
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = temp*a(k,k)
                           ELSE
                               temp = temp*dconjg(a(k,k))
                           END IF
                       END IF
                       IF (temp.NE.one) THEN
                           DO 270 i = 1,m
                               b(i,k) = temp*b(i,k)
   270                     CONTINUE
                       END IF
   280             CONTINUE
               ELSE
                   DO 320 k = n,1,-1
                       DO 300 j = k + 1,n
                           IF (a(j,k).NE.zero) THEN
                               IF (noconj) THEN
                                   temp = alpha*a(j,k)
                               ELSE
                                   temp = alpha*dconjg(a(j,k))
                               END IF
                               DO 290 i = 1,m
                                   b(i,j) = b(i,j) + temp*b(i,k)
   290                         CONTINUE
                           END IF
   300                 CONTINUE
                       temp = alpha
                       IF (nounit) THEN
                           IF (noconj) THEN
                               temp = temp*a(k,k)
                           ELSE
                               temp = temp*dconjg(a(k,k))
                           END IF
                       END IF
                       IF (temp.NE.one) THEN
                           DO 310 i = 1,m
                               b(i,k) = temp*b(i,k)
   310                     CONTINUE
                       END IF
   320             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRMM .
 !
  end subroutine ztrmm

  !  =====================================================================
  subroutine zlarft( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
  !
  !  -- LAPACK auxiliary routine (version 3.6.0) --
  !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !     November 2015
  !
  !     .. Scalar Arguments ..
        CHARACTER          DIRECT, STOREV
        INTEGER            K, LDT, LDV, N
  !     ..
  !     .. Array Arguments ..
        COMPLEX*16         T( ldt, * ), TAU( * ), V( ldv, * )
  !     ..
  !
  !  =====================================================================
  !
  !     .. Parameters ..
        COMPLEX*16         ONE, ZERO
        parameter( one = ( 1.0d0, 0.0d0 ), &
                           zero = ( 0.0d0, 0.0d0 ) )
  !     ..
  !     .. Local Scalars ..
        INTEGER            I, J, PREVLASTV, LASTV
  !     ..
  !     .. External Subroutines ..
  !      EXTERNAL           zgemv, zlacgv, ztrmv, zgemm
  !     ..
  !     .. External Functions ..
  !      LOGICAL            LSAME
  !      EXTERNAL           lsame
  !     ..
  !     .. Executable Statements ..
  !
  !     Quick return if possible
  !
        IF( n.EQ.0 ) &
           RETURN
  !
        IF( lsame( direct, 'F' ) ) THEN
           prevlastv = n
           DO i = 1, k
              prevlastv = max( prevlastv, i )
              IF( tau( i ).EQ.zero ) THEN
  !
  !              H(i)  =  I
  !
                 DO j = 1, i
                    t( j, i ) = zero
                 END DO
              ELSE
  !
  !              general case
  !
                 IF( lsame( storev, 'C' ) ) THEN
  !                 Skip any trailing zeros.
                    DO lastv = n, i+1, -1
                       IF( v( lastv, i ).NE.zero ) EXIT
                    END DO
                    DO j = 1, i-1
                       t( j, i ) = -tau( i ) * conjg( v( i , j ) )
                    END DO
                    j = min( lastv, prevlastv )
  !
  !                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
  !
                    CALL zgemv( 'Conjugate transpose', j-i, i-1, &
                                -tau( i ), v( i+1, 1 ), ldv, &
                                v( i+1, i ), 1, one, t( 1, i ), 1 )
                 ELSE
  !                 Skip any trailing zeros.
                    DO lastv = n, i+1, -1
                       IF( v( i, lastv ).NE.zero ) EXIT
                    END DO
                    DO j = 1, i-1
                       t( j, i ) = -tau( i ) * v( j , i )
                    END DO
                    j = min( lastv, prevlastv )
  !
  !                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
  !
                    CALL zgemm( 'N', 'C', i-1, 1, j-i, -tau( i ), &
                                v( 1, i+1 ), ldv, v( i, i+1 ), ldv, &
                                one, t( 1, i ), ldt )
                 END IF
  !
  !              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
  !
                 CALL ztrmv( 'Upper', 'No transpose', 'Non-unit', i-1, t, &
                             ldt, t( 1, i ), 1 )
                 t( i, i ) = tau( i )
                 IF( i.GT.1 ) THEN
                    prevlastv = max( prevlastv, lastv )
                 ELSE
                    prevlastv = lastv
                 END IF
               END IF
           END DO
        ELSE
           prevlastv = 1
           DO i = k, 1, -1
              IF( tau( i ).EQ.zero ) THEN
  !
  !              H(i)  =  I
  !
                 DO j = i, k
                    t( j, i ) = zero
                 END DO
              ELSE
  !
  !              general case
  !
                 IF( i.LT.k ) THEN
                    IF( lsame( storev, 'C' ) ) THEN
  !                    Skip any leading zeros.
                       DO lastv = 1, i-1
                          IF( v( lastv, i ).NE.zero ) EXIT
                       END DO
                       DO j = i+1, k
                          t( j, i ) = -tau( i ) * conjg( v( n-k+i , j ) )
                       END DO
                       j = max( lastv, prevlastv )
  !
  !                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
  !
                       CALL zgemv( 'Conjugate transpose', n-k+i-j, k-i, &
                                   -tau( i ), v( j, i+1 ), ldv, v( j, i ), &
                                   1, one, t( i+1, i ), 1 )
                    ELSE
  !                    Skip any leading zeros.
                       DO lastv = 1, i-1
                          IF( v( i, lastv ).NE.zero ) EXIT
                       END DO
                       DO j = i+1, k
                          t( j, i ) = -tau( i ) * v( j, n-k+i )
                       END DO
                       j = max( lastv, prevlastv )
  !
  !                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
  !
                       CALL zgemm( 'N', 'C', k-i, 1, n-k+i-j, -tau( i ), &
                                   v( i+1, j ), ldv, v( i, j ), ldv, &
                                   one, t( i+1, i ), ldt )
                    END IF
  !
  !                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
  !
                    CALL ztrmv( 'Lower', 'No transpose', 'Non-unit', k-i, &
                                t( i+1, i+1 ), ldt, t( i+1, i ), 1 )
                    IF( i.GT.1 ) THEN
                       prevlastv = min( prevlastv, lastv )
                    ELSE
                       prevlastv = lastv
                    END IF
                 END IF
                 t( i, i ) = tau( i )
              END IF
           END DO
        END IF
        RETURN
  !
  !     End of ZLARFT
  !
  end subroutine zlarft

 !  =====================================================================
  subroutine ztrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
 !
 !  -- Reference BLAS level2 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER INCX,LDA,N
       CHARACTER DIAG,TRANS,UPLO
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),X(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,IX,J,JX,KX
       LOGICAL NOCONJ,NOUNIT
 !     ..
 !     .. External Functions ..
 !      LOGICAL LSAME
 !      EXTERNAL lsame
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
           info = 1
       ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND. &
                .NOT.lsame(trans,'C')) THEN
           info = 2
       ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
           info = 3
       ELSE IF (n.LT.0) THEN
           info = 4
       ELSE IF (lda.LT.max(1,n)) THEN
           info = 6
       ELSE IF (incx.EQ.0) THEN
           info = 8
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZTRMV ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF (n.EQ.0) RETURN
 !
       noconj = lsame(trans,'T')
       nounit = lsame(diag,'N')
 !
 !     Set up the start point in X if the increment is not unity. This
 !     will be  ( N - 1 )*INCX  too small for descending loops.
 !
       IF (incx.LE.0) THEN
           kx = 1 - (n-1)*incx
       ELSE IF (incx.NE.1) THEN
           kx = 1
       END IF
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
       IF (lsame(trans,'N')) THEN
 !
 !        Form  x := A*x.
 !
           IF (lsame(uplo,'U')) THEN
               IF (incx.EQ.1) THEN
                   DO 20 j = 1,n
                       IF (x(j).NE.zero) THEN
                           temp = x(j)
                           DO 10 i = 1,j - 1
                               x(i) = x(i) + temp*a(i,j)
    10                     CONTINUE
                           IF (nounit) x(j) = x(j)*a(j,j)
                       END IF
    20             CONTINUE
               ELSE
                   jx = kx
                   DO 40 j = 1,n
                       IF (x(jx).NE.zero) THEN
                           temp = x(jx)
                           ix = kx
                           DO 30 i = 1,j - 1
                               x(ix) = x(ix) + temp*a(i,j)
                               ix = ix + incx
    30                     CONTINUE
                           IF (nounit) x(jx) = x(jx)*a(j,j)
                       END IF
                       jx = jx + incx
    40             CONTINUE
               END IF
           ELSE
               IF (incx.EQ.1) THEN
                   DO 60 j = n,1,-1
                       IF (x(j).NE.zero) THEN
                           temp = x(j)
                           DO 50 i = n,j + 1,-1
                               x(i) = x(i) + temp*a(i,j)
    50                     CONTINUE
                           IF (nounit) x(j) = x(j)*a(j,j)
                       END IF
    60             CONTINUE
               ELSE
                   kx = kx + (n-1)*incx
                   jx = kx
                   DO 80 j = n,1,-1
                       IF (x(jx).NE.zero) THEN
                           temp = x(jx)
                           ix = kx
                           DO 70 i = n,j + 1,-1
                               x(ix) = x(ix) + temp*a(i,j)
                               ix = ix - incx
    70                     CONTINUE
                           IF (nounit) x(jx) = x(jx)*a(j,j)
                       END IF
                       jx = jx - incx
    80             CONTINUE
               END IF
           END IF
       ELSE
 !
 !        Form  x := A**T*x  or  x := A**H*x.
 !
           IF (lsame(uplo,'U')) THEN
               IF (incx.EQ.1) THEN
                   DO 110 j = n,1,-1
                       temp = x(j)
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 90 i = j - 1,1,-1
                               temp = temp + a(i,j)*x(i)
    90                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 100 i = j - 1,1,-1
                               temp = temp + dconjg(a(i,j))*x(i)
   100                     CONTINUE
                       END IF
                       x(j) = temp
   110             CONTINUE
               ELSE
                   jx = kx + (n-1)*incx
                   DO 140 j = n,1,-1
                       temp = x(jx)
                       ix = jx
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 120 i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + a(i,j)*x(ix)
   120                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 130 i = j - 1,1,-1
                               ix = ix - incx
                               temp = temp + dconjg(a(i,j))*x(ix)
   130                     CONTINUE
                       END IF
                       x(jx) = temp
                       jx = jx - incx
   140             CONTINUE
               END IF
           ELSE
               IF (incx.EQ.1) THEN
                   DO 170 j = 1,n
                       temp = x(j)
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 150 i = j + 1,n
                               temp = temp + a(i,j)*x(i)
   150                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 160 i = j + 1,n
                               temp = temp + dconjg(a(i,j))*x(i)
   160                     CONTINUE
                       END IF
                       x(j) = temp
   170             CONTINUE
               ELSE
                   jx = kx
                   DO 200 j = 1,n
                       temp = x(jx)
                       ix = jx
                       IF (noconj) THEN
                           IF (nounit) temp = temp*a(j,j)
                           DO 180 i = j + 1,n
                               ix = ix + incx
                               temp = temp + a(i,j)*x(ix)
   180                     CONTINUE
                       ELSE
                           IF (nounit) temp = temp*dconjg(a(j,j))
                           DO 190 i = j + 1,n
                               ix = ix + incx
                               temp = temp + dconjg(a(i,j))*x(ix)
   190                     CONTINUE
                       END IF
                       x(jx) = temp
                       jx = jx + incx
   200             CONTINUE
               END IF
           END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRMV .
 !
  end subroutine ztrmv

 !  =====================================================================
  subroutine zung2l( M, N, K, A, LDA, TAU, WORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, K, LDA, M, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d0, 0.0d0 ), &
                          zero = ( 0.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, II, J, L
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zlarf, zscal
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
          info = -2
       ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -5
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZUNG2L', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.LE.0 ) &
          RETURN
 !
 !     Initialise columns 1:n-k to columns of the unit matrix
 !
       DO 20 j = 1, n - k
          DO 10 l = 1, m
             a( l, j ) = zero
    10    CONTINUE
          a( m-n+j, j ) = one
    20 CONTINUE
 !
       DO 40 i = 1, k
          ii = n - k + i
 !
 !        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
 !
          a( m-n+ii, ii ) = one
          CALL zlarf( 'Left', m-n+ii, ii-1, a( 1, ii ), 1, tau( i ), a, &
                      lda, work )
          CALL zscal( m-n+ii-1, -tau( i ), a( 1, ii ), 1 )
          a( m-n+ii, ii ) = one - tau( i )
 !
 !        Set A(m-k+i+1:m,n-k+i) to zero
 !
          DO 30 l = m - n + ii + 1, m
             a( l, ii ) = zero
    30    CONTINUE
    40 CONTINUE
       RETURN
 !
 !     End of ZUNG2L
 !
  end subroutine zung2l

 !  =====================================================================
 subroutine zlarf( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          SIDE
       INTEGER            INCV, LDC, M, N
       COMPLEX*16         TAU
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         C( ldc, * ), V( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d0, 0.0d0 ), &
                          zero = ( 0.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            APPLYLEFT
       INTEGER            I, LASTV, LASTC
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           zgemv, zgerc
 !     ..
 !     .. External Functions ..
 !      LOGICAL            LSAME
 !      INTEGER            ILAZLR, ILAZLC
 !      EXTERNAL           lsame, ilazlr, ilazlc
 !     ..
 !     .. Executable Statements ..
 !
       applyleft = lsame( side, 'L' )
       lastv = 0
       lastc = 0
       IF( tau.NE.zero ) THEN
 !     Set up variables for scanning V.  LASTV begins pointing to the end
 !     of V.
          IF( applyleft ) THEN
             lastv = m
          ELSE
             lastv = n
          END IF
          IF( incv.GT.0 ) THEN
             i = 1 + (lastv-1) * incv
          ELSE
             i = 1
          END IF
 !     Look for the last non-zero row in V.
          DO WHILE( lastv.GT.0 .AND. v( i ).EQ.zero )
             lastv = lastv - 1
             i = i - incv
          END DO
          IF( applyleft ) THEN
 !     Scan for the last non-zero column in C(1:lastv,:).
             lastc = ilazlc(lastv, n, c, ldc)
          ELSE
 !     Scan for the last non-zero row in C(:,1:lastv).
             lastc = ilazlr(m, lastv, c, ldc)
          END IF
       END IF
 !     Note that lastc.eq.0 renders the BLAS operations null; no special
 !     case is needed at this level.
       IF( applyleft ) THEN
 !
 !        Form  H * C
 !
          IF( lastv.GT.0 ) THEN
 !
 !           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
 !
             CALL zgemv( 'Conjugate transpose', lastv, lastc, one, &
                  c, ldc, v, incv, zero, work, 1 )
 !
 !           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
 !
             CALL zgerc( lastv, lastc, -tau, v, incv, work, 1, c, ldc )
          END IF
       ELSE
 !
 !        Form  C * H
 !
          IF( lastv.GT.0 ) THEN
 !
 !           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
 !
             CALL zgemv( 'No transpose', lastc, lastv, one, c, ldc, &
                  v, incv, zero, work, 1 )
 !
 !           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
 !
             CALL zgerc( lastc, lastv, -tau, work, 1, v, incv, c, ldc )
          END IF
       END IF
       RETURN
 !
 !     End of ZLARF
 !
  end subroutine zlarf

 !  =====================================================================
  subroutine zgerc(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
 !
 !  -- Reference BLAS level2 routine (version 3.4.0) --
 !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       COMPLEX*16 ALPHA
       INTEGER INCX,INCY,LDA,M,N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16 A(lda,*),X(*),Y(*)
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16 ZERO
       parameter(zero= (0.0d0,0.0d0))
 !     ..
 !     .. Local Scalars ..
       COMPLEX*16 TEMP
       INTEGER I,INFO,IX,J,JY,KX
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL xerbla
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC dconjg,max
 !     ..
 !
 !     Test the input parameters.
 !
       info = 0
       IF (m.LT.0) THEN
           info = 1
       ELSE IF (n.LT.0) THEN
           info = 2
       ELSE IF (incx.EQ.0) THEN
           info = 5
       ELSE IF (incy.EQ.0) THEN
           info = 7
       ELSE IF (lda.LT.max(1,m)) THEN
           info = 9
       END IF
       IF (info.NE.0) THEN
           CALL xerbla('ZGERC ',info)
           RETURN
       END IF
 !
 !     Quick return if possible.
 !
       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
 !
 !     Start the operations. In this version the elements of A are
 !     accessed sequentially with one pass through A.
 !
       IF (incy.GT.0) THEN
           jy = 1
       ELSE
           jy = 1 - (n-1)*incy
       END IF
       IF (incx.EQ.1) THEN
           DO 20 j = 1,n
               IF (y(jy).NE.zero) THEN
                   temp = alpha*dconjg(y(jy))
                   DO 10 i = 1,m
                       a(i,j) = a(i,j) + x(i)*temp
    10             CONTINUE
               END IF
               jy = jy + incy
    20     CONTINUE
       ELSE
           IF (incx.GT.0) THEN
               kx = 1
           ELSE
               kx = 1 - (m-1)*incx
           END IF
           DO 40 j = 1,n
               IF (y(jy).NE.zero) THEN
                   temp = alpha*dconjg(y(jy))
                   ix = kx
                   DO 30 i = 1,m
                       a(i,j) = a(i,j) + x(ix)*temp
                       ix = ix + incx
    30             CONTINUE
               END IF
               jy = jy + incy
    40     CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZGERC .
 !
  end subroutine zgerc

 !  =====================================================================
  integer function ilazlr( M, N, A, LDA )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       INTEGER            M, N, LDA
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16       ZERO
       parameter( zero = (0.0d0, 0.0d0) )
 !     ..
 !     .. Local Scalars ..
       INTEGER I, J
 !     ..
 !     .. Executable Statements ..
 !
 !     Quick test for the common case where one corner is non-zero.
       IF( m.EQ.0 ) THEN
          ilazlr = m
       ELSE IF( a(m, 1).NE.zero .OR. a(m, n).NE.zero ) THEN
          ilazlr = m
       ELSE
 !     Scan up each column tracking the last zero row seen.
          ilazlr = 0
          DO j = 1, n
             i=m
             DO WHILE((a(max(i,1),j).EQ.zero).AND.(i.GE.1))
                i=i-1
             ENDDO
             ilazlr = max( ilazlr, i )
          END DO
       END IF
       RETURN
   end function ilazlr

 !  =====================================================================
  integer function ilazlc( M, N, A, LDA )
 !
 !  -- LAPACK auxiliary routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       INTEGER            M, N, LDA
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16       ZERO
       parameter( zero = (0.0d0, 0.0d0) )
 !     ..
 !     .. Local Scalars ..
       INTEGER I
 !     ..
 !     .. Executable Statements ..
 !
 !     Quick test for the common case where one corner is non-zero.
       IF( n.EQ.0 ) THEN
          ilazlc = n
       ELSE IF( a(1, n).NE.zero .OR. a(m, n).NE.zero ) THEN
          ilazlc = n
       ELSE
 !     Now scan each column from the end, returning with the first non-zero.
          DO ilazlc = n, 1, -1
             DO i = 1, m
                IF( a(i, ilazlc).NE.zero ) RETURN
             END DO
          END DO
       END IF
       RETURN
  end function ilazlc

 !  =====================================================================
  subroutine zungqr( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, K, LDA, LWORK, M, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ZERO
       parameter( zero = ( 0.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            LQUERY
       INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, &
                          lwkopt, nb, nbmin, nx
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zlarfb, zlarft, zung2r
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. External Functions ..
 !      INTEGER            ILAENV
 !      EXTERNAL           ilaenv
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
       nb = ilaenv( 1, 'ZUNGQR', ' ', m, n, k, -1 )
       lwkopt = max( 1, n )*nb
       work( 1 ) = lwkopt
       lquery = ( lwork.EQ.-1 )
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
          info = -2
       ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -5
       ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
          info = -8
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZUNGQR', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.LE.0 ) THEN
          work( 1 ) = 1
          RETURN
       END IF
 !
       nbmin = 2
       nx = 0
       iws = n
       IF( nb.GT.1 .AND. nb.LT.k ) THEN
 !
 !        Determine when to cross over from blocked to unblocked code.
 !
          nx = max( 0, ilaenv( 3, 'ZUNGQR', ' ', m, n, k, -1 ) )
          IF( nx.LT.k ) THEN
 !
 !           Determine if workspace is large enough for blocked code.
 !
             ldwork = n
             iws = ldwork*nb
             IF( lwork.LT.iws ) THEN
 !
 !              Not enough workspace to use optimal NB:  reduce NB and
 !              determine the minimum value of NB.
 !
                nb = lwork / ldwork
                nbmin = max( 2, ilaenv( 2, 'ZUNGQR', ' ', m, n, k, -1 ) )
             END IF
          END IF
       END IF
 !
       IF( nb.GE.nbmin .AND. nb.LT.k .AND. nx.LT.k ) THEN
 !
 !        Use blocked code after the last block.
 !        The first kk columns are handled by the block method.
 !
          ki = ( ( k-nx-1 ) / nb )*nb
          kk = min( k, ki+nb )
 !
 !        Set A(1:kk,kk+1:n) to zero.
 !
          DO 20 j = kk + 1, n
             DO 10 i = 1, kk
                a( i, j ) = zero
    10       CONTINUE
    20    CONTINUE
       ELSE
          kk = 0
       END IF
 !
 !     Use unblocked code for the last or only block.
 !
       IF( kk.LT.n ) &
          CALL zung2r( m-kk, n-kk, k-kk, a( kk+1, kk+1 ), lda, &
                       tau( kk+1 ), work, iinfo )
 !
       IF( kk.GT.0 ) THEN
 !
 !        Use blocked code
 !
          DO 50 i = ki + 1, 1, -nb
             ib = min( nb, k-i+1 )
             IF( i+ib.LE.n ) THEN
 !
 !              Form the triangular factor of the block reflector
 !              H = H(i) H(i+1) . . . H(i+ib-1)
 !
                CALL zlarft( 'Forward', 'Columnwise', m-i+1, ib, &
                             a( i, i ), lda, tau( i ), work, ldwork )
 !
 !              Apply H to A(i:m,i+ib:n) from the left
 !
                CALL zlarfb( 'Left', 'No transpose', 'Forward', &
                             'Columnwise', m-i+1, n-i-ib+1, ib, &
                             a( i, i ), lda, work, ldwork, a( i, i+ib ), &
                             lda, work( ib+1 ), ldwork )
             END IF
 !
 !           Apply H to rows i:m of current block
 !
             CALL zung2r( m-i+1, ib, ib, a( i, i ), lda, tau( i ), work, &
                          iinfo )
 !
 !           Set rows 1:i-1 of current block to zero
 !
             DO 40 j = i, i + ib - 1
                DO 30 l = 1, i - 1
                   a( l, j ) = zero
    30          CONTINUE
    40       CONTINUE
    50    CONTINUE
       END IF
 !
       work( 1 ) = iws
       RETURN
 !
 !     End of ZUNGQR
 !
  end subroutine zungqr

 !  =====================================================================
  subroutine zung2r( M, N, K, A, LDA, TAU, WORK, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       INTEGER            INFO, K, LDA, M, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d0, 0.0d0 ), &
                          zero = ( 0.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       INTEGER            I, J, L
 !     ..
 !     .. External Subroutines ..
 !      EXTERNAL           xerbla, zlarf, zscal
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input arguments
 !
       info = 0
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
          info = -2
       ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -5
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZUNG2R', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.LE.0 ) &
          RETURN
 !
 !     Initialise columns k+1:n to columns of the unit matrix
 !
       DO 20 j = k + 1, n
          DO 10 l = 1, m
             a( l, j ) = zero
    10    CONTINUE
          a( j, j ) = one
    20 CONTINUE
 !
       DO 40 i = k, 1, -1
 !
 !        Apply H(i) to A(i:m,i:n) from the left
 !
          IF( i.LT.n ) THEN
             a( i, i ) = one
             CALL zlarf( 'Left', m-i+1, n-i, a( i, i ), 1, tau( i ), &
                         a( i, i+1 ), lda, work )
          END IF
          IF ( i.LT.m ) &
             CALL zscal( m-i, -tau( i ), a( i+1, i ), 1 )
          a( i, i ) = one - tau( i )
 !
 !        Set A(1:i-1,i) to zero
 !
          DO 30 l = 1, i - 1
             a( l, i ) = zero
    30    CONTINUE
    40 CONTINUE
       RETURN
 !
 !     End of ZUNG2R
 !
  end subroutine zung2r

 !  =====================================================================
  subroutine ztrtri( UPLO, DIAG, N, A, LDA, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.0) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     November 2011
 !
 !     .. Scalar Arguments ..
       CHARACTER          DIAG, UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE, ZERO
       parameter( one = ( 1.0d0, 0.0d0 ), &
                          zero = ( 0.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            NOUNIT, UPPER
       INTEGER            J, JB, NB, NN
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       INTEGER            ILAENV
       EXTERNAL           lsame, ilaenv
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           xerbla, ztrmm, ztrsm, ztrti2
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max, min
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       info = 0
       upper = lsame( uplo, 'U' )
       nounit = lsame( diag, 'N' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
          info = -2
       ELSE IF( n.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZTRTRI', -info )
          RETURN
       END IF
 !
 !     Quick return if possible
 !
       IF( n.EQ.0 ) &
          RETURN
 !
 !     Check for singularity if non-unit.
 !
       IF( nounit ) THEN
          DO 10 info = 1, n
             IF( a( info, info ).EQ.zero ) &
                RETURN
    10    CONTINUE
          info = 0
       END IF
 !
 !     Determine the block size for this environment.
 !
       nb = ilaenv( 1, 'ZTRTRI', uplo // diag, n, -1, -1, -1 )
       IF( nb.LE.1 .OR. nb.GE.n ) THEN
 !
 !        Use unblocked code
 !
          CALL ztrti2( uplo, diag, n, a, lda, info )
       ELSE
 !
 !        Use blocked code
 !
          IF( upper ) THEN
 !
 !           Compute inverse of upper triangular matrix
 !
             DO 20 j = 1, n, nb
                jb = min( nb, n-j+1 )
 !
 !              Compute rows 1:j-1 of current block column
 !
                CALL ztrmm( 'Left', 'Upper', 'No transpose', diag, j-1, &
                            jb, one, a, lda, a( 1, j ), lda )
                CALL ztrsm( 'Right', 'Upper', 'No transpose', diag, j-1, &
                            jb, -one, a( j, j ), lda, a( 1, j ), lda )
 !
 !              Compute inverse of current diagonal block
 !
                CALL ztrti2( 'Upper', diag, jb, a( j, j ), lda, info )
    20       CONTINUE
          ELSE
 !
 !           Compute inverse of lower triangular matrix
 !
             nn = ( ( n-1 ) / nb )*nb + 1
             DO 30 j = nn, 1, -nb
                jb = min( nb, n-j+1 )
                IF( j+jb.LE.n ) THEN
 !
 !                 Compute rows j+jb:n of current block column
 !
                   CALL ztrmm( 'Left', 'Lower', 'No transpose', diag, &
                               n-j-jb+1, jb, one, a( j+jb, j+jb ), lda, &
                               a( j+jb, j ), lda )
                   CALL ztrsm( 'Right', 'Lower', 'No transpose', diag, &
                               n-j-jb+1, jb, -one, a( j, j ), lda, &
                               a( j+jb, j ), lda )
                END IF
 !
 !              Compute inverse of current diagonal block
 !
                CALL ztrti2( 'Lower', diag, jb, a( j, j ), lda, info )
    30       CONTINUE
          END IF
       END IF
 !
       RETURN
 !
 !     End of ZTRTRI
 !
  end subroutine ztrtri

 !  =====================================================================
  subroutine ztrti2( UPLO, DIAG, N, A, LDA, INFO )
 !
 !  -- LAPACK computational routine (version 3.4.2) --
 !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !     September 2012
 !
 !     .. Scalar Arguments ..
       CHARACTER          DIAG, UPLO
       INTEGER            INFO, LDA, N
 !     ..
 !     .. Array Arguments ..
       COMPLEX*16         A( lda, * )
 !     ..
 !
 !  =====================================================================
 !
 !     .. Parameters ..
       COMPLEX*16         ONE
       parameter( one = ( 1.0d0, 0.0d0 ) )
 !     ..
 !     .. Local Scalars ..
       LOGICAL            NOUNIT, UPPER
       INTEGER            J
       COMPLEX*16         AJJ
 !     ..
 !     .. External Functions ..
       LOGICAL            LSAME
       EXTERNAL           lsame
 !     ..
 !     .. External Subroutines ..
       EXTERNAL           xerbla, zscal, ztrmv
 !     ..
 !     .. Intrinsic Functions ..
       INTRINSIC          max
 !     ..
 !     .. Executable Statements ..
 !
 !     Test the input parameters.
 !
       info = 0
       upper = lsame( uplo, 'U' )
       nounit = lsame( diag, 'N' )
       IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
          info = -1
       ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
          info = -2
       ELSE IF( n.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, n ) ) THEN
          info = -5
       END IF
       IF( info.NE.0 ) THEN
          CALL xerbla( 'ZTRTI2', -info )
          RETURN
       END IF
 !
       IF( upper ) THEN
 !
 !        Compute inverse of upper triangular matrix.
 !
          DO 10 j = 1, n
             IF( nounit ) THEN
                a( j, j ) = one / a( j, j )
                ajj = -a( j, j )
             ELSE
                ajj = -one
             END IF
 !
 !           Compute elements 1:j-1 of j-th column.
 !
             CALL ztrmv( 'Upper', 'No transpose', diag, j-1, a, lda, &
                         a( 1, j ), 1 )
             CALL zscal( j-1, ajj, a( 1, j ), 1 )
    10    CONTINUE
       ELSE
 !
 !        Compute inverse of lower triangular matrix.
 !
          DO 20 j = n, 1, -1
             IF( nounit ) THEN
                a( j, j ) = one / a( j, j )
                ajj = -a( j, j )
             ELSE
                ajj = -one
             END IF
             IF( j.LT.n ) THEN
 !
 !              Compute elements j+1:n of j-th column.
 !
                CALL ztrmv( 'Lower', 'No transpose', diag, n-j, &
                            a( j+1, j+1 ), lda, a( j+1, j ), 1 )
                CALL zscal( n-j, ajj, a( j+1, j ), 1 )
             END IF
    20    CONTINUE
       END IF
 !
       RETURN
 !
 !     End of ZTRTI2
 !
  end subroutine ztrti2

  end module MatFunc
