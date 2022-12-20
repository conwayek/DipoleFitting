*> \brief <b> DGELSS solves overdetermined or underdetermined systems for GE matrices</b>
 *
 *  =========== DOCUMENTATION ===========
 *
 * Online html documentation available at
 *            http://www.netlib.org/lapack/explore-html/
 *
 *> \htmlonly
 *> Download DGELSS + dependencies
 *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgelss.f">
 *> [TGZ]</a>
 *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgelss.f">
 *> [ZIP]</a>
 *> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgelss.f">
 *> [TXT]</a>
 *> \endhtmlonly
 *
 *  Definition:
 *  ===========
 *
 *       SUBROUTINE DGELSS( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
 *                          WORK, LWORK, INFO )
 *
 *       .. Scalar Arguments ..
 *       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
 *       DOUBLE PRECISION   RCOND
 *       ..
 *       .. Array Arguments ..
 *       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
 *       ..
 *
 *
 *> \par Purpose:
 *  =============
 *>
 *> \verbatim
 *>
 *> DGELSS computes the minimum norm solution to a real linear least
 *> squares problem:
 *>
 *> Minimize 2-norm(| b - A*x |).
 *>
 *> using the singular value decomposition (SVD) of A. A is an M-by-N
 *> matrix which may be rank-deficient.
 *>
 *> Several right hand side vectors b and solution vectors x can be
 *> handled in a single call; they are stored as the columns of the
 *> M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
 *> X.
 *>
 *> The effective rank of A is determined by treating as zero those
 *> singular values which are less than RCOND times the largest singular
 *> value.
 *> \endverbatim
 *
 *  Arguments:
 *  ==========
 *
 *> \param[in] M
 *> \verbatim
 *>          M is INTEGER
 *>          The number of rows of the matrix A. M >= 0.
 *> \endverbatim
 *>
 *> \param[in] N
 *> \verbatim
 *>          N is INTEGER
 *>          The number of columns of the matrix A. N >= 0.
 *> \endverbatim
 *>
 *> \param[in] NRHS
 *> \verbatim
 *>          NRHS is INTEGER
 *>          The number of right hand sides, i.e., the number of columns
 *>          of the matrices B and X. NRHS >= 0.
 *> \endverbatim
 *>
 *> \param[in,out] A
 *> \verbatim
 *>          A is DOUBLE PRECISION array, dimension (LDA,N)
 *>          On entry, the M-by-N matrix A.
 *>          On exit, the first min(m,n) rows of A are overwritten with
 *>          its right singular vectors, stored rowwise.
 *> \endverbatim
 *>
 *> \param[in] LDA
 *> \verbatim
 *>          LDA is INTEGER
 *>          The leading dimension of the array A.  LDA >= max(1,M).
 *> \endverbatim
 *>
 *> \param[in,out] B
 *> \verbatim
 *>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
 *>          On entry, the M-by-NRHS right hand side matrix B.
 *>          On exit, B is overwritten by the N-by-NRHS solution
 *>          matrix X.  If m >= n and RANK = n, the residual
 *>          sum-of-squares for the solution in the i-th column is given
 *>          by the sum of squares of elements n+1:m in that column.
 *> \endverbatim
 *>
 *> \param[in] LDB
 *> \verbatim
 *>          LDB is INTEGER
 *>          The leading dimension of the array B. LDB >= max(1,max(M,N)).
 *> \endverbatim
 *>
 *> \param[out] S
 *> \verbatim
 *>          S is DOUBLE PRECISION array, dimension (min(M,N))
 *>          The singular values of A in decreasing order.
 *>          The condition number of A in the 2-norm = S(1)/S(min(m,n)).
 *> \endverbatim
 *>
 *> \param[in] RCOND
 *> \verbatim
 *>          RCOND is DOUBLE PRECISION
 *>          RCOND is used to determine the effective rank of A.
 *>          Singular values S(i) <= RCOND*S(1) are treated as zero.
 *>          If RCOND < 0, machine precision is used instead.
 *> \endverbatim
 *>
 *> \param[out] RANK
 *> \verbatim
 *>          RANK is INTEGER
 *>          The effective rank of A, i.e., the number of singular values
 *>          which are greater than RCOND*S(1).
 *> \endverbatim
 *>
 *> \param[out] WORK
 *> \verbatim
 *>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
 *>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *> \endverbatim
 *>
 *> \param[in] LWORK
 *> \verbatim
 *>          LWORK is INTEGER
 *>          The dimension of the array WORK. LWORK >= 1, and also:
 *>          LWORK >= 3*min(M,N) + max( 2*min(M,N), max(M,N), NRHS )
 *>          For good performance, LWORK should generally be larger.
 *>
 *>          If LWORK = -1, then a workspace query is assumed; the routine
 *>          only calculates the optimal size of the WORK array, returns
 *>          this value as the first entry of the WORK array, and no error
 *>          message related to LWORK is issued by XERBLA.
 *> \endverbatim
 *>
 *> \param[out] INFO
 *> \verbatim
 *>          INFO is INTEGER
 *>          = 0:  successful exit
 *>          < 0:  if INFO = -i, the i-th argument had an illegal value.
 *>          > 0:  the algorithm for computing the SVD failed to converge;
 *>                if INFO = i, i off-diagonal elements of an intermediate
 *>                bidiagonal form did not converge to zero.
 *> \endverbatim
 *
 *  Authors:
 *  ========
 *
 *> \author Univ. of Tennessee
 *> \author Univ. of California Berkeley
 *> \author Univ. of Colorado Denver
 *> \author NAG Ltd.
 *
 *> \ingroup doubleGEsolve
 *
 *  =====================================================================
       SUBROUTINE dgelss( M, N, NRHS, A, LDA, B, LDB, S, RCOND, RANK,
      $                   WORK, LWORK, INFO )
 *
 *  -- LAPACK driver routine --
 *  -- LAPACK is a software package provided by Univ. of Tennessee,    --
 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 *
 *     .. Scalar Arguments ..
       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS, RANK
       DOUBLE PRECISION   RCOND
 *     ..
 *     .. Array Arguments ..
       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), S( * ), WORK( * )
 *     ..
 *
 *  =====================================================================
 *
 *     .. Parameters ..
       DOUBLE PRECISION   ZERO, ONE
       parameter( zero = 0.0d+0, one = 1.0d+0 )
 *     ..
 *     .. Local Scalars ..
       LOGICAL            LQUERY
       INTEGER            BDSPAC, BL, CHUNK, I, IASCL, IBSCL, IE, IL,
      $                   itau, itaup, itauq, iwork, ldwork, maxmn,
      $                   maxwrk, minmn, minwrk, mm, mnthr
       INTEGER            LWORK_DGEQRF, LWORK_DORMQR, LWORK_DGEBRD,
      $                   lwork_dormbr, lwork_dorgbr, lwork_dormlq,
      $                   lwork_dgelqf
       DOUBLE PRECISION   ANRM, BIGNUM, BNRM, EPS, SFMIN, SMLNUM, THR
 *     ..
 *     .. Local Arrays ..
       DOUBLE PRECISION   DUM( 1 )
 *     ..
 *     .. External Subroutines ..
       EXTERNAL           dbdsqr, dcopy, dgebrd, dgelqf, dgemm, dgemv,
      $                   dgeqrf, dlabad, dlacpy, dlascl, dlaset, dorgbr,
      $                   dormbr, dormlq, dormqr, drscl, xerbla
 *     ..
 *     .. External Functions ..
       INTEGER            ILAENV
       DOUBLE PRECISION   DLAMCH, DLANGE
       EXTERNAL           ilaenv, dlamch, dlange
 *     ..
 *     .. Intrinsic Functions ..
       INTRINSIC          max, min
 *     ..
 *     .. Executable Statements ..
 *
 *     Test the input arguments
 *
       info = 0
       minmn = min( m, n )
       maxmn = max( m, n )
       lquery = ( lwork.EQ.-1 )
       IF( m.LT.0 ) THEN
          info = -1
       ELSE IF( n.LT.0 ) THEN
          info = -2
       ELSE IF( nrhs.LT.0 ) THEN
          info = -3
       ELSE IF( lda.LT.max( 1, m ) ) THEN
          info = -5
       ELSE IF( ldb.LT.max( 1, maxmn ) ) THEN
          info = -7
       END IF
 *
 *     Compute workspace
 *      (Note: Comments in the code beginning "Workspace:" describe the
 *       minimal amount of workspace needed at that point in the code,
 *       as well as the preferred amount for good performance.
 *       NB refers to the optimal block size for the immediately
 *       following subroutine, as returned by ILAENV.)
 *
       IF( info.EQ.0 ) THEN
          minwrk = 1
          maxwrk = 1
          IF( minmn.GT.0 ) THEN
             mm = m
             mnthr = ilaenv( 6, 'DGELSS', ' ', m, n, nrhs, -1 )
             IF( m.GE.n .AND. m.GE.mnthr ) THEN
 *
 *              Path 1a - overdetermined, with many more rows than
 *                        columns
 *
 *              Compute space needed for DGEQRF
                CALL dgeqrf( m, n, a, lda, dum(1), dum(1), -1, info )
                lwork_dgeqrf=dum(1)
 *              Compute space needed for DORMQR
                CALL dormqr( 'L', 'T', m, nrhs, n, a, lda, dum(1), b,
      $                   ldb, dum(1), -1, info )
                lwork_dormqr=dum(1)
                mm = n
                maxwrk = max( maxwrk, n + lwork_dgeqrf )
                maxwrk = max( maxwrk, n + lwork_dormqr )
             END IF
             IF( m.GE.n ) THEN
 *
 *              Path 1 - overdetermined or exactly determined
 *
 *              Compute workspace needed for DBDSQR
 *
                bdspac = max( 1, 5*n )
 *              Compute space needed for DGEBRD
                CALL dgebrd( mm, n, a, lda, s, dum(1), dum(1),
      $                      dum(1), dum(1), -1, info )
                lwork_dgebrd=dum(1)
 *              Compute space needed for DORMBR
                CALL dormbr( 'Q', 'L', 'T', mm, nrhs, n, a, lda, dum(1),
      $                b, ldb, dum(1), -1, info )
                lwork_dormbr=dum(1)
 *              Compute space needed for DORGBR
                CALL dorgbr( 'P', n, n, n, a, lda, dum(1),
      $                   dum(1), -1, info )
                lwork_dorgbr=dum(1)
 *              Compute total workspace needed
                maxwrk = max( maxwrk, 3*n + lwork_dgebrd )
                maxwrk = max( maxwrk, 3*n + lwork_dormbr )
                maxwrk = max( maxwrk, 3*n + lwork_dorgbr )
                maxwrk = max( maxwrk, bdspac )
                maxwrk = max( maxwrk, n*nrhs )
                minwrk = max( 3*n + mm, 3*n + nrhs, bdspac )
                maxwrk = max( minwrk, maxwrk )
             END IF
             IF( n.GT.m ) THEN
 *
 *              Compute workspace needed for DBDSQR
 *
                bdspac = max( 1, 5*m )
                minwrk = max( 3*m+nrhs, 3*m+n, bdspac )
                IF( n.GE.mnthr ) THEN
 *
 *                 Path 2a - underdetermined, with many more columns
 *                 than rows
 *
 *                 Compute space needed for DGELQF
                   CALL dgelqf( m, n, a, lda, dum(1), dum(1),
      $                -1, info )
                   lwork_dgelqf=dum(1)
 *                 Compute space needed for DGEBRD
                   CALL dgebrd( m, m, a, lda, s, dum(1), dum(1),
      $                      dum(1), dum(1), -1, info )
                   lwork_dgebrd=dum(1)
 *                 Compute space needed for DORMBR
                   CALL dormbr( 'Q', 'L', 'T', m, nrhs, n, a, lda,
      $                dum(1), b, ldb, dum(1), -1, info )
                   lwork_dormbr=dum(1)
 *                 Compute space needed for DORGBR
                   CALL dorgbr( 'P', m, m, m, a, lda, dum(1),
      $                   dum(1), -1, info )
                   lwork_dorgbr=dum(1)
 *                 Compute space needed for DORMLQ
                   CALL dormlq( 'L', 'T', n, nrhs, m, a, lda, dum(1),
      $                 b, ldb, dum(1), -1, info )
                   lwork_dormlq=dum(1)
 *                 Compute total workspace needed
                   maxwrk = m + lwork_dgelqf
                   maxwrk = max( maxwrk, m*m + 4*m + lwork_dgebrd )
                   maxwrk = max( maxwrk, m*m + 4*m + lwork_dormbr )
                   maxwrk = max( maxwrk, m*m + 4*m + lwork_dorgbr )
                   maxwrk = max( maxwrk, m*m + m + bdspac )
                   IF( nrhs.GT.1 ) THEN
                      maxwrk = max( maxwrk, m*m + m + m*nrhs )
                   ELSE
                      maxwrk = max( maxwrk, m*m + 2*m )
                   END IF
                   maxwrk = max( maxwrk, m + lwork_dormlq )
                ELSE
 *
 *                 Path 2 - underdetermined
 *
 *                 Compute space needed for DGEBRD
                   CALL dgebrd( m, n, a, lda, s, dum(1), dum(1),
      $                      dum(1), dum(1), -1, info )
                   lwork_dgebrd=dum(1)
 *                 Compute space needed for DORMBR
                   CALL dormbr( 'Q', 'L', 'T', m, nrhs, m, a, lda,
      $                dum(1), b, ldb, dum(1), -1, info )
                   lwork_dormbr=dum(1)
 *                 Compute space needed for DORGBR
                   CALL dorgbr( 'P', m, n, m, a, lda, dum(1),
      $                   dum(1), -1, info )
                   lwork_dorgbr=dum(1)
                   maxwrk = 3*m + lwork_dgebrd
                   maxwrk = max( maxwrk, 3*m + lwork_dormbr )
                   maxwrk = max( maxwrk, 3*m + lwork_dorgbr )
                   maxwrk = max( maxwrk, bdspac )
                   maxwrk = max( maxwrk, n*nrhs )
                END IF
             END IF
             maxwrk = max( minwrk, maxwrk )
          END IF
          work( 1 ) = maxwrk
 *
          IF( lwork.LT.minwrk .AND. .NOT.lquery )
      $      info = -12
       END IF
 *
       IF( info.NE.0 ) THEN
          CALL xerbla( 'DGELSS', -info )
          RETURN
       ELSE IF( lquery ) THEN
          RETURN
       END IF
 *
 *     Quick return if possible
 *
       IF( m.EQ.0 .OR. n.EQ.0 ) THEN
          rank = 0
          RETURN
       END IF
 *
 *     Get machine parameters
 *
       eps = dlamch( 'P' )
       sfmin = dlamch( 'S' )
       smlnum = sfmin / eps
       bignum = one / smlnum
       CALL dlabad( smlnum, bignum )
 *
 *     Scale A if max element outside range [SMLNUM,BIGNUM]
 *
       anrm = dlange( 'M', m, n, a, lda, work )
       iascl = 0
       IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
 *
 *        Scale matrix norm up to SMLNUM
 *
          CALL dlascl( 'G', 0, 0, anrm, smlnum, m, n, a, lda, info )
          iascl = 1
       ELSE IF( anrm.GT.bignum ) THEN
 *
 *        Scale matrix norm down to BIGNUM
 *
          CALL dlascl( 'G', 0, 0, anrm, bignum, m, n, a, lda, info )
          iascl = 2
       ELSE IF( anrm.EQ.zero ) THEN
 *
 *        Matrix all zero. Return zero solution.
 *
          CALL dlaset( 'F', max( m, n ), nrhs, zero, zero, b, ldb )
          CALL dlaset( 'F', minmn, 1, zero, zero, s, minmn )
          rank = 0
          GO TO 70
       END IF
 *
 *     Scale B if max element outside range [SMLNUM,BIGNUM]
 *
       bnrm = dlange( 'M', m, nrhs, b, ldb, work )
       ibscl = 0
       IF( bnrm.GT.zero .AND. bnrm.LT.smlnum ) THEN
 *
 *        Scale matrix norm up to SMLNUM
 *
          CALL dlascl( 'G', 0, 0, bnrm, smlnum, m, nrhs, b, ldb, info )
          ibscl = 1
       ELSE IF( bnrm.GT.bignum ) THEN
 *
 *        Scale matrix norm down to BIGNUM
 *
          CALL dlascl( 'G', 0, 0, bnrm, bignum, m, nrhs, b, ldb, info )
          ibscl = 2
       END IF
 *
 *     Overdetermined case
 *
       IF( m.GE.n ) THEN
 *
 *        Path 1 - overdetermined or exactly determined
 *
          mm = m
          IF( m.GE.mnthr ) THEN
 *
 *           Path 1a - overdetermined, with many more rows than columns
 *
             mm = n
             itau = 1
             iwork = itau + n
 *
 *           Compute A=Q*R
 *           (Workspace: need 2*N, prefer N+N*NB)
 *
             CALL dgeqrf( m, n, a, lda, work( itau ), work( iwork ),
      $                   lwork-iwork+1, info )
 *
 *           Multiply B by transpose(Q)
 *           (Workspace: need N+NRHS, prefer N+NRHS*NB)
 *
             CALL dormqr( 'L', 'T', m, nrhs, n, a, lda, work( itau ), b,
      $                   ldb, work( iwork ), lwork-iwork+1, info )
 *
 *           Zero out below R
 *
             IF( n.GT.1 )
      $         CALL dlaset( 'L', n-1, n-1, zero, zero, a( 2, 1 ), lda )
          END IF
 *
          ie = 1
          itauq = ie + n
          itaup = itauq + n
          iwork = itaup + n
 *
 *        Bidiagonalize R in A
 *        (Workspace: need 3*N+MM, prefer 3*N+(MM+N)*NB)
 *
          CALL dgebrd( mm, n, a, lda, s, work( ie ), work( itauq ),
      $                work( itaup ), work( iwork ), lwork-iwork+1,
      $                info )
 *
 *        Multiply B by transpose of left bidiagonalizing vectors of R
 *        (Workspace: need 3*N+NRHS, prefer 3*N+NRHS*NB)
 *
          CALL dormbr( 'Q', 'L', 'T', mm, nrhs, n, a, lda, work( itauq ),
      $                b, ldb, work( iwork ), lwork-iwork+1, info )
 *
 *        Generate right bidiagonalizing vectors of R in A
 *        (Workspace: need 4*N-1, prefer 3*N+(N-1)*NB)
 *
          CALL dorgbr( 'P', n, n, n, a, lda, work( itaup ),
      $                work( iwork ), lwork-iwork+1, info )
          iwork = ie + n
 *
 *        Perform bidiagonal QR iteration
 *          multiply B by transpose of left singular vectors
 *          compute right singular vectors in A
 *        (Workspace: need BDSPAC)
 *
          CALL dbdsqr( 'U', n, n, 0, nrhs, s, work( ie ), a, lda, dum,
      $                1, b, ldb, work( iwork ), info )
          IF( info.NE.0 )
      $      GO TO 70
 *
 *        Multiply B by reciprocals of singular values
 *
          thr = max( rcond*s( 1 ), sfmin )
          IF( rcond.LT.zero )
      $      thr = max( eps*s( 1 ), sfmin )
          rank = 0
          DO 10 i = 1, n
             IF( s( i ).GT.thr ) THEN
                CALL drscl( nrhs, s( i ), b( i, 1 ), ldb )
                rank = rank + 1
             ELSE
                CALL dlaset( 'F', 1, nrhs, zero, zero, b( i, 1 ), ldb )
             END IF
    10    CONTINUE
 *
 *        Multiply B by right singular vectors
 *        (Workspace: need N, prefer N*NRHS)
 *
          IF( lwork.GE.ldb*nrhs .AND. nrhs.GT.1 ) THEN
             CALL dgemm( 'T', 'N', n, nrhs, n, one, a, lda, b, ldb, zero,
      $                  work, ldb )
             CALL dlacpy( 'G', n, nrhs, work, ldb, b, ldb )
          ELSE IF( nrhs.GT.1 ) THEN
             chunk = lwork / n
             DO 20 i = 1, nrhs, chunk
                bl = min( nrhs-i+1, chunk )
                CALL dgemm( 'T', 'N', n, bl, n, one, a, lda, b( 1, i ),
      $                     ldb, zero, work, n )
                CALL dlacpy( 'G', n, bl, work, n, b( 1, i ), ldb )
    20       CONTINUE
          ELSE
             CALL dgemv( 'T', n, n, one, a, lda, b, 1, zero, work, 1 )
             CALL dcopy( n, work, 1, b, 1 )
          END IF
 *
       ELSE IF( n.GE.mnthr .AND. lwork.GE.4*m+m*m+
      $         max( m, 2*m-4, nrhs, n-3*m ) ) THEN
 *
 *        Path 2a - underdetermined, with many more columns than rows
 *        and sufficient workspace for an efficient algorithm
 *
          ldwork = m
          IF( lwork.GE.max( 4*m+m*lda+max( m, 2*m-4, nrhs, n-3*m ),
      $       m*lda+m+m*nrhs ) )ldwork = lda
          itau = 1
          iwork = m + 1
 *
 *        Compute A=L*Q
 *        (Workspace: need 2*M, prefer M+M*NB)
 *
          CALL dgelqf( m, n, a, lda, work( itau ), work( iwork ),
      $                lwork-iwork+1, info )
          il = iwork
 *
 *        Copy L to WORK(IL), zeroing out above it
 *
          CALL dlacpy( 'L', m, m, a, lda, work( il ), ldwork )
          CALL dlaset( 'U', m-1, m-1, zero, zero, work( il+ldwork ),
      $                ldwork )
          ie = il + ldwork*m
          itauq = ie + m
          itaup = itauq + m
          iwork = itaup + m
 *
 *        Bidiagonalize L in WORK(IL)
 *        (Workspace: need M*M+5*M, prefer M*M+4*M+2*M*NB)
 *
          CALL dgebrd( m, m, work( il ), ldwork, s, work( ie ),
      $                work( itauq ), work( itaup ), work( iwork ),
      $                lwork-iwork+1, info )
 *
 *        Multiply B by transpose of left bidiagonalizing vectors of L
 *        (Workspace: need M*M+4*M+NRHS, prefer M*M+4*M+NRHS*NB)
 *
          CALL dormbr( 'Q', 'L', 'T', m, nrhs, m, work( il ), ldwork,
      $                work( itauq ), b, ldb, work( iwork ),
      $                lwork-iwork+1, info )
 *
 *        Generate right bidiagonalizing vectors of R in WORK(IL)
 *        (Workspace: need M*M+5*M-1, prefer M*M+4*M+(M-1)*NB)
 *
          CALL dorgbr( 'P', m, m, m, work( il ), ldwork, work( itaup ),
      $                work( iwork ), lwork-iwork+1, info )
          iwork = ie + m
 *
 *        Perform bidiagonal QR iteration,
 *           computing right singular vectors of L in WORK(IL) and
 *           multiplying B by transpose of left singular vectors
 *        (Workspace: need M*M+M+BDSPAC)
 *
          CALL dbdsqr( 'U', m, m, 0, nrhs, s, work( ie ), work( il ),
      $                ldwork, a, lda, b, ldb, work( iwork ), info )
          IF( info.NE.0 )
      $      GO TO 70
 *
 *        Multiply B by reciprocals of singular values
 *
          thr = max( rcond*s( 1 ), sfmin )
          IF( rcond.LT.zero )
      $      thr = max( eps*s( 1 ), sfmin )
          rank = 0
          DO 30 i = 1, m
             IF( s( i ).GT.thr ) THEN
                CALL drscl( nrhs, s( i ), b( i, 1 ), ldb )
                rank = rank + 1
             ELSE
                CALL dlaset( 'F', 1, nrhs, zero, zero, b( i, 1 ), ldb )
             END IF
    30    CONTINUE
          iwork = ie
 *
 *        Multiply B by right singular vectors of L in WORK(IL)
 *        (Workspace: need M*M+2*M, prefer M*M+M+M*NRHS)
 *
          IF( lwork.GE.ldb*nrhs+iwork-1 .AND. nrhs.GT.1 ) THEN
             CALL dgemm( 'T', 'N', m, nrhs, m, one, work( il ), ldwork,
      $                  b, ldb, zero, work( iwork ), ldb )
             CALL dlacpy( 'G', m, nrhs, work( iwork ), ldb, b, ldb )
          ELSE IF( nrhs.GT.1 ) THEN
             chunk = ( lwork-iwork+1 ) / m
             DO 40 i = 1, nrhs, chunk
                bl = min( nrhs-i+1, chunk )
                CALL dgemm( 'T', 'N', m, bl, m, one, work( il ), ldwork,
      $                     b( 1, i ), ldb, zero, work( iwork ), m )
                CALL dlacpy( 'G', m, bl, work( iwork ), m, b( 1, i ),
      $                      ldb )
    40       CONTINUE
          ELSE
             CALL dgemv( 'T', m, m, one, work( il ), ldwork, b( 1, 1 ),
      $                  1, zero, work( iwork ), 1 )
             CALL dcopy( m, work( iwork ), 1, b( 1, 1 ), 1 )
          END IF
 *
 *        Zero out below first M rows of B
 *
          CALL dlaset( 'F', n-m, nrhs, zero, zero, b( m+1, 1 ), ldb )
          iwork = itau + m
 *
 *        Multiply transpose(Q) by B
 *        (Workspace: need M+NRHS, prefer M+NRHS*NB)
 *
          CALL dormlq( 'L', 'T', n, nrhs, m, a, lda, work( itau ), b,
      $                ldb, work( iwork ), lwork-iwork+1, info )
 *
       ELSE
 *
 *        Path 2 - remaining underdetermined cases
 *
          ie = 1
          itauq = ie + m
          itaup = itauq + m
          iwork = itaup + m
 *
 *        Bidiagonalize A
 *        (Workspace: need 3*M+N, prefer 3*M+(M+N)*NB)
 *
          CALL dgebrd( m, n, a, lda, s, work( ie ), work( itauq ),
      $                work( itaup ), work( iwork ), lwork-iwork+1,
      $                info )
 *
 *        Multiply B by transpose of left bidiagonalizing vectors
 *        (Workspace: need 3*M+NRHS, prefer 3*M+NRHS*NB)
 *
          CALL dormbr( 'Q', 'L', 'T', m, nrhs, n, a, lda, work( itauq ),
      $                b, ldb, work( iwork ), lwork-iwork+1, info )
 *
 *        Generate right bidiagonalizing vectors in A
 *        (Workspace: need 4*M, prefer 3*M+M*NB)
 *
          CALL dorgbr( 'P', m, n, m, a, lda, work( itaup ),
      $                work( iwork ), lwork-iwork+1, info )
          iwork = ie + m
 *
 *        Perform bidiagonal QR iteration,
 *           computing right singular vectors of A in A and
 *           multiplying B by transpose of left singular vectors
 *        (Workspace: need BDSPAC)
 *
          CALL dbdsqr( 'L', m, n, 0, nrhs, s, work( ie ), a, lda, dum,
      $                1, b, ldb, work( iwork ), info )
          IF( info.NE.0 )
      $      GO TO 70
 *
 *        Multiply B by reciprocals of singular values
 *
          thr = max( rcond*s( 1 ), sfmin )
          IF( rcond.LT.zero )
      $      thr = max( eps*s( 1 ), sfmin )
          rank = 0
          DO 50 i = 1, m
             IF( s( i ).GT.thr ) THEN
                CALL drscl( nrhs, s( i ), b( i, 1 ), ldb )
                rank = rank + 1
             ELSE
                CALL dlaset( 'F', 1, nrhs, zero, zero, b( i, 1 ), ldb )
             END IF
    50    CONTINUE
 *
 *        Multiply B by right singular vectors of A
 *        (Workspace: need N, prefer N*NRHS)
 *
          IF( lwork.GE.ldb*nrhs .AND. nrhs.GT.1 ) THEN
             CALL dgemm( 'T', 'N', n, nrhs, m, one, a, lda, b, ldb, zero,
      $                  work, ldb )
             CALL dlacpy( 'F', n, nrhs, work, ldb, b, ldb )
          ELSE IF( nrhs.GT.1 ) THEN
             chunk = lwork / n
             DO 60 i = 1, nrhs, chunk
                bl = min( nrhs-i+1, chunk )
                CALL dgemm( 'T', 'N', n, bl, m, one, a, lda, b( 1, i ),
      $                     ldb, zero, work, n )
                CALL dlacpy( 'F', n, bl, work, n, b( 1, i ), ldb )
    60       CONTINUE
          ELSE
             CALL dgemv( 'T', m, n, one, a, lda, b, 1, zero, work, 1 )
             CALL dcopy( n, work, 1, b, 1 )
          END IF
       END IF
 *
 *     Undo scaling
 *
       IF( iascl.EQ.1 ) THEN
          CALL dlascl( 'G', 0, 0, anrm, smlnum, n, nrhs, b, ldb, info )
          CALL dlascl( 'G', 0, 0, smlnum, anrm, minmn, 1, s, minmn,
      $                info )
       ELSE IF( iascl.EQ.2 ) THEN
          CALL dlascl( 'G', 0, 0, anrm, bignum, n, nrhs, b, ldb, info )
          CALL dlascl( 'G', 0, 0, bignum, anrm, minmn, 1, s, minmn,
      $                info )
       END IF
       IF( ibscl.EQ.1 ) THEN
          CALL dlascl( 'G', 0, 0, smlnum, bnrm, n, nrhs, b, ldb, info )
       ELSE IF( ibscl.EQ.2 ) THEN
          CALL dlascl( 'G', 0, 0, bignum, bnrm, n, nrhs, b, ldb, info )
       END IF
 *
    70 CONTINUE
       work( 1 ) = maxwrk
       RETURN
 *
 *     End of DGELSS
 *
       END
