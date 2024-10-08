! *****************************************************************************
! Heteroskedastic robust supLM Statistics
! Simone Giannerini + Greta Goracci 2022
! *****************************************************************************
! *****************************************************************************
SUBROUTINE ARvsTARh(x,n,d,p,a,b,neff,testv,xth,coef,s2,eps)
    USE TARMA_MOD
    use ISO_C_BINDING
    USE m_mrgrnk
    IMPLICIT NONE
    ! *****************************************************************************
    ! ARvsTARh Heteroskedastic robust supLM Statistics
    ! SG+GG 2022
    ! *****************************************************************************
    INTEGER,intent(IN):: p,d,a,b,n,neff
    integer :: i,lags(p+1),ind0,k,nr
    INTEGER::indt(neff)
    REAL(C_DOUBLE),intent(IN):: x(n)
    REAL(C_DOUBLE),intent(OUT):: testv(b-a+1,2),xth(neff),coef(p+1),s2,eps(neff)
    REAL(C_DOUBLE):: Xlag(neff,p+1),&
    dphik(neff,(p+1)),dpsik(neff,(p+1)),Ir(neff),scoreht(p+1,neff),M11is(p+1,neff),&
    epsdpsik(neff,(p+1))
    REAL(C_DOUBLE),ALLOCATABLE:: xreg(:,:),coef2(:),epst(:,:)
    REAL(C_DOUBLE):: tol,r,M11((p+1),(p+1)),M11i((p+1),(p+1)),&
    M21((p+1),(p+1)), Mi((p+1),(p+1)), score(p+1),Vi((p+1),(p+1)),&
    scoreh(p+1),cde(p+1,1),cdem(p+1,p+1),scd(p+1)!,M12((p+1),(p+1)), M22((p+1),(p+1))

    lags = [(i, i=1,p,1), 0]
    lags(p+1) = d
    nr = b-a+1
    call EMBED(x,lags,xreg)   ! X_{t} ...  X_{t-p},  X_{t-d} # (p+2) cols
!    neff = SIZE(xreg, DIM=1)
    k    = MAX(p,d)
    tol  = 1e-7
    CALL lsfit(xreg(:,2:(p+1)),xreg(:,1:1),tol,coef2,epst,1)
    eps(:) = epst(:,1)
    coef   = coef2
    s2     = SUM(eps**2)/neff
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    dphik = -Xlag                     ! neff * (p+1)
    M11   = matmul(transpose(dphik),dphik)
    CALL INVERSE(M11,p+1,M11i)
    testv(:,:) = 0
    !! computes the quantities that depend upon the threshold
    xth    = xreg(:,p+2)
    CALL mrgrnk(xth, indt) ! equivalent to R order
    ind0   = indt(a)
    r      = xth(ind0) ! first element of the threshold range
    Ir     = 0
    WHERE(xth <= r) Ir = 1 ! indicator function
    dpsik    = -Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)) ! neff * (p+1)
    epsdpsik = dpsik*SPREAD(eps,DIM=2,NCOPIES=(p+1))
!    M22      = matmul(transpose(dpsik),dpsik)
    M21      = matmul(transpose(dpsik),dphik)
!    M12      = transpose(M21)
    M11is    = matmul(M11i,transpose(dphik*SPREAD(eps,DIM=2,NCOPIES=(p+1))))
    score    = -SUM(epsdpsik,DIM=1)! (p+1) colSums
    scoreht  = -transpose(epsdpsik) + matmul(M21,M11is)
    scoreh   = SUM(scoreht,DIM=2)! rowSums (p+1)
    CALL INVERSE(matmul(scoreht,transpose(scoreht)),p+1,Vi)
!    CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
    CALL INVERSE(M21 - matmul(M21,matmul(M11i,M21)),p+1,Mi)
    testv(1,1)   = dot_product(score,matmul(Mi,score))
    testv(1,2)   = dot_product(scoreh,matmul(Vi,scoreh))

! start cycle over the threshold grid
    do i=2,nr,1 ! cycle over the threshold grid
        ind0    = indt(a+i-1)
        cde(:,1)= - Xlag(ind0,:) ! element of the derivative of eps that updates
        cdem    = matmul(cde,transpose(cde)) ! update
        scd     = cde(:,1)*eps(ind0)
        score   = score - scd
        scoreht(:,ind0) = scoreht(:,ind0) - scd
        scoreht = scoreht + matmul(cdem,M11is)
        M21     = M21 + cdem
!        M22     = M22 + cdem
!        M12     = transpose(M21)
        scoreh  = SUM(scoreht,DIM=2)! rowSums (p+1)
        CALL INVERSE(matmul(scoreht,transpose(scoreht)),p+1,Vi)
!        CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
        CALL INVERSE(M21 - matmul(M21,matmul(M11i,M21)),p+1,Mi)
        testv(i,1)   = dot_product(score,matmul(Mi,score))
        testv(i,2)   = dot_product(scoreh,matmul(Vi,scoreh))
    enddo
    testv(:,1) = testv(:,1)/s2
 END SUBROUTINE ARvsTARh
! *****************************************************************************
SUBROUTINE ARvsTAR_HB(x,n,d,p,a,b,neff,nrep,test,testb)
    USE TARMA_MOD
    use ISO_C_BINDING
    USE m_mrgrnk
    IMPLICIT NONE
    ! *****************************************************************************
    ! ARvsTARh Heteroskedastic robust supLM Statistics
    ! BOOTSTRAP version of Hansen(1996)
    ! *****************************************************************************
    INTEGER,intent(IN):: n,d,p,a,b,neff,nrep
    integer :: i,lags(p+1),ind0,k,nr
    INTEGER::indt(neff)
    REAL(C_DOUBLE),intent(IN):: x(n)
    REAL(C_DOUBLE),intent(OUT):: test(2),testb(nrep,2)
    REAL(C_DOUBLE):: Xlag(neff,p+1),testv(b-a+1,2),testvb(b-a+1,nrep,2),&
    dphik(neff,(p+1)),dpsik(neff,(p+1)),Ir(neff),scoreht(p+1,neff),M11is(p+1,neff),&
    epsdpsik(neff,(p+1)),xth(neff),coef(p+1),s2,eps(neff),scoret(p+1,neff)
    REAL(C_DOUBLE),ALLOCATABLE:: xreg(:,:),coef2(:),epst(:,:)
    REAL(C_DOUBLE):: tol,r,M11((p+1),(p+1)),M11i((p+1),(p+1)),&
    M21((p+1),(p+1)), Mi((p+1),(p+1)), score(p+1),Vi((p+1),(p+1)),&
    scoreh(p+1),cde(p+1,1),cdem(p+1,p+1),scd(p+1),&!,M12((p+1),(p+1)), M22((p+1),(p+1))
    eta(neff,nrep),dum(neff,nrep),scorehb(p+1,nrep),scoreb(p+1,nrep)

    call randnorm(dum,neff*nrep)
    eta = RESHAPE(dum,SHAPE=[neff,nrep])
    lags = [(i, i=1,p,1), 0]
    lags(p+1) = d
    nr = b-a+1
    call EMBED(x,lags,xreg)   ! X_{t} ...  X_{t-p},  X_{t-d} # (p+2) cols
!    neff = SIZE(xreg, DIM=1)
    k    = MAX(p,d)
    tol  = 1e-7
    CALL lsfit(xreg(:,2:(p+1)),xreg(:,1:1),tol,coef2,epst,1)
    eps(:) = epst(:,1)
    coef   = coef2
    s2     = SUM(eps**2)/neff
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    dphik = -Xlag                     ! neff * (p+1)
    M11   = matmul(transpose(dphik),dphik)
    CALL INVERSE(M11,p+1,M11i)
    testv(:,:) = 0
    !! computes the quantities that depend upon the threshold
    xth    = xreg(:,p+2)
    CALL mrgrnk(xth, indt) ! equivalent to R order
    ind0   = indt(a)
    r      = xth(ind0) ! first element of the threshold range
    Ir     = 0
    WHERE(xth <= r) Ir = 1 ! indicator function
    dpsik    = -Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)) ! neff * (p+1)
    epsdpsik = dpsik*SPREAD(eps,DIM=2,NCOPIES=(p+1))  ! neff * (p+1)
!    M22      = matmul(transpose(dpsik),dpsik)
    M21      = matmul(transpose(dpsik),dphik)
!    M12      = transpose(M21)
    M11is    = matmul(M11i,transpose(dphik*SPREAD(eps,DIM=2,NCOPIES=(p+1))))
    score    = -SUM(epsdpsik,DIM=1)! (p+1) colSums
    scoret   = -transpose(epsdpsik)
    scoreht  = -transpose(epsdpsik) + matmul(M21,M11is)
    scoreh   = SUM(scoreht,DIM=2)! rowSums (p+1)
    CALL INVERSE(matmul(scoreht,transpose(scoreht)),p+1,Vi)
!    CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
    CALL INVERSE(M21 - matmul(M21,matmul(M11i,M21)),p+1,Mi)
    testv(1,1)   = dot_product(score,matmul(Mi,score))
    testv(1,2)   = dot_product(scoreh,matmul(Vi,scoreh))
    scoreb       = matmul(scoret,eta)
    scorehb      = matmul(scoreht,eta)
    testvb(1,:,1)= SUM(scorehb*matmul(Mi,scorehb),DIM=1)
    testvb(1,:,2)= SUM(scorehb*matmul(Vi,scorehb),DIM=1)
! start cycle over the threshold grid
    do i=2,nr,1 ! cycle over the threshold grid
        ind0    = indt(a+i-1)
        cde(:,1)= - Xlag(ind0,:) ! element of the derivative of eps that updates
        cdem    = matmul(cde,transpose(cde)) ! update
        scd     = cde(:,1)*eps(ind0)
        score   = score - scd
        scoret(:,ind0)  = scoret(:,ind0) - scd
        scoreht(:,ind0) = scoreht(:,ind0) - scd
        scoreht = scoreht + matmul(cdem,M11is)
        scoreb  = matmul(scoret,eta)
        scorehb = matmul(scoreht,eta)
        scoreh   = SUM(scoreht,DIM=2)! rowSums (p+1)
        M21     = M21 + cdem
!        M22     = M22 + cdem
!        M12     = transpose(M21)
        CALL INVERSE(matmul(scoreht,transpose(scoreht)),p+1,Vi) ! robust
!        CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
        CALL INVERSE(M21 - matmul(M21,matmul(M11i,M21)),p+1,Mi) ! standard
        testv(i,1)   = dot_product(score,matmul(Mi,score))
        testv(i,2)   = dot_product(scoreh,matmul(Vi,scoreh))
        testvb(i,:,1)= SUM(scorehb*matmul(Mi,scorehb),DIM=1)
        testvb(i,:,2)= SUM(scorehb*matmul(Vi,scorehb),DIM=1)
    enddo
    test       = MAXVAL(testv, DIM = 1)  ! max for each col
    testb      = MAXVAL(testvb, DIM = 1) ! max for each col
    test(1)    = test(1)/s2
    testb(:,1) = testb(:,1)/s2
 END SUBROUTINE ARvsTAR_HB
! *****************************************************************************

