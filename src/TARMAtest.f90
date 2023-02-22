! *****************************************************************************
!! Set of Fortran subroutines/functions for TARMA testing
!! companion to the tseriesTARMA package
!! Simone Giannerini
!! 2021 - 
! *****************************************************************************

! *****************************************************************************

SUBROUTINE TARMAGARCH(x,eps,h,n,d,trange,nr,p,ma,q,aa,m,bb,s,testv)
    USE TARMA_MOD
    IMPLICIT NONE
 ! *****************************************************************************
 ! ARMA-GARCH vs TARMA-GARCH  supLM STATISTIC: test on both the AR and MA parameters
 ! SG 2021
 ! *****************************************************************************
    INTEGER,intent(in):: n,nr,p,d,q,m,s
    REAL(KIND=REAL64),intent(in) :: x(n),eps(n),trange(nr),ma(q),aa(m),bb(s),h(n)
    REAL(KIND=REAL64),intent(out):: testv(nr)
    integer :: i,ii,k,neff,arlags(p+m+2),malags(q+m+2),glags(s+2)
    REAL(KIND=REAL64) :: r,M11((p+q+1),(p+q+1)),M11i((p+q+1),(p+q+1)), M12((p+q+1),(p+q+1)), &
    M21((p+q+1),(p+q+1)), M22((p+q+1),(p+q+1)), Mi((p+q+1),(p+q+1)), score(p+q+1), &
    scoreH0(p+q+1),scoreH1(p+q+1)
    REAL(KIND=REAL64), ALLOCATABLE :: xreg(:,:),epsreg(:,:),hreg(:,:),epst(:), Xlag(:,:),elag(:,:),xth(:), &
    Xlagh(:,:),elagh(:,:),dedphi(:,:),dedth(:,:),dpsiar(:,:),dpsima(:,:),dedH0(:,:),ht(:), &
    dhdH0(:,:),dhdHdum(:,:),ones(:),Xdum(:,:),edum(:,:),hlag(:,:),eti(:), &
    dedphih(:,:),dedthh(:,:),dhdH02(:,:),dedH02(:,:),dhdH12(:,:),dedH12(:,:),dedH1(:,:),dhdH1(:,:)
    !,dhda0(:),dhdb(:,:),dhdai(:,:)
    INTEGER, ALLOCATABLE :: Ir(:)
    k    = MAX(p,d,q)
    arlags = [(i, i=1,p+m,1), 0,0]
    arlags(p+m+1) = d
    arlags(p+m+2) = k+m
    malags = [(i, i=1,q+m,1), 0,0]
    malags(q+m+1) = d
    malags(q+m+2) = k+m
    glags = [(i, i=1,s,1), 0,0]
    glags(s+1) = d
    glags(s+2) = k+m

    call EMBED(x,arlags,xreg)     ! X_{t}   ...     X_{t-m-p},   X_{t-d},   X_{t-k-m} # (p+m+3) cols
    call EMBED(eps,malags,epsreg) ! eps_{t} ... eps_{t-m-q}, eps_{t-d}, eps_{t-k-m} # (q+m+3) cols
    call EMBED(h,glags,hreg)      ! h_{t}   ... h_{t-s}, h_{t-d}, h_{t-k-m} # (s+3) cols
    neff = SIZE(xreg, DIM=1)

    ALLOCATE(epst(neff),xth(neff),Ir(neff),Xlag(neff,(p+1)),Xlagh(neff,(p+m+1)), &
    elag(neff,q), elagh(neff,q+m), hlag(neff,s), dedphi(neff,(p+1)),ht(neff), dedth(neff,q),&
    dhdH0(neff,(p+q+1)),dhdHdum(neff,(p+q+1)),dhdH02(neff,(p+q+1)), ones(neff),&
    dhdH12(neff,(p+q+1)),& !dhda0(neff),dhdai(neff,m),dhdb(neff,s)
    Xdum(neff,(p+1)),edum(neff,q),dedphih(neff,(p+1)),dedthh(neff,q),eti(neff), &
    dedH0(neff,(p+q+1)),dedH02(neff,(p+q+1)),dedH1(neff,(p+q+1)),dedH12(neff,(p+q+1)),&
    dhdH1(neff,(p+q+1)),dpsiar(neff,(p+1)),dpsima(neff,q))
    ones(:)            = 1
    Xlag(:,1)          = 1
    Xlag(:,2:(p+1))    = xreg(:,2:(p+1))   ! # 1, X_t-1, X_t-2 ... X_t-p
    Xlagh(:,1)         = 1
    Xlagh(:,2:(p+m+1)) = xreg(:,2:(p+m+1)) ! # 1, X_t-1, X_t-2 ... X_t-m-p
    xth(:)             = xreg(:,(p+m+2))   ! # X_t-d threshold variable
    epst               = epsreg(:,1)
    elag(:,:)          = epsreg(:,2:(q+1))   ! # eps_t-1, eps_t-2 ... eps_t-q
    elagh(:,:)         = epsreg(:,2:(q+m+1)) ! # eps_t-1, eps_t-2 ... eps_t-(m+q)
    ht(:)              = hreg(:,1)
    hlag(:,:)          = hreg(:,2:(s+1))   ! # h_t-1, ... h_t-s

    CALL rfilterm(-Xlag,ma,q,neff,(p+1),dedphi) !# d(eps_t)/d(phi)
    CALL rfilterm(elag,ma,q,neff,q,dedth)       !# d(eps_t)/d(theta)
    testv = 0
    dedH0(:,1:(p+1))       = dedphi(:,:)
    dedH0(:,(p+2):(p+q+1)) = dedth(:,:) !# partial derivatives of e_t w.r.t the parameters under H0
    dhdHdum =  0
!    CALL rfilter(ones,bb,s,neff,dhda0)      ! d(h_t)/d(a_0)
!    CALL rfilterm(hlag,bb,s,neff,s,dhdb)    ! d(h_t)/d(b_j) j = 1,...,s
!    dhdai(:,:) = 0
    do i=1,m,1
        Xdum(:,1)       = 1
        Xdum(:,2:(p+1)) = Xlagh(:,(2+i):(p+i+1))           ! 1, X_{t-i-1} ... X_{t-i-p}
        edum = elagh(:,(1+i):(q+i))                      ! eps_{t-i-1} ... eps_{t-i-q}
        eti  = elagh(:,i)                                ! eps_{t-i}
        CALL rfilterm(-Xdum,ma,q,neff,(p+1),dedphih)     ! d(eps_t-i)/d(phi)
        CALL rfilterm(edum,ma,q,neff,q,dedthh)           ! d(eps_t-i)/d(theta)
        dhdHdum = dhdHdum + 2*RESHAPE([dedphih,dedthh],[neff,(p+q+1)])*aa(i)*SPREAD(eti, DIM=2, NCOPIES=p+q+1) ! (neff * (p+q+1))
!        CALL rfilter(eti**2,bb,s,neff,dhdai(:,i))        ! d(h_t)/d(a_i)
    enddo
    CALL rfilterm(dhdHdum,bb,s,neff,(p+q+1),dhdH0)     ! # (neff * (p+q+1)) d(h_t)/d(Psi_1) partial derivatives w.r.t the ARMA parameters
    dedH02  = dedH0/SPREAD(sqrt(ht), DIM=2, NCOPIES=p+q+1)
    dhdH02  = dhdH0/SPREAD(ht, DIM=2, NCOPIES=p+q+1)
    M11     = (matmul(transpose(dedH02),dedH02)+0.5*matmul(transpose(dhdH02),dhdH02))/neff
    CALL INVERSE(M11,p+q+1,M11i)
    scoreH0  = SUM(-dedH0*SPREAD(epst/ht,DIM=2,NCOPIES=(p+q+1)) + &  ! score vector (p+q+1)
               0.5*dhdH0/SPREAD(ht,DIM=2,NCOPIES=(p+q+1))*(SPREAD(epst**2/ht,DIM=2,NCOPIES=(p+q+1))-1),DIM=1)/sqrt(dble(neff))
    do i=1,nr,1 ! cycle over the threshold grid
        r  = trange(i)
        Ir = 0
        WHERE(xth <= r) Ir = 1 ! indicator function
        CALL rfilterm(-Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)),ma,q,neff,(p+1),dpsiar) ! d(eps_t)/d(varphi)   AR part
        CALL rfilterm(elag*SPREAD(Ir, DIM=2, NCOPIES=q),ma,q,neff,q,dpsima)          ! d(eps_t)/d(vartheta) MA part
        dedH1(:,1:(p+1)) = dpsiar(:,:)
        dedH1(:,(p+2):(p+q+1)) = dpsima(:,:) ! partial derivatives of eps_t w.r.t the parameters under H1 (Psi_2)
        dhdHdum =  0
        do ii=1,m,1
            Xdum(:,1)       = 1
            Xdum(:,2:(p+1)) = Xlagh(:,(2+ii):(p+ii+1))       ! 1, X_{t-i-1} ... X_{t-i-p}
            edum = elagh(:,(1+ii):(q+ii))                      ! eps_{t-i-1} ... eps_{t-i-q}
            eti  = elagh(:,ii)                                ! eps_{t-i}
            CALL rfilterm(-Xdum*SPREAD(Ir, DIM=2, NCOPIES=(p+1)),ma,q,neff,(p+1),dedphih) ! d(eps_t-i)/d(varphi)
            CALL rfilterm(edum*SPREAD(Ir, DIM=2, NCOPIES=q),ma,q,neff,q,dedthh)           ! d(eps_t-i)/d(vartheta)
            dhdHdum = dhdHdum + 2*RESHAPE([dedphih,dedthh],[neff,(p+q+1)])*aa(ii)*SPREAD(eti, DIM=2, NCOPIES=p+q+1) ! (neff * (p+q+1))
        enddo
        CALL rfilterm(dhdHdum,bb,s,neff,(p+q+1),dhdH1)     ! # (neff * (p+q+1)) d(h_t)/d(Psi_2) partial derivatives w.r.t the tested parameters
        dedH12  = dedH1/SPREAD(sqrt(ht), DIM=2, NCOPIES=p+q+1)
        dhdH12  = dhdH1/SPREAD(ht, DIM=2, NCOPIES=p+q+1)
        M22     = (matmul(transpose(dedH12),dedH12)+0.5*matmul(transpose(dhdH12),dhdH12))/neff
        M21     = (matmul(transpose(dedH12),dedH02)+0.5*matmul(transpose(dhdH12),dhdH02))/neff
        M12     = transpose(M21)
        scoreH1 = SUM(-dedH1*SPREAD(epst/ht,DIM=2,NCOPIES=(p+q+1)) + &  ! score vector (p+q+1)
               0.5*dhdH1/SPREAD(ht,DIM=2,NCOPIES=(p+q+1))*(SPREAD(epst**2/ht,DIM=2,NCOPIES=(p+q+1))-1),DIM=1)/sqrt(dble(neff))
        score    = scoreH1 - matmul(matmul(M21,M11i),scoreH0)
        CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+q+1,Mi)
        testv(i)   = dot_product(score,matmul(Mi,score))
    enddo
 END SUBROUTINE TARMAGARCH


! *****************************************************************************

! *****************************************************************************

SUBROUTINE ARvsTAR(x,eps,n,d,trange,nr,s2,p,testv)
    USE TARMA_MOD
    IMPLICIT NONE
! *****************************************************************************
! ARvsTAR supLM STATISTICS
! SG 2020
! *****************************************************************************
    INTEGER,intent(in):: n,nr,p,d
    REAL(KIND=REAL64),intent(in) :: x(n),eps(n),trange(nr),s2
    REAL(KIND=REAL64),intent(out):: testv(nr)
    integer :: i,k,neff,lags(p+1)
    REAL(KIND=REAL64) :: r,M11((p+1),(p+1)),M11i((p+1),(p+1)), M12((p+1),(p+1)), M21((p+1),(p+1)),&
    M22((p+1),(p+1)), Mi((p+1),(p+1)), score(p+1)
    REAL(KIND=REAL64), ALLOCATABLE :: xreg(:,:),epst(:),Xlag(:,:),xth(:),dphik(:,:),dpsik(:,:)
    INTEGER, ALLOCATABLE :: Ir(:)
    lags = [(i, i=1,p,1), 0]
    lags(p+1) = d
    call EMBED(x,lags,xreg)   ! X_{t} ...  X_{t-p},  X_{t-d} # (p+2) cols
    neff = SIZE(xreg, DIM=1)
    k    = MAX(p,d)
    ALLOCATE(epst(neff),xth(neff),Ir(neff),Xlag(neff,(p+1)),dphik(neff,(p+1)),dpsik(neff,(p+1)))
    epst = eps((k+1):n)
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    xth  = xreg(:,(p+2))              ! # X_t-d threshold variable

    dphik = -Xlag                     ! neff * (p+1)
    testv = 0
    M11   = matmul(transpose(dphik),dphik)/(s2*neff)
    CALL INVERSE(M11,p+1,M11i)
    do i=1,nr,1 ! cycle over the threshold grid
        r  = trange(i)
        Ir = 0
        WHERE(xth <= r) Ir = 1 ! indicator function
        dpsik      = -Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)) ! neff * (p+1)
        score      = -SUM(dpsik*SPREAD(epst, DIM=2, NCOPIES=(p+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+1)
        M22        = matmul(transpose(dpsik),dpsik)/(s2*neff)
        M21        = matmul(transpose(dpsik),dphik)/(s2*neff)
        M12        = transpose(M21)
        CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
!        testv(i)   = RESHAPE(matmul(transpose(score),matmul(Mi,score)),[1])
        testv(i)   = dot_product(score,matmul(Mi,score))
    enddo
 END SUBROUTINE ARvsTAR
! *****************************************************************************

SUBROUTINE ARvsTARboot(x,eps,n,d,trange,nr,s2,p,testb,B)
    USE TARMA_MOD
    IMPLICIT NONE
! *********************************************************************
! ARvsTAR supLM STATISTICS - bootstrap version (Hansen Econometrica 96)
! SG 2020
! *********************************************************************
    INTEGER,intent(in):: n,nr,p,d,B
    REAL(KIND=REAL64),intent(in) :: x(n),eps(n),trange(nr),s2
    REAL(KIND=REAL64),intent(out):: testb(B)
    integer :: i,j,k,neff,lags(p+1)
    REAL(KIND=REAL64) :: r,M11((p+1),(p+1)),M11i((p+1),(p+1)), M12((p+1),(p+1)), M21((p+1),(p+1)),&
    M22((p+1),(p+1)), Mi((p+1),(p+1)), scorepsi(p+1),scorephi(p+1),score(p+1),testv(nr)
    REAL(KIND=REAL64), ALLOCATABLE :: xreg(:,:),epst(:),Xlag(:,:),xth(:),dphik(:,:),dpsik(:,:)
    REAL(KIND=REAL64), ALLOCATABLE :: Ir(:),dum(:),eta(:,:)!,eta(:,:,:)
    lags = [(i, i=1,p,1), 0]
    lags(p+1) = d
    call EMBED(x,lags,xreg)   ! X_{t} ...  X_{t-p},  X_{t-d} # (p+2) cols
    neff = SIZE(xreg, DIM=1)
    k    = MAX(p,d)
    ALLOCATE(epst(neff),xth(neff),dum(neff*B),eta(neff,B),Ir(neff)) !dum(neff*nr*B),eta(neff,nr,B)
    ALLOCATE(Xlag(neff,(p+1)),dphik(neff,(p+1)),dpsik(neff,(p+1)))
    epst = eps((k+1):n)
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    xth  = xreg(:,(p+2))              ! # X_t-d threshold variable
    dphik = -Xlag                     ! neff * (p+1)
    testv = 0
    M11   = matmul(transpose(dphik),dphik)/(s2*neff)
    CALL INVERSE(M11,p+1,M11i)
! per generare i numeri casuali in una botta sola ma in R a volte è più lento
    call randnorm(dum,neff*B)
    eta = RESHAPE(dum ,[neff,B])
    eta(:,1) = 1  ! sets the first column to 1 and get the original statistic
!    DEALLOCATE(dum)
    do j=1,B,1 ! bootstrap cycle
!        call randnorm(eta,neff)
        do i=1,nr,1 ! cycle over the threshold grid
            r  = trange(i)
            Ir = 0
            WHERE(xth <= r) Ir = 1 ! indicator function
            dpsik      = -Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)) ! neff * (p+1)
            scorepsi  = -SUM(dpsik*SPREAD(epst*eta(:,j), DIM=2, NCOPIES=(p+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+1)
            scorephi  = -SUM(dphik*SPREAD(epst*eta(:,j), DIM=2, NCOPIES=(p+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+1)
            M22        = matmul(transpose(dpsik),dpsik)/(s2*neff)
            M21        = matmul(transpose(dpsik),dphik)/(s2*neff)
            M12        = transpose(M21)
            CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
            score      = scorepsi - matmul(M21,matmul(M11i,scorephi))
            testv(i)   = dot_product(score,matmul(Mi,score))
        enddo
        testb(j) = MAXVAL(testv)
    enddo
END SUBROUTINE ARvsTARboot
! *****************************************************************************

SUBROUTINE ARMAvsTARMA(x,eps,n,d,trange,nr,s2,p,ma,q,testv)
    USE TARMA_MOD
    IMPLICIT NONE
! *****************************************************************************
! ARMAvsTARMA supLM STATISTICS
! SG 2020
! *****************************************************************************
    INTEGER,intent(in):: n,nr,p,d,q
    REAL(KIND=REAL64),intent(in) :: x(n),eps(n),trange(nr),s2,ma(q)
    REAL(KIND=REAL64),intent(out):: testv(nr)
    integer :: i,k,neff,lags(p+1)
    REAL(KIND=REAL64) :: r,M11((p+1),(p+1)),M11i((p+1),(p+1)), M12((p+1),(p+1)), M21((p+1),(p+1)),&
    M22((p+1),(p+1)), Mi((p+1),(p+1)), score(p+1)
    REAL(KIND=REAL64), ALLOCATABLE :: xreg(:,:),epst(:),Xlag(:,:),xth(:),dphik(:,:),dpsik(:,:)
    INTEGER, ALLOCATABLE :: Ir(:)
    lags = [(i, i=1,p,1), 0]
    lags(p+1) = d
    call EMBED(x,lags,xreg)   ! X_{t} ...  X_{t-p},  X_{t-d} # (p+2) cols
    neff = SIZE(xreg, DIM=1)
    k    = MAX(p,d)
    ALLOCATE(epst(neff),xth(neff),Ir(neff),Xlag(neff,(p+1)),dphik(neff,(p+1)),dpsik(neff,(p+1)))
    epst = eps((k+1):n)
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    xth  = xreg(:,(p+2))              ! # X_t-d threshold variable
    CALL rfilterm(-Xlag,ma,q,neff,(p+1),dphik)
    testv = 0
    M11   = matmul(transpose(dphik),dphik)/(s2*neff)
    CALL INVERSE(M11,p+1,M11i)
    do i=1,nr,1 ! cycle over the threshold grid
        r  = trange(i)
        Ir = 0
        WHERE(xth <= r) Ir = 1 ! indicator function
        CALL rfilterm(-Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)),ma,q,neff,(p+1),dpsik)
        score      = -SUM(dpsik*SPREAD(epst, DIM=2, NCOPIES=(p+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+1)
        M22        = matmul(transpose(dpsik),dpsik)/(s2*neff)
        M21        = matmul(transpose(dpsik),dphik)/(s2*neff)
        M12        = transpose(M21)
        CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
        testv(i)   = dot_product(score,matmul(Mi,score))
    enddo
 END SUBROUTINE ARMAvsTARMA

! *****************************************************************************

SUBROUTINE ARMAvsTARMAboot(x,eps,n,d,trange,nr,s2,p,ma,q,testb,B)
    USE TARMA_MOD
    IMPLICIT NONE
! *********************************************************************
! ARMAvsTARMA supLM STATISTICS - bootstrap version (Hansen Econometrica 96)
! SG 2020
! *********************************************************************
    INTEGER,intent(in) :: n,nr,p,d,B,q
    REAL(KIND=REAL64),intent(in) :: x(n),eps(n),trange(nr),s2,ma(q)
    REAL(KIND=REAL64),intent(out):: testb(B)
    integer :: i,j,k,neff,lags(p+1)
    REAL(KIND=REAL64) :: r,M11((p+1),(p+1)),M11i((p+1),(p+1)), M12((p+1),(p+1)), M21((p+1),(p+1)),&
    M22((p+1),(p+1)), Mi((p+1),(p+1)), scorepsi(p+1),scorephi(p+1),score(p+1),testv(nr)
    REAL(KIND=REAL64), ALLOCATABLE :: xreg(:,:),epst(:),Xlag(:,:),xth(:),dphik(:,:),dpsik(:,:)
    REAL(KIND=REAL64), ALLOCATABLE :: Ir(:),dum(:),eta(:,:)!,eta(:,:,:)
    lags = [(i, i=1,p,1), 0]
    lags(p+1) = d
    call EMBED(x,lags,xreg)   ! X_{t} ...  X_{t-p},  X_{t-d} # (p+2) cols
    neff = SIZE(xreg, DIM=1)
    k    = MAX(p,d)
    ALLOCATE(epst(neff),xth(neff),dum(neff*B),eta(neff,B),Ir(neff)) !dum(neff*nr*B),eta(neff,nr,B)
    ALLOCATE(Xlag(neff,(p+1)),dphik(neff,(p+1)),dpsik(neff,(p+1)))
    epst = eps((k+1):n)
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    xth  = xreg(:,(p+2))              ! # X_t-d threshold variable
    CALL rfilterm(-Xlag,ma,q,neff,(p+1),dphik)
    testv = 0
    M11   = matmul(transpose(dphik),dphik)/(s2*neff)
    CALL INVERSE(M11,p+1,M11i)
! per generare i numeri casuali in una botta sola ma in R a volte è più lento
    call randnorm(dum,neff*B)
    eta = RESHAPE(dum ,[neff,B])
    eta(:,1) = 1  ! sets the first column to 1 and gets the original statistic
!    DEALLOCATE(dum)
    do j=1,B,1 ! bootstrap cycle
!        call randnorm(eta,neff)
        do i=1,nr,1 ! cycle over the threshold grid
            r  = trange(i)
            Ir = 0
            WHERE(xth <= r) Ir = 1 ! indicator function
            CALL rfilterm(-Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)),ma,q,neff,(p+1),dpsik) ! neff * (p+1)
            scorepsi  = -SUM(dpsik*SPREAD(epst*eta(:,j), DIM=2, NCOPIES=(p+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+1)
            scorephi  = -SUM(dphik*SPREAD(epst*eta(:,j), DIM=2, NCOPIES=(p+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+1)
            M22        = matmul(transpose(dpsik),dpsik)/(s2*neff)
            M21        = matmul(transpose(dpsik),dphik)/(s2*neff)
            M12        = transpose(M21)
            CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+1,Mi)
            score      = scorepsi - matmul(M21,matmul(M11i,scorephi))
            testv(i)   = dot_product(score,matmul(Mi,score))
        enddo
        testb(j) = MAXVAL(testv)
    enddo
END SUBROUTINE ARMAvsTARMAboot
! *****************************************************************************

SUBROUTINE ARMAvsTARMAg(x,eps,n,d,trange,nr,s2,p,ma,q,testv)
    USE TARMA_MOD
    IMPLICIT NONE
! *****************************************************************************
! ARMAvsTARMAg supLM STATISTICS: test on both the AR and MA parameters
! SG 2021
! *****************************************************************************
    INTEGER,intent(in):: n,nr,p,d,q
    REAL(KIND=REAL64),intent(in) :: x(n),eps(n),trange(nr),s2,ma(q)
    REAL(KIND=REAL64),intent(out):: testv(nr)
    integer :: i,k,neff,arlags(p+2),malags(q+2)
    REAL(KIND=REAL64) :: r,M11((p+q+1),(p+q+1)),M11i((p+q+1),(p+q+1)), M12((p+q+1),(p+q+1)), M21((p+q+1),(p+q+1)),&
    M22((p+q+1),(p+q+1)), Mi((p+q+1),(p+q+1)), score(p+q+1)
    REAL(KIND=REAL64), ALLOCATABLE :: xreg(:,:),epsreg(:,:),epst(:),Xlag(:,:),epslag(:,:),xth(:),&
    dphik(:,:),dths(:,:),dpsiar(:,:),dpsima(:,:),dedH0(:,:),dedH1(:,:)
    INTEGER, ALLOCATABLE :: Ir(:)
    k    = MAX(p,d,q)
    neff = n-k
    arlags = [(i, i=1,p,1), 0,0]
    arlags(p+1) = d
    arlags(p+2) = k
    malags = [(i, i=1,q,1), 0,0]
    malags(q+1) = d
    malags(q+2) = k
    call EMBED(x,arlags,xreg)     ! X_{t} ...     X_{t-p},   X_{t-d},   X_{t-k} # (p+3) cols
    call EMBED(eps,malags,epsreg) ! eps_{t} ... eps_{t-q}, eps_{t-d}, eps_{t-k} # (q+3) cols
    ALLOCATE(epst(neff),xth(neff),Ir(neff),Xlag(neff,(p+1)),epslag(neff,q),dphik(neff,(p+1)),&
    dths(neff,q),dpsiar(neff,(p+1)),dpsima(neff,q),dedH0(neff,(p+q+1)),dedH1(neff,(p+q+1)))
    epst = epsreg(:,1)
    epslag(:,:) = epsreg(:,2:(q+1)) ! # eps_t-1, eps_t-2 ... eps_t-q
    Xlag(:,1) = 1
    Xlag(:,2:(p+1)) = xreg(:,2:(p+1)) ! # 1, X_t-1, X_t-2 ... X_t-p
    xth  = xreg(:,(p+2))              ! # X_t-d threshold variable
    CALL rfilterm(-Xlag,ma,q,neff,(p+1),dphik) !# d(eps_t)/d(phi_k)
    CALL rfilterm(epslag,ma,q,neff,q,dths) !# d(eps_t)/d(theta_s)
    testv = 0
    dedH0(:,1:(p+1))       = dphik(:,:)
    dedH0(:,(p+2):(p+q+1)) = dths(:,:)
    M11   = matmul(transpose(dedH0),dedH0)/(s2*neff)
    CALL INVERSE(M11,p+q+1,M11i)
    do i=1,nr,1 ! cycle over the threshold grid
        r  = trange(i)
        Ir = 0
        WHERE(xth <= r) Ir = 1 ! indicator function
        CALL rfilterm(-Xlag*SPREAD(Ir, DIM=2, NCOPIES=(p+1)),ma,q,neff,(p+1),dpsiar)  !# d(eps_t)/d(Psi_k) k=0,1,..,p
        CALL rfilterm(epslag*SPREAD(Ir, DIM=2, NCOPIES=q),ma,q,neff,q,dpsima)         !# d(eps_t)/d(Psi_s) k=1,2,..,q
        dedH1(:,1:(p+1)) = dpsiar(:,:)
        dedH1(:,(p+2):(p+q+1)) = dpsima(:,:)
        score      = -SUM(dedH1*SPREAD(epst, DIM=2, NCOPIES=(p+q+1)),DIM=1)/(sqrt(dble(neff))*s2) ! (p+q+1)
        M22        = matmul(transpose(dedH1),dedH1)/(s2*neff)
        M21        = matmul(transpose(dedH1),dedH0)/(s2*neff)
        M12        = transpose(M21)
        CALL INVERSE(M22 - matmul(M21,matmul(M11i,M12)),p+q+1,Mi)
        testv(i)   = dot_product(score,matmul(Mi,score))
    enddo
 END SUBROUTINE ARMAvsTARMAg

! *****************************************************************************
