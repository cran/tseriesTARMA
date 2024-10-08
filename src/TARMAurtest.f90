! *****************************************************************************
  SUBROUTINE IMAvsTARMA(x,eps,n,trange,nr,s2,ma,testv)
  USE TARMA_MOD
  use ISO_C_BINDING
  IMPLICIT NONE
  ! *****************************************************************************
    ! IMAvsTARMA supLM STATISTICS
    ! SG 2021
    ! *****************************************************************************
    INTEGER,intent(in):: n,nr
REAL(C_DOUBLE),intent(in) :: x(n),eps(n),trange(nr),s2,ma(1)
  REAL(C_DOUBLE),intent(out):: testv(nr)
  integer :: i,neff
REAL(C_DOUBLE) :: r,M11(1,1), M21(2,1),M22(2,2),&
  score(2),I22(nr,4),I21(nr,2),M1(nr,4),M(nr,4),Mi(nr,4),sc(nr,2),detM(nr)
  REAL(C_DOUBLE), ALLOCATABLE :: epst(:),Xlag(:,:),xth(:),dphik(:,:),dpsik(:,:)
  INTEGER, ALLOCATABLE :: Ir(:)
  neff = n-1
ALLOCATE(epst(neff),xth(neff),Ir(neff),Xlag(neff,2),dphik(neff,1),dpsik(neff,2))
  epst      = eps(2:n)
  Xlag(:,1) = 1
Xlag(:,2) = x(1:neff) ! # X_t-1
xth  = x(1:neff)      ! # X_t-d threshold variable
CALL rfilter(-Xlag(:,1),ma,1,neff,dphik)
  !    call dblepr ("dphik", -1, dphik, neff)
  testv = 0
M11   = matmul(transpose(dphik),dphik)/(s2*neff)
  !    call dblepr ("M11", -1, M11, 1)
  do i=1,nr,1 ! cycle over the threshold grid
  r  = trange(i)
    Ir = 0
  WHERE(xth <= r) Ir = 1 ! indicator function
  CALL rfilterm(-Xlag*SPREAD(Ir, DIM=2, NCOPIES=2),ma,1,neff,2,dpsik)
    score      = -SUM(dpsik*SPREAD(epst, DIM=2, NCOPIES=2),DIM=1)! (p+1)
    score(1)   = score(1)/(sqrt(dble(neff))*s2)
    score(2)   = score(2)/(dble(neff)*s2)
    sc(i,:)    = score
  !        call dblepr ("score", -1, score, 2)
    M22        = matmul(transpose(dpsik),dpsik)/s2
  M22(2,2)   = M22(2,2)/(dble(neff)**2)
    M22(1,1)   = M22(1,1)/(dble(neff))
    M22(2,1)   = M22(2,1)/(dble(neff)**(1.5))
    M22(1,2)   = M22(1,2)/(dble(neff)**(1.5))
    I22(i,1)   = M22(1,1)
    I22(i,2)   = M22(2,1)
    I22(i,3)   = M22(1,2)
    I22(i,4)   = M22(2,2)
    M21        = matmul(transpose(dpsik),dphik)/(s2*neff)
    M21(2,1)   = M21(2,1)/sqrt(dble(neff))
    I21(i,1)   = M21(1,1)
    I21(i,2)   = M21(2,1)
    enddo
  M1(:,1) = I21(:,1)*I21(:,1)/M11(1,1)
    M1(:,2) = I21(:,1)*I21(:,2)/M11(1,1)
    M1(:,3) = I21(:,2)*I21(:,1)/M11(1,1)
    M1(:,4) = I21(:,2)*I21(:,2)/M11(1,1)
    M       = (I22 - M1)
    detM    = M(:,1)*M(:,4)-M(:,2)*M(:,3)
    Mi(:,1) =  M(:,4)/detM
  Mi(:,2) = -M(:,2)/detM
  Mi(:,3) = -M(:,3)/detM
  Mi(:,4) =  M(:,1)/detM
  testv = (sc(:,1)*(sc(:,1)*Mi(:,1) + sc(:,2)*Mi(:,2)) + sc(:,2)*(sc(:,1)*Mi(:,3) + sc(:,2)*Mi(:,4)))
    !   call dblepr ("s2", -1, s2, 1)
    END SUBROUTINE IMAvsTARMA
    
    ! *****************************************************************************
      