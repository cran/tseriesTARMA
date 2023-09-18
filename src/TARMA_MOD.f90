! *****************************************************************************
!! Set of Fortran utility subroutines/functions companion to the
!! tseriesTARMA package
!! Simone Giannerini
!! Last modified Sep 2023
! *****************************************************************************

MODULE TARMA_MOD
   use ISO_C_BINDING
   USE ISO_Fortran_env
   implicit none

    REAL(C_DOUBLE),PARAMETER :: M_PI         = 3.141592653589793238462643383280
    REAL(C_DOUBLE),PARAMETER :: M_1_SQRT_2PI = 0.398942280401432677939946059934 ! 1/(sqrt(2*pi))

interface ROW
    procedure ROW_REAL, ROW_LOGI
end interface
interface COL
    procedure COL_REAL, COL_LOGI
end interface
interface
    subroutine GetRNGstate() bind(C,name='GetRNGstate')
    end subroutine GetRNGstate

    subroutine PutRNGstate() bind(C,name='PutRNGstate')
    end subroutine PutRNGstate

    function norm_rand() bind(C,name='norm_rand')
     import
     implicit none
     real(C_DOUBLE) norm_rand
    end function norm_rand

    function unif_rand() bind(C,name='unif_rand')
     import
     implicit none
     real(C_DOUBLE) unif_rand
    end function unif_rand

    subroutine rqsortc(x,i,j) bind(C,name='R_qsort')
     import
     implicit none
     INTEGER(C_INT),INTENT(IN):: i,j
     REAL(C_DOUBLE),INTENT(INOUT):: x(j)
    end subroutine rqsortc

    subroutine rsortc(x,n) bind(C,name='R_isort')
     import
     implicit none
     INTEGER(C_INT),value:: n
     INTEGER(C_INT),INTENT(INOUT):: x(n)
    end subroutine rsortc

    function sumc(x,n) bind(C,name='sumc')
     use iso_c_binding
     implicit none
     integer(c_int), value :: n
     REAL(C_DOUBLE),INTENT(IN):: x(n)
     REAL(C_DOUBLE):: sumc
    end function sumc

end interface
    public GetRNGstate, PutRNGstate
    public norm_rand, unif_rand
!    public dnorm
!    public gammafn
    public rqsortc, PERM, rsortc, EMBED
    public :: ROW, COL
    private :: ROW_REAL, ROW_LOGI, COL_REAL, COL_LOGI

CONTAINS
! *****************************************************************************
SUBROUTINE EMBED(X,lags,OUT)
	! time delay embedding
	! X[t], X[t-i(1)], ... ,X[t - i(lmax)]
    !
    INTEGER, INTENT(IN)  :: lags(:)
    REAL(C_DOUBLE), INTENT(IN)  :: X(:)
    REAL(C_DOUBLE), ALLOCATABLE,INTENT(OUT) :: OUT(:,:)
    INTEGER :: i,n,d,nlag,lmax
    INTEGER, ALLOCATABLE :: lags2(:)
    nlag = size(lags)+1 ! dimension of the state vector
    lmax = MAXVAL(lags) ! max lag
    n    = size(X,1)
    ALLOCATE(OUT((n-lmax),nlag),lags2(nlag))
    lags2(1) = 0
    lags2(2:nlag) = lags
    do i=1,nlag
       d        = lags2(i)
        OUT(:,i) = X((lmax+1-d):(n-d))
    enddo
END SUBROUTINE EMBED
! *****************************************************************************
FUNCTION BOOT(N,size)
! gives a s sample with replacement of length size from the first N integers
    INTEGER,INTENT(IN) :: N,size
    INTEGER :: BOOT(size)
    REAL(C_DOUBLE) :: u(size)
    CALL randunif(u,size)
    BOOT = int(u*N)+1   ! random number between fra 1 e N
END FUNCTION BOOT

!! *****************************************************************************
FUNCTION ROW_LOGI(X)
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: X(:,:)
    INTEGER,ALLOCATABLE:: IND(:),ROW_LOGI(:,:)
    INTEGER:: nr,nc,i
    nr  = SIZE(X,DIM=1) ! nrows
    nc  = SIZE(X,DIM=2) ! ncols
    Allocate(IND(nr))
    IND = [(i,i=1,nr)]
    ROW_LOGI = SPREAD(IND,DIM= 2,NCOPIES= nc)
END FUNCTION ROW_LOGI

!! ****************************************************************************

FUNCTION ROW_REAL(X)
    IMPLICIT NONE
    REAL(C_DOUBLE),INTENT(IN):: X(:,:)
    INTEGER,ALLOCATABLE:: IND(:),ROW_REAL(:,:)
    INTEGER:: nr,nc,i
    nr  = SIZE(X,DIM=1) ! nrows
    nc  = SIZE(X,DIM=2) ! ncols
    Allocate(IND(nr))
    IND = [(i,i=1,nr)]
    ROW_REAL = SPREAD(IND,DIM= 2,NCOPIES= nc)
END FUNCTION ROW_REAL

!! ****************************************************************************

FUNCTION COL_REAL(X)
    IMPLICIT NONE
    REAL(C_DOUBLE),INTENT(IN):: X(:,:)
    INTEGER,ALLOCATABLE:: IND(:),COL_REAL(:,:)
    INTEGER:: nr,nc,i
    nr  = SIZE(X,DIM=1) ! nrows
    nc  = SIZE(X,DIM=2) ! ncols
    Allocate(IND(nc))
    IND = [(i,i=1,nc)]
    COL_REAL = SPREAD(IND,DIM= 1,NCOPIES= nr)
END FUNCTION COL_REAL

!! ****************************************************************************

FUNCTION COL_LOGI(X)
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: X(:,:)
    INTEGER,ALLOCATABLE:: IND(:),COL_LOGI(:,:)
    INTEGER:: nr,nc,i
    nr  = SIZE(X,DIM=1) ! nrows
    nc  = SIZE(X,DIM=2) ! ncols
    Allocate(IND(nc))
    IND = [(i,i=1,nc)]
    COL_LOGI = SPREAD(IND,DIM= 1,NCOPIES= nr)
END FUNCTION COL_LOGI

!! ****************************************************************************

!    function dnorm(x,mu,sigma,give_log) bind(C,name='dnorm')
!     import
!     implicit none
!     real(C_DOUBLE):: dnorm,x,mu,sigma
!     integer(C_INT):: give_log
!    end function dnorm
!! ****************************************************************************

!    function gammafn(x) bind(C,name='gammafn')
!     import
!     implicit none
!     real(C_DOUBLE):: gammafn,x
!    end function gammafn
!! ****************************************************************************

!FUNCTION WHICHM(X)
!! Equivalent to the function which in R USES F2018 findloc
!!
!    IMPLICIT NONE
!    LOGICAL,INTENT(IN):: X(:)
!    INTEGER,ALLOCATABLE:: IND(:),WHICHM(:)
!    INTEGER:: n,i
!    n = size(X)
!    ALLOCATE(IND(n))
!    IND = 0
!    WHERE(X) IND = 1
!    WHICHM = findloc(IND,1)
!END FUNCTION WHICHM

!! ****************************************************************************
    FUNCTION PERM(N)
        ! GIVES A RANDOM PERMUTATION OF THE FIRST N INTEGERS
    IMPLICIT NONE
        INTEGER,INTENT(IN) :: N
        INTEGER :: PERM(N)
        INTEGER :: B(N)
        REAL(C_DOUBLE) :: u(N)
        INTEGER :: i,j,k

        B   =  (/(i,i=1,N)/)
        PERM = 0  ! vettore degli n numeri estratti
        CALL randunif(u,N)
        do i = N,1,-1
            j = int(u(i)*i)+1   ! numero casuale compreso fra 1 e i
            PERM(i) = B(j)
            k = B(j)
            B(j)=B(i)
            B(i)=k
        enddo
    END FUNCTION PERM
!! ****************************************************************************

!! ****************************************************************************
    FUNCTION BOOTR(N,size)
! sampling without replacement
! gives a random sample without replacement of length size from the first N integers
! generalizes the function PERM
    use ISO_C_BINDING
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: N,size
    INTEGER :: BOOTR(size)
    INTEGER :: B(N)
    REAL(C_DOUBLE) :: u(N)
    INTEGER :: i,j,k

    B     =  (/(i,i=1,N)/)
    BOOTR = 0  ! vettore degli n numeri estratti
    CALL randunif(u,N)
    do i = N,(N-size+1),-1
        j = int(u(i)*i)+1   ! numero casuale compreso fra 1 e i
        BOOTR((N-i+1)) = B(j)
        k = B(j)
        B(j)=B(i)
        B(i)=k
    enddo
    END FUNCTION BOOTR
!! ****************************************************************************

! *****************************************************************************
FUNCTION WHICH(X)
! Equivalent to the function which in R
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: X(:)
    INTEGER,ALLOCATABLE:: IND(:),WHICH(:)
    INTEGER:: n,i
    n = size(X)
    ALLOCATE(IND(n))
    IND = (/ (i,i=1,n) /)
    WHICH = pack(IND, mask = X)
END FUNCTION WHICH
!! ****************************************************************************
!! ****************************************************************************
FUNCTION WHICHM(X)
! Equivalent to the function which(, arr.ind=TRUE) in R
       use ISO_C_BINDING
    IMPLICIT NONE
    LOGICAL,INTENT(IN):: X(:,:)
    INTEGER,ALLOCATABLE:: WHICHM(:,:),rg(:),cg(:),rows(:),cols(:)
    LOGICAL,ALLOCATABLE:: ind(:)
    INTEGER:: n,nr,nc,ng
    nr  = SIZE(X,DIM=1) ! nrows
    nc  = SIZE(X,DIM=2) ! ncols
    n   = SIZE(X) ! size as a vector
    ALLOCATE(ind(n))
    ind  = RESHAPE(X ,[n])
    rows = RESHAPE(ROW(X),[n])
    cols = RESHAPE(COL(X),[n])
    rg   = pack(rows,mask=ind)
    cg   = pack(cols,mask=ind)
    ng   = size(rg)
    WHICHM(:,:) = RESHAPE([rg,cg],[ng,2])
END FUNCTION WHICHM
!! ****************************************************************************

subroutine randnorm(x,n)
! generates n standard normal random numbers
   use ISO_C_BINDING
   implicit none
   INTEGER,INTENT(IN):: n
   real(C_DOUBLE),intent(out):: x(n)
   integer :: i
   call GetRNGstate();
    do i=1,n
        x(i) = norm_rand()
    end do
   call PutRNGstate();
end subroutine randnorm

!! *****************************************************************************

subroutine randunif(x,n)
! generates n continuous uniform random numbers in [0,1]
   use ISO_C_BINDING
   implicit none
   INTEGER,INTENT(IN):: n
   real(C_DOUBLE),intent(out):: x(n)
   integer :: i
   call GetRNGstate();
    do i=1,n
        x(i) = unif_rand()
    end do
   call PutRNGstate();
end subroutine randunif

! *****************************************************************************

SUBROUTINE randunifd(x,n,k)

! generates n discrete uniform random numbers in 1:k

    use ISO_C_BINDING
    IMPLICIT NONE
    INTEGER,INTENT(IN):: n,k
    INTEGER,INTENT(OUT):: x(n)
    REAL(C_DOUBLE):: y(n)
    CALL randunif(y,n)
    x = INT(y*k)+1
END SUBROUTINE randunifd

!! *****************************************************************************
!! *****************************************************************************
SUBROUTINE rsort(x,n)
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
     INTEGER(C_INT),INTENT(IN):: n
     INTEGER(C_INT),INTENT(INOUT):: x(n)
    call rsortc(x,n)
END SUBROUTINE rsort
!! DOES NOT WORK !!
! *****************************************************************************
SUBROUTINE qsort(x,i,j)
! sorts the vector x(i:j)
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
     INTEGER(C_INT),INTENT(IN):: i,j
     REAL(C_DOUBLE),INTENT(INOUT):: x(j)
    call qsort3(x,i,j)
END SUBROUTINE qsort

! *****************************************************************************

!! ****************************************************************************
SUBROUTINE permute(x,n)
! random permutation of a vector x
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
     INTEGER(C_INT),INTENT(IN):: n
     REAL(C_DOUBLE),INTENT(INOUT):: x(n)
     INTEGER:: IND(n)
     IND = PERM(n)
     x   = x(IND)
END SUBROUTINE permute
!! ****************************************************************************

SUBROUTINE samplenr(x,n,out,size)
! sampling without replacement
! INPUT:  vector x of length n
! OUTPUT: random vector out of length size
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
     INTEGER(C_INT),INTENT(IN):: n,size
     REAL(C_DOUBLE),INTENT(IN):: x(n)
     REAL(C_DOUBLE),INTENT(OUT):: out(size)
     INTEGER:: ind(size)
     ind = BOOTR(n,size)
     out   = x(ind)
END SUBROUTINE samplenr
!! ****************************************************************************

!SUBROUTINE fminc(f,ax,bx,tol,res)
!    use, intrinsic :: iso_c_binding
!     REAL(C_DOUBLE),INTENT(IN):: ax,bx,tol
!     REAL(C_DOUBLE) :: f
!       INTERFACE
!         FUNCTION f(x) bind(C)
!           use ISO_C_BINDING
!           REAL(C_DOUBLE) :: x
!         END FUNCTION f
!       END INTERFACE
!     REAL(C_DOUBLE),INTENT(OUT):: res
!
!    res = fmin(f,ax,bx,tol)
!END SUBROUTINE fminc
!! ****************************************************************************


SUBROUTINE Per(A,N,OUT)
	! Calcola il quadrato di un vettore (REAL(C_DOUBLE)), cioe' A*A'
	! Dato un vettore di N*1 restituisce una matrice quadrata N*N
    use ISO_C_BINDING
    integer, INTENT(IN) :: N
    REAL(C_DOUBLE), INTENT(IN) :: A(N)
    REAL(C_DOUBLE), INTENT(OUT) :: OUT(N,N)  ! must be square
    INTEGER :: i
    do i=1,N
        OUT(:,i)=A(i)*A
    enddo
END SUBROUTINE Per
! -------------------------------------------------------

SUBROUTINE KGAUSS(x,d,K)
    use ISO_C_BINDING
    IMPLICIT NONE
!! ********************************
!! Multiplicative Gaussian Kernel
!! ********************************
    INTEGER,INTENT(IN)  :: d
    REAL(C_DOUBLE),INTENT(IN)  :: x(d)
    REAL(C_DOUBLE),INTENT(OUT) :: K
    K = (M_1_SQRT_2PI**d)*EXP(-0.5*(SUM(x**2)))
END SUBROUTINE KGAUSS
! -------------------------------------------------------

    SUBROUTINE INVERSE(A,N,AINV)
!! ***************************************
!! MATRIX INVERSION
!! A    : N by N REAL64 matrix
!! Ainv : Inverse of A
!! WRAPPER FOR THE LAPACK SUBROUTINE DGESV
!! ***************************************
    use ISO_C_BINDING
    IMPLICIT NONE
	
    INTEGER,INTENT(IN):: N
    REAL(C_DOUBLE),INTENT(IN)::  A(N,N)
    REAL(C_DOUBLE),INTENT(INOUT)::AINV(N,N)
    INTEGER INFO,I
    REAL(C_DOUBLE)::IPIV(N),ADUM(N,N)
	
    ADUM = A
    AINV = 0
    FORALL (I=1:N)  AINV(I, I) = 1
	
    CALL DGESV(N, N, ADUM, N, IPIV, AINV, N, INFO)
END  SUBROUTINE INVERSE

! *****************************************************************************

subroutine rfilterm(x,f,p,n,k,y)
! equivalent to filter(method=recursive) in R
! Matrix version applies the filter to each column of the matrix
    use ISO_C_BINDING
   implicit none
   INTEGER,INTENT(IN):: p,n,k
   REAL(C_DOUBLE),intent(IN):: x(n,k),f(p)
   REAL(C_DOUBLE),intent(OUT):: y(n,k)
   integer :: i
   REAL(C_DOUBLE) :: f2(p,k)
   REAL(C_DOUBLE) :: y2(-p:n,k)
   y2 = 0
   f2 = SPREAD(f, DIM=2, NCOPIES=k) ! p*k
    do i=1,n
        y2(i,:) = x(i,:) + SUM(f2*y2((i-1):(i-p):(-1),:),DIM=1)
    end do
    y = y2(1:n,:)
end subroutine rfilterm

! *****************************************************************************

subroutine rfilter(x,f,p,n,y)
! equivalent to filter(method=recursive) in R for input vector x
    use ISO_C_BINDING
   implicit none
   INTEGER,INTENT(IN):: p,n
   REAL(C_DOUBLE),intent(IN):: x(n),f(p)
   REAL(C_DOUBLE),intent(OUT):: y(n)
   integer :: i
   REAL(C_DOUBLE):: y2(-p:n)
   y2 = 0
!    do i =1,p
!        call dblepr("f2", -1, f2(i,:),p)
!    end do
    do i=1,n
        y2(i) = x(i) + dot_product(y2((i-1):(i-p):(-1)),f)
!        call dblepr("y2", -1, SUM(y2(((i+p-1):i:-1),:),DIM=1),k)
    end do
    y = y2(1:n)
!    y = x
end subroutine rfilter

! *****************************************************************************
SUBROUTINE lsfit(x,y,tol,b,rsd,itc)
!! --------------------------------------------------------
!! Wrapper for the (LINPACK) subroutine DQRLS
!! SG 2022
!! --------------------------------------------------------
!    USE TARMA_MOD
    use ISO_C_BINDING
 !   USE ISO_Fortran_env
    IMPLICIT NONE
    INTEGER,INTENT(IN):: itc
    REAL(C_DOUBLE),INTENT(IN):: x(:,:),y(:,:),tol
    REAL(C_DOUBLE),ALLOCATABLE,INTENT(OUT):: b(:),rsd(:,:)
    REAL(C_DOUBLE),ALLOCATABLE:: qty(:,:),qraux(:),work(:),x2(:,:)
    INTEGER,ALLOCATABLE:: jpvt(:)
    INTEGER:: n,p,k,ny
    REAL(C_DOUBLE):: tol2
    n  = size(x,1) ! nrows
    p  = size(x,2) ! ncols
    ny = size(y,2) ! ncols
! !   jpvt = [(i, i=1,p,1)]
   if(itc==1) then
     p = p+1
     ALLOCATE(x2(n,p))
     x2(:,1) = 1
     x2(:,2:p) = x(:,:)
   else
    ALLOCATE(x2(n,p))
    x2 = x
   end if
    ALLOCATE(b(p),rsd(n,ny),qty(n,ny),qraux(p),work(2*p),jpvt(p))
    tol2 = epsilon(tol)/maxval(abs(x))
    CALL dqrls(x2,n,p,y,ny,tol2,b,rsd,qty,k,jpvt,qraux,work)
END  SUBROUTINE lsfit

END MODULE TARMA_MOD

!! *****************************************************************************
