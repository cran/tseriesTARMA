! *****************************************************************************
!! Set of Fortran subroutines/functions for TARMA modelling
!! companion to the tseriesTARMA package
!! Simone Giannerini
!! Last modified August 2021

! *****************************************************************************
! *****************************************************************************
SUBROUTINE tarmaLS(x,nx,th,k,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,L)
!! Least Squares objective function for a TARMA model
    USE ISO_Fortran_env
    use ISO_C_BINDING
    IMPLICIT NONE
    INTEGER,intent(in):: nx,k,p1,p2,q1,q2,Ir(nx),tlag1(p1-1),tlag2(p2-1),&
    mlag1(q1),mlag2(q2)
    REAL(C_DOUBLE),intent(in) :: x(nx),th(p1+p2+q1+q2)
    REAL(C_DOUBLE),intent(out):: L
    integer :: i
    REAL(C_DOUBLE) :: eps(nx),h,phi1(p1),phi2(p2),th1(q1),th2(q2)
    eps  = 0
    L    = 0
    h    = 1
    phi1 = th(1:p1)
    phi2 = th((p1+1):(p1+p2))
    th1  = th((p1+p2+1):(p1+p2+q1))
    th2  = th((p1+p2+q1+1):(p1+p2+q1+q2))
    do i=(k+1),nx,1
        eps(i) = x(i)-(dot_product(phi1,[h,x(i-tlag1)])+dot_product(th1,eps(i-mlag1)))*Ir(i)&
        - (dot_product(phi2,[h,x(i-tlag2)])+dot_product(th2,eps(i-mlag2)))*(1-Ir(i))
        L = L + eps(i)**2
    enddo
 END SUBROUTINE tarmaLS
! *****************************************************************************
!SUBROUTINE tarmaLSW(x,nx,th,k,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,alpha,L)
!!! Robust Least Squares objective function for a TARMA model
!    USE ISO_Fortran_env
!    IMPLICIT NONE
!    INTEGER,intent(in):: nx,k,p1,p2,q1,q2,Ir(nx),tlag1(p1-1),tlag2(p2-1),&
!    mlag1(q1),mlag2(q2)
!    REAL(C_DOUBLE),intent(in) :: x(nx),th(p1+p2+q1+q2),alpha
!    REAL(C_DOUBLE),intent(out):: L
!    integer :: i
!    REAL(C_DOUBLE) :: eps(nx),h,phi1(p1),phi2(p2),th1(q1),th2(q2)
!    eps  = 0
!    L    = 0
!    h    = 1
!    phi1 = th(1:p1)
!    phi2 = th((p1+1):(p1+p2))
!    th1  = th((p1+p2+1):(p1+p2+q1))
!    th2  = th((p1+p2+q1+1):(p1+p2+q1+q2))
!    do i=(k+1),nx,1
!        eps(i) = x(i)-(dot_product(phi1,[h,x(i-tlag1)])+dot_product(th1,eps(i-mlag1)))*Ir(i)&
!        - (dot_product(phi2,[h,x(i-tlag2)])+dot_product(th2,eps(i-mlag2)))*(1-Ir(i))
!        L = L + exp(-(alpha*eps(i)**2))
!    enddo
! END SUBROUTINE tarmaLSW
! *****************************************************************************
SUBROUTINE tarmaLSW(x,nx,th,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,wt,indg,ng,L)
!! Robust Least Squares objective function for a TARMA model
    USE ISO_Fortran_env
    use ISO_C_BINDING
    IMPLICIT NONE
    INTEGER,intent(in):: nx,p1,p2,q1,q2,Ir(nx),tlag1(p1-1),tlag2(p2-1),&
    mlag1(q1),mlag2(q2),ng,indg(ng)
    REAL(C_DOUBLE),intent(in) :: x(nx),wt(nx),th(p1+p2+q1+q2)
    REAL(C_DOUBLE),intent(out):: L
    integer :: i,j
    REAL(C_DOUBLE) :: eps(nx),h,phi1(p1),phi2(p2),th1(q1),th2(q2)
    eps  = 0
    L    = 0
    h    = 1
    phi1 = th(1:p1)
    phi2 = th((p1+1):(p1+p2))
    th1  = th((p1+p2+1):(p1+p2+q1))
    th2  = th((p1+p2+q1+1):(p1+p2+q1+q2))
    do j=1,ng,1
        i = indg(j)
        eps(i) = x(i)-(dot_product(phi1,[h,x(i-tlag1)])+dot_product(th1,eps(i-mlag1)))*Ir(i)&
        - (dot_product(phi2,[h,x(i-tlag2)])+dot_product(th2,eps(i-mlag2)))*(1-Ir(i))
        L = L + wt(i)*eps(i)**2
    enddo
 END SUBROUTINE tarmaLSW
! *****************************************************************************
SUBROUTINE tarmaDLS(x,nx,th,k,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,DL)
!! Gradient of the Least Squares objective function for a TARMA model
    USE ISO_Fortran_env
    use ISO_C_BINDING
    IMPLICIT NONE
    INTEGER,intent(in):: nx,k,p1,p2,q1,q2,Ir(nx),tlag1(p1-1),tlag2(p2-1),&
        mlag1(q1),mlag2(q2)
    REAL(C_DOUBLE),intent(in) :: x(nx),th(p1+p2+q1+q2)
    REAL(C_DOUBLE),intent(out):: DL(p1+p2+q1+q2)
    integer :: i
    REAL(C_DOUBLE) :: h,eps(nx),phi1(p1),phi2(p2),th1(q1),th2(q2),&
    deps(nx,p1+p2+q1+q2)
    h     = 1.0
    deps  = 0.0
    eps   = 0.0
    DL    = 0.0
    phi1  = th(1:p1)
    phi2  = th((p1+1):(p1+p2))
    th1   = th((p1+p2+1):(p1+p2+q1))
    th2   = th((p1+p2+q1+1):(p1+p2+q1+q2))
    do i=(k+1),nx,1
        eps(i) = x(i)-(dot_product(phi1,[h,x(i-tlag1)])+dot_product(th1,eps(i-mlag1)))*Ir(i)&
        - (dot_product(phi2,[h,x(i-tlag2)])+dot_product(th2,eps(i-mlag2)))*(1-Ir(i))
        deps(i,1) =  -Ir(i)*(1 + dot_product(th1,deps(i-mlag1,1)))&
        -(1-Ir(i))*dot_product(th2,deps(i-mlag2,1))
        deps(i,2:p1) = -Ir(i)*(x(i-tlag1) + matmul(th1,deps(i-mlag1,2:p1)))&
        -(1-Ir(i))*(matmul(th2,deps(i-mlag2,2:p1)))
        deps(i,(p1+1)) = -Ir(i)*(dot_product(th1,deps(i-mlag1,(p1+1))))&
        -(1-Ir(i))*(1+dot_product(th2,deps(i-mlag2,(p1+1))))
        deps(i,(p1+2):(p1+p2)) = -Ir(i)*(matmul(th1,deps(i-mlag1,(p1+2):(p1+p2))))&
        -(1-Ir(i))*(x(i-tlag2)+ matmul(th2,deps(i-mlag2,(p1+2):(p1+p2))))
        deps(i,(p1+p2+1):(p1+p2+q1)) = -Ir(i)*(eps(i-mlag1)+matmul(th1,deps(i-mlag1,(p1+p2+1):(p1+p2+q1)))) &
        -(1-Ir(i))*(matmul(th2,deps(i-mlag2,(p1+p2+1):(p1+p2+q1))))
        deps(i,(p1+p2+q1+1):(p1+p2+q1+q2)) = -Ir(i)*(matmul(th1,deps(i-mlag1,(p1+p2+q1+1):(p1+p2+q1+q2)))) &
        -(1-Ir(i))*(eps(i-mlag2)  +matmul(th2,deps(i-mlag2,(p1+p2+q1+1):(p1+p2+q1+q2))))
    enddo
    DL = 2*SUM(SPREAD(eps, DIM=2, NCOPIES=(p1+p2+q1+q2))*deps,DIM=1)
END SUBROUTINE tarmaDLS
! *****************************************************************************
SUBROUTINE tarmaDLSW(x,nx,th,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,wt,indg,ng,DL)
!! Gradient of the robust Least Squares objective function for a TARMA model
    USE ISO_Fortran_env
    use ISO_C_BINDING
    IMPLICIT NONE
    INTEGER,intent(in):: nx,p1,p2,q1,q2,Ir(nx),tlag1(p1-1),tlag2(p2-1),&
        mlag1(q1),mlag2(q2),ng,indg(ng)
    REAL(C_DOUBLE),intent(in) :: x(nx),wt(nx),th(p1+p2+q1+q2)
    REAL(C_DOUBLE),intent(out):: DL(p1+p2+q1+q2)
    integer :: i,j
    REAL(C_DOUBLE) :: h,eps(nx),phi1(p1),phi2(p2),th1(q1),th2(q2),&
    deps(nx,p1+p2+q1+q2),epsm(nx,p1+p2+q1+q2)
    h     = 1.0
    deps  = 0.0
    eps   = 0.0
    DL    = 0.0
    phi1  = th(1:p1)
    phi2  = th((p1+1):(p1+p2))
    th1   = th((p1+p2+1):(p1+p2+q1))
    th2   = th((p1+p2+q1+1):(p1+p2+q1+q2))
    do j=1,ng,1
        i = indg(j)
        eps(i) = x(i)-(dot_product(phi1,[h,x(i-tlag1)])+dot_product(th1,eps(i-mlag1)))*Ir(i)&
        - (dot_product(phi2,[h,x(i-tlag2)])+dot_product(th2,eps(i-mlag2)))*(1-Ir(i))
        deps(i,1) =  -Ir(i)*(1 + dot_product(th1,deps(i-mlag1,1)))&
        -(1-Ir(i))*dot_product(th2,deps(i-mlag2,1))
        deps(i,2:p1) = -Ir(i)*(x(i-tlag1) + matmul(th1,deps(i-mlag1,2:p1)))&
        -(1-Ir(i))*(matmul(th2,deps(i-mlag2,2:p1)))
        deps(i,(p1+1)) = -Ir(i)*(dot_product(th1,deps(i-mlag1,(p1+1))))&
        -(1-Ir(i))*(1+dot_product(th2,deps(i-mlag2,(p1+1))))
        deps(i,(p1+2):(p1+p2)) = -Ir(i)*(matmul(th1,deps(i-mlag1,(p1+2):(p1+p2))))&
        -(1-Ir(i))*(x(i-tlag2)+ matmul(th2,deps(i-mlag2,(p1+2):(p1+p2))))
        deps(i,(p1+p2+1):(p1+p2+q1)) = -Ir(i)*(eps(i-mlag1)+matmul(th1,deps(i-mlag1,(p1+p2+1):(p1+p2+q1)))) &
        -(1-Ir(i))*(matmul(th2,deps(i-mlag2,(p1+p2+1):(p1+p2+q1))))
        deps(i,(p1+p2+q1+1):(p1+p2+q1+q2)) = -Ir(i)*(matmul(th1,deps(i-mlag1,(p1+p2+q1+1):(p1+p2+q1+q2)))) &
        -(1-Ir(i))*(eps(i-mlag2)  +matmul(th2,deps(i-mlag2,(p1+p2+q1+1):(p1+p2+q1+q2))))
    enddo
    epsm = SPREAD(eps, DIM=2, NCOPIES=(p1+p2+q1+q2))
    DL   = SUM(2*SPREAD(wt, DIM=2, NCOPIES=(p1+p2+q1+q2))*epsm*deps,DIM=1)
END SUBROUTINE tarmaDLSW
! *****************************************************************************
!! *****************************************************************************
!SUBROUTINE tarmaDLSW(x,nx,th,k,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,alpha,DL)
!!! Gradient of the robust Least Squares objective function for a TARMA model
!    USE ISO_Fortran_env
!    use ISO_C_BINDING
!    IMPLICIT NONE
!    INTEGER,intent(in):: nx,k,p1,p2,q1,q2,Ir(nx),tlag1(p1-1),tlag2(p2-1),&
!        mlag1(q1),mlag2(q2)
!    REAL(C_DOUBLE),intent(in) :: x(nx),th(p1+p2+q1+q2),alpha
!    REAL(C_DOUBLE),intent(out):: DL(p1+p2+q1+q2)
!    integer :: i
!    REAL(C_DOUBLE) :: h,eps(nx),phi1(p1),phi2(p2),th1(q1),th2(q2),&
!    deps(nx,p1+p2+q1+q2),epsm(nx,p1+p2+q1+q2)
!    h     = 1.0
!    deps  = 0.0
!    eps   = 0.0
!    DL    = 0.0
!    phi1  = th(1:p1)
!    phi2  = th((p1+1):(p1+p2))
!    th1   = th((p1+p2+1):(p1+p2+q1))
!    th2   = th((p1+p2+q1+1):(p1+p2+q1+q2))
!    do i=(k+1),nx,1
!        eps(i) = x(i)-(dot_product(phi1,[h,x(i-tlag1)])+dot_product(th1,eps(i-mlag1)))*Ir(i)&
!        - (dot_product(phi2,[h,x(i-tlag2)])+dot_product(th2,eps(i-mlag2)))*(1-Ir(i))
!        deps(i,1) =  -Ir(i)*(1 + dot_product(th1,deps(i-mlag1,1)))&
!        -(1-Ir(i))*dot_product(th2,deps(i-mlag2,1))
!        deps(i,2:p1) = -Ir(i)*(x(i-tlag1) + matmul(th1,deps(i-mlag1,2:p1)))&
!        -(1-Ir(i))*(matmul(th2,deps(i-mlag2,2:p1)))
!        deps(i,(p1+1)) = -Ir(i)*(dot_product(th1,deps(i-mlag1,(p1+1))))&
!        -(1-Ir(i))*(1+dot_product(th2,deps(i-mlag2,(p1+1))))
!        deps(i,(p1+2):(p1+p2)) = -Ir(i)*(matmul(th1,deps(i-mlag1,(p1+2):(p1+p2))))&
!        -(1-Ir(i))*(x(i-tlag2)+ matmul(th2,deps(i-mlag2,(p1+2):(p1+p2))))
!        deps(i,(p1+p2+1):(p1+p2+q1)) = -Ir(i)*(eps(i-mlag1)+matmul(th1,deps(i-mlag1,(p1+p2+1):(p1+p2+q1)))) &
!        -(1-Ir(i))*(matmul(th2,deps(i-mlag2,(p1+p2+1):(p1+p2+q1))))
!        deps(i,(p1+p2+q1+1):(p1+p2+q1+q2)) = -Ir(i)*(matmul(th1,deps(i-mlag1,(p1+p2+q1+1):(p1+p2+q1+q2)))) &
!        -(1-Ir(i))*(eps(i-mlag2)  +matmul(th2,deps(i-mlag2,(p1+p2+q1+1):(p1+p2+q1+q2))))
!    enddo
!    epsm = SPREAD(eps, DIM=2, NCOPIES=(p1+p2+q1+q2))
!    DL   = SUM(exp(-alpha*epsm**2)*(-2*alpha*epsm*deps),DIM=1)
!END SUBROUTINE tarmaDLSW
!! *****************************************************************************
