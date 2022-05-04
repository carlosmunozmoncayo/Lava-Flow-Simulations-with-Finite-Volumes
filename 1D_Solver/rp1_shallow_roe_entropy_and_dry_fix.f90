! =====================================================
subroutine rp1(maxmx,num_eqn,num_waves,num_aux,num_ghost,num_cells, &
               ql,qr,auxl,auxr,wave,s,amdq,apdq)

! =====================================================

! Roe-solver for the 1D shallow water equations with entropy and dry states fix
! We add an extra equation for the temperature of the fluid.

! waves: 3
! equations: 3

! Conserved quantities:
!       1 depth
!       2 x_momentum
!       3 T_momentum

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.

! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!   f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.

! Note that the i'th Riemann problem has left state qr(:,i-1)
!                                    and right state ql(:,i)
! From the basic clawpack routines, this routine is called with ql = qr


! This Riemann solver differs from the original clawpack Riemann solver
! for the interleaved indices

    implicit none

    integer, intent(in) :: maxmx, num_eqn, num_waves, num_aux, num_ghost, &
                            num_cells
    real(kind=8), intent(in), dimension(num_eqn,1-num_ghost:maxmx+num_ghost) :: ql, qr
    real(kind=8), intent(in), dimension(num_aux,1-num_ghost:maxmx+num_ghost) :: auxl, auxr
    real(kind=8), intent(out) :: s(num_waves, 1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(out) :: wave(num_eqn, num_waves, 1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(out), dimension(num_eqn,1-num_ghost:maxmx+num_ghost) :: amdq,apdq

!   # Gravity constant set in the shallow1D.py or setprob.f90 file
    real(kind=8):: grav
    common /cparam/ grav

!   # Roe averages quantities of each interface
    !parameter (maxm2 = 1800)
    real(kind=8) :: h,u,T,a
    real(kind=8), dimension(3) :: delta 
    real(kind=8) :: a1,a2,a3,a4
    real(kind=8) :: hl,Tl,sl,hr,Tr,sr,ul, ur
    real(kind=8) :: hsqrtl, hsqrtr, hsq2
    integer :: m, mw, i
!   # Note that notation for u and v reflects assumption that the
!   # Riemann problems are in the x-direction with u in the normal
!   # direciton and v in the orthogonal direcion, but with the above
!   # definitions of mu and mv the routine also works with ixy=2
!   # and returns, for example, f0 as the Godunov flux g0 for the
!   # Riemann problems u_t + g(u)_y = 0 in the y-direction.


!   # Compute the Roe-averaged variables needed in the Roe solver.

    do 10 i=2-num_ghost,num_cells+num_ghost
        h = (qr(1,i-1)+ql(1,i))*0.50d0
        hsqrtl = dsqrt(qr(1,i-1))
        hsqrtr = dsqrt(ql(1,i))
        hsq2 = hsqrtl + hsqrtr
        u = (qr(2,i-1)/hsqrtl + ql(2,i)/hsqrtr) / hsq2
        T = (qr(3,i-1)/hsqrtl + ql(3,i)/hsqrtr) / hsq2
        a = dsqrt(grav*h)


!   # now split the jump in q at each interface into waves

        delta(1) = ql(1,i) - qr(1,i-1)
        delta(2) = ql(2,i) - qr(2,i-1)
        !delta(3) = ql(mv,i) - qr(mv,i-1)
        delta(3) = ql(3,i) - qr(3,i-1)
!   # find a1 thru a4, the coefficients of the 4 eigenvectors:

        a1 = ((u+a)*delta(1) - delta(2))*(0.50d0/a)
        a2 = -T*delta(1) + delta(3)
        !a3 = -T(i)*delta(4) + delta(4)
        a3 = ((-u+a)*delta(1) + delta(2))*(0.50d0/a)
    
!   # Compute the waves and speeds, apply entropy fix if needed as in (Monthe, 1999).
!   # Note that there is no need to check for transonic rarefactions
!   # in the 2nd and 3rd wave because they are linearly degenerate
    
        wave(1,1,i) = a1
        wave(2,1,i) = a1*(u-a)
        !wave(mv,1,i) = a1*v(i)
        wave(3,1,i) = a1*T
        s(1,i) = u-a
        !left state of 1-wave  Ql
        hl = qr(1,i-1)  
        ul = qr(2,i-1)/hl
        sl = ul-dsqrt(hl*grav)
        !right state of 1-wave  Ql+W1
        hr = hl+wave(1,1,i)
        ur = (qr(2,i-1)+wave(2,1,i))/hr
        sr = ur-dsqrt(hr*grav)
        !Check for a transonic 1-rarefaction
        if (sl <0 .AND. sr >0) then
            s(1,i)=sr*(s(1,i)-sl)/(sr-sl)
        END IF

        wave(1,2,i) = 0.0d0
        wave(2,2,i) = 0.0d0
        wave(3,2,i) = a2
        !wave(4,2,i) = 0.0d0
        s(2,i) = u 

    
        wave(1,3,i) = a3
        wave(2,3,i) = a3*(u+a)
        !wave(,4,i) = a4*v(i)
        wave(3,3,i)= a3*T
        s(3,i) = u+a
        !right state of 3-wave Qr
        hr = ql(1,i)
        ur = ql(2,i)/hr
        sr = ur+dsqrt(hr*grav)
        !left state of 3-wave Qr-W3
        hl = hr-wave(1,3,i)
        ul = (ql(2,i)-wave(2,3,i))/hl
        sl = ul+dsqrt(hl*grav)
        !Checking for a transonic 4-rarefaction
        if (sl <0 .AND. sr >0) then
            s(1,i)=sr*(s(1,i)-sl)/(sr-sl)
        END IF

    10 END DO


!    # compute flux differences amdq and apdq.
!    ---------------------------------------

!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    do 100 m=1,3
        do 30 i=2-num_ghost,num_cells+num_ghost
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,num_waves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            90 END DO
        30 END DO
    100 END DO
    
    return
    end subroutine rp1


