module newton_solver_mod
  implicit none
  private
  public :: newton_batch, set_params, print_params

  double precision :: nur, Bparam, Tc, Tenv, Tr, cH, cK, Wparam, Eparam !Unwrapping params
  double precision :: dry_tolerance, grav

  save nur, Bparam, Tc, Tenv, Tr, cH, cK, Wparam, Eparam
  save dry_tolerance, grav
  
contains

  pure double precision function l2norm(vec) result(n)
    double precision, intent(in) :: vec(:)
    n = sqrt(sum(vec*vec))
  end function l2norm

  subroutine set_params(params_in)
    implicit none
    double precision, intent(in) :: params_in(:)
    nur = params_in(1)
    Bparam = params_in(2)
    Tc = params_in(3)
    Tenv = params_in(4)
    Tr = params_in(5)
    cH = params_in(6)
    cK = params_in(7)
    Wparam = params_in(8)
    Eparam = params_in(9)
    dry_tolerance = params_in(10)
    grav = params_in(11)
  end subroutine set_params

  subroutine print_params()
    print *, "Parameters in Newton solver module:"
    print *, "nur = ", nur
    print *, "Bparam = ", Bparam
    print *, "Tc = ", Tc
    print *, "Tenv = ", Tenv
    print *, "Tr = ", Tr
    print *, "cH = ", cH
    print *, "cK = ", cK
    print *, "Wparam = ", Wparam
    print *, "Eparam = ", Eparam
    print *, "dry_tolerance = ", dry_tolerance
    print *, "grav = ", grav
  end subroutine print_params

  ! f(q, qn) -> fq
  subroutine f(q, qn,  dt, fq)
    !This routines expects primitive variables, i.e. q = [h, u, T]
    !Dry tolerance control expected to be handled outside this routine
    implicit none
    double precision, intent(in)  :: q(:), qn(:), dt
    double precision, intent(out) :: fq(:)

    double precision ::  Hcap, Kcap
    double precision ::  h, u, T

    h = qn(1)  
    u = q(2)
    T = q(3)
  
    Hcap = cH/h
    Kcap = cK/h


    fq(1) = 0.d0 ! No equation for eta
    fq(2) = (u-qn(2))*h/dt + (3*nur/h**2.d0)*exp(-Bparam*(T-Tr))
    fq(3) = (T-qn(3))*h/dt + Eparam*(T**4.d0-Tenv**4.d0)+ &
            Wparam*(T-Tenv) +Hcap*(T-Tc) - Kcap*(u**2.d0)*exp(-Bparam*(T-Tr))
  end subroutine f

  ! Jf(x, a, params) -> J
  subroutine jf(q, qn, dt, Jout)
    implicit none
    double precision, intent(in)  :: q(:), qn(:), dt
    double precision, intent(out) :: Jout(:,:)

    double precision ::  Hcap, Kcap
    double precision ::  h, u, T
    Jout = 0.d0
    h = qn(1)
    u = q(2)
    T = q(3)
    Jout(2,2) = h/dt
    Jout(2,3) = -3*Bparam*nur*exp(-Bparam*(T-Tr))/h**2.d0
    Jout(3,2) = -2*Kcap*u*exp(-Bparam*(T-Tr))
    Jout(3,3) = h/dt + Hcap +Bparam*Kcap*u**2.d0*exp(-Bparam*(T-Tr)) + &
                4*Eparam*T**3.d0 + Wparam
  end subroutine jf

  subroutine newton_batch(m, nx, qn, dt, qnp1, tol_f, tol_step, maxit, status, iters)
    integer, intent(in) :: m, nx, maxit
    double precision, intent(in)  :: qn(m, nx), dt
    double precision, intent(in)  :: tol_f, tol_step
    double precision, intent(out) :: qnp1(m, nx)
    integer, intent(out) :: status(nx), iters(nx)

    double precision :: x(m), fx(m), J(m, m), dx(m)
    integer :: i, k, info
    integer :: ipiv(m)

    interface ! LAPACK routine for solving linear systems
      subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
        integer, intent(in) :: n, nrhs, lda, ldb
        integer, intent(out) :: ipiv(*)
        double precision, intent(inout) :: a(lda,*), b(ldb,*)
        integer, intent(inout) :: info
      end subroutine dgesv
    end interface

    do k = 1, nx
      !Initial guess
      x = qn(:, k)
      status(k) = 1
      do i = 1, maxit
        !Newton iterations to solve f(q^{n+1}(:,k)) = 0
        call f(x, qn(:, k), dt, fx)
        !Residual check for convergence
        if (l2norm(fx) < tol_f) then
          status(k) = 0
          exit
        end if
        !Compute Jf(x_i) (Jacobian of f)
        call jf(x, qn(:, k), dt, J)
        !Solve the linear system Jf(x_i)dx = -f(x_i),
        !where dx=(x_{i+1}-x_i) (RHS is overwritten)
        dx = -fx
        call dgesv(m, 1, J, m, ipiv, dx, m, info)
        if (info /= 0) then
          status(k) = 2
          exit
        end if
        !Update x_{i+1} = x_i + dx
        x = x + dx
        !Self-convergence check
        if (l2norm(dx) < tol_step) then
          status(k) = 0
          exit
        end if
      end do
      iters(k) = i
      qnp1(:, k) = x
    end do
  end subroutine newton_batch

end module newton_solver_mod