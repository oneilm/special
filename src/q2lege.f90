!
! (c) 2015 Mike O'Neil
!
! this code containes routines for evaluate the half-order legendre
! polynomials of the second-kind, and scalings to compute the
! axisymmetric laplace kernels
!




subroutine torlaps_eps(eps, src, targ, nterms, lvals, vals)
  implicit real *8 (a-h,o-z)
  real *8 :: src(2), targ(2), vals(0:lvals)

  !
  ! input:
  !   eps - the decay required in the last term
  !   src - source in r,z coordinates
  !   targ - the target in r,z coordinates
  !   lvals - the length of vals, vals(0:lvals)
  !
  ! ouput:
  !   nterms - the number of terms returned, vals(0:nterms)
  !   vals - the values of
  !         \frac{1}{2 \pi} \int_0^{2\pi} exp(-i*mode*\theta)
  !         \frac{1}{4 \pi R} d\theta
  !



  return
end subroutine torlaps_eps





subroutine torlap1_kernel(par0, src, targ, mode, work, &
    val, grad, hess)
  implicit real *8 (a-h,o-z)
  real *8 :: src(2), targ(2), pars2(1), grad(2), hess(2,2)
  real *8, parameter :: done=1, half=0.5d0
  real *8 :: vals(0:10000)
  
  !
  ! this returns the value of a single mode of the axisymmetric
  ! Laplace Greens's function:
  !
  !   val = \frac{1}{2 \pi} \int_0^{2\pi} exp(-i*mode*\theta)
  !         frac{1}{4 \pi R} d\theta
  !
  ! input:
  !   work - work array, must be sufficiently long, say 30000
  !
  ! output:
  !

  lwork = 9900
  i1 = 1
  i2 = i1 + lwork + 10
  i3 = i2 + lwork + 10
  lvals = lwork/4

  
!!!!call torlaps(ier, nterms, src, targ, vals, lvals, ntop)

  lvals = 10000
  ier = 0
  call torlaps(ier, mode, src, targ, vals, lvals, ntop)
  call prinf('after torlaps, ntop = *', ntop, 1)

  call prin2('all laps = *', vals, ntop+1)

  !!call torlaps_all(ier, nterms, src, targ, work(i1), lvals, &
  !!    work(i2), work(i3), ntop)

  if (ier .ne. 0) then
    call prinf('ier = *', ier, 1)
    stop
  endif

  val = vals(abs(mode))
  return

  
  ival = mode+1
  igrad = i2 + (mode)*2
  ihess = i3 + mode*4

  val = work(mode+1)
  !!grad(1) = work(igrad)
  !!grad(2) = work(igrad+1)

  !!hess(1,1) = work(ihess)
  !!hess(2,1) = work(ihess+1)
  !!hess(1,2) = work(ihess+2)
  !!hess(2,2) = work(ihess+3)

  return
end subroutine torlap1_kernel





subroutine torlaps_kernel(par0, src, targ, modemax, pars2, &
    vals, grads, hesses)
  implicit real *8 (a-h,o-z)
  real *8 :: src(2), targ(2), pars2(1)
  real *8 :: vals(0:1), grad(2,0:1), hesses(2,2,0:1)
  real *8, parameter :: done=1, half=0.5d0

  !
  ! this returns the value of a single mode of the axisymmetric
  ! Laplace Greens's function:
  !
  !   val = \frac{1}{2 \pi} \int_0^{2\pi} exp(-i*mode*\theta)
  !         frac{1}{4 \pi R} d\theta
  !
  ! input:
  !   work - work array, must be sufficiently long, say 30000
  !
  ! output:
  !   vals - the values
  !
  !
  !   NOTE: the above arrays must be say, 10000 long each -
  !   just in case
  !

!!!!call torlaps(ier, nterms, src, targ, vals, lvals, ntop)
  !!call torlaps_all(ier, nterms, src, targ, work(i1), lvals, &
  !!    work(i2), work(i3), ntop)

  if (ier .ne. 0) then
    call prinf('ier = *', ier, 1)
    stop
  endif

  ival = mode+1
  igrad = i2 + (mode)*2
  ihess = i3 + mode*4

  !val = work(mode+1)
  !grad(1) = work(igrad)
  !grad(2) = work(igrad+1)

  !hess(1,1) = work(ihess)
  !hess(2,1) = work(ihess+1)
  !hess(1,2) = work(ihess+2)
  !hess(2,2) = work(ihess+3)

  return
end subroutine torlaps_kernel





subroutine torlaps(ier, nterms, src, targ, vals, lvals, ntop)
  implicit real *8 (a-h,o-z)
  real *8 :: vals(0:lvals), src(2), targ(2)
  real *8, parameter :: done=1, half=0.5d0

  !
  ! this is basically a wrapper for q2leges that will return
  ! the normalized version which are the fourier modes of 1/r.
  ! The integral computed is:
  !
  !     a_n = \frac{1}{2 \pi} \int_0^{2\pi} exp(-i*n*\theta)
  !         frac{1}{4 \pi R} d\theta
  !
  ! Just a reminder, this is an even integral, so a_{-n} = a_n
  !
  ! input:
  !   nterms - the number of terms to compute, 0 through nterms
  !       are returned
  !   src - the source in r,z coordinates
  !   targ - the target in r,z coordiates
  !   lvals - the length of vals
  !
  ! output:
  !   vals -
  !   ntop -
  !

  x = targ(1)**2 + src(1)**2 + (targ(2) - src(2))**2
  x = x/2/src(1)/targ(1)
  xminus = (targ(2) - src(2))**2 + (targ(1) - src(1))**2
  xminus = xminus/2/src(1)/targ(1)
  call q2leges(ier, nterms, x, xminus, vals, lvals, ntop)

  call prin2('after q2leges, vals = *', vals, nterms+1)

  eps = 1.0d-16
  call q2leges_eps(ier, eps, x, xminus, vals, lvals, nterms)
  call prin2('after q2leges_eps, vals = ', vals, nterms+1)
  stop
  
  pi = 4*atan(done)
  fac = 4*pi*pi*sqrt(src(1)*targ(1))
  do i = 0,ntop
    vals(i) = vals(i)/fac
  enddo

  return
end subroutine torlaps





subroutine q2leges_eps(ier, eps, x, xminus, vals, lvals, nterms)
  implicit real *8 (a-h,o-z)
  real *8 :: vals(0:lvals)
  real *8, parameter :: done = 1

  !
  ! evaluate all of the functions
  !
  !   Q_{-1/2}, Q_{1/2}, ...
  !
  ! that are larger than eps to a *relative* precision of eps.
  ! the only difference between this routine and q2leges is memory
  ! management, and output data.
  !
  ! input:
  !   eps - absolute size of last term desired
  !   x - the argument
  !   xminus - should be set to 1-x, evaluated as accurately as possible
  !   lvals - the length of vals, assuming vals(0:lvals)
  !
  ! output:
  !   ier - error code, ier=4 means lvals is too small
  !   vals - the values of the legendre functions
  !   nterms - vals(0:nterms) will be larger than eps
  !
  
  ier = 0
  eps2 = eps**2/100

  call q2lege01(x, xminus, q0, q1)

  !
  ! run the forward recurrence now until 1/eps2 blowup is seen
  !

  vprev = 0
  v = 1
  upbound = 1/eps2
  half = done/2

  do i=1,lvals
    vnext = (2*i*x*v - (i-half)*vprev)/(i+half)
    if (abs(vnext) .ge. upbound) then
      nterms = i+1
      goto 1001
    endif
    vprev = v
    v = vnext
  enddo

  !
  ! if here, bomb and return
  !
  ier = 4
  call prin2('lvals not large enough in q2leges_eps!*', done,0)
  call prinf('lvals = *', lvals, 1)
  call prin2('vnext = *', vnext, 1)
  stop

 1001 continue
  
  !
  ! now run the recursion down
  !
  if (nterms .lt. 10) then
    nterms = 10
  endif

  v = 1
  vnext = 0

  vals(nterms) = 1
  vals(nterms+1) = 0
  
  do i = nterms,1,-1
    vals(i-1) = (2*i*x*vals(i) - (i+half)*vals(i+1))/(i-half)
  enddo

  ratio = q0/vals(0)
  do i = 0,nterms
    vals(i) = vals(i)*ratio
  enddo

  do i = 0,nterms
    if (abs(vals(i) .lt. eps)) then
      nterms = i
      return
    end if
  end do
  
  return
end subroutine q2leges_eps
  




subroutine q2leges(ier, nterms, x, xminus, vals, lvals, ntop)
  implicit real *8 (a-h,o-z)
  real *8 :: vals(0:lvals)
  real *8, parameter :: done=1

  !
  ! evaluate the first nterms+1 legendre functions
  ! this routine uses Miller's algorithm: recurse up, recurse
  ! down, and normalize. It will bomb if len is not large
  ! enough...
  !
  ! input:
  !   nterms - return values of Q_{n-1/2} from 0 through nterms
  !   x - the argument
  !   xminus - should be set to 1-x, evaluated as accurately as possible
  !   lvals - the total length of vals, should be sufficiently
  !       larger than nterms to handle an upward recursion
  !
  ! output:
  !   ier - error code
  !   vals - the values of Q_{n-1/2}
  !   ntop - the top value needed in the recursion
  !

  if (nterms .lt. 0) return
  ier = 0
  
  call q2lege01(x, xminus, q0, q1)

  !if (x .le. 0.001d0) then
  !  stop
  !end if

  if (nterms .eq. 0) then
    vals(0) = q0
    return
  endif

  if (nmax .eq. 1) then
    vals(0) = q0
    vals(1) = q1
    return
  endif

  !
  ! run the recurrence now, starting from nterms to ensure
  ! relative precision in all values
  !
  vprev = 0
  v = 1
  upbound = 1.0d20
  half=done/2

  do i=nterms,lvals-1
    vnext = (2*i*x*v - (i-half)*vprev)/(i+half)
    if (abs(vnext) .ge. upbound) then
      ntop = i+1
      exit
    endif
    vprev = v
    v = vnext
  enddo

  !
  ! now run the recursion down
  !
  if (ntop .lt. nterms) then
    ntop = nterms+5
  endif

  v = 1
  vnext = 0

  do i = ntop,nterms,-1
    vprev = (2*i*x*v - (i+half)*vnext)/(i-half)
    vnext = v
    v = vprev
  enddo

  vals(nterms) = vnext
  vals(nterms-1) = v
  do i = nterms-1,1,-1
    vals(i-1) = (2*i*x*vals(i) - (i+half)*vals(i+1))/(i-half)
  enddo

  ratio = q0/vals(0)
  do i = 0,nterms
    vals(i) = vals(i)*ratio
  enddo

  return
end subroutine q2leges





subroutine q2leges_all(ier, nterms, x, xminus, vals, ders, &
    der2s, lvals, ntop)
  implicit real *8 (a-h,o-z)
  real *8 :: vals(0:lvals), ders(0:lvals), der2s(0:lvals)
  real *8, parameter :: done=1

  !
  ! evaluate the first nterms+1 legendre functions
  !
  !   Q_{-1/2}, Q_{1/2}, ..., Q_{nterms-1/2}
  !
  ! and their derivatives...
  !
  ! input:
  !   nterms - return values of Q_{n-1/2} from 0 through nterms
  !   x - the argument
  !   xminus - should be set to 1-x, evaluated as accurately as possible
  !   lvals - the total length of vals, should be sufficiently
  !       larger than nterms to handle an upward recursion
  !
  ! output:
  !   ier - error code
  !   vals - the values of Q_{n-1/2}
  !   ders - the values of Q_{n-1/2}
  !   der2s - the values of Q_{n-1/2}
  !   ntop - the top value needed in the recursion
  !

  if (nterms .lt. 0) return

  call q2lege01_all(x, xminus, q0, q1, q0der, q1der, &
      q0der2, q1der2)

  if (nterms .eq. 0) then
    vals(0) = q0
    ders(0) = q0der
    der2s(0) = q0der2
    ier = 0
    return
  endif

  if (nmax .eq. 1) then
    vals(0) = q0
    ders(0) = q0der
    der2s(0) = q0der2
    vals(1) = q1
    ders(1) = q1der
    der2s(1) = q1der2
    ier = 0
    return
  endif

  !
  ! run the recurrence now, starting from nterms to ensure
  ! relative precision in all values
  !

  vals(nterms-1) = 0
  vals(nterms) = 1
  ders(nterms-1) = 0
  ders(nterms) = 1
  der2s(nterms-1) = 0
  der2s(nterms) = 1
  upbound = 1.0d20

  half=done/2

  do i=nterms,lvals-1
    vals(i+1) = (2*i*x*vals(i) - (i-half)*vals(i-1))/(i+half)
    ders(i+1) = (2*i*(vals(i) + x*ders(i)) &
        - (i-half)*ders(i-1))/(i+half)
    der2s(i+1) = (2*i*(2*ders(i) + x*der2s(i)) &
        - (i-half)*der2s(i-1))/(i+half)
    if (abs(vals(i+1)) .ge. upbound) then
      if (abs(der7) .ge. upbound) then
        if (abs(der27) .ge. upbound) then
          ntop = i+1
          exit
        endif
      endif
    endif
    der0 = der1
    der1 = der7
    der20 = der21
    der21 = der27
  enddo

  !
  ! now run the recursion down
  !
  vals(ntop) = 1
  vals(ntop+1) = 0

  ders(ntop) = 1
  ders(ntop+1) = 0

  der2s(ntop) = 1
  der2s(ntop+1) = 0

  do i = ntop,1,-1
    vals(i-1) = (2*i*x*vals(i) - (i+half)*vals(i+1))/(i-half)
    ders(i-1) = (2*i*(vals(i) + x*ders(i)) - &
        (i+half)*ders(i+1))/(i-half)
    der2s(i-1) = (2*i*(2*ders(i) + x*der2s(i)) - &
        (i+half)*der2s(i+1))/(i-half)
  enddo

  ratio = q0/vals(0)
  ratio1 = q0der/ders(0)
  ratio2 = q0der2/der2s(0)

  do i = 0,ntop
    vals(i) = vals(i)*ratio
    ders(i) = ders(i)*ratio1
    der2s(i) = der2s(i)*ratio2
  enddo

  return
end subroutine q2leges_all





subroutine q2lege01(x, xminus, val0, val1)
  implicit double precision (a-h,o-z)
  !
  ! this subroutine returns the value of the first two
  ! half-order legendre functions of the
  ! second kind, q_-1/2 and q_1/2. this routine requires that
  ! the evaluation point x be larger than 1, i.e. in the toroidal
  ! regime of the associated functions
  !
  ! input:
  !
  !   x - evaluation point, must be larger than 1 or else you
  !       will get garbage
  !
  ! output:
  !
  !   val0 - value of Q_{-1/2}(x)
  !   val1 - value of Q_{1/2}(x)
  !

  xarg = 2/(1+x)
  rk = sqrt(xarg)
  rkp = sqrt(xminus/(x+1))

  ifder = 0
  call ellipticke(rk, rkp, valk, vale, ifder, dk, de)
  val0 = rk*valk
  val1 = x*rk*valk - sqrt(2*(x+1))*vale

  return
end subroutine q2lege01





subroutine q2lege01_all(x, xminus, val0, val1, der0, der1, &
    der20, der21)
  implicit double precision (a-h,o-z)
  !
  ! this subroutine returns the value and derivative
  ! of the first two half-order legendre functions of the
  ! second kind, Q_-1/2 and Q_1/2. this routine requires that
  ! the evaluation point x be larger than 1, i.e. in the toroidal
  ! regime of the associated functions
  !
  ! input:
  !   x - evaluation point, must be larger than 1
  !   xminus - the value of x-1
  !
  ! output:
  !   val0 - value of Q_{-1/2}(x)
  !   val1 - value of Q_{1/2}(x)
  !   der0 - value of Q'_{-1/2}(x)
  !   der1 - value of Q'_{1/2}(x)
  !   der20 - value of Q''_{-1/2}(x)
  !   der21 - value of Q''_{1/2}(x)
  !

  xarg = 2/(1+x)
  rk = sqrt(xarg)
  rkp = sqrt(xminus/(x+1))

  ifder = 1
  call ellipticke(rk, rkp, valk, vale, ifder, dk, de)


  ! xarg=2/(1+x)
  ! xarg2=sqrt(xarg)
  ! xarg2p=sqrt(2-xarg2**2)
  ! ifder=1
  ! call ellipke(xarg2,xarg2p,valk,vale,ifder,fkd,fed)

  val0 = rk*valk
  val1 = x*rk*valk - sqrt(2*(x+1))*vale

  der0 = (val1 - x*val0)/2/(x+1)/xminus
  der1 = (-val0 + x*val1)/2/(x+1)/xminus

  der20 = (2*(x+1)*xminus*(der1-val0-x*der0) &
      - 4*x*(val1 - x*val0))/4/(x+1)**2/xminus**2
  der21 = -(2*(x+1)*xminus*(der0-val1-x*der1) &
      - 4*x*(val0 - x*val1))/4/(x+1)**2/xminus**2

  return
end subroutine q2lege01_all
