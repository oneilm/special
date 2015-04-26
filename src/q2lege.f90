        !
        ! (c) 2015 Mike O'Neil
        !
        ! this code containes routines for evaluate the half-order legendre
        ! polynomials of the second-kind
        !




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

        pi = 4*atan(done)
        fac = 4*pi*pi*sqrt(src(1)*targ(1))
        do i = 0,ntop
          vals(i) = vals(i)/fac
        enddo
        
        return
      end subroutine torlaps




      subroutine q2leges(ier, nterms, x, xminus, vals, lvals, ntop)
        implicit real *8 (a-h,o-z)
        real *8 :: vals(0:lvals)
        real *8, parameter :: done=1

        !
        ! evaluate the first nterms+1 legendre functions
        !
        !   Q_{-1/2}, Q_{1/2}, ..., Q_{nterms-1/2}
        !
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

        call q2lege01(x, xminus, q0, q1)

        !if (x .le. 0.001d0) then
        !  stop
        !end if
        
        if (nterms .eq. 0) then
          vals(0) = q0
          ier = 0
          return
        endif
        
        if (nmax .eq. 1) then
          vals(0) = q0
          vals(1) = q1
          ier = 0
          return
        endif

        !
        ! run the recurrence now, starting from nterms to ensure
        ! relative precision in all values
        !

        vals(nterms-1) = 0
        vals(nterms) = 1
        upbound = 1.0d20
        
        half=done/2
        
        do i=nterms,lvals-1
          vals(i+1) = (2*i*x*vals(i) - (i-half)*vals(i-1))/(i+half)
          if (abs(vals(i+1)) .ge. upbound) then
            ntop = i+1
            exit
          endif
        enddo

        !
        ! now run the recursion down
        !
        vals(ntop) = 1
        vals(ntop+1) = 0

        do i = ntop,1,-1
          vals(i-1) = (2*i*x*vals(i) - (i+half)*vals(i+1))/(i-half)
        enddo

        !!!call prin2('after bckwd recursion, vals = *', vals, ntop+1)
        ratio = q0/vals(0)

        do i = 0,ntop
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

        call q2lege01(x, xminus, q0, q1)

        !if (x .le. 0.001d0) then
        !  stop
        !end if
        
        if (nterms .eq. 0) then
          vals(0) = q0
          stop
          ier = 0
          return
        endif
        
        if (nmax .eq. 1) then
          vals(0) = q0
          stop
          vals(1) = q1
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

        !!!call prin2('after bckwd recursion, vals = *', vals, ntop+1)
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

        der20 = -(2*(x+1)*xminus*(der1-val0-x*der0) &
            - 4*x*(val1 - x*val0))/4/(x+1)**2/xminus**2
        der21 = -1

        return
      end subroutine q2lege01_all
