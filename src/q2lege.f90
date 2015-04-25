        !
        ! (c) 2015 Mike O'Neil
        !
        ! this code containes routines for evaluate the half-order legendre
        ! polynomials of the second-kind
        !

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

        call q2lege01(x, q0, q1)

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

        !!!call prin2('upbound = *', upbound, 1)
        !!!call prin2('after fwd recursion, vals = *', vals, ntop+1)

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





      subroutine q2lege01(x, val0, val1)
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
        
        xarg=2/(1+x)
        xarg2=sqrt(xarg)
        call elliptic_ke(xarg2,valk,vale)

        val0=sqrt(xarg)*valk
        val1=x*sqrt(xarg)*valk-sqrt(2*(x+1))*vale        

        return
      end subroutine q2lege01
