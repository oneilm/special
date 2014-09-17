c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c      This is the end of the testing code and the beginning of the
c      code for the evaluation of the Legendre functions Q_i
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
        subroutine qlege01(x, q0, q1)
        implicit real *8 (a-h,o-z)
c
c       this routine evaluates Q_0 and Q_1 for |x| .neq. 1
c
        done = 1
        d = log( abs((done+x)/(done-x)) )
        q0 = d/2
        q1 = x/2*d-1
c
        return
        end
c
c
c
c
c
        subroutine zqlege01(z, q0, q1)
        implicit real *8 (a-h,o-z)
        complex *16 z, q0, q1, d
c
c       this routine evaluates Q_0 and Q_1 for |x| .neq. 1
c
        done = 1
        d = log((done+z)/(done-z))
        q0 = d/2
        q1 = x/2*d-1
c
        return
        end
c
c
c
c
c
        subroutine zqneval(z, n, qfuns)
        implicit real *8 (a-h,o-z)
        complex *16 z, qfuns(0:n), d
c
c       This subroutine evaluates Q_0(z), ..., Q_n(z) for arbitrary
c       complex values of z. The upward recurrence is used everywhere, as
c       Q_n blows up for values of z off the real-line. If you want to
c       evaluate Q_n on the real-line, use the routine qneval - this routine
c       will NOT return the principal value of Q_n for z with arbitrarily
c       small imaginary parts.
c
c       construct Q_0(z) and Q_1(z)
c 
        done=1
        i0=0
        d=log( (done+z)/(z-1) ) 
c 
        qfuns(i0) = d/2
        qfuns(1) = z/2*d-1
c 
c       recurse up
c 
        do i=1,n-1
          qfuns(i+1)=( (2*i+1)*z*qfuns(i)-i*qfuns(i-1) ) /(i+1)
        enddo
c 
        return
        end
c
c
c
c
c
        subroutine qneval(x, n, qfuns)
        implicit real *8 (a-h,o-z)
        real *8 qfuns(0:1)
c 
c        for user-specified REAL x and integer n, this subroutine
c        evaluates the Legendre functions Q_0, Q_1, . . . ,Q_n at
c        the point x anywhere on the real line (except |x|=1).
c
c        If |x| > 1, what is actually returned is the average of
c        Q(x+i\epsilon) and Q(x-i\epsilon) - which is real. Off the real-line,
c        the imaginary party of Q_n explodes exponentially
c
c                        Input parameters:
c 
c  x - the point where the functions Q_i will be evaluated
c  n - maximum order of Q to be evaluated
c 
c                        Output parameters:
c 
c  qfuns - the values of the Legendre functions at the
c      point x (n+1 of them things) - this actually needs to be
c      something like n+3 because of up/down recurrences
c
c  NOTE: the indexing of qfuns begins at 0 !!!
c
c  NOTE 2: you will lose digits if x is extremely close to +/- 1 because
c     of floating point errors in subtraction - this cannot be avoided.
c     See the special purpose routine ... for this case.
c 
c        . . . if x is inside the interval [-1,1], or sufficiently
c              close to 1 or -1, evaluate the Legenedre functions
c              via the simple recursion
c 
        if(dabs(x) .lt. 1) goto 1100
c 
        delta=dabs(x)-1
        b=1+delta+sqrt(2*delta+delta**2)
c 
        if( (n+1)*log(b) .gt. 2.3) goto 2000
c 
 1100 continue
c 
c        construct Q_0(x),Q_1(x)
c 
        done=1
        i0=0
        d=log( abs((done+x)/(done-x)) )
c 
        qfuns(i0)=d/2
        qfuns(1)=x/2*d-1
c 
c       recurse up
c 
        do i=1,n-1
          qfuns(i+1)=( (2*i+1)*x*qfuns(i)-i*qfuns(i-1) ) /(i+1)
        enddo
c 
        return
c 
 2000 continue
c 
c       the point x is outside the interval [-1,1].
c       determine how far we will have to start to get it by
c       the combination of recursing up and scaling
c 
        eps=1.0d-20
c 
        delta=dabs(x)-1
        b=1+delta+sqrt(2*delta+delta**2)
        nn=-log(eps)/log(b)+1
c 
cccc        call prinf('nn as obtained in qninte is*',nn,1)
c 
c       use recursion down coupled with scaling to get the
c       Legendre functions Q_i
c 
        done=1
        i0=0
        d=log( abs((done+x)/(done-x)) )
c 
        q0=d/2
c 
        do i=i0,n
          qfuns(i)=0
        enddo

c 
c       if nn is greater than n+1, do the preliminary recursion
c       to get the starting values for qfuns(n+1),qfuns(n)
c 
        fip1=0
        fi  =1
c 
        if(nn .le. n+1) goto 2150
c 
        do i=nn,n+1,-1
          fim1=( (2*i+1)*x*fi-(i+1)*fip1 ) / i
          fip1=fi
          fi=fim1
        enddo
c 
        nn=n
 2150 continue
c 
        qfuns(nn+1)=fip1
        qfuns(nn)=fi
c 
c       recurse down starting with i=n
c 
        do i=nn,1,-1
          qfuns(i-1)=( (2*i+1)*x*qfuns(i)-(i+1)*qfuns(i+1) ) / i
        enddo
c 
c        scale the values of the recursed (unscaled) Q_i to
c        obtain the correct one
c 
        rat=q0/qfuns(i0)
cccc        call prin2('and rat=*',rat,1)
c 
        do i=i0,n+1
          qfuns(i)=rat*qfuns(i)
        enddo
c 
        return
        end
  
  
