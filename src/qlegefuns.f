c       This is a library file for Cauchy transforms of Legendre polys
c
c       Copyright (c) 2014  Vladimir Rokhlin
c                           Michael O'Neil
c                           oneil@cims.nyu.edu
c
c       This program is free software: you can redistribute it and/or modify
c       it under the terms of the GNU General Public License as published by
c       the Free Software Foundation, either version 2 of the License, or
c       (at your option) any later version.
c
c       This program is distributed in the hope that it will be useful,
c       but WITHOUT ANY WARRANTY; without even the implied warranty of
c       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c       GNU General Public License for more details.
c
c       You should have received a copy of the GNU General Public License
c       along with this program.  If not, see <http://www.gnu.org/licenses/>.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c      This is the end of the testing code and the beginning of the
c      code for the evaluation of the Legendre functions Q_n and other
c      routines...
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c




        subroutine hilbert_legendre_ab(a, b, z, k, vals, pot)
        implicit real *8 (a-h,o-z)
        complex *16 ima, z, vals(k), pot, z2
        data ima/(0,1)/
c
c       this routine computes the integral:
c
c          pot = \int_a^b  f(x)/(z - x) dx
c
c       where f(x) is a (complex) function sampled at k legendre 
c       nodes on [a,b]. The target z can be arbitrary, as the
c       Hilbert transform is evaluated using the following P_n, Q_n
c       identity
c
c          Q_n = \frac{1}{2} \int_{-1}^1 P_n / (z - x) dx
c
c       all of the following calculations are done assuming that we
c       define the branch cut in Q_n along the interval [-1,1] so that
c       there is a jump in the Hilbert transform across [-1,1].
c
c       input:
c         z - target point
c         k - number of legendre nodes at which f is evaluated
c         vals - the values of f at the k legendre nodes, complex valued
c
c       output:
c         pot - the Hilbert transform of f on [a,b]
c
c       NOTE: this routine should more or less be as accurate as you
c       have resolved the function f on [a,b] - meaning that if the
c       Legendre expansion of f has decayed to 10^{-12}, then the
c       integral is accurate to 10^{-12} (or better}
c
c
c       first: compute the equivalent target relative to [-1,1]
c
        done = 1
        z2 = 2*(z-a)/(b-a) - done
cccc        call prin2('z2 = *', z2, 2)
        call hilbert_legendre(z2, k, vals, pot)
c
        return
        end
c
c
c
c
c
        subroutine hilbert_legendre(z, k, vals, pot)
        implicit real *8 (a-h,o-z)
        complex *16 ima, z, vals(k), pot
c
        data ima/(0,1)/
        real *8 xs(1000), ys(1000), xnodes(1000)
        real *8 u(1000000), v(1000000), whts(1000)
        real *8 coefs_real(0:1000), coefs_imag(0:1000)
        complex *16 coefs(0:1000), qfuns(0:1000), zterms(0:1000)
c
c       this routine computes the integral:
c
c          pot = \int_{-1}^1  f(x)/(z - x) dx
c
c       where f(x) is a (complex) function sampled at k legendre 
c       nodes on [-1,1]. The target z can be arbitrary, as the
c       Hilbert transform is evaluated using the following P_n, Q_n
c       identity
c
c          Q_n = \frac{1}{2} \int_{-1}^1 P_n / (z - x) dx
c
c       all of the following calculations are done assuming that we
c       define the branch cut in Q_n along the interval [-1,1] so that
c       there is a jump in the Hilbert transform across [-1,1].
c
c       input:
c         z - target point
c         k - number of legendre nodes at which f is evaluated
c         vals - the values of f at the k legendre nodes, complex valued
c
c       output:
c         pot - the Hilbert transform of f
c
c
c       get the expansion coefficients of vals
c
        if (k .gt. 900) then
          call prinf('k = *', k, 1)
          call prinf('bomb!! k is too large*', ima, 0)
          stop
        endif
c
        itype = 2
        call legeexps(itype, k, xnodes, u, v, whts)
c        
        do i = 1,k
          xs(i) = vals(i)
          ys(i) = -ima*vals(i)
        enddo
c
        call matvec(k, u, xs, coefs_real)
        call matvec(k, u, ys, coefs_imag)
c
        do i = 0,k-1
          coefs(i) = coefs_real(i) + ima*coefs_imag(i)
        enddo

c
c       evaluate the first k-1 Q_n's
c
        call zqneval(z, k-1, qfuns)
c
cccc        call prin2('coefs = *', coefs, 2*nterms)
cccc        call prin2('qfuns = *', qfuns, 2*nterms)

c
c       and evaluate the expansion
c
        pot = 0
        do i = 0,k-1
          pot = pot + 2*coefs(i)*qfuns(i)
cccc          zterms(i) = 2*coefs(i)*qfuns(i)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine matvec(n, a, x, y)
        implicit real *8 (a-h,o-z)
        real *8 a(n,n), x(n), y(n)
c
        do i = 1,n
          d = 0
          do j = 1,n
            d = d + a(i,j)*x(j)
          enddo
          y(i) = d
        enddo
c
        return
        end
c
c
c
c
c
        subroutine qlege01(x, q0, q1)
        implicit real *8 (a-h,o-z)
c
c       this routine evaluates Q_0 and Q_1 on R for |x| .neq. 1
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
c       this routine evaluates Q_0 and Q_1 for complex argument
c       *off* of the real line - this assumes the branch cut is
c       on [-1,1]
c
c       input:
c         z - a complex number
c
c       output:
c         q0 - the value of Q_0(z)  
c         q1 - the value of Q_1(z)  
c
        done = 1
        d = log((done+z)/(z - done))
        q0 = d/2
        q1 = z/2*d-1
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
        complex *16 z, qfuns(0:n), d, q0, q1, qnext, qn
        complex *16 ratio
c
c       Evaluate Q_0 .. Q_n at the point z in the complex plane.
c       Taking the branch cut along [-1,1] means that Q_n decays
c       as z \to \infty. We must use an up, down, and normalize
c       recursion.
c
c       input:
c         z - the complex target point
c         n - max order of Q_n to evaluate
c
c       ouput:
c         qfuns - returns Q_0(z), ..., Q_n(z) with the INDEX 
c             STARTING AT ZERO
c
c       first, construct Q_0(z) and Q_1(z)
c
c
        call zqlege01(z, q0, q1)
        qfuns(0) = q0
        if (n .eq. 0) return
c
        qfuns(1) = q1
        if (n .eq. 1) return
c
c       recurse up until things have blown up sufficiently
c 
        qlarge = 1.0d16/abs(q0)
        qlarge = qlarge**2
c
        nmax = 10000
        ntop = 0
c
        do i=1,nmax
          qn = ( (2*i+1)*z*q1-i*q0 ) /(i+1) 
          q0 = q1
          q1 = qn
          if (abs(qn) .gt. qlarge) then
            ntop = i+1
            goto 1300
          endif 
        enddo
c 
 1300 continue

c
c       if it did not blowup, just recurse up
c
        if (ntop .eq. 0) then
          call zqneval_up(z, n, qfuns)
          return
        endif

c
c       otherwise recurse down
c
        ntop = ntop+5
cccc        call prinf('ntop = *', ntop, 1)
c
        qfuns(ntop) = 1
        qfuns(ntop+1) = 0
c
        do i = ntop,1,-1
          qfuns(i-1) = (2*i+1)*z*qfuns(i)/i - (i+1)*qfuns(i+1)/i
        enddo
c
c       and normalize
c
        call zqlege01(z, q0, q1)
        ratio = q0/qfuns(0)
c
        do i = 0,ntop
          qfuns(i) = qfuns(i)*ratio
        enddo

        return
        end
        
c
c
c
c
c
        subroutine zqneval_up(z, n, qfuns)
        implicit real *8 (a-h,o-z)
        complex *16 z, qfuns(0:n), d, q0, q1
c
c       Evaluate Q_0 .. Q_n at the point z in the complex plane
c       ONLY USING AN UPWARD RECURRENCE!!! (see zqneval).
c       Branch cut is assumed to be on [-1,1].
c
c       input:
c         z - the complex target point
c         n - max order of Q_n to evaluate
c
c       ouput:
c         qfuns - returns Q_0(z), ..., Q_n(z) with the INDEX 
c             STARTING AT ZERO
c
c       NOTE: 
c       If you want to evaluate Q_n on the real-line, use the 
c       routine qneval - this routine will NOT return the principal 
c       value of Q_n for z with arbitrarily small imaginary parts.
c
c       NOTE 2: The branch cut is taken to be on [-1,1], consistent with
c       the identity:
c
c               Q_n(z) = 1/2 \int_{-1}^1 P_n(x)/(z - x) dx
c
c       first, construct Q_0(z) and Q_1(z)
c 
        call zqlege01(z, q0, q1)
        qfuns(0) = q0
        if (n .eq. 0) return
c
        qfuns(1) = q1
        if (n .eq. 1) return
c
c       recurse up until 
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
c     x - the point where the functions Q_i will be evaluated
c     n - maximum order of Q to be evaluated
c 
c                        Output parameters:
c 
c     qfuns - the values of the Legendre functions at the
c          point x (n+1 of them things) - this actually needs to be
c          something like n+3 because of up/down recurrences
c
c     NOTE: the indexing of qfuns begins at 0 !!!
c
c     NOTE 2: you will lose digits if x is extremely close to +/- 1 because
c         of floating point errors in subtraction - this cannot be avoided.
c 
c        . . . if x is inside the interval [-1,1], or sufficiently
c              close to 1 or -1, evaluate the Legenedre functions
c              via the simple recursion
c
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
