c       This is a driver file for Cauchy transforms of Legendre polys
c
c       Copyright (c) 2014 Michael O'Neil
c                          oneil@cims.nyu.edu
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
c
c
        implicit real *8 (a-h,o-z)
        real *8 qfuns(100000),qfuns3(100000)
        real *8 xs(10000), whts(10000)
        real *8 errs(10000)
        complex *16 z, zfuns(0:100000), ima, zfuns2(0:100000)
        complex *16 zint, zint1, q0, q1, vals(100000), pot
        complex *16 cvals(100000), cvals2(10000), z2
c 

        call prini(6,13)
        ima = (0,1)
        done = 1
        pi = 4*atan(done)
C 
C       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n'
        READ *, n
        CALL PRInf('n=*',n,1 )
C 
        x = .5d0
        call prin2('x = *', x, 1)
        call prin2('and x-1=*',x-1,1)
c 
c        evaluate Q_n (x) via qlegfuns
c 
        call qlegfuns(x,n,qfuns)
        call prin2('qfuns after qlegfuns is*',qfuns,n+1)
c 
c        now, evaluate the functions via the qneval
c 
        call qneval(x, n, qfuns3)
        call prin2('and qfuns3 via qneval*',qfuns3,n+1)
c 
        do i = 1,n+1
          errs(i) = (qfuns(i) - qfuns3(i))/qfuns3(i)
        enddo
c
        call prin2('relative diffs = *', errs, n+1)


        call prinf2('from qneval:*', x, 0)
        do i = 1,n+1
          write(6,*) 'i = ', i-1, 'val = ', qfuns3(i)
          write(13,*) 'i = ', i-1, 'val = ', qfuns3(i)
        enddo


c
c       test the complex routine for Q_0 and Q_1
c
        print *
        print *

        small = 1.0d-1
        small = 1.0d0
        small = 1.0d-6
        z = .5d0 - ima*small
        call prin2('z = *', z, 2)
        call zqlege01(z, q0, q1)
c
        call prin2('from zqlege01, q0 = *', q0, 2)
        call prin2('from zqlege01, q1 = *', q1, 2)

c
c       and compare with the actual hilbert transform of P_0 and P_1
c
        k = 32 
        ifwhts = 1
        call legewhts(k, xs, whts, ifwhts)
c
        zint = 0
        zint1 = 0
        do i = 1,k
          zint = zint + whts(i)/(-xs(i) + z)/2
          zint1 = zint1 + whts(i)*xs(i)/(-xs(i) + z)/2
        enddo
c
        call prin2('Cauchy integral of 1 = *', zint, 2)
        call prin2('Cauchy integral of x = *', zint1, 2)
c
        call prin2('error in Q_0 = *', q0-zint, 2)
        call prin2('error in Q_1 = *', q1-zint1, 2)

c
c       try the cauchy integral of P_m
c
        print *
        print *
        m = 16 
        call prinf('degree of P_m, m = *', m, 1)
        call prin2('target = *', z, 2)
c
        zint = 0
        do i = 1,k
          call legepol(xs(i), m, pol, der)
          zint = zint + whts(i)*pol/(-xs(i) + z)/2
        enddo
c
        call prin2('Cauchy integral of P_m = *', zint, 2)
c
        call zqneval(z, m, zfuns)
        call prin2('Q_m(z) = *', zfuns(m), 2)
        call prin2('diff = *', zint - zfuns(m), 2)
c
c       find the ellipse on which you can just evaluate the
c       Cauchy transform
c
        ntheta = 200
        sc = 2.0d0
        a = 2d0
        b = 1d0
        a = a*sc
        b = b*sc
c
        do j = 1,ntheta
          theta = 2*pi*(j-1)/ntheta
          z2 = a*cos(theta) + ima*b*sin(theta)
          zint = 0
          do i = 1,k
            call legepol(xs(i), m, pol, der)
            zint = zint + whts(i)*pol/(z2 - xs(i))/2
          enddo
          cvals(j) = zint
          call zqneval(z2, m, zfuns)
          cvals2(j) = zfuns(m)
        enddo
c
        do j = 1,ntheta
          errs(j) = cvals(j) - cvals2(j)
        enddo
c
cccc        call prin2('cvals = *', cvals, 2*ntheta)
        call prin2('from recursion, cvals2 = *', cvals2, 2*ntheta)
cccc        call prin2('errs on ellipse = *', errs, ntheta)

        errmax = -1
        do i = 1,ntheta
          if (abs(errs(i)) .gt. errmax) errmax = abs(errs(i))
        enddo
        call prin2('maximum error on the ellipse = *', errmax, 1)
c
        call zqneval_up(z, m, zfuns)
        call prin2('zfuns_up = *', zfuns, 2*m+2)

        call zqneval(z, m, zfuns)
        call prin2('zfuns = *', zfuns, 2*m+2)
        write(6,*) 'm = ', m, 'q_m = ', zfuns(m)

        stop


c
c       test the routine for computing the Cauchy transform
c       of a Legendre expansion
c
        do i = 1,k
          vals(i) = cos(pi*xs(i)) + ima*sin(pi*xs(i))
          vals(i) = exp(ima*14.5d0*xs(i))
        enddo
c
        z = 1.25d0 + ima*.5d0
        call hilbert_legendre(z, k, vals, pot)
        call prin2('from hilbert_legendre, pot = *', pot, 2)
c
        zint = 0
        do i = 1,k
          zint = zint + whts(i)*vals(i)/(z - xs(i))
        enddo
c
        call prin2('discretely, pot = *', zint, 2)
        call prin2('difference = *', pot - zint, 2)
        
c
c       test the piecewise routine
c
        print *
        print *
        print *
c
        a = .2d0
        b = .8d0
        h = b-a

c        call prin2('before, xs = *', xs, k)
c        call prin2('before, whts = *', whts, k)

        do i = 1,k
          xs(i) = (xs(i)+1)/2*h + a
          whts(i) = whts(i)/2*h
        enddo
c
c        call prin2('after, xs = *', xs, k)
c        call prin2('after, whts = *', whts, k)
c        stop

        do i = 1,k
          vals(i) = cos(pi*xs(i)) + ima*sin(pi*xs(i))
cc          vals(i) = exp(ima*14.5d0*xs(i))
cc          vals(i) = 1
        enddo
c
        z = 1.25d0 + ima*.5d0
        call hilbert_legendre_ab(a, b, z, k, vals, pot)
        call prin2('from hilbert_legendre_ab, pot = *', pot, 2)
c
        zint = 0
        do i = 1,k
          zint = zint + whts(i)*vals(i)/(z - xs(i))
        enddo
c
        call prin2('discretely, pot = *', zint, 2)
        call prin2('difference = *', pot - zint, 2)


        stop
        end
c 
c 
c
c
c
        subroutine prinzfull(n, zs)
        implicit real *8 (a-h,o-z)
        complex *16 zs(n)
c
        do i = 1,n
          x = zs(i)
          y = imag(zs(i))
          write(6,1200) x, y 
        enddo
c
 1200 FORMAT(e23.16, ' + ', e23.16, '*I')
c
        return
        end
c
c       
c
c
c
        subroutine qlegfuns(x,n,qfuns)
        implicit real *8 (a-h,o-z)
        dimension qfuns(0:n)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c 
c        construct Q_0(x),Q_1(x)
c 
        done=1
        i0=0
        delta=1.0d-20
c 
        d=log( (done+x+ima*delta)/(done-x+ima*delta) )
c 
        qfuns(i0)=d/2
        qfuns(1)=x/2*d-1
c 
c       recurse up
c 
        do 1200 i=1,n-1
c 
        qfuns(i+1)=( (2*i+1)*x*qfuns(i)-i*qfuns(i-1) ) /(i+1)
 1200 continue
c 
        return
        end
c 
