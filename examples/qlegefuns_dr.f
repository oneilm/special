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
        implicit real *8 (a-h,o-z)
        real *8 qfuns(100000),qfuns3(100000)
        real *8 xs(10000), whts(10000)
        complex *16 z, zfuns(0:100000), ima, zfuns2(0:100000)
        complex *16 zint, q0, q1, vals(100000), pot
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
c 
        call prin2('qfuns after qlegfuns is*',qfuns,n+1)
c 
c        now, evaluate the functions via the qneval
c 
        call qneval(x, n, qfuns3)
c 
        call prin2('and qfuns3 via qneval*',qfuns3,n+1)
c 
        do i = 1,n+1
ccc          qfuns(i) = (qfuns(i) - qfuns3(i))/qfuns3(i)
        enddo
c
cccc        call prin2('relative diffs = *', qfuns, n+1)


        call prinf2('from qneval:*', x, 0)
        do i = 1,n+1
          write(6,*) 'i = ', i-1, 'val = ', qfuns3(i)
          write(13,*) 'i = ', i-1, 'val = ', qfuns3(i)
        enddo


c
c       test the complex routine, which just runs the forward
c       recursion
c
c       NOTE: several calculations have been checked against
c       40 digit calculations run in Mathematica and reliably agree
c       to 13 or 14 digits, unless x is very close to 1 or -1
c
        eps = 1.0d-13
        z = .5d0 - ima*eps
        call prin2('z = *', z, 2)
        call zqneval(z, n, zfuns)
c
        call prinf2('from zqneval:*', x, 0)
        do i = 0,n
          write(6,*) 'i = ', i-1, 'val = ', zfuns(i)
          write(13,*) 'i = ', i-1, 'val = ', zfuns(i)
        enddo
c

c
c       compute the hilbert transform of P_n and compare to see if we're
c       getting the right answer
c
        z = 3.5+ima/10
        k = 50
        ifwhts = 1
        call legewhts(k, xs, whts, ifwhts)
        call prin2('roots = *', xs, k)
        call prin2('weights = *', whts, k)
c
        zint = 0
        do i = 1,k
          zint = zint + whts(i)/(-xs(i) + z)/2
        enddo
c
        call prin2('Cauchy integral of 1 = *', zint, 2)
c
        call zqlege01(z, q0, q1)
        call prin2('Q_0 = *', q0, 2)

        call prin2('diff = *', zint - q0, 2)

c
c       try the cauchy integral of P_m
c
        print *
        print *
        m = 5
        call prinf('degree of P_m, m = *', m, 1)
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
