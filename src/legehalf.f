c       This is a library file for half-order legendre functions.
c
c       Copyright (c) 2014  Michael O'Neil
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       this is the end of the debugging code and the beginning
c       of the half order legendre functions code
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       the following subroutines are user-callable:
c
c         torlaps - evaluate a bunch of fourier modes of
c             1/(4*pi*r)
c
c         torlap0 - evaluate the 0th order fourier mode for 3D laplace
c
c         torlap1 - evaluate the 1st order fourier mode for 3D laplace
c
c         p2lege01 - evaluate P_{-1/2} and P_{1/2}, the first two half-degree
c           legendre functions of the first kind
c
c         p2lege01der - evaluate the derivatives of P_{-1/2} and P_{1/2}
c
c         p2leges - evaluate the first nmax+1 functions P_{-1/2} through
c             P_{nmax-1/2}
c
c         p2legeders - evaluate the first nmax+1 functions P_{-1/2} through
c             P_{nmax-1/2} along with their derivatives
c
c         q2lege01 - evaluate Q_{-1/2} and Q_{1/2}, the first two half-degree
c           legendre functions of the second kind
c
c         q2lege01der - evaluate Q_{-1/2} and Q_{1/2} along with
c             their derivatives
c
c         q2legeratio - computes the ratio Q_{n-1/2}(x)/Q_{n-3/2}(x)
c             using a continued fraction expansion
c
c         q2leges - evaluate the first nmax+1 Q_{n-1/2} using
c             an up-down recursion, slower but accurate
c
c         pq2leges - evaluates the first nmax+1 half-degree legendre
c           functions of the first and second kinds of order 0
c
c         pq2legeders - evaluates the first nmax+1 half-degree legendre
c           functions of the first and second kinds of order 0 along
c           with their derivatives
c
c         ellipke - finally, a robust evaluator for elliptic integrals
c             which should provide near machine precision for
c             all values in [0,1), pending you provide k and k' (see
c             abramowitz & stegun)
c
c         elliptic_ke - tabulated elliptic integral routine, good
c             on the interal [0,1023/1024), will lose digits near
c             k=1 because of cancellation (i.e. the number 1-k is 
c             computed for use in a series expansion)
c
c         ke_eva - relies on elliptic_ke, but allows for more accurate
c             computation near the singularity at k=1
c
c       note - in all the routines that follow, we follow the convention
c           that the coordinate system is cylindrical, with the target
c           located at r,z,theta and the source at r0,z0,0 - only the
c           difference in theta matters, and the axisymmetric kernels
c           integrate it out.
c
c       note 2 - the only routine that has not been relatively
c           robustly tested is q2legeratio, which uses a continued
c           fraction expansion to evaluate the ratio of two
c           successive q2 functions
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
        subroutine q2leges(eps, x, nmax, qvals)
        implicit real *8 (a-h, o-z)
        real *8 qvals(0:1)
c
        done = 1
        blowup = 1/eps*10
c
c       run the recursion up until it blows, then run it down
c
        call q2lege01(x, val0, val1)
c
        ntop = 1
        half=done/2
c
        do i = 1,100000
          val =2*i*x*val1/(i+half)-(i-half)*val0/(i+half)
          val0 = val1
          val1 = val
          if (abs(val) .gt. blowup) ntop = i+2
        enddo
c
c       now run the recursion down
c

c
        return
        end
c
c
c
c
c
        subroutine torlaps(r,z,r0,z0,nmax,vals,ifder,drs,dzs,
     1      dr0s,dz0s)
        implicit real *8 (a-h,o-z)
        real *8 vals(0:1),drs(0:1),dzs(0:1),dr0s(0:1),dz0s(0:1)
        real *8, allocatable :: pvals(:),qders(:)
c
c       this routine calculates the first nmax+1 laplace
c       fourier modes, i.e. the functions
c
c           1/(2*pi) 1/sqrt(r*r0) Q_{n-1/2}(\xi)
c
c       for n=0, ..., nmax. this function is equivalent to
c       calculating the integral
c
c           \int_0^{2*pi} 1/(4*pi*R) * exp(-i*n*\theta) d\theta
c
c       input:
c
c         r,z - cylindrical coordinates of the target point
c         r0,z0 - cylindrical coordinates of the source point
c         nmax - the highest fourier mode to evaluate, modes
c             0 through nmax will be returned (note that the
c             arrays are dimensioned to begin at 0!!)
c         ifder - switch to evaluate partial derivatives as well
c
c       output:
c
c         vals - values of the fourier modes
c         drs,dzs - partial derivatives at the target
c         dr0s,dz0s - partial derivatives at the source
c
c
        allocate(pvals(0:nmax+5))
        allocate(qders(0:nmax+5))
        nmax2=nmax+2
c
        done=1
        h2=done/2
        pi=4*atan(done)
c
        if (nmax .eq. 0) then
c
            call torlap0(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
            vals(0)=val
c
            if (ifder .ne. 1) return
c
            drs(0)=dr
            dzs(0)=dz
            dr0s(0)=dr0
            dz0s(0)=dz0
            return
c
        endif
c
        if (nmax .eq. 1) then
c
            call torlap0(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
            vals(0)=val
            call torlap1(r,z,r0,z0,val1,ifder,dr1,dz1,dr01,dz01)
            vals(1)=val1
c
            if (ifder .ne. 1) return
c
            drs(0)=dr
            dzs(0)=dz
            dr0s(0)=dr0
            dz0s(0)=dz0
            drs(1)=dr1
            dzs(1)=dz1
            dr0s(1)=dr01
            dz0s(1)=dz01
            return
c
        endif

c
c       nmax > 1, so run the recursion now
c       . . . first evaluate p_{-1/2} and p_{1/2} accurately as possible
c
        zdiff=z-z0
        xi=(r*r+r0*r0+zdiff**2)/2/r/r0
        rmu=sqrt(2/(xi+1))
c
        xi1=((r-r0)**2+zdiff**2)/2/r/r0
        x0=sqrt(xi1/(xi+1))
        x0p=rmu
        call ellipke(x0,x0p,valk0,vale0,ifder,dk0,de0)
        p0=2/pi*rmu*valk0
c
        ximinus=(r**2-r0**2)**2+2*r**2*zdiff**2+2*r0**2*zdiff**2
     1      +zdiff**4
        ximinus=ximinus/4/r**2/r0**2
        x0=sqrt(ximinus)
c
        x1=sqrt(2*x0/(xi+x0))
        x1p=1/(xi+x0)
        call ellipke(x1,x1p,valk1,vale1,ifder,dk1,de1)
        p1=2/pi*sqrt(xi+x0)*vale1
c
        pvals(0)=p0
        pvals(1)=p1        

c
c       . . . run the recurrences now 
c
        do 1800 i=1,nmax2-1
        coef0=2*i/(i+h2)
        coef1=(i-h2)/(i+h2)
        pvals(i+1)=coef0*xi*pvals(i)-coef1*pvals(i-1)
 1800 continue
c

        call q2legeratio(xi,nmax2,ratio,err)
        call prin2('inside torlaps, error in ratio=*',err,1)
c
        d=(nmax2-h2)*(pvals(nmax2)-ratio*pvals(nmax2-1))
        vals(nmax2-1)=1/d
        vals(nmax2)=vals(nmax2-1)*ratio

c
c       and run the downward recurrence for the q's
c       
        do 2200 i=nmax2-1,1,-1
        vals(i-1)=2*i*xi*vals(i)/(i-h2)-(i+h2)*vals(i+1)/(i-h2)
 2200 continue
c
        if (ifder .ne. 1) then
            scale=done/(2*pi)/sqrt(r*r0)
            do 2400 i=0,nmax
            vals(i)=vals(i)*scale
 2400 continue
            return
        endif

c
c       now evaluate the derivatives if need be...
c
        i=nmax2
        qders(i)=(2*i-done)/2/ximinus*(xi*vals(i)-vals(i-1))
c
        i=nmax2-1
        qders(i)=(2*i-done)/2/ximinus*(xi*vals(i)-vals(i-1))
c
        do 3200 i=nmax2-1,1,-1
        coef0=2*i/(i-h2)
        coef1=(i+h2)/(i-h2)
        qders(i-1)=coef0*(vals(i)+xi*qders(i))-coef1*qders(i+1)
 3200 continue

c
c       now scale the vals and qders array to generate 
c       all the partial derivatives
c
        scale=done/(2*pi)/sqrt(r*r0)
        scale2r=-scale/2/r
        scale2r0=-scale/2/r0
c
        dxidr=(r**2-r0**2-zdiff**2)/2/r**2/r0
        dxidr0=(r0**2-r**2-zdiff**2)/2/r/r0**2
        dxidz=2*zdiff/2/r/r0
        dxidz0=-dxidz
c
        do 3600 i=0,nmax
        drs(i)=scale2r*vals(i)+scale*qders(i)*dxidr
        dzs(i)=scale*qders(i)*dxidz
        dr0s(i)=scale2r0*vals(i)+scale*qders(i)*dxidr0
        dz0s(i)=scale*qders(i)*dxidz0
        vals(i)=scale*vals(i)
 3600 continue
c
        return
        end
c
c
c
c
c
        subroutine torlap0(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
        implicit real *8 (a-h,o-z)
c
c       this routine evaluates the 0th order laplace fourier mode,
c       give by the formula
c
c           1/sqrt(4 pi^2 r r0) * Q_{-1/2}(xi)
c
c       this routine does NOT rely on q2lege01 because it separately
c       computes the elliptic integrals, possible in terms of 
c       sqrt(1-xi**2) which can be computed accurately knowing
c       r, r0, z, z0
c
c       input:
c
c         r,z - cylindrical coordinates of the target
c         r0,z0 - cylindrical coordinates of the source
c         ifder - switch determining whether to calculate derivatives,
c             this should be set to 1 for derivative calculation
c
c       output:
c
c         val - the value of the 0th order laplace fourier mode
c         dr,dz - partial derivatives d/dr and d/dz
c         dr0,dz0 - partial derivatives d/dr and d/dz0
c         
c
        done=1
        two=2
        pi=4*atan(done)
c
        zdiff=z-z0
        xi=(r*r+r0*r0+zdiff**2)/2/r/r0
c
        x=sqrt(2/(1+xi))
        xp=((r-r0)**2+(z-z0)**2)/2/r/r0/(1+xi)
        xp=sqrt(xp)
        call ellipke(x,xp,valk,vale,ifder,dk,de)
c
        q0=x*valk
        val=1/sqrt(r*r0)*q0/2/pi
        if (ifder .ne. 1) return

c
c       and the partial derivatives . . . we are assuming that
c       valk and vale are computed as accurately as possible from ellipke
c
        fac1=1/two/pi/sqrt(r*r0)
        dxidr=(r**2-r0**2-zdiff**2)/2/r**2/r0
        dxidr0=(r0**2-r**2-zdiff**2)/2/r/r0**2
        dxidz=2*zdiff/2/r/r0
c
        q1=xi*x*valk-sqrt(2*(xi+1))*vale
        ximinus=(r**2-r0**2)**2+2*r**2*zdiff**2+2*r0**2*zdiff**2
     1      +zdiff**4
        ximinus=ximinus/4/r**2/r0**2
        der0=(q1-xi*q0)/2/ximinus
c
        dr=fac1*(-1/r/2*q0+der0*dxidr)
        dr0=fac1*(-1/r0/2*q0+der0*dxidr0)
        dz=fac1*der0*dxidz
        dz0=-dz
c
        return
        end
c
c
c
c
c
        subroutine torlap1(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
        implicit real *8 (a-h,o-z)
c
c       this routine evaluates the 1st order laplace fourier mode
c       give by the formula
c
c           1/sqrt(4 pi^2 r r0) * Q_{1/2}(xi)
c
c       this routine does NOT rely on q2lege01 because it separately
c       computes the elliptic integrals, possibly in terms of 
c       sqrt(1-xi**2) which can be computed accurately knowing
c       r, r0, z, z0
c
c       input:
c
c         r,z - cylindrical coordinates of the target
c         r0,z0 - cylindrical coordinates of the source
c         ifder - switch determining whether to calculate dr and dz,
c             this should be set to 1 for derivative calculation
c
c       output:
c
c         val - the value of the 1st order laplace fourier mode
c         dr,dz - partial derivatives d/dr, d/dz
c         dr0,dz0 - partial derivatives d/dr0, and d/dz0
c         
c
        done=1
        two=2
        pi=4*atan(done)
c
        zdiff=z-z0
        xi=(r*r+r0*r0+zdiff**2)/2/r/r0
c
        x=sqrt(2/(1+xi))
        xp=((r-r0)**2+(z-z0)**2)/2/r/r0/(1+xi)
        xp=sqrt(xp)
        call ellipke(x,xp,valk,vale,ifder,dk,de)
c
        q1=xi*x*valk-sqrt(2*(xi+1))*vale
        val=1/sqrt(r*r0)*q1/2/pi
        if (ifder .ne. 1) return

c
c       and the partial derivatives . . . we are assuming that
c       valk and vale are computed as accurately as possible from ellipke
c
        fac1=1/two/pi/sqrt(r*r0)
        dxidr=(r**2-r0**2-zdiff**2)/2/r**2/r0
        dxidr0=(r0**2-r**2-zdiff**2)/2/r/r0**2
        dxidz=2*zdiff/2/r/r0
c
        q0=x*valk
        ximinus=(r**2-r0**2)**2+2*r**2*zdiff**2+2*r0**2*zdiff**2
     1      +zdiff**4
        ximinus=ximinus/4/r**2/r0**2
        der1=(-q0+xi*q1)/2/ximinus
c
        dr=fac1*(-1/r/2*q1+der1*dxidr)
        dr0=fac1*(-1/r0/2*q1+der1*dxidr0)
        dz=fac1*der1*dxidz
        dz0=-dz
c
        return
        end     
c
c
c
c
c
        subroutine ellipke(rk,rkp,valk,vale,ifder,dk,de)
        implicit real *8 (a-h,o-z)
        real *8 coefsk(13),coefse(13),pi
        data pi/3.1415926535897931d0/
        data coefsk /0.25000000000000000d0,
     2      0.28125000000000000d0,
     4      0.29296875000000000d0,
     6      0.29907226562500000d0,
     8      0.30281066894531250d0,
     *      0.30533409118652344d0,
     2      0.30715155601501465d0,
     4      0.30852276831865311d0,
     6      0.30959402793087065d0,
     8      0.31045401134178974d0,
     *      0.31115958864029380d0,
     2      0.31174890604302163d0,
     4      0.31224850364885981d0/
c
        data coefse/0.25000000000000000d0,
     2      9.37500000000000000d-2,
     4      5.85937500000000000d-2,
     6      4.27246093750000000d-2,
     8      3.36456298828125000d-2,
     *      2.77576446533203125d-2,
     2      2.36270427703857422d-2,
     4      2.05681845545768738d-2,
     6      1.82114134076982737d-2,
     8      1.63396848074626178d-2,
     *      1.48171232685854193d-2,
     2      1.35543002627400710d-2,
     4      1.24899401459543924d-2/
c
c       this routine returns, to full double precision, the
c       values of K, E, dK/dk and dE/dk - it asks that the user
c       provide both k AND k' = sqrt(1-k**2) to full double
c       precision
c
c       input:
c
c         rk - the point at which to evaluate the elliptic integrals
c         rkp - this should be equal to sqrt(1-rk**2), the routine
c             WILL bomb if it is not
c
c       output:
c
c         valk, vale - values of K and E
c         dk, de - derivatives of K and E
c
c
        done=1
        call ke_eva(rk,rkp,valk,vale)
        if (ifder .ne. 1) return

c
c       . . . and then compute the derivatives if desired
c
        thresh=.2d0
        if (rk .gt. thresh) then
            dk=vale/rkp**2/rk-valk/rk
            de=(vale-valk)/rk
            return
        endif

c
c       use a series expansion
c
        dk=0
        de=0
        do 1400 i=1,13
        i2=2*i-1
        dk=dk+coefsk(i)*rk**i2
        de=de+coefse(i)*rk**i2
 1400 continue
c
        dk=dk*pi
        de=-de*pi
c
        return
        end
c
c
c
c
c
        subroutine ke_eva(rk,rkp,valk,vale)
        implicit real *8 (a-h,o-z)
        real *8 cpe(6)
        data cpe /0.5000000000000000d+00,
     2            0.1083333333333333d+01,
     3            0.1200000000000000d+01,
     4            0.1251190476190476D+01,
     5            0.1280158730158730D+01,
     6            0.1298845598845599D+01/
c
c       computes the complete elliptic integrals E and K
c       to full double precision for any value of rk in [0,1)
c
c       K=int_0^pi/2 {1/(1-rk^2 sin(t)^2)}
c
c       E=int_0^pi/2 {1-rk^2 sin(t)^2}
c
c       The code uses Chebyshev interpolation for values of rk
c       in the interval [0 , 1023/1024] and a power series expansion
c       for values of rk in the interval (1023/1024 , 1)
c
c       input:
c         rk - the elliptic integral parameter above
c         rkp - the user is instructed to supply this number CORRECTLY
c           if rk is anywhere near 1, rkp = sqrt(1-rk**2)
c
c       output:
c
c         valk - the value of K(rk)
c         vale - the value of E(rk)
c 
c
        done=1.0d0
        b=1023*done/1024
c

        call elliptic_ke(rk,valk,vale)
        if (rk .lt. b) return

c
c       if here, then use rkp which the user should have supplied
c       CORRECTLY as the value sqrt(1-rk**2)
c

        rlog=log(4/rkp)

        coef2=done/2/2*(rlog-1)
        coef4=(done*3/8)**2*(rlog-1-done/6)
        coef6=(done*15/48)**2*(rlog-1-done/6-done/15)
        coef8=(done*105/384)**2*(rlog-1-done/6-done/15-done/28)
        coef10=(done*945/3840)**2*(rlog-1-done/6-done/15-done/28
     1      -done/45)
        coef12=(done*10395/46080)**2*(rlog-1-done/6-done/15-done/28
     1      -done/45-done/66)
        coef14=(done*135135/645120)**2*(rlog-1-done/6-done/15
     1      -done/28-done/45-done/66-done/91)
c
        valk=rlog+coef2*rkp**2
     1      +coef4*rkp**4
     2      +coef6*rkp**6
     3      +coef8*rkp**8
     3      +coef10*rkp**10
     3      +coef12*rkp**12
     3      +coef14*rkp**14

c
c       now calculate E
c

        xp=rkp
        xplog=-log(xp/4)
c
        fe=1
        fe=fe+(xplog-cpe(1))/2*xp**2
        fe=fe+3*(xplog-cpe(2))/16*xp**4
c
        fe=fe+45*(xplog-cpe(3))/384*xp**6
        fe=fe+1575*(xplog-cpe(4))/18432*xp**8
c
        fe=fe+99225*(xplog-cpe(5))/1474560*xp**10
        fe=fe+9823275*(xplog-cpe(6))/176947200*xp**12
c
        vale=fe
c
        return
        end
c
c
c
c
c
        subroutine p2lege01(x,p0,p1)
        implicit real *8 (a-h,o-z)
c
c       this subroutine returns the value of the first two
c       half-order legendre functions of the 
c       first kind, P_-1/2 and P_1/2. this routine requires that
c       the evaluation point x be larger than 1, i.e. in the toroidal
c       regime of the associated functions
c
c       input:
c
c         x - evaluation point, must be larger than 1
c
c       output:
c
c         p0 - value of P_{-1/2}(x)
c         p1 - value of P_{1/2}(x)
c
c
        done=1
        pi=4*atan(done)
c
        if (abs(x) .lt. 1) then
          call prin2('x is less than 1!!! x=*',x,1)
          stop
        endif
c
        x0=sqrt((x-1)/(x+1))
        call elliptic_ke(x0,fk,fe)
        p0=2/pi*sqrt(2/(x+1))*fk
c
        x0=sqrt(x*x-1)
        x1=sqrt(2*x0/(x+x0))
        call elliptic_ke(x1,fk,fe)
        p1=2/pi*sqrt(x+x0)*fe
c
        return
        end
c
c
c
c
c
        subroutine p2lege01der(x,p0,p1)
        implicit real *8 (a-h,o-z)
c
c       this subroutine returns the derivative of the first two
c       half-order legendre functions of the 
c       first kind, P'_-1/2 and P'_1/2. this routine requires that
c       the evaluation point x be larger than 1, i.e. in the toroidal
c       regime of the associated functions
c
c       input:
c
c         x - evaluation point, must be larger than 1
c
c       output:
c
c         p0 - value of P'_{-1/2}(x)
c         p1 - value of P'_{1/2}(x)
c
c
        done=1
        pi=4*atan(done)
c
        if (abs(x) .lt. 1) then
          call prin2('x is less than 1!!! x=*',x,1)
          stop
        endif
c
        xi=sqrt((x-1)/(x+1))
        ifder=1
        xip=sqrt(1-xi*xi)
        call ellipke(xi,xip,fk,fe,ifder,dk,de)
c
        p0=-2/pi*(sqrt((x+1)/2)/(x+1)**2*fk
     1      -sqrt(2/(x+1))*dk/(x+1)**1.5d0/sqrt(x-1))
c
        x0=sqrt(x*x-1)
        x1=sqrt(2*x0/(x+x0))
        x1p=sqrt(1-x1*x1)
        call ellipke(x1,x1p,fk,fe,ifder,dk,de)
c
        p1=1/pi*(1+x/x0)/sqrt(x+x0)*fe+
     1      2/pi*sqrt(x+x0)*de*sqrt(done/2)*(x/(x0*(x+x0))-
     2      x0*(1+x/x0)/(x+x0)**2)/sqrt(x0/(x+x0))
c
        return
        end
c
c
c
c
c
        subroutine pq2legeders(x,nmax,pvals,qvals,pders,qders)
        implicit real *8 (a-h,o-z)
        real *8 pvals(0:1),qvals(0:1),pders(0:1),qders(0:1)
c
c       evaluate the values of P_(n-1/2) and Q_(n-1/2) 
c       of the first and second kinds along with their derivates - 
c       the first nmax+3 elements are used in qvals and pvals because
c       of the need for derivatives
c
        done=1
        h2=done/2
c
c       . . . first run the recurrence up for the p's
c
        nmax2=nmax+2
        call p2legeders(x,nmax2,pvals,pders)

c
c       use ratio of the q's and the wronskian to calculate the high
c       degree q's
c
        call q2legeratio(x,nmax2,ratio,err)
        call prin2('inside pq2legeders, error in ratio=*',err,1)
c
        d=(nmax2-h2)*(pvals(nmax2)-ratio*pvals(nmax2-1))
        qvals(nmax2-1)=1/d
        qvals(nmax2)=qvals(nmax2-1)*ratio
c
c       and run the downward recurrence for the q's
c
        i=nmax2
        qders(i)=(2*i-done)/2/(x*x-1)*(x*qvals(i)-qvals(i-1))
c
        i=nmax2-1
        qders(i)=(2*i-done)/2/(x*x-1)*(x*qvals(i)-qvals(i-1))
c
        do 1600 i=nmax2-1,1,-1
        coef0=2*i/(i-h2)
        coef1=(i+h2)/(i-h2)
        qvals(i-1)=coef0*x*qvals(i)-coef1*qvals(i+1)
        qders(i-1)=coef0*(qvals(i)+x*qders(i))-coef1*qders(i+1)
 1600 continue
c
        return
        end
c
c
c
c
c
        subroutine pq2leges(x,nmax,pvals,qvals)
        implicit real *8 (a-h,o-z)
        real *8 pvals(0:1),qvals(0:1)
c
c       evaluate the first nmax+1 half-degree legendre functions
c       of the first and second kinds
c
        done=1
        h2=done/2
c
c       . . . first run the recurrence up for the p's
c
        call p2leges(x,nmax,pvals)

c
c       use ratio of the q's and the wronskian to calculate the high
c       degree q's
c
        call q2legeratio(x,nmax,ratio,err)
        call prin2('inside pq2leges, error in ratio=*',err,1)
c
        d=(nmax-h2)*(pvals(nmax)-ratio*pvals(nmax-1))
        qvals(nmax-1)=1/d
        qvals(nmax)=qvals(nmax-1)*ratio

c
c       and run the downward recurrence for the q's
c
        do 1600 i=nmax-1,1,-1
        qvals(i-1)=2*i*x*qvals(i)/(i-h2)-(i+h2)*qvals(i+1)/(i-h2)
 1600 continue

c
        return
        end
c
c
c
c
c
        subroutine p2leges(x,nmax,vals)
        implicit real *8 (a-h,o-z)
        real *8 vals(0:1)
c
c       this routine returns the first nmax+1 (!!!) half-degree
c       legendre functions of the first kind (and zeroth order).
c       the P's are the maximal solution to the recurrence, so
c       they WILL blow up!
c
c       NOTE: the vals array starts at 0 in this routine....and upon
c       output, the first nmax+1 values are occupied
c
c       input:
c
c         x - the point at which to evaluate the polynomials, we
c           are assuming that x >= 1 !!!
c         nmax - the highest degree to evaluate
c
c       output:
c
c         vals - the values P_{-1/2}(x), ..., P_{nmax-1/2}(x)
c
c
        done=1
        pi=4*atan(done)
c
c       return if nmax is less than 0
c
        if (nmax .lt. 0) return
c
        call p2lege01(x,p0,p1)
c
        if (nmax .ge. 0) vals(0)=p0
        if (nmax .ge. 1) vals(1)=p1
        if (nmax .lt. 2) return

c
c       . . . run the recurrence now
c
        half=done/2
c
        do 1800 i=1,nmax-1
        vals(i+1)=2*i*x*vals(i)/(i+half)-(i-half)*vals(i-1)/(i+half)
 1800 continue
c
        return
        end
c
c
c
c
c
        subroutine p2legeders(x,nmax,vals,ders)
        implicit real *8 (a-h,o-z)
        real *8 vals(0:1),ders(0:1)
c
c       this routine returns the value and derivative of the first 
c       nmax+1 (!!!) half-degree
c       legendre functions of the first kind (and zeroth order).
c       the P's are the maximal solution to the recurrence, so
c       they WILL blow up!
c
c       NOTE: the vals array starts at 0 in this routine....and upon
c       output, the first nmax+1 values are occupied
c
c       input:
c
c         x - the point at which to evaluate the polynomials, we
c           are assuming that x >= 1 !!!
c         nmax - the highest degree to evaluate
c
c       output:
c
c         vals - the values P_{-1/2}(x), ..., P_{nmax-1/2}(x)
c         ders - the values P'_{-1/2}(x), ..., P'_{nmax-1/2}(x)
c
c
        done=1
        pi=4*atan(done)
c
c       return if nmax is less than 0
c
        if (nmax .lt. 0) return
c
        call p2lege01(x,p0,p1)
        call p2lege01der(x,der0,der1)
c
        if (nmax .ge. 0) then
            vals(0)=p0
            ders(0)=der0
        endif
c
        if (nmax .ge. 1) then
            vals(1)=p1
            ders(1)=der1
        endif
c
        if (nmax .lt. 2) return

c
c       . . . run the recurrence now
c
        half=done/2
c
        do 1800 i=1,nmax-1
        coef0=2*i/(i+half)
        coef1=(i-half)/(i+half)
        vals(i+1)=coef0*x*vals(i)-coef1*vals(i-1)
        ders(i+1)=coef0*(vals(i)+x*ders(i))-coef1*ders(i-1)
 1800 continue
c
        return
        end
c
c
c
c
c
        subroutine q2legeratio(x,n,val,err)
        implicit real *8 (a-h,o-z)
c
c       compute the ratio Q_{n-1/2}(x)/Q_{n-1-1/2}(x) using a
c       continued fraction expansion. will use modified lentz's
c       algorithm, see numerical recipes.
c
c       NOTE: this routine is using in converting an upward p_{n-1/2}
c         recurrence into a downward q_{n-1/2} recurrence. usually loss
c         in precision in q_{n-1/2} is caused by loss of precision in
c         this routine - which is why 'err' is returned, which is the
c         relative error obtained in this routine. swap out this routine
c         if you are not happy.
c
c       input:
c
c         x - the point at which to compute the ratio, larger than 1
c         n - the order of the ratio to compute, Q_{n-1/2}/Q_{n-3/2}
c
c       output:
c
c         val - the value of the ratio
c         err - the relative error obtained in the ratio
c
c
        done=1
        half=done/2
        v=n-half
        err=1
c
        tiny=1.0d-30
        f0=tiny
        c0=f0
        d0=0
c
        ifdone=0
        thresh=1.0d-13
c
        maxiter=100000
        do 2000 i=1,maxiter
c
        if (i .eq. 1) a1=1
        if (i .ne. 1) a1=-(1+1/(v+i-2))
        b1=(2+1/(v+i-1))*x
c
        c1=b1+a1/c0
        if (c1 .eq. 0) c1=tiny
c
        d1=b1+a1*d0
        if (d1 .eq. 0) d1=tiny
        d1=1/d1

        t1=c1*d1
        f1=f0*t1
        err=abs(f1-f0)/abs(f1)
        if (err .lt. thresh) ifdone=ifdone+1
        if (ifdone .eq. 10) goto 2100
cccc        call prin2('err=*',err,1)
c
        a0=a1
        b0=b1
        c0=c1
        d0=d1
        f0=f1
c
 2000 continue
c
        call prinf('cont frac did not converge!!*',a0,0)
        call prin2('requested precision=*',thresh,1)
        call prin2('achieved relative precision=*',err,1)
        call prinf('maxiter=*',maxiter,1)
c
 2100 continue
c
        val=f1
c
        return
        end
c
c
c
c
c
        subroutine q2lege01(x,val0,val1)
        implicit real *8 (a-h,o-z)
c
c       this subroutine returns the value of the first two
c       half-order legendre functions of the 
c       second kind, q_-1/2 and q_1/2. this routine requires that
c       the evaluation point x be larger than 1, i.e. in the toroidal
c       regime of the associated functions
c
c       input:
c
c         x - evaluation point, must be larger than 1
c
c       output:
c
c         val0 - value of Q_{-1/2}(x)
c         val1 - value of Q_{1/2}(x)
c
c
        done=1
        if (abs(x) .lt. 1) then
          call prin2('x is less than 1!!! x=*',x,1)
          stop
        endif
c
        xarg=2/(1+x)
        xarg2=sqrt(xarg)
        call elliptic_ke(xarg2,valk,vale)
c
        val0=sqrt(xarg)*valk
        val1=x*sqrt(xarg)*valk-sqrt(2*(x+1))*vale        
c
        return
        end
c
c
c
c
c
        subroutine q2lege01der(x,val0,val1,der0,der1)
        implicit real *8 (a-h,o-z)
c
c       this subroutine returns the value and derivative 
c       of the first two half-order legendre functions of the 
c       second kind, Q_-1/2 and Q_1/2. this routine requires that
c       the evaluation point x be larger than 1, i.e. in the toroidal
c       regime of the associated functions
c
c       input:
c
c         x - evaluation point, must be larger than 1
c
c       output:
c
c         val0 - value of Q_{-1/2}(x)
c         val1 - value of Q_{1/2}(x)
c         der0 - value of Q'_{-1/2}(x)
c         der1 - value of Q'_{1/2}(x)
c
c
        done=1
        if (abs(x) .lt. 1) then
          call prin2('x is less than 1!!! x=*',x,1)
          stop
        endif
c
        xarg=2/(1+x)
        xarg2=sqrt(xarg)
        xarg2p=sqrt(2-xarg2**2)
        ifder=1
        call ellipke(xarg2,xarg2p,valk,vale,ifder,fkd,fed)
c
        val0=sqrt(xarg)*valk
        val1=x*sqrt(xarg)*valk-sqrt(2*(x+1))*vale
c
        der0=(val1-x*val0)/2/(x*x-1)
        der1=(-val0+x*val1)/2/(x*x-1)
c
        return
        end
c
c
c
c
c
        subroutine elliptic_ke(rk,fk,fe)
        implicit real *8 (a-h,o-z)
        real *8 ck1(21),ce1(19),ck2(21),ce2(19),
     1      ck3(21),ce3(18),ck4(21),ce4(18),ck5(21),ce5(18),
     2      ck6(21),ce6(17),ck7(21),ce7(17),ck8(21),ce8(17),
     3      ck9(21),ce9(16),ck10(21),ce10(16),cpk(6),cpe(6)
c
c       computes the complete elliptic integrals E and K
c       to (almost) full double precision for any value of rk in [0,1)
c
c       K=int_0^pi/2 {1/(1-rk^2 sin(t)^2)}
c
c       E=int_0^pi/2 {1-rk^2 sin(t)^2}
c
c       The code uses Chebyshev interpolation for values of rk
c       in the interval [0 , 1023/1024] and a power series expansion
c       for values of rk in the interval (1023/1024 , 1)
c
c
        data ck1 /0.1612037992010708D+01,
     2            0.5625155780525147D-01,
     3            0.1601139981954603D-01,
     4            0.1199136013244665D-02,
     5            0.2196984906892129D-03,
     6            0.2572714003455883D-04,
     7            0.4157373351133303D-05,
     8            0.5790779632254359D-06,
     9            0.9091363437458903D-07,
     *            0.1363119008094121D-07,
     1            0.2141749241065508D-08,
     2            0.3326026853745003D-09,
     3            0.5269913445724438D-10,
     4            0.8343497599664959D-11,
     5            0.1334096507318333D-11,
     6            0.2138095607981829D-12,
     7            0.3446303376046879D-13,
     8            0.5570800755507768D-14,
     9            0.9039540708624828D-15,
     *            0.1470717395525351D-15,
     1            0.2399650200550059D-16/
c
        data ce1 /0.1532599271134926D+01,
     2           -0.5131294197315643D-01,
     3           -0.1341338599523514D-01,
     4           -0.3497194527665565D-03,
     5           -0.5599972340446773D-04,
     6           -0.4328696740208485D-05,
     7           -0.6078542460845124D-06,
     8           -0.6737801968566317D-07,
     9           -0.9321801192985317D-08,
     *           -0.1202983594275060D-08,
     1           -0.1695974623341116D-09,
     2           -0.2356411622366691D-10,
     3           -0.3404799129275913D-11,
     4           -0.4930604596817409D-12,
     5           -0.7284503841104897D-13,
     6           -0.1083134044853210D-13,
     7           -0.1629868279897924D-14,
     8           -0.2469004892495648D-15,
     9           -0.3770730824485522D-16/
c
        data ck2 /0.1784515365808022D+01,
     2            0.1110391126936834D+00,
     3            0.1364432924323473D-01,
     4            0.1551430978212123D-02,
     5            0.2061645208821321D-03,
     6            0.2854685434399510D-04,
     7            0.4114351694566976D-05,
     8            0.6080806922347191D-06,
     9            0.9163457174139017D-07,
     *            0.1401506372433434D-07,
     1            0.2169035193158199D-08,
     2            0.3389332582667541D-09,
     3            0.5338600216206609D-10,
     4            0.8465677951543258D-11,
     5            0.1350180228861333D-11,
     6            0.2164105515699072D-12,
     7            0.3483753773175031D-13,
     8            0.5629549294560753D-14,
     9            0.9127878988896837D-15,
     *            0.1484496168760237D-15,
     1            0.2420850830386146D-16/
c
        data ce2 /0.1398413751717964D+01,
     2           -0.7422936211166669D-01,
     3           -0.5421759273405274D-02,
     4           -0.2631240618414422D-03,
     5           -0.2453671874081236D-04,
     6           -0.2527875092124950D-05,
     7           -0.2923226226972259D-06,
     8           -0.3598949017837802D-07,
     9           -0.4648855205888318D-08,
     *           -0.6220270969051253D-09,
     1           -0.8555558050476965D-10,
     2           -0.1202963669580172D-10,
     3           -0.1722225378078289D-11,
     4           -0.2502977134142151D-12,
     5           -0.3684266869970842D-13,
     6           -0.5482576315224507D-14,
     7           -0.8236206323351144D-15,
     8           -0.1247573715880566D-15,
     9           -0.1903616725601511D-16/
c
        data ck3 /0.2034149071283014D+01,
     2            0.1356022774700329D+00,
     3            0.1387500556122127D-01,
     4            0.1617321640157409D-02,
     5            0.2107737133804371D-03,
     6            0.2911430602275378D-04,
     7            0.4179500647510380D-05,
     8            0.6163356466601956D-06,
     9            0.9271234239707845D-07,
     *            0.1416079838902866D-07,
     1            0.2189229190726057D-08,
     2            0.3417899101676672D-09,
     3            0.5379697280593070D-10,
     4            0.8525647942900036D-11,
     5            0.1359037585909114D-11,
     6            0.2177324451166015D-12,
     7            0.3503661639036443D-13,
     8            0.5659770355633076D-14,
     9            0.9174080327089753D-15,
     *            0.1491603832189214D-15,
     1            0.2431847102727987D-16/
c
        data ce3 /0.1262369746058800D+01,
     2           -0.5853793826087978D-01,
     3           -0.2565750034777385D-02,
     4           -0.1414885053628810D-03,
     5           -0.1247567895174542D-04,
     6           -0.1293285404561683D-05,
     7           -0.1485818978102447D-06,
     8           -0.1825658743108401D-07,
     9           -0.2353322622754271D-08,
     *           -0.3144218151774452D-09,
     1           -0.4319546323184765D-10,
     2           -0.6067793096753708D-11,
     3           -0.8680171663028279D-12,
     4           -0.1260696693391248D-12,
     5           -0.1854654360441891D-13,
     6           -0.2758599606590789D-14,
     7           -0.4142385935821100D-15,
     8           -0.6272362991886799D-16/
c
        data ck4 /0.2322517890327813D+01,
     2            0.1497771936934633D+00,
     3            0.1420871176099745D-01,
     4            0.1648767109859717D-02,
     5            0.2135681326487684D-03,
     6            0.2941497621320073D-04,
     7            0.4214584912456625D-05,
     8            0.6206930235118721D-06,
     9            0.9327813627971957D-07,
     *            0.1423679797280886D-07,
     1            0.2199712624674734D-08,
     2            0.3432672664426349D-09,
     3            0.5400885908900624D-10,
     4            0.8556487427945484D-11,
     5            0.1363582560636546D-11,
     6            0.2184094828648336D-12,
     7            0.3513841423667947D-13,
     8            0.5675201890871851D-14,
     9            0.9197642323772305D-15,
     *            0.1495224619151855D-15,
     1            0.2437443275449716D-16/
c
        data ce4 /0.1162632765594922D+01,
     2           -0.3968514093321079D-01,
     3           -0.1274454526186430D-02,
     4           -0.7292538649443840D-04,
     5           -0.6334370425994085D-05,
     6           -0.6545703045020360D-06,
     7           -0.7500687896271960D-07,
     8           -0.9201011581786124D-08,
     9           -0.1184631017571241D-08,
     *           -0.1581363895290874D-09,
     1           -0.2171010889659148D-10,
     2           -0.3048035842538741D-11,
     3           -0.4358398443605881D-12,
     4           -0.6327777959519529D-13,
     5           -0.9306161708398539D-14,
     6           -0.1383830169918344D-14,
     7           -0.2077527271603695D-15,
     8           -0.3145154390456146D-16/
c
        data ck5 /0.2633917629481799D+01,
     2            0.1584759400179597D+00,
     3            0.1443514026431293D-01,
     4            0.1665522234518389D-02,
     5            0.2150652679062538D-03,
     6            0.2957241223839394D-04,
     7            0.4232790212329976D-05,
     8            0.6229390127335453D-06,
     9            0.9356843033280709D-07,
     *            0.1427565803194010D-07,
     1            0.2205058922051074D-08,
     2            0.3440191144859930D-09,
     3            0.5411650912985906D-10,
     4            0.8572133744391795D-11,
     5            0.1365885735836113D-11,
     6            0.2187522304292906D-12,
     7            0.3518990460567593D-13,
     8            0.5683001463921423D-14,
     9            0.9209543438643490D-15,
     *            0.1497052404492962D-15,
     1            0.2440266773335938D-16/
c
        data ce5 /0.1097142142624543D+01,
     2           -0.2505893163011859D-01,
     3           -0.6420581919102552D-03,
     4           -0.3705473933475876D-04,
     5           -0.3196441594943055D-05,
     6           -0.3294742130443460D-06,
     7           -0.3769756457765553D-07,
     8           -0.4619887829895383D-08,
     9           -0.5944169309149156D-09,
     *           -0.7931031492752457D-10,
     1           -0.1088428217542038D-10,
     2           -0.1527678668420280D-11,
     3           -0.2183923324934515D-12,
     4           -0.3170137764181948D-13,
     5           -0.4661519777952799D-14,
     6           -0.6930750948306479D-15,
     7           -0.1040382692332133D-15,
     8           -0.1574867191709941D-16/
c
        data ck6 /0.2959489686284203D+01,
     2            0.1638289698634456D+00,
     3            0.1456756980554160D-01,
     4            0.1674337912268060D-02,
     5            0.2158413820446733D-03,
     6            0.2965315429421077D-04,
     7            0.4242073368329189D-05,
     8            0.6240800612970612D-06,
     9            0.9371553720586108D-07,
     *            0.1429531405902026D-07,
     1            0.2207759368837520D-08,
     2            0.3443984591499076D-09,
     3            0.5417077573536188D-10,
     4            0.8580015330053540D-11,
     5            0.1367045215889798D-11,
     6            0.2189246893203585D-12,
     7            0.3521580126730418D-13,
     8            0.5686922678509032D-14,
     9            0.9215524664009104D-15,
     *            0.1497970731966793D-15,
     1            0.2441684994727412D-16/
c
        data ce6 /0.1056546065213113D+01,
     2           -0.1516138877268048D-01,
     3           -0.3241597221726031D-03,
     4           -0.1869271621900106D-04,
     5           -0.1606307414704859D-05,
     6           -0.1653184718313071D-06,
     7           -0.1889952044223614D-07,
     8           -0.2314960279285094D-08,
     9           -0.2977489070258700D-09,
     *           -0.3971715105904223D-10,
     1           -0.5449600454313243D-11,
     2           -0.7647713012424263D-12,
     3           -0.1093161997436108D-12,
     4           -0.1586654118552050D-13,
     5           -0.2332897956714794D-14,
     6           -0.3468311253221947D-15,
     7           -0.5206007217477787D-16/
c
        data ck7 /0.3293731730309479D+01,
     2            0.1670721829402436D+00,
     3            0.1464012727025866D-01,
     4            0.1678885869220422D-02,
     5            0.2162369764309903D-03,
     6            0.2969406684212659D-04,
     7            0.4246762308022841D-05,
     8            0.6246552729437665D-06,
     9            0.9378959596426587D-07,
     *            0.1430520006287684D-07,
     1            0.2209116572094867D-08,
     2            0.3445890037466461D-09,
     3            0.5419802140294816D-10,
     4            0.8583970964673986D-11,
     5            0.1367626957926976D-11,
     6            0.2190111936895563D-12,
     7            0.3522878793774243D-13,
     8            0.5688888704650069D-14,
     9            0.9218523026321839D-15,
     *            0.1498431015614956D-15,
     1            0.2442395738940912D-16/
c
        data ce7 /0.1032286146823760D+01,
     2           -0.8908618289978243D-02,
     3           -0.1633949907294598D-03,
     4           -0.9391613202465229D-05,
     5           -0.8052923511978042D-06,
     6           -0.8280952483488412D-07,
     7           -0.9462748159600979D-08,
     8           -0.1158757161840068D-08,
     9           -0.1490115581163899D-09,
     *           -0.1987426547028807D-10,
     1           -0.2726684878582413D-11,
     2           -0.3826208117400699D-12,
     3           -0.5468835688059231D-13,
     4           -0.7937263258206990D-14,
     5           -0.1166987123272577D-14,
     6           -0.1734894139379159D-15,
     7           -0.2604032928871931D-16/
c
        data ck8 /0.3633184545027836D+01,
     2            0.1689983893947868D+00,
     3            0.1467845650912843D-01,
     4            0.1681200288243994D-02,
     5            0.2164367692202118D-03,
     6            0.2971466355202746D-04,
     7            0.4249118927149408D-05,
     8            0.6249440739772014D-06,
     9            0.9382675370703423D-07,
     *            0.1431015775771155D-07,
     1            0.2209796938894353D-08,
     2            0.3446844963241054D-09,
     3            0.5421167258159797D-10,
     4            0.8585952522604433D-11,
     5            0.1367918333149771D-11,
     6            0.2190545150711659D-12,
     7            0.3523529091787843D-13,
     8            0.5689873080309319D-14,
     9            0.9220024155400447D-15,
     *            0.1498661438841441D-15,
     1            0.2442751521382831D-16/
c
        data ce8 /0.1018158757102843D+01,
     2           -0.5122914154214850D-02,
     3           -0.8217011671331295D-04,
     4           -0.4707869559416428D-05,
     5           -0.4031977474265149D-06,
     6           -0.4144300872755771D-07,
     7           -0.4734661073448387D-08,
     8           -0.5797008146600572D-09,
     9           -0.7454031195163281D-10,
     *           -0.9941079743185240D-11,
     1           -0.1363816128661815D-11,
     2           -0.1913694720565795D-12,
     3           -0.2735177424781746D-13,
     4           -0.3969633542645326D-14,
     5           -0.5836285569310349D-15,
     6           -0.8676322752589878D-16,
     7           -0.1302274537145517D-16/
c
        data ck9 /0.3975707886238759D+01,
     2            0.1701205730484671D+00,
     3            0.1469826397805721D-01,
     4            0.1682368536798688D-02,
     5            0.2165371812613249D-03,
     6            0.2972499770330323D-04,
     7            0.4250300314821449D-05,
     8            0.6250887763340238D-06,
     9            0.9384536493122885D-07,
     *            0.1431264030337606D-07,
     1            0.2210137566164818D-08,
     2            0.3447322979641566D-09,
     3            0.5421850528933919D-10,
     4            0.8586944240548844D-11,
     5            0.1368064147280727D-11,
     6            0.2190761931201500D-12,
     7            0.3523854482675549D-13,
     8            0.5690365609780252D-14,
     9            0.9220775208246445D-15,
     *            0.1498776720984849D-15,
     1            0.2442929515433969D-16/
c
        data ce9 /0.1010090220020678D+01,
     2           -0.2897362069082847D-02,
     3           -0.4124113031059618D-04,
     4           -0.2357080798395684D-05,
     5           -0.2017391026655395D-06,
     6           -0.2073115441101683D-07,
     7           -0.2368157432138485D-08,
     8           -0.2899313455529224D-09,
     9           -0.3727882178379426D-10,
     *           -0.4971529726722949D-11,
     1           -0.6820268057370082D-12,
     2           -0.9569953736769597D-13,
     3           -0.1367779005067876D-13,
     4           -0.1985067721322444D-14,
     5           -0.2918480850270836D-15,
     6           -0.4338625115994152D-16/
c
        data ck10 /0.4320007671151427D+01,
     2            0.1707633032997668D+00,
     3            0.1470836334087705D-01,
     4            0.1682955573557723D-02,
     5            0.2165875185445303D-03,
     6            0.2973017381825033D-04,
     7            0.4250891783435439D-05,
     8            0.6251612033599909D-06,
     9            0.9385467866478125D-07,
     *            0.1431388250393704D-07,
     1            0.2210307991093370D-08,
     2            0.3447562126575236D-09,
     3            0.5422192342687092D-10,
     4            0.8587440334747216D-11,
     5            0.1368137086035028D-11,
     6            0.2190870364915949D-12,
     7            0.3524017238684874D-13,
     8            0.5690611960050915D-14,
     9            0.9221150856912632D-15,
     *            0.1498834379711634D-15,
     1            0.2443018538198001D-16/
c
        data ce10 /0.1005551466563614D+01,
     2           -0.1617149642459013D-02,
     3           -0.2066939246131722D-04,
     4           -0.1179348158538705D-05,
     5           -0.1009049202576981D-06,
     6           -0.1036800107568071D-07,
     7           -0.1184286093029554D-08,
     8           -0.1449859552941835D-09,
     9           -0.1864158144856706D-10,
     *           -0.2486012717444621D-11,
     1           -0.3410431283095187D-12,
     2           -0.4785347338709464D-13,
     3           -0.6839371255821170D-14,
     4           -0.9925966573161736D-15,
     5           -0.1459325013935453D-15,
     6           -0.2169428583899399D-16/
c
        data cpk /0.2772588722239781d+01,
     2           -0.3862943611198906d+00,
     3            0.5393397569993164d-01,
     4            0.8079550343381188d-02,
     5           -0.1621012950014027d-01,
     6            0.1347996438028915d-01/
c
        data cpe /0.5000000000000000d+00,
     2            0.1083333333333333d+01,
     3            0.1200000000000000d+01,
     4            0.1251190476190476D+01,
     5            0.1280158730158730D+01,
     6            0.1298845598845599D+01/
c
c
c       if rk is in [0 , 1/2]...
c
        done=1
        a=0
        b=done/2
        if (rk .gt. b) goto 2000
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=19
        call chebseval(tk,fk,ck1,mk-1)
        call chebseval(tk,fe,ce1,me-1)
        return
c
 2000   continue
c
c       if rk is in [1/2 , 3/4]...
c
        a=done/2
        b=3*done/4
        if (rk .gt. b) goto 2400
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=19
        call chebseval(tk,fk,ck2,mk-1)
        call chebseval(tk,fe,ce2,me-1)
        return
c
 2400   continue
c
c       if rk is in [3/4 , 7/8]...
c
        a=3*done/4
        b=7*done/8
        if (rk .gt. b) goto 2800
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=18
        call chebseval(tk,fk,ck3,mk-1)
        call chebseval(tk,fe,ce3,me-1)
        return
c
 2800   continue
c
c       if rk is in [7/8 , 15/16]...
c
        a=7*done/8
        b=15*done/16
        if (rk .gt. b) goto 3200
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=18
        call chebseval(tk,fk,ck4,mk-1)
        call chebseval(tk,fe,ce4,me-1)
        return
c
 3200   continue
c
c       if rk is in [15/16 , 31/32]...
c
        a=15*done/16
        b=31*done/32
        if (rk .gt. b) goto 3600
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=18
        call chebseval(tk,fk,ck5,mk-1)
        call chebseval(tk,fe,ce5,me-1)
        return
c
 3600   continue
c
c       if rk is in [31/32 , 63/64]...
c
        a=31*done/32
        b=63*done/64
        if (rk .gt. b) goto 4000
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=17
        call chebseval(tk,fk,ck6,mk-1)
        call chebseval(tk,fe,ce6,me-1)
        return
c
 4000   continue
c
c       if rk is in [63/64 , 127/128]...
c
        a=63*done/64
        b=127*done/128
        if (rk .gt. b) goto 4200
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=17
        call chebseval(tk,fk,ck7,mk-1)
        call chebseval(tk,fe,ce7,me-1)
        return
c
 4200   continue
c
c       if rk is in [127/128 , 255/256]...
c
        a=127*done/128
        b=255*done/256
        if (rk .gt. b) goto 4600
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=17
        call chebseval(tk,fk,ck8,mk-1)
        call chebseval(tk,fe,ce8,me-1)
        return
c
 4600   continue
c
c       if rk is in [255/256 , 511/512]...
c
        a=255*done/256
        b=511*done/512
        if (rk .gt. b) goto 4800
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=16
        call chebseval(tk,fk,ck9,mk-1)
        call chebseval(tk,fe,ce9,me-1)
        return
c
 4800   continue
c
c       if rk is in [511/512 , 1023/1024]...
c
        a=511*done/512
        b=1023*done/1024
        if (rk .gt. b) goto 5200
c
        tk=2*(rk-a)/(b-a)-1
c
        mk=21
        me=16
        call chebseval(tk,fk,ck10,mk-1)
        call chebseval(tk,fe,ce10,me-1)
        return
c
 5200   continue
c
c       rk is in [1023/1024 , 1)  - use taylor series about rk=1 for fk
c
        x=rk-1
c
c       first series for fk...
c
        dd1=(1-x/2+5*x**2/16-7*x**3/32+169*x**4/1024-269*x**5/2048)
c
        dd2=-dlog(done*2)-dlog(dabs(x))-x/2+x**2/8-x**3/24+x**4/64
        dd2=dd2-x**5/160
c
        dd3=cpk(1)+cpk(2)*x+cpk(3)*x**2+cpk(4)*x**3+cpk(5)*x**4
        dd3=dd3+cpk(6)*x**5
c
        fk=(dd1*dd2+dd3)/2
c
c       now series for fe...
c
        xp=dsqrt(1.0d0-rk**2)
        xplog=-dlog(xp/4)
c
        fe=1
        fe=fe+(xplog-cpe(1))/2*xp**2
        fe=fe+3*(xplog-cpe(2))/16*xp**4
c
        fe=fe+45*(xplog-cpe(3))/384*xp**6
        fe=fe+1575*(xplog-cpe(4))/18432*xp**8
c
        fe=fe+99225*(xplog-cpe(5))/1474560*xp**10
        fe=fe+9823275*(xplog-cpe(6))/176947200*xp**12
c
        return
        end
c
c
c
c
c
      subroutine chebseval(x,val,texp,n)
      implicit real *8 (a-h,o-z)
      dimension texp(1)
c 
c     this subroutine computes the value of a chebychev
c     expansion with coefficients texp at point x in interval [-1,1]
c 
c     input:
c 
c       x - evaluation point
c       texp - expansion coefficients
c       n - order of expansion
c
c       NOTE: n is the order of the expansion, which is
c         one less than the number of terms in the expansion!!
c
c     output:
c 
c       val - computed value
c
c
        done=1
        pjm2=1
        pjm1=x
c 
        val=texp(1)*pjm2+texp(2)*pjm1
c 
        do j = 2,n
           pj= 2*x*pjm1-pjm2
           val=val+texp(j+1)*pj
           pjm2=pjm1
           pjm1=pj
       enddo
c 
        return
        end
