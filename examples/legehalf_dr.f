c
c       Copyright (C) March 27, 2013:
c
c               Mike O'Neil
c               oneil@cims.nyu.eu
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        implicit real *8 (a-h,o-z)
        dimension pols(0:10000)
        complex *16 ima,cdz,cval,cval2,zk,cval3
        real *8 tarray(2),result,vals(100000),
     1      pvals(100000),qvals(100000),pders(100000),qders(100000),
     2      drs(100000),dr0s(100000),dzs(100000),dz0s(100000)
c
c
        call prini(6,13)
c
        done=1
        pi=4*atan(done)
        ima=(0,1)
c
        print *, 'enter nmax:'
        read *, nmax
        call prinf('nmax=*',nmax,1)
c
        x=1.5d0
        call prin2('x=*',x,1)
c
        call q2lege01(x,q0,q1)
        call prin2('q0=*',q0,1)
        call prin2('q1=*',q1,1)
c
        call p2leges(x,nmax,vals)
        call prin2('vals=*',vals,nmax+1)
c
        call q2legeratio(x,nmax,ratio,err)
        call prin2('and q ratio=*',ratio,1)
        call prin2('error in q ratio=*',err,1)
c
        call pq2leges(x,nmax,pvals,qvals)        
        call prin2('pvals=*',pvals,nmax+1)
        call prin2('qvals=*',qvals,nmax+1)
c
        print *
        print *
        call prin2('error in q0=*',q0-qvals(1),1)
        call prin2('error in q1=*',q1-qvals(2),1)
c
c       and the derivative routines
c
        print *
        print *
c
        xi=.6d0
        h=1.0d-5
        call elliptic_ke(xi-h,fk1,fe1)
        call elliptic_ke(xi+h,fk2,fe2)
        fkder=(fk2-fk1)/2/h
        feder=(fe2-fe1)/2/h
c
        call prin2('from finite diff, der of fk is *',fkder,1)
        call prin2('from finite diff, der of fe is *',feder,1)
c
        ifder=1
        xip=sqrt(1-xi**2)
        call ellipke(xi,xip,d0,d1,ifder,fkd,fed)
        call prin2('from formula, der of fk is *',fkd,1)
        call prin2('error is *',fkder-fkd,1)
        call prin2('from formula, der of fe is *',fed,1)
        call prin2('error is *',feder-fed,1)
c
        print *
        print *
        call p2lege01der(x,p0der,p1der)
        call prin2('via formula, p0der=*',p0der,1)
        call prin2('via formula, p1der=*',p1der,1)

        h=1.0d-6
        call p2lege01(x-h,p01,p11)
        call p2lege01(x+h,p02,p12)
        p0d=(p02-p01)/2/h
        p1d=(p12-p11)/2/h
c
        call prin2('and via finite diff, p0d=*',p0d,1)
        call prin2('error is*',p0der-p0d,1)
        call prin2('and via finite diff, p1d=*',p1d,1)
        call prin2('error is*',p1der-p1d,1)
c
        call pq2legeders(x,nmax,pvals,qvals,pders,qders)
        call prin2('pders=*',pders,nmax+1)
        call prin2('qders=*',qders,nmax+1)
c
        call q2lege01der(x,val0,val1,der0,der1)
        call prin2('diff in der of q0=*',qders(1)-der0,1)
        call prin2('diff in der of q1=*',qders(2)-der1,1)

c
c       calculate the 0th and 1st order laplace fourier modes
c
        r=1
        r0=.5d0
        z=1
        z0=-.1d0
        xi=(r*r+r0*r0+(z-z0)**2)/2/r/r0
c
        ifder=1
        call torlap0(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
c
        npts=1000
        h=2*pi/npts
        th0=0
        th1=0
c
        do 2000 i=1,npts
c
           theta=h*(i-1)
           x=r*cos(theta)
           y=r*sin(theta)
c
           theta0=0
           x0=r0*cos(theta0)
           y0=r0*sin(theta0)
c
           rrr=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
           th0=th0+1/rrr/4/pi*h
           th1=th1+1/rrr/4/pi*h*cos(theta)
 2000 continue
c
        print *
        print *
        print *
        call prin2('tor zero from integrating a ring=*',th0,1)
        call prin2('tor zero from subroutine=*',val,1)
        call prin2('ratio=*',th0/val,1)
        call prin2('error=*',th0-val,1)
c
        stop
c
        h=1.0d-6
        call torlap0(r-h,z,r0,z0,val1,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r+h,z,r0,z0,val2,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r,z-h,r0,z0,val3,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r,z+h,r0,z0,val4,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r,z,r0-h,z0,val5,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r,z,r0+h,z0,val6,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r,z,r0,z0-h,val7,ifder,dr1,dz1,dr01,dz01)
        call torlap0(r,z,r0,z0+h,val8,ifder,dr1,dz1,dr01,dz01)
c
        print *
        print *
        print *

        call torlap0(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)

        drest=(val2-val1)/2/h
        dzest=(val4-val3)/2/h
        dr0est=(val6-val5)/2/h
        dz0est=(val8-val7)/2/h
c
        err1=dr-drest
        err2=dz-dzest
        err3=dr0-dr0est
        err4=dz0-dz0est
c
        call prin2('error in r derivative=*',err1,1)
        call prin2('error in z derivative=*',err2,1)
        call prin2('error in r0 derivative=*',err3,1)
        call prin2('error in z0 derivative=*',err4,1)
c
c       check the 1st order laplace fourier mode
c
        print *
        print *
        call prin2('from integration, th1=*',th1,1)
c
        call torlap1(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
        call prin2('from subroutine, val=*',val,1)
        call prin2('diff=*',val-th1,1)
c
        print *
        h=1.0d-6
        call torlap1(r-h,z,r0,z0,val1,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r+h,z,r0,z0,val2,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r,z-h,r0,z0,val3,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r,z+h,r0,z0,val4,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r,z,r0-h,z0,val5,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r,z,r0+h,z0,val6,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r,z,r0,z0-h,val7,ifder,dr1,dz1,dr01,dz01)
        call torlap1(r,z,r0,z0+h,val8,ifder,dr1,dz1,dr01,dz01)
c
        drest=(val2-val1)/2/h
        dzest=(val4-val3)/2/h
        dr0est=(val6-val5)/2/h
        dz0est=(val8-val7)/2/h
c
        err1=dr-drest
        err2=dz-dzest
        err3=dr0-dr0est
        err4=dz0-dz0est
c
        call prin2('error in r derivative=*',err1,1)
        call prin2('error in z derivative=*',err2,1)
        call prin2('error in r0 derivative=*',err3,1)
        call prin2('error in z0 derivative=*',err4,1)

c
c       test the big long recursion...
c
        call torlaps(r,z,r0,z0,nmax,vals,ifder,drs,dzs,
     1      dr0s,dz0s)
c
        print *
        print *
        call prin2('torlaps=*',vals,nmax+1)
c
        call torlap0(r,z,r0,z0,val,ifder,dr,dz,dr0,dz0)
        call torlap1(r,z,r0,z0,val1,ifder,dr1,dz1,dr01,dz01)
c
        call prin2('relative error in tl0=*',(val-vals(1))/val,1)
        call prin2('relative error in tl1=*',(val1-vals(2))/val1,1)
c
        call prin2('rel err dr for order 0=*',(drs(1)-dr)/dr,1)
        call prin2('rel err dz for order 0=*',(dzs(1)-dz)/dz,1)
        call prin2('rel err dr0 for order 0=*',(dr0s(1)-dr0)/dr0,1)
        call prin2('rel err dz0 for order 0=*',(dz0s(1)-dz0)/dz0,1)
c
        call prin2('rel err dr for order 1=*',(drs(2)-dr1)/dr1,1)
        call prin2('rel err dz for order 1=*',(dzs(2)-dz1)/dz1,1)
        call prin2('rel err dr0 for order 1=*',(dr0s(2)-dr01)/dr01,1)
        call prin2('rel err dz0 for order 1=*',(dz0s(2)-dz01)/dz01,1)
c

        stop
        end
