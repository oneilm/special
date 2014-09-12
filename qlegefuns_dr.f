        implicit real *8 (a-h,o-z)
        dimension qfuns(100000),qfuns3(100000)
        complex *16 z, zfuns(100000), ima
c 

        call prini(6,13)
        ima = (0,1)
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
        eps = 1.0d-13
        z = 1.5d0 + ima*eps
        call prin2('z = *', z, 2)
        call zqneval(z, n, zfuns)
c
        call prinf2('from zqneval:*', x, 0)
        do i = 1,n+1
          write(6,*) 'i = ', i-1, 'val = ', zfuns(i)
          write(13,*) 'i = ', i-1, 'val = ', zfuns(i)
        enddo


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
