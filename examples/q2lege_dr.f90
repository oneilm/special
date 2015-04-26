        !
        ! Copyright (C) April 6, 2015:
        !
        !               Mike O'Neil
        !               oneil@cims.nyu.eu
        !

      program q2legetest
        implicit double precision (a-h,o-z)
        real *8 :: vals(0:1000000), src(10), targ(10)
        real *8 :: vals2(0:1000000), errs(0:1000000)
        real *8 :: ders(0:1000000), der2s(0:1000000)
        complex *16 :: ima

        call prini(6,13)

        done=1
        pi=4*atan(done)
        ima=(0,1)

        print *, 'enter nterms:'
        read *, nterms
        call prinf('nterms = *', nterms, 1)

        x=1.5d0
        call prin2('x=*',x,1)
        write(6,*) 'x = ', x

        !
        ! calculate Q_-1/2 and Q_1/1
        !

        xminus = x-1
        call q2lege01(x, xminus, val0, val1)
        write(6,*) 'Q_-1/2 = ', val0
        write(6,*) 'Q_1/2 =  ', val1
        write(13,*) 'Q_-1/2 = ', val0
        write(13,*) 'Q_1/2 =  ', val1

        !
        ! test the derivative routine
        !
        call q2lege01_all(x, xminus, w0, w1, der0, der1, &
            der20, der21)

        print *
        call prin2('from q2lege_all, val0 = *', w0, 1)
        call prin2('from q2lege_all, val1 = *', w1, 1)
        call prin2('err in val0 = *', val0-w0, 1)
        call prin2('err in val1 = *', val1-w1, 1)

        h = .0001d0
        y = x+h
        yminus = y-1
        call q2lege01(y, yminus, b0, b1)
        y = x-h
        yminus = y-1
        call q2lege01(y, yminus, a0, a1)

        u0 = (b0-a0)/2/h
        u1 = (b1-a1)/2/h

        print *
        call prin2('from finite difference, der0 = *', u0, 1)
        call prin2('from finite difference, der1 = *', u1, 1)
        call prin2('from q2lege_all, der0 = *', der0, 1)
        call prin2('from q2lege_all, der1 = *', der1, 1)
        call prin2('err in der0 = *', u0-der0, 1)
        call prin2('err in der1 = *', u1-der1, 1)

        y = x+h
        yminus = y-1
        call q2lege01_all(y, yminus, w0, w1, a0, a1, &
            z0, z1)

        y = x-h
        yminus = y-1
        call q2lege01_all(y, yminus, w0, w1, b0, b1, &
            z0, z1)

        u0 = (b0-a0)/2/h
        u1 = (b1-a1)/2/h
        print *
        call prin2('from finite difference, der20 = *', u0, 1)
        call prin2('from finite difference, der21 = *', u1, 1)
        call prin2('from q2lege_all, der20 = *', der20, 1)
        call prin2('from q2lege_all, der21 = *', der21, 1)
        call prin2('err in der20 = *', u0-der20, 1)
        call prin2('err in der21 = *', u1-der21, 1)





        
        stop

        
        
        !
        ! calculate the first nterms
        !
        xminus = x - done
        lvals = 1000000
        call q2leges(ier, nterms, x, xminus, vals, lvals, ntop)
        call prinf('after q2leges, ier = *', ier, 1)
        call prinf('after q2leges, ntop = *', ntop, 1)
        call prin2('q-halves are = *', vals, nterms+1)

        do i = 0,nterms
          write(6,*), 'i = ', i, 'q_i-1/2 = ', vals(i)
          write(13,*), 'i = ', i, 'q_i-1/2 = ', vals(i)
        enddo


        !
        ! calculate the first nterms along with first and second derivatives
        !
        call q2leges_all(ier, nterms, x, xminus, vals, ders, &
          der2s, lvals, ntop)
        call prinf('after q2leges_all, ier = *', ier, 1)
        call prinf('after q2leges_all, ntop = *', ntop, 1)
        call prin2('q-halves are = *', vals, nterms+1)
        call prin2('q-halves ders = *', ders, nterms+1)
        call prin2('q-halves der2s = *', der2s, nterms+1)

        stop

        
        !
        ! call the normalized routine to compute the Fourier modes of
        ! 1/r
        !
        src(1) = 1.0d0
        src(2) = 1.0d0

        targ(1) = 1.5d0
        targ(2) = 1.5d0
        
        call torlaps(ier, nterms, src, targ, vals, lvals, ntop)
        call prinf('after torlaps, ier = *', ier, 1)
        call prinf('after torlaps, ntop = *', ntop, 1)
        call prin2('torlaps are = *', vals, nterms+1)

        !
        ! do the direct integration now
        !
        n = 1000
        h = 2*pi/n
        call prin2('h = *', h, 1)
        
        do i = 0,nterms
          vals2(i) = 0
        
          do j = 1,n
            theta = (j-1)*h
            r = src(1)**2 + targ(1)**2 - 2*src(1)*targ(1)*cos(theta) &
                + (src(2) - targ(2))**2
            r = sqrt(r)
            vals2(i) = vals2(i) + h*exp(-ima*i*theta)/(4*pi*r)
          enddo

          vals2(i) = vals2(i)/2/pi
        enddo

        call prin2('from direct integration, vals2 = *', vals2, &
            nterms+1)

        do i = 0,nterms
          errs(i) = vals2(i) - vals(i)
          vals2(i) = vals2(i)/vals(i)
        enddo        

        call prin2('errs = *', errs, nterms+1)
        call prin2('ratios = *', vals2, nterms+1)
        
        stop

        
      end program q2legetest
      
