        !
        ! Copyright (C) April 6, 2015:
        !
        !               Mike O'Neil
        !               oneil@cims.nyu.eu
        !

      program q2legetest
        implicit double precision (a-h,o-z)
        real *8 :: vals(0:1000000)

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

        call q2lege01(x, val0, val1)
        write(6,*) 'Q_-1/2 = ', val0
        write(6,*) 'Q_1/2 =  ', val1
        write(13,*) 'Q_-1/2 = ', val0
        write(13,*) 'Q_1/2 =  ', val1

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
        ! call the normalized routine to compute the Fourier modes of
        ! 1/r
        !
        call torlaps(ier, nterms, src, targ, vals, lvals, ntop)
        call prinf('after torlaps, ier = *', ier, 1)


        
        
        stop
      end program q2legetest
      
