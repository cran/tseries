      subroutine muin2ser (information, ndata, freq, x, y,
     +                     intx, inty, s, q, q_unsort,
     +                     indices_x, indices_y, position, 
     +                     maxbit, confidence, trace)


c------------------------------------------------------------------------
c
c     Compute the mutual information I(S,Q) from two time series (x,y).
c
c     Theory and algorithm are described in:
c
c     A.M. Fraser & H.L. Swinney, "Independent Coordinates for Strange
c     Coordinates from Mutual Information", Phys.Rev. 33A, 1134 (1986).
c
c     Original version of the program by Dr.Katharina Krischer,
c     <krischer@prince.princeton.edu>, <krischer@anika.rz-berlin.mpg.de>
c
c     Adaptation and complete revision by Th.-M. Kruel,
c     <kruel@phys-chemie.uni-wuerzburg.dbp.de>
c
c     Version 1.21     June 17, 1992
c
c
c     Changes for R version, A. Trapletti, 26.3.99
c
c
c     Important variables and constants:
c
c     x:          the original time series
c     npt:        the number of data in the time series
c     intx, inty: the time series in integer representation
c                 rescaled to the interval [0 <= intx < 2^maxbit]
c     maxbit:     maximum resolution of the time series in integer
c                 representation in bit.
c     confidence: confidence level for chi^2-test upon uniform
c                 distribution
c
c
c
c     calling scheme:
c
c     main ->   rescale
c          \->  quicksort
c           \-> muin ->    zaehle
c                    \->   unterteile ->  begin
c                     \               \-> zaehle
c                      \-> leastsquare
c
c------------------------------------------------------------------------



*      implicit none
      parameter   (maxut=25)
      integer     intx(1), inty(1), intmax,
     +            s(1), q_unsort(1), q(1),
     +            indices_x(1), indices_y(1), position(1),
     +            npt, ndata, ichi,
     +            maxbit, maxprec, nabbruch
      real*8      x(1), y(1), freq, information
      real*8      confidence, alpha(15), sqchi3(15), sqchi15(15),
     +            chi3alpha, chi15alpha, cdiff
      logical*4   trace

c--- the following data are taken from Tab.1.1.2.10, p.21,
c    Bronstein & Semendjajew, Taschenbuch der Mathematik,
c    25.Aufl., Teubner, Leipzig, 1991
      data   alpha    /0.99,   0.98,   0.95,   0.90,   0.80,
     +                 0.70,   0.60,   0.50,   0.40,   0.30,
     +                 0.20,   0.10,   0.05,   0.02,   0.01/
      data   sqchi3   /0.0383, 0.0617, 0.117,  0.195,  0.335,
     +                 0.475,  0.629,  0.789,  0.985,  1.223,
     +                 1.547,  2.10,   2.60,   3.27,   3.77/
      data   sqchi15  /0.347,  0.400,  0.487,  0.567,  0.687,
     +                 0.780,  0.870,  0.953,  1.049,  1.153,
     +                 1.287,  1.487,  1.667,  1.887,  2.040/
      common /bounds/ npt, maxprec
      common /chi/    chi3alpha, chi15alpha


c---  get the limits for the chi^2-Test given the confidence level
c     (also: significance niveau) 'alpha':
      if (((confidence .gt. 1.0) .or. (confidence .lt. 0.0))
     +      .and. trace) then
         write (0,*) 'WARNING: Confidence level illegal; '//
     +               'set to default value of 20%'
         ichi = 11
         goto 15
      endif

      cdiff = abs(confidence-alpha(1))
      ichi = 1
      do 10 i=2,15
         if (abs(confidence-alpha(i)) .le. cdiff) then
            cdiff = abs(confidence-alpha(i))
            ichi = i
         endif
 10   continue
      if ((confidence .ne. alpha(ichi)) .and. trace) then
         write (0,'(a,i3,a)') 'NOTE: Confidence level has been '//
     +          'rounded to ', nint(100*alpha(ichi)), '%'
      endif

 15   chi3alpha  = sqchi3(ichi)
      chi15alpha = sqchi15(ichi)

     
      npt = ndata

c---  rescale the data:
      maxprec = maxbit
      call rescale (x, y, intx, inty, intmax, trace)

c---  sort the time series 'intx':
      call quicksort (intx, intmax, position, s, indices_x)

c---  arrange the variable 'y' in that way, in which 'x' is sorted
c     in ascending order:
      call quicksort (inty, intmax, position, q_unsort, indices_y)
      do 30 i=1,npt
         q(i) = q_unsort(indices_x(i))
 30   continue


c---  write the head of output:
      if (trace) then
         write (*,'(i5,a)')   maxprec,     '     ! max.bit'
         write (*,'(f5.2,a)') alpha(ichi), '     ! confidence level'
      endif

c---  compute the mutual information:
      call muin (q, information, nabbruch)
      if (information .lt. 0.0) then
         information = 0.0
      endif
      if (trace) then
         write (*,'(a,g14.7,a)') 'information = ', information, ' bit'
      endif
      if (trace .and. (nabbruch .gt. 0)) then
         write (0,'(a,i2,a)')
     +         'WARNING: There exists substructure '//
     +         'beyond the maximum resolution of ',
     +         maxbit, ' bit.'
         write (0,'(a,g14.7,a,i6,a)')
     +      '         information =', information,
     +      '; substructure in', nabbruch, ' elements.'
      endif

      return
      end                                           

c------------------------------------------------------------------------
c------------------------------------------------------------------------

      subroutine rescale (x, y, intx, inty, intmax, trace)

c------------------------------------------------------------------------
c     Rescales the time series in 'x' to an integer representation 'xint'
c     in the interval [0,xintmax[
c------------------------------------------------------------------------


*      implicit none
      real*8      x(1), xmax, xmin, xdiff, xdiffmin, xextent,
     +            y(1), ymax, ymin, ydiff, ydiffmin, yextent,
     +            xfactor, x_precision, intx_precision,
     +            yfactor, y_precision, inty_precision,
     +            roundadd
      integer     intx(1), inty(1), intxmax, intymax, intmax,
     +            npt, maxprec, len
      character   runform*10
      logical     trace

      common /bounds/ npt, maxprec
      data   roundadd /1.e-6/


      xmax = x(1)
      xmin = x(1)
      ymax = y(1)
      ymin = y(1)
      xdiffmin = 1.e+38
      ydiffmin = 1.e+38
      do 100 i=2,npt
         xmin = min (x(i),xmin)
         xmax = max (x(i),xmax)
         ymin = min (y(i),ymin)
         ymax = max (y(i),ymax)
         xdiff = abs (x(i)-x(i-1))
         ydiff = abs (y(i)-y(i-1))
         if (xdiff .ne. 0.0) xdiffmin = min (xdiffmin,xdiff)
         if (ydiff .ne. 0.0) ydiffmin = min (ydiffmin,ydiff)
 100  continue
      xextent = xmax - xmin
      yextent = ymax - ymin

c---  compute the precision of the data [in bit]:
      x_precision = log (xextent/xdiffmin) / log(2.0)
      y_precision = log (yextent/ydiffmin) / log(2.0)
      if (trace) then
         write (0,'(a,f6.3,a)') 'The resolution of the 1st time '//
     +                          'series is ', x_precision, ' bit.'
         write (0,'(a,f6.3,a)') 'The resolution of the 2nd time '//
     +                          'series is ', y_precision, ' bit.'
      endif

c---  limit the maximum possible resolution:
      intx_precision = min (dble(maxprec), x_precision)
      inty_precision = min (dble(maxprec), y_precision)
      if (trace .and. (intx_precision .lt. x_precision)) then
         write (0,'(a,i2,a)')
     +      'The precision of the 1st time series has been '//
     +      'limited to ', maxprec, ' bit.'
      endif
      if (trace .and. (inty_precision .lt. y_precision)) then
         write (0,'(a,i2,a)')
     +      'The precision of the 2nd time series has been '//
     +      'limited to ', maxprec, ' bit.'
      endif

c---  compute the maximum necessary integer value for the rescaled
c     variables:
      intxmax = int ((2**intx_precision) + roundadd) + 1
      intymax = int ((2**inty_precision) + roundadd) + 1

c---  rescaling factor:
      xfactor = dble(intxmax-1)/xextent
      yfactor = dble(intymax-1)/yextent

      do 200 i=1,npt
         intx(i) = int ((x(i) - xmin) * xfactor  +  roundadd) + 1
         inty(i) = int ((y(i) - ymin) * yfactor  +  roundadd) + 1
 200  continue

      intxmax = intxmax - xmin + 1
      intymax = intymax - ymin + 1
      intmax = max(intxmax, intymax)    ! it should be intxmax = intymax
      intmax = min(2**maxprec+1,intmax)
      if (trace) then
         len = int(log(float(intmax))/log(10.0))+1
         write (runform,'(a,i1,a)') '(a,i', len, ',a)'
         write (0,runform) 'The original time series were rescaled '//
     +                     'to the interval [1,', intmax,'].'
         write (0,*)
      endif

      return
      end
                      
c------------------------------------------------------------------------
                                          
      subroutine quicksort (ix, ixmax, position, ix_bin, indices)

c------------------------------------------------------------------------
c     Sort the values of the array 'ix' in ascending order such that the 
c     respective positions of the values are stored in 'indices'.
c
c     'ix_bin' contains the unsorted array 'ix' in integer 
c     representation [1...npt].
c     'ix_bin(i)' thereby specifies the location of the i-th element of x
c     in an array containing the values of 'ix' in ascending order.
c     Thus it is a projection of 'ix' to the interval [1...npt]
c     with equal probability distribution.
c     'indices(i)' gives the location of the i-th element of the sorted
c     field within the array 'ix'.
c
c     Example:      ix      = {6,8,7,5}
c                => ix_bin  = {2,4,3,1}
c                => indices = {4,1,3,2}
c------------------------------------------------------------------------

*      implicit none
      integer     position(1), indices(1), npt,
     +            ix_bin(1), ix(1), ixmax, k, maxprec

      common /bounds/ npt, maxprec

      do 5 i=1,ixmax
         position(i) = 0
 5    continue

      do 10 i=1,npt
         indices(i) = i
         ix_bin(i) = position(ix(i))   ! if x(i) is known beforehand,
                                       !  store its position
         position(ix(i)) = i           ! store position of x(i)
 10   continue


      i = 0
      do 30 j=1,ixmax         ! for all values of 'ix+1'
         k = position(j)      ! get the stored position
 20      if (k .ne. 0) then   ! this value is known already
            i = i+1           ! increment counter
            indices(i) = k    ! store position
            k = ix_bin(k)     ! check, if this position has been
            goto 20           !  present before
         endif
 30   continue

      do 40 j=1,npt                 ! for all points
         ix_bin(indices(j)) = j     ! the position 'indices(j)' holds
 40   continue                      !  the j-th largest value of 'ix'
    

      return
      end

c------------------------------------------------------------------------

      subroutine muin (q, information, nabbruch)

c------------------------------------------------------------------------
c
c     Compute the mutual information from (s,q)
c
c     Important constants and variables:
c
c     m:          depth of partitioning
c     maxut:      maximum depth of partitioning
c     num:        storage for the tree of probabilites in the
c                 squares R_m(K_m)
c     g(m):       size of the squares R_m(K_m)
c     startx, starty:   starting coordinates of the lower left corner of 
c                       the square R_m(K_m)
c     k(m):       index of R_m(K_{m-1},k(m)) in level m;
c                 may have values between 0 and 3
c     teilinfo:   partial contribution to F[R_0(K_0)] from the actual
c                 structure
c     substruct:  flag indicating the presence of substructure below the
c                 actual level
c------------------------------------------------------------------------

*      implicit none
      parameter   (maxut=25)            ! max. number of levels
      integer     npt, num, n4,
     +            m, startx, starty,
     +            g(maxut), k(maxut),
     +            q(1), maxprec
      real*8      information, teilinfo, f, xhelp, dlog2
      logical     debug, substruct, g0
      parameter   (debug=.false.)       ! for debugging only
      dimension   n4(4), num(0:20)                                          

      common /bounds/   npt, maxprec
      common /geometry/ startx, starty, g
      common /indices/  k

      external    dlog2
                                                      

c---  partition the axes into units of 2^m:
c     (works optimal, when 'npt' is a power of 2)
      do 10 m=1,maxut                                                
         k(m) = 0 
         g(m) = npt/(2**m)
         if (g(m) .eq. 1) goto 20     ! no further intervals
 10   continue

 20   m = 1
 
c---  partition the total square into 4 parts and count the points:
      startx = 1
      starty = 1
      call zaehle (q, 1, n4)
      if (debug) write (*,*)
     +   'Results of the first partition: ',(n4(i),i=1,4)
                                                 
c---  partition the total square (level 0) into 16 parts and count the
c     points:
      call unterteile (q, 1, num)
      if (debug) write (*,30) (num(i),i=5,20)
 30   format ('Results of the second partition: ',/,(3x,i7))

c---  test the total area upon substructure:
      g0 = .true.      ! test square G0
      call leastsquare (m, g0, num, f, substruct)
      if (.not. substruct) then    ! no substructure at all
         information = 0.
         return
      endif
      
      g0 = .false.
      teilinfo = f
      
      if (debug) write (*,*) 'G0 has substructure. ',
     +                       'teilinfo = ',teilinfo

c---  systematically test refined partitions upon the presence of
c     substructure.
c     start with the lower left square R_1(0)
      m = 1
      k(m) = 0
      nabbruch = 0
      
 50   call unterteile (q, m, num)   ! partition the square R_m(K_m) into
                                    ! 16 parts and count the points
      call leastsquare (m, g0, num, f, substruct)   ! test the structure
      if (substruct) then
         if ((m+1) .gt. (maxprec-3)) then
           xhelp = dble(num(0))
           f = xhelp * dlog2(xhelp)
           nabbruch = nabbruch+1
         else                                                               
           teilinfo = teilinfo + f  ! add up the actual contribution
           m = m + 1                ! step one level deeper..
           goto 50                  ! ..and continue partitioning
         endif
      endif

  60  teilinfo = teilinfo + f
      if (debug) then 
         write (*,*) 'No further substructure for ',
     +               'm=', m, ' k(m)=', k(m)
         write (*,*) 'teilinfo = ', teilinfo
      endif
      
c---  There has been no further substructure for this square,
c     continue with the next square of the level:
  70  k(m) = k(m) + 1
      if (k(m) .le. 3) goto 50   ! not yet finished with this level

      k(m) = 0                   ! reset the index of the finished
                                 ! sublevel
      m = m - 1                  ! step one level higher
      if (m .gt. 0) goto 70      ! 0-level not yet reached => continue
      
c---  end of partitions.

      information = (1./dble(npt)) * teilinfo - dlog2(dble(npt))
      
      return
      end


c------------------------------------------------------------------------

      subroutine zaehle (q, m ,n)

c------------------------------------------------------------------------
c     Partition a square of size 2*g(m) with its lower left corner
c     located at 'startx, starty' (the lower left corner of the total 
c     square has the coordinates (1,1)) into 4 subsquares and count the
c     number of points within these subsquares
c
c     Scheme of squares:
c
c
c                  2*g(m) -> +---------------+
c                            |       |       |
c                            |   2   |   3   |
c                            |       |       |
c                    g(m) -> |-------+-------|
c                            |       |       |
c                            |   0   |   1   |
c                            |       |       |
c                  starty -> +---------------+
c                            ^       ^       ^
c                       startx       g(m)    2*g(m)
c
c
c------------------------------------------------------------------------
                                                            
*      implicit none
      parameter   (maxut=25)            ! max. number of levels
      integer     q(1), m, npt, 
     +            startx, starty, g(maxut), n(4),
     +            k(maxut), maxprec

      common /bounds/   npt, maxprec
      common /geometry/ startx, starty, g
      common /indices/  k


      
      do 5 i=1,4
         n(i) = 0       ! preset the sums to zero
 5    continue 

c---  search the lower left (0) upper left (2) square:
      do 10 j = startx, startx+g(m)-1
         if ((q(j) .ge. starty) .and. (q(j) .lt. starty+2*g(m))) then
c---        the point if within one of the squares (0,2)
            if (q(j) .lt. (starty+g(m))) then
               n(1) = n(1) + 1            ! square (0)
            else
               n(3) = n(3) + 1            ! square (2)
            endif
         endif
 10   continue

c---  search the lower right (1) upper right (3) square:
      do 20 j = startx+g(m), startx + 2*g(m) -1
         if ((q(j) .ge. starty) .and. (q(j) .lt. starty+2*g(m))) then
c---        the point if within one of the squares (1,3)
            if (q(j) .lt. (starty+g(m))) then
               n(2) = n(2) + 1            ! square (1)
            else
               n(4) = n(4) + 1            ! square (3)
            endif
         endif
 20   continue

      return 
      end

c------------------------------------------------------------------------

      subroutine leastsquare (m, g0, num, f, substruct)
   
c------------------------------------------------------------------------
c     Test the square R_m(K_m) upon the presence of substructure
c     according to equations (20) u.(21) in
c     A.M.Fraser and H.L.Swinney, Phys.Rev. 33A(2), 1134 (1986).
c------------------------------------------------------------------------

*      implicit none
      real*8      a, b, lb4
      parameter   (a=16./9., b=256./225., ! for chi^2-test
     +             lb4=2.0)               ! log_2 (4.0)
      integer     m, num(0:20), npt, ipunkt, num1, maxprec
      real*8      f, chi3, chi15, num1vt, num1inv, num1sz, dlog2,
     +            chi3alpha, chi15alpha
      logical     g0, substruct

      common /bounds/ npt, maxprec
      common /chi/    chi3alpha, chi15alpha

      external dlog2
      
    
      if (g0) then      ! special case G0 (total square)
         num1 = npt     ! all points are in the total square (trivial)
      else
         num1 = num(0)         ! the number of points within this level
         if (num1 .eq. 0) then      ! no points at all
            substruct = .false.
            f = 0.
            return
         else
            m = m+1                 ! look at a level deeper
         end if 
      endif
      
      num1vt = dble(num1)/4.
      num1inv = 1./dble(num1)
      num1sz = dble(num1)/16.

      ipunkt = 1
      chi3 = 0.
      do 10 i=0,3
         chi3 = chi3 + (dble(num(ipunkt+i)) - num1vt)**2   !  (Gl.21)
 10   continue
      chi3 = a * num1inv * chi3
      
      if (chi3 .ge. chi3alpha) then ! test according to eqn.21 failed
         substruct = .true.         ! => there is a substructure
         m = m-1                    ! restore the level
         f = dble(num1) * lb4       ! F(R_m) = N(R_m) * log(4)       (Gl.20b)
         return
      endif

      m = m+1                 ! look at two levels deeper
      ipunkt = 5
      chi15 = 0
      do 20 i=0,15
         chi15 = chi15 + (dble(num(ipunkt+i)) - num1sz)**2  ! (Gl.22)
 20   continue
      chi15 = b * num1inv * chi15
      
      if (chi15 .ge. chi15alpha) then ! test according to eqn.22 failed
         substruct = .true.           ! => there is a substructure
         f = dble(num1) * lb4         ! F(R_m) = N(R_m) * log(4)       (Gl.20b)
      else 
         substruct = .false.
         f = dble(num1) * dlog2(dble(num1)) ! F(R_m)=N(R_m)*log[N(R_m)] (Gl.20a)
      endif
      
      m = m-2                 ! restore the level
      
      return
      end

c------------------------------------------------------------------------

      subroutine begin (m)

c------------------------------------------------------------------------
c     Compute the coordinates of the lower left corner of the square
c     R_m(k1,k2,..,km)
c
c     scheme:
c
c                +-------+
c          km =  | 2 | 3 |
c                |---+---|
c                | 0 | 1 |
c                +-------+
c
c------------------------------------------------------------------------

*      implicit none
      parameter   (maxut=25)            ! max. number of levels
      integer     m, startx, starty, g(maxut),
     +            k(maxut)

      common /geometry/ startx, starty, g
      common /indices/  k


      startx = 1
      starty = 1

      do 10 i=1,m-1
         startx = mod(k(i),2) * g(i) + startx
c---        addend positive for k(i)=1 or k(i)=3 (right squares),
c           else zero (left squares)
         starty = (k(i)/2) * g(i) + starty         
c---        addend positive for k(i)=2 or k(i)=3 (upper squares),
c           else zero (lower squares)
 10   continue

      return
      end

c------------------------------------------------------------------------

      subroutine unterteile (q, m, num)
                                                            
c------------------------------------------------------------------------
c     Partition the square R_m(K_m) into 16 subsquares and count the
c     points therein
c------------------------------------------------------------------------

*      implicit none
      parameter   (maxut=25)            ! max. number of levels
      integer     ksave, k(maxut), ipunkt,
     +            m, q(1), n4(4), num(0:20)

      common /indices/  k

    
      ksave = k(m+1)      ! save the actual index in level m+1

c---  count the points within this level:
      call begin (m)
      call zaehle (q, m, n4)
c---  store the total number of points of the actual square in num(0):
       num(0) = n4(k(m)+1)

c---  count the points one level deeper:
      call begin (m+1)           ! set the starting values of
                                 !  R_{m+1}(K_{m+1}) within the
                                 !  sublevel m+1
      call zaehle (q, m+1, n4)
      do 5 i=1,4
         num(i) = n4(i)
 5    continue
      
c---  count the points two levels deeper:
      do 10 i=0,3         ! run through all 4 subsquares
         k(m+1) = i       !  of the level m+1
         call begin (m+2)           ! set the starting values of
                                    !  R_{m+1}(K_{m+1},k(m+2)) within
                                    !  the sublevel m+2
         call zaehle (q, m+2, n4)   ! partition the square of level m+1
                                    !  into 4 further subsquares of 
                                    !  level m+2 and count the points
         ipunkt=4*(i+1)
         do 15 j=1,4
            num(ipunkt+j) = n4(j)
 15      continue
 10   continue
 
      k(m+1) = ksave      ! restore the index of level m+1

      return
      end
                                 
c------------------------------------------------------------------------

      real*8 function dlog2 (arg)
      real*8  arg

      dlog2 = dlog(arg) / dlog(2.0d0)

      return
      end
