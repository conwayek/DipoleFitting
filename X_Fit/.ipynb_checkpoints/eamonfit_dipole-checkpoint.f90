
  module fit_module


    implicit none
    !
    integer, parameter    :: verbose  = 3 ! Verbosity level    
    !
    public fitting
    !
    private
    EXTERNAL DGELSS
    !
    integer, parameter :: enermax=6000    ! maximal number of calculated energies for a given J and a given symmetry 
    !integer, parameter :: obsmax=5000     ! maximal number of experimental energies
    integer, parameter :: pointsmax=7000  ! maximal number of potential data points 
    !
    integer, parameter :: f_inp = 5, f_out = 6   ! read/write units 
    !
    double precision,save :: pi                  ! the pi number 
    !
    integer,parameter          ::  findif = 3    ! Dinite difference differentiation with 2 or 3 points
    !
    double precision,parameter :: stadev_best=1e-06,stab_best=1e-12  ! best standard error and srability
                                                 ! the fit will finish when these reached. 
    !
    double precision,parameter :: fitfactordeltax=0.1   ! parameter for the finite differncies differention
    !
    character(len=70) :: deriv_type = 'hellman  '  ! hellman or direct - how we calculate derivatives. 
                                                 ! direct means the finite diferentiation of energies wrt parameters.
    !
    integer,parameter :: nofititer=6             ! the Jacobi matrix can be evaluated only each nofititer+1 step. For all 
                                           ! iterations in between the the same jacobi matrix is used
    integer,parameter :: f_en =206               ! All computed energies (not only those that have obs. counterparts) will be written into here.
    integer,parameter :: f_pot=207               ! The current potential parameters will be written here. 
    integer,parameter :: f_potpoints =208        ! Computed at the current iteration step potential energy points are printed here. 
    integer,parameter :: f_res=209               ! it is where we write the fitting results 
    integer,parameter :: rms_plot=99
    !

    !
   
    contains
     !
     subroutine fitting 


      integer           :: parmax                  ! total number of potential parameters 


      integer           :: alloc,info, ierror      ! error state variables 

      integer           :: npts,pot_npts,nused      ! number of data points: current, all obs. energies, all pot. points, and actuall used.
      !
      integer           :: nrow,ncol               ! loop-variables
      integer           :: numpar                  ! number of varying parameters 
      integer           :: ndigits                 ! number of rounded digits in the parameters output. 
      integer           :: itmax                   ! Number of fitting iterations.
      !
      ! objects to store the pot. parameters, geometries, values of poten. function.  
      !
      double precision,allocatable :: parold(:),pot_values(:),potparam(:),energy(:) !,local()
double precision, ALLOCATABLE :: local(:,:)

      character(len=8),allocatable :: nampar(:)    ! parameter names 
      
      integer,allocatable :: ivar(:)               ! switch for the parameters fit: ivar(i) = 0 for the parameters that do not vary.

      !

      !
      ! objects to store the energies 
double precision,allocatable :: eps(:),abinitio_en(:)


      !
      ! objects to needed for the fitting procedure: jacobi matrix, derivatives, matrices to solve the 
      ! the ststem of linear equations, standard errors and so on. 
      double precision,allocatable :: rjacob(:,:),al(:,:),bl(:),dx(:),ai(:,:),sterr(:)

      double precision,allocatable :: wt_bit(:),wtall(:),wtnormal(:) ! weight factors - only for the energies and total.
      !
      !
      double precision :: wtsum,fit_factor,r1,r2,theta,ssq,rms,sum_sterr,conf_int
      double precision :: stadev_old,stability,stadev,tempx,deltax,v,potright,potleft
      double precision :: ssq1,ssq2,rms1,rms2
      logical          :: still_run
      integer          :: iener1,iener2,j,meval,i,l
      integer          :: iJ,ipar,isym,jsym,iener,irow,icolumn,ncards,icard
      !
      ! different character objects 
      !
      character(len=68)  :: fmt1,fmt2   ! to have nice parameters output   
      character(len=2)   :: fmt0,ch2    
      character(len=1)   :: ch1
      character(len=6)   :: fit_type    ! to switch between fitting methods. 
      character(len=7)   :: ch_tmp      ! parameter names 
      character(len=80)  :: char_tmp    ! to read a current line in the input file.
      character(len=80)  :: char_pnt
      character(len=200) :: char_job_   
      !
      ! parameters needed for the dgelss lapack routine: 
      !
      double precision, allocatable :: Tsing(:,:)
      integer                       :: rank0,jlistmax,neval_t(4),lwork
      double precision,allocatable  :: wspace(:)
      !
      integer                     :: itime0,itime2,itime,irate2,imax2 ! time control variables 
      !
      integer                     :: fititer   ! iteration variable
      !
      logical                     :: ifopen,isys
      !
      logical                     ::  SYSTEMQQ  ! system function for  calling extyernal programs, 
                                                ! platform dependent.  

       !
       ! START PROGRAM HERE 
       !
       if (verbose>=2) write(f_out,"(/'The fitting procedure starts...'/)")
    
       !
       open(f_res,file='res.fit',status='replace')
       !

       !
       open(f_en,file='en.fit',status='replace')
       !
       !
       pi = 4.0d0 * atan2(1.0d0,1.0d0)

       !
       read(f_inp,*) itmax
   
       !
       read(f_inp,"(a6)") fit_type
       !
    
       !
       write(f_out,"('itmax = ',i5)") itmax
       write(f_out,"('fit_type = ',a6)") trim(fit_type)
       !
       ! skip a line 
       read(f_inp,"(a80)") char_tmp
       !
       ! Number of the potential parameters 
       read(f_inp,*) parmax
       !
       ! Allocate all arrays that relate to the potential parameters 
       !
       allocate (parold(parmax),potparam(parmax),nampar(parmax),ivar(parmax),&
                 al(parmax,parmax),ai(parmax,parmax),bl(parmax),dx(parmax),sterr(parmax),&
                 Tsing(parmax,parmax),stat=alloc)
       if (alloc/=0) then
         write (f_out,"(' Error ',I4.0,' initializing potparam objects')") alloc
         stop 'potparam - alloc'
       end if
       !
       do i=1,parmax  
          !
          read(f_inp,"(a80)") char_tmp
          read(char_tmp(1:10),"(a10)") nampar(i)
          read(char_tmp(11:80),*) ivar(i),potparam(i)
      
          !
       end do

 pot_npts=0
 char_pnt(1:3)="" 
do while(char_pnt(1:3) /= "---")

pot_npts=pot_npts+1
read      (f_inp,"(a80)") char_pnt 

end do
pot_npts=pot_npts-1
npts=pot_npts

allocate (wtall(npts),wt_bit(npts),wtnormal(npts),&
local(1:3,npts),pot_values(npts),energy(npts),abinitio_en(npts),stat=alloc)

rewind(f_inp)



read(f_inp,*) itmax
read(f_inp,"(a6)") fit_type
read(f_inp,"(a80)") char_tmp
read(f_inp,*) parmax
do i=1,parmax  !parmax defined abovw
          read(f_inp,"(a80)") char_tmp
          read(char_tmp(1:10),"(a10)") nampar(i)
          read(char_tmp(11:80),*) ivar(i),potparam(i)
end do
     
       !
       do i=1,pot_npts
         !
        read(f_inp,*) local(1:3,i),pot_values(i),wtall(i),abinitio_en(i)
        wtnormal(i) = wtall(i)
         !
       enddo 
       !
       close(f_inp,status="keep")
       !
       ! write the potential parameters into the file, where  
       ! we store the current values 
       !
       open(f_pot,file='pot.fit',status='replace')
       !
       write(f_pot,*) parmax
       !
       do i=1,parmax
          write(f_pot,"(i8,e24.12,4x,a8)") ivar(i),potparam(i),nampar(i)
       end do
       !
       close(f_pot)

       !
          
       wtsum = sum(wtall(1:npts))
       wtall(1:npts) = wtall(1:npts)/wtsum
       
       ! 

       !
       nused=0
 !      wtsum=0.0d0
       !
       wt_bit = 0 
       !
       do i=1,pot_npts
         if (wtall(i) .gt. 0.0d0) then 
           nused=nused+1
           !
           wt_bit(i) = 1.0d0
           !
         endif
       enddo
       !
       write(f_out,"('Number of data points used in the fit: ',i9)") nused
       !
       ! Allocate objects, that will be used for the fitting procedure:

       !
       allocate (rjacob(npts,parmax),eps(npts),stat=alloc)
       if (alloc/=0) then
           write (f_out,"(' Error ',I4.0,' initializing fitting objects')") alloc
           stop 'derj0 - alloc'
       end if
       !
       !
       ! The last object to allocate - the lapack related work array
       !
       lwork = 50*parmax
       !
       allocate (wspace(lwork),stat=alloc)
       if (alloc/=0) then 
        write(f_out,"('wspace - out of memory')")  
        stop 'wspace - out of memory'
       endif
       !
       ! The fitting loop is about to start. 
       !
       ! fititer will count the iterations. 
       !
       fititer = 0
       !
       ! The initial parameters to be remembered. 
       !
       parold=potparam

  stadev_old = 1e10
  stability = 1e10
  stadev    = 1e10


       !
       still_run = .true.
       outer_loop: do while (still_run)  
          !

          !
          numpar  = 0
          !
          ! Count the actually varying number of parameters:
          !
          numpar  = 0
          do i=1,parmax
            if (ivar(i) .gt. 0) numpar=numpar+1
          enddo 
          !
          rjacob = 0 

          ! The loop starts here. 
          !
          
          do while( fititer.le.itmax)
            !
            fititer = fititer + 1
            !
            ! If fit_factor is set zero - there would be no fit to the experimental energies. 
            !
            !
            write(f_out,"(/'Iteration = ',i8)") fititer

            !
            ! Here the potential energy section starts. 
            !
            do nrow=1,pot_npts
              !
     
              r1   =local(1,nrow)
              r2   =local(2,nrow)
              theta=local(3,nrow)*pi/1.8d+02

        !
              ! Call the potential function. It has to be included into the "poten" subroutine
              !
              call poten(potparam,v,r1,r2,theta)


              !
              eps(nrow) = pot_values(nrow)-v


     

              !
              if (itmax.ge.1) then !.and.(mod(fititer-1,nofititer+1).eq.0)) then
                !
                ncol=0
                numpar = 0
                do  i=1,parmax
                  if (ivar(i) .ne. 0) then
                     !
                     ncol=ncol+1
                     tempx=potparam(i)

                     deltax=fitfactordeltax*abs(tempx)
                     if (deltax .le. 1e-15) deltax=1e-5
      
                     potparam(i)=tempx+deltax
                    
                     call poten(potparam,potright,r1,r2,theta)
                     !
                     potparam(i)=tempx-deltax
                     !
                     call poten(potparam,potleft,r1,r2,theta)
                     !
                     potparam(i)=tempx
                     rjacob(nrow,ncol)=(potright-potleft)/(2*deltax) 
                    
                     !
                  endif
                enddo ! --- ncol
                !
                numpar = ncol
                !
              endif     
              !
            enddo  ! ---  nrow
            !
            ! ssq  - weighted rms**2, rms  - root mean square deviation. 
            !
            ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
            rms=dsqrt(sum(eps(1:npts)*eps(1:npts))/npts)
            !
            ! Prepare the linear system a x = b as in the Newton fitting approach.  
            !
            if (itmax.ge.1) then
               !----- form the a and b matrix ------c
               ! form A matrix 
               do irow=1,numpar       
                 do icolumn=1,irow 
   
                   al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wtall(1:npts))

                   al(icolumn,irow)=al(irow,icolumn)
                  

                 enddo
               enddo
               !
               ! form B matrix 
               do irow=1,numpar      
                 bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wtall(1:npts))
               enddo   


               !
               ! Two types of the linear solver are availible: 
               ! 1. linur (integrated into the program, from Ulenikov Oleg)
               ! 2. dgelss - Lapack routine (recommended).
               !
               select case (trim(fit_type)) 
               !
               case default
                 !
                 write (f_out,"('fit_type ',a,' unknown')") trim(fit_type)
                 stop 'fit_type unknown'
                 !
               case('linur') 
                 !
                 call linur(numpar,numpar,al(1:numpar,1:numpar),bl(1:numpar),dx(1:numpar),ierror)
                 !
                 ! In case of dependent parameters  "linur" reports an error = ierror, 
                 ! which is a number of the dependent parameter. We can remove this paramter 
                 ! from the fit and set its value to zero. And start the iteration again. 
                 !
                 if (ierror.ne.0) then 
                   ncol=0
                   do i=1,parmax
                      if (ivar(i) .ne. 0) then
                        ncol=ncol+1
                        if  ( ncol.eq.ierror ) then 
                            ivar(i) = 0
                            potparam(i) = parold(i)
                            write(f_out,"(I4.0,'-th is out - ',a8)") i,nampar(i)
                        endif 
                      endif 
                   enddo 
                   cycle outer_loop    
                 endif 
                  !
               case ('dgelss')
                 !
                 ai(1:numpar,1:numpar) = al(1:numpar,1:numpar) 
               call dgelss(numpar,numpar,1,ai(1:numpar,1:numpar),numpar,bl(1:numpar),numpar,Tsing,-1.D-40,RANK0,wspace,lwork,info)

                 if (info/=0) then
                   write(f_out,"('dgelss:error',I4.0)") info
                   stop 'dgelss'
                 endif
                 !
                 dx(1:numpar) = bl(1:numpar)
              

               end select 
               !
               !----- update the parameter values ------!
               !
               
               ncol=0
               do i=1,parmax
                if (ivar(i) > 0) then
                     ncol=ncol+1
                     potparam(i)=potparam(i)+dx(ncol)
                    ! potparam(i)=dx(ncol)
                  endif
               enddo
          
               !
               ! write the potential parameters into the file, where  
               ! we store the current values 
               !
              200 continue

               open(f_pot,file='pot.fit',status='replace')
               !
               write(f_pot,*) parmax
               !
               do  i=1,parmax
                  write(f_pot,"(i8,e24.12,4x,a8)") ivar(i),potparam(i),nampar(i)
               end do
               !
               close(f_pot)
               !
               ! Estimate standard deviation error. 
               !
               if ( nused.ne.numpar ) then 
                 stadev=dsqrt(ssq/float(nused-numpar))
               else 
                 stadev=dsqrt(ssq/nused)
               endif
               !
               ! Estimate the standard errors for each parameter using 
               ! the inverse matrix of a. 
               !
               call invmat(al,ai,numpar,parmax)
               !
               sum_sterr=0.d0
               ncol = 0 
               do i=1,parmax
                  if (ivar(i) > 0) then
                      ncol=ncol+1
                     if (nused.eq.numpar) then  
                        sterr(ncol)=0
                     else
                        sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
                        sum_sterr=sum_sterr+abs(sterr(ncol)/potparam(i))
                     endif
                   endif
               enddo    
               !
               sum_sterr=sum_sterr/numpar 
               !
               ! This is how we define stability of the fit:
               ! as a relative change of stadev comparing with the step before. 
               !
               stability=abs( (stadev-stadev_old)/stadev )
               stadev_old=stadev
               !
            else
               !
               stadev=dsqrt(ssq/nused)
               !
            endif
            !
            ! Print the updated parameters. 
            !
            write(f_out,"(/'Potential paramters:')")
            !
            do i=1,parmax
              write (f_out,"(a8,4x,i2,e22.14)") nampar(i),ivar(i),potparam(i)
            enddo
            !
            ! Print the potential energy points into a separate unit. 
                 !
            inquire(f_potpoints,opened=ifopen)
            if ( ifopen ) then
               rewind(f_potpoints)
            else
               open (f_potpoints,file='abinitio.fit',status='replace' )
            endif
            !
            write(f_potpoints,"(1h1,5x,12('*'),' ab initio points ',  & 
                 12('*')// & 
                 4x,'r1',5x,'r2',5x,'theta',7x,'ab initio PES',3x, & 
                 'cal.PES',3x,'a-c',3x,'weight',3x,'energy'/)")
            !
            do nrow=1,pot_npts
              !
              r1   =local(1,nrow)
              r2   =local(2,nrow)
              theta=local(3,nrow)
              
               v = pot_values(nrow)-eps(nrow)
               write (f_potpoints,"(f6.4,' ',f6.4,' ',f6.2,' ',3(f15.8,' '),' ',e20.6,' ',f15.8)") r1,r2,theta, &
                      pot_values(nrow),v, & 
                      eps(nrow),wtnormal(nrow),abinitio_en(nrow)
              !
            enddo
            !
            ! Output some staqtistics and results 
            !
            !  only if we are fitting:  
            !
            if (itmax.ne.0) then
           
              !
              !still_run = .false.
              ! 
              write (f_out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
              l = 0 
              do i=1,parmax
                if (ivar(i) .ne. 0) then
                   l=l+1
                   ndigits = 0
                   conf_int = sterr(l)
                   do while (conf_int.le.10.0)
                     ndigits = ndigits +1 
                     conf_int = conf_int*10
                   enddo
                   !
                   !if (ndigits<1) ndigits = 12
                   !
                   if (conf_int>1e8) conf_int = 0 
                   !
                   write(fmt0,"(i2)") ndigits
                   fmt0 = adjustl(fmt0)
                   fmt1 = '(a8,i4,2x,f22.'//fmt0//'    ''    ('',i14,'')'')'
                   write (f_out,FMT=fmt1) nampar(i),ivar(i),potparam(i),nint(conf_int)
                else 
                   ndigits =2
                   if (potparam(i).ne.0.0) ndigits = 8

                   write(fmt0,"(i2)") ndigits
                   fmt0 = adjustl(fmt0)
                   fmt1 = '(a8,I4,2x,f22.'//fmt0//')'
                   write (f_out,FMT=fmt1) nampar(i),ivar(i),potparam(i)
                endif
              enddo  ! --- i
              !
            endif
!
!
!
!
!
!
!
!

          rms=dsqrt(sum(eps(1:npts)*eps(1:npts))/npts) ! NEW RMS!!!!!!!!!!!!!!


 !         stadev=dsqrt( sum( eps(1:npts)*eps(1:npts)*wtall(1:npts) )/npts )

          stadev=dsqrt( sum( eps(1:npts)*eps(1:npts)*wtnormal(1:npts) )/wtsum )



            !
            still_run = .false.
            !
            rewind(f_res)
            !
                    !
                        !
            ! Print out the ssq for the rovib. energies and pot. data points separetely:
            !
            ssq1 = 0 ; ssq2 = 0 
            !
            !
            !
            !wtsum = sum(wt_bit(1:npts))
            !
           ! if (wtsum/=0) ssq2 = sqrt( sum(eps(1:npts)**2*dble(wt_bit(1:npts)))/wtsum )

            rms2=sqrt(sum(eps(1:npts)**2)/pot_npts)
            rms1=sqrt(sum(eps(1:npts)**2)/pot_npts)
            !
            write (f_out,6550) fititer,nused,numpar,stadev,rms,stability
            write (f_res,6553) fititer,nused,numpar,stadev,rms1,rms2,stability
            !
            !
          enddo  ! --- fititer
          !
       enddo outer_loop

       !
       ! Print out the final information: paramteres and eneries.
       !
       rewind(f_res)
       !
       write (f_res,"(/36('-')/' Parameter | W |      Value        |',/36('-')/)")
       !
       do i=1,parmax
          write (f_res,"(A8,1X,I4,2X,e18.8)") nampar(i),ivar(i),potparam(i)
          write(33,'(g22.12)') potparam(i)
       enddo
       !

       !
       !
       do nrow=1,pot_npts
         !
         r1   =local(1,nrow)
         r2   =local(2,nrow)
         theta=local(3,nrow)
         !
         v = pot_values(nrow)-eps(nrow)
         write (f_res,"(f6.4,' ',f6.4,' ',f6.2,' ',3(f14.5,' '),' ',e1.6,' ',f14.5)") r1,r2,theta, &
                 pot_values(nrow),v, & 
                 eps(nrow),wtnormal(nrow),abinitio_en(nrow)

       open(rms_plot,file='rms.plot',status='unknown',action='write')
       write(rms_plot,"(f15.8,',',f15.8,',', f20.10,',',f15.8)") eps(nrow),abs(eps(nrow)), pot_values(nrow),abinitio_en(nrow)
         !
       enddo
       
       !
       ssq1 = 0 ; ssq2 = 0 
       !
       !
       !
       wtsum = sum(wt_bit(1:npts))
       !
       if (wtsum/=0) ssq2 = sqrt( sum(eps(1:npts)**2*dble(wt_bit(1:npts)))/wtsum )
       !
       write (f_res,6551) fititer,nused,numpar,stadev,ssq1,ssq2
       write (f_res,6550) fititer,nused,numpar,stadev,rms,stability


  

6550   format(/3X,67('-')/'   |  Iter  | Points | Params |   Weight RMS    |',&
       '    rms      | Stability |'/&
       3X,67('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',/3X,67('-')/)


6551   format(/3X,67('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq1     |   ssq2    |'/&
       3X,67('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',/3X,67('-')/)

6552   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    ssq_ener |   ssq_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,80('-')/)

6553   format(/3X,80('-')/'   |  Iter  | Points | Params |   Deviat    |',&
       '    rms_ener |   rms_pot | Stability |'/&
       3X,80('-')/,&
       '   | ',I6,' | ',I6,' | ',I6,' |  ',E12.5,' | ',E12.5,'  |  ',&
            E10.3,' |',E10.3,' |',/3X,80('-')/)

       !




  end subroutine fitting
  !



subroutine diag(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)

double precision,intent(in)   :: a(1:m,1:n),b(1:m),rcond
double precision,intent(out)  :: s(1:m,1:n),work(1:m)
integer,intent(in)            :: lwork,m,n,nrhs,lda,ldb
integer,intent(out)           :: info,rank




 call DGELSS(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)

return

end subroutine diag



  !
  subroutine linur(dimen,npar,coeff,constant,solution,error)

  integer,intent(in)  :: dimen,npar
  integer,intent(out) :: error 
  double precision,intent(in)  :: coeff(npar,npar),constant(npar)
  double precision,intent(out) :: solution(npar)
  double precision          :: a0(npar,npar)
  double precision          :: c
  integer                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(f_out,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(f_out,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1                                                                           
         c=0.0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0.0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0.0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine linur


!------------------------------------------!
  subroutine invmat(al,ai,dimen,npar)
  integer,intent(in)           :: npar,dimen
  double precision,intent(in)  :: al(npar,npar)
  double precision,intent(out) :: ai(npar,npar)
  double precision             :: h(npar),p,q
  integer                      :: i1,i2,k,i,j,k8,k9
      

    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine invmat
!------------------------------------------!



  end module fit_module
  !
  program driver

    use fit_module
    call fitting()
  end program driver


