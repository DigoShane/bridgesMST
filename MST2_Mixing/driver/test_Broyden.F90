!This is to test out Broyden Mixing Scheme
!!1. The main objective here is to test out how much error we generate if we did the
!!   broyden Mix calculation with a vector of size vlen withn vlen-MAxBasis zeros padded
!!   as opposed to just doing the calculation with MaxBasis elements.
!!1b. IS running the Broyden calculation for specific inputs based on results form crashed broyden calc.
!!    This also requires that we just run the Broyden iter for iter >1. Thus iterloop 2,MaxNumITer
!!2. We want only one proc to perform broyden mixing and then at the end it will BCAST
!!   its result to the other PEs. !!NB!! 2 has been skipped half way.
!!3. here we store the values in chunks, Basis< MaxBasis. We store 1:basis, then
!!   21:20+basis and then 41:40+Basis, the rest are 0, we wanted to test if it messes with anything.
!!4. Here we want to test out the case where we perform BroydenMixing with 10x1 array.
!!   We remove all the communication between the PEs, and want to see if all the PEs
!!   are returning the same answer. In order to check that the calculation done by
!!   each PE independently is the same as that being done in parallel, which we know is correct,
!!   we will run the code twice, once wil the communication turned off and then with communication.
!!   They are labelled 4a and 4b respectively. To make the calculation easier, we set the following
!!   Basis=MaxBAsis=vlen=8 and run it on 4 PEs. The no. of elements should be divisible by #PEs.
!!   The outputs are printed on different files with labels.

      program broyden
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : zero, one, two, three, ten2m2, ten2m3,ten2m8, &
                               ten2m9, czero, cone, PI4, THIRD, sqrtm1

    use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MyPE, GlobalSum, &
                         getCommunicator
   use GroupCommModule, only : initGroupComm, endGroupComm
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
   use GroupCommModule, only : syncAllPEsInGroup
   use GroupCommModule, only : createProcGrid, createBoxGroupFromGrid, &
                               getProcGridID, GlobalSumInGroup

   implicit none

  !!xDx!! Looping variables
   integer (kind=IntKind) :: ii, jj, j, k, i, typ
   integer(kind=IntKind) :: idt, id, idq, iter_broy!!xDx!!vlen=1031
  !!xDx!! Parameters
!!3.!!   integer (kind=IntKind) :: Basis=10, MaxBasis=50, vlen=1031, MaxRPMIter=10
!!4.!!   integer (kind=IntKind) :: Basis=8, MaxBasis=8, vlen=8, MaxRPMIter=5 
  integer (kind=IntKind) :: Basis=10,MaxBasis=10, vlen=10, MaxRPMIter=10
  !!xDx!! Iter variables
   integer (kind=IntKind) :: lastIter, nn, lastm1, NumRPMIter_broy=10, NumBroydenIter=10
!!xDx!! Mixing Parameters 
real(kind=RealKind) :: alpha=0.1, beta=1-0.1
!!xDx!! Used as pointer to original variable and pointer to unstable subspace variable.
   real(kind=RealKind), allocatable :: pvect_old(:), pvect_new(:)
   real(kind=RealKind), allocatable :: pdf(:), pf(:)
   real(kind=RealKind), allocatable :: pvold(:)
   real(kind=RealKind), allocatable :: pu(:,:), pw(:), pvt(:,:)!! MaxRPMIter=10
   real(kind=RealKind), allocatable :: pn(:), po(:)
   !!xDx!! Variables related to Broyden part
   real (kind=RealKind) :: w0, wtmp, dfnorm, fnorm, fac2, fac1
   real (kind=RealKind) :: aij, cmj, gmi
 !!xDx!! Message passing
   real (kind=RealKind) :: msgbuf(1:2)
   real (kind=RealKind), allocatable :: msgbuf1(:) 
   real (kind=RealKind), allocatable :: msgbuf2(:,:) 
   real(kind=RealKind), pointer :: pad(:), pbd(:)
   integer (kind=IntKind) :: INFO
   integer (kind=IntKind) :: dim(1), box(1), color, key, my_rank, group_size
   integer (kind=IntKind) :: itmp, ictxt1, ictxt2, grid, group
 !
   integer (kind=IntKind), allocatable :: ipiv(:)
!
   real(kind=RealKind), allocatable:: a_r(:,:), b_r(:,:)
   real(kind=RealKind), allocatable:: d_r(:,:), cm_r(:)
!
   complex(kind=CmplxKind), allocatable:: a_c(:,:), b_c(:,:)
   complex(kind=CmplxKind), allocatable:: d_c(:,:), cm_c(:)
!
  integer(kind=IntKind) :: GroupID

  real (kind=RealKind) :: rms
!
  character (len=50) :: PEStr

!! MPI initialization
   call initMPP()

   call initGroupComm()

   dim(1) = NumPEs
   call createProcGrid(1,'1-d Grid',dim)
   box(1)=NumPEs/2
!!??!!  call createBoxGroupFromGrid(grid,box,'Half')
!!??!!  GroupID = getGroupID('Half') 
   GroupID = getGroupID('1-d Grid') 


!! Opening file for output
   write(PEStr, '(I2)') MyPE+10
   open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )

                      
TYPITER: do typ=1,2   
if (typ .eq. 2) then  
!!4.!!MaxBasis=MaxBasis/NumPEs
!!4.!!Basis=MaxBasis
vlen=MaxBasis         
endif                 

allocate(pn(MaxBasis))
allocate(po(MaxBasis))

!!Allocation
allocate(pvold(vlen))
allocate(pvect_old(vlen))
allocate(pvect_new(vlen))
allocate(pdf(vlen))
allocate(pf(vlen))
allocate(pu(vlen,MaxRPMIter))
allocate(pw(MaxRPMIter))
allocate(pvt(vlen,MaxRPMIter))
      allocate( a_r( NumBroydenIter,NumBroydenIter ) )
      allocate( b_r( NumBroydenIter,NumBroydenIter ) )
      allocate( d_r( NumBroydenIter,NumBroydenIter ) )
      allocate( cm_r( NumBroydenIter ) )
      a_r = ZERO; b_r = ZERO; d_r = ZERO; cm_r = ZERO


!!xDx!! Allocating input -----------------------------------------------------------------------
   pvect_old=ZERO
   pvect_new=ZERO
!!3.!!   do i=1,Basis         
!!3.!!   pvect_old(i)=2*i
!!3.!!   pvect_new(i)=0.5*i
!!3.!!   pvect_old(20+i)=2*i   
!!3.!!   pvect_new(20+i)=0.5*i 
!!3.!!   pvect_old(40+i)=2*i   
!!3.!!   pvect_new(40+i)=0.5*i 
!!3.!!   enddo

!!4.!!   do i=1,Basis        
!!4.!!     if (typ .eq. 1) then
!!4.!!        pvect_old(i)=2*i
!!4.!!        pvect_new(i)=0.5*i
!!4.!!     else
!!4.!!        pvect_old(i)=2*(i+MyPE*Basis)
!!4.!!        pvect_new(i)=0.5*(i+MYPE*Basis)
!!4.!!     endif
!!4.!!   enddo

!!1.!!   do i=1,Basis        
!!1.!!        pvect_old(i)=2*i
!!1.!!        pvect_new(i)=0.5*i
!!1.!!   enddO
        pvect_old= (/ 0.0000721253, -0.0000392751, -0.0000392751,  &
                   -0.0000392751, -0.0000392751, -0.0000392750, -0.0000392750,&
                     -0.0000392750, -0.0000392750, -0.0000392749/)
        pvect_new= (/ 116.0718025336, -57.9967237465, -57.9966842149,  &
                   -57.9966442063, -57.9966037150, -57.9965627351, -57.9965212607,&
                  -57.9964792860, -57.9964368047, -57.9963938110/)


!!3.!!	if ( MyPE .eq. 0 ) then
!!3.!!	print *, 'printing out the input'            
!!3.!!	do i=1,MaxBasis                              
!!3.!!	print *, 'pvect_old(i),i',pvect_old(i),i     
!!3.!!	print *, 'pvect_new(i),i',pvect_new(i),i     
!!3.!!	enddo                                        
!!3.!!	endif
!!3.!!	call syncAllPEs()

!!Testing for 4.!!
!!Testing for 1.!!	call syncAllPEs()
!!Testing for 1.!!	print *, 'printing out the input for typ',typ            
!!Testing for 1.!!	do ii=0,NumPEs-1
!!Testing for 1.!!	   if ( MyPE .eq. ii ) then
!!Testing for 1.!!	      do i=1,MaxBasis                              
!!Testing for 1.!!	         print *, 'pvect_old(i),i,MyPE',pvect_old(i),i,MyPE     
!!Testing for 1.!!	         print *, 'pvect_new(i),i,MyPE',pvect_new(i),i,MyPE     
!!Testing for 1.!!	      enddo                                        
!!Testing for 1.!!	   endif
!!Testing for 1.!!	enddo
!!Testing for 1.!!	call syncAllPEs()
!!Testing for 4.!!



ITER:   do iter_broy=2,MaxRPMIter!!1a.!!1,MaxRPMIter, is for the rest.

!!1.!!if (typ .eq. 1 .and. iter_broy .eq. 10) then          
   pn(1:MaxBasis)=pvect_new(1:MaxBasis)               
   po(1:MaxBasis)=pvect_old(1:MaxBasis)               
!!1.!!elseif (typ .eq. 2 .and. iter_broy .eq. 10) then      
!!1.!!   pn(1:MaxBasis)=pn(1:MaxBasis)-pvect_new(1:MaxBasis)
!!1.!!   po(1:MaxBasis)=po(1:MaxBasis)-pvect_old(1:MaxBasis)
!!1.!!endif                                                 

   if ( iter_broy .eq. 1) then

!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
            do k = 1,vlen
               pf(k) = pvect_new(k) - pvect_old(k) !!cDt!!
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k) !!cDt!!
            enddo
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k) !!cDt!!
            enddo
!

   elseif ( iter_broy .gt. 1 ) then

!     ===============================================================
!     iter_count > 1: this is where the non-linear mixing is done
!     ===============================================================
      lastIter = lastIter + 1      ! update pointers
      lastm1   = lastIter - 1
!
!     ===============================================================
!     set current length of broyden cycle
!     ===============================================================
      if (  iter_broy > NumRPMIter_broy ) then 
         nn = NumRPMIter_broy
      else
         nn = lastIter
      endif
!
         w0= ten2m2
            dfnorm = zero
            fnorm = zero
               do k = 1,vlen !!cDt!!x2
                  pdf(k) = pvect_new(k) - pvect_old(k) - pf(k)!! dR^(iter)=R^(iter)-R^(iter-1)
                  pf(k)  = pvect_new(k) - pvect_old(k) !!R^(iter)
               enddo

!              =========================================================
!              find: fnorm  := |f|
!                    dfnorm := |df|
!              =========================================================
               do k = 1,vlen  !!xDx!! These are sums over all the atoms in the proc. 
                  dfnorm = dfnorm + pdf(k)*pdf(k)!!cDt!!
                  fnorm  = fnorm  + pf(k)*pf(k)!!cDt!!
               enddo
!
!!Testing 4.!!
!!Testing 4.!!!!Testing for 4.!!
!!Testing 4.!!	if ( (typ .eq. 1 .and. MyPE .eq. 0) .or. (typ .eq. 2)  )then
!!Testing 4.!!	print *, 'df^2,typ,MYPE,iter',dfnorm,typ,MYPE,iter_broy
!!Testing 4.!!	print *, ' f^2,typ,MYPE,iter', fnorm,typ,MYPE,iter_broy
!!Testing 4.!!	endif
!!Testing 4.!!!!Testing for 4.!!
!!Testing 4.!!
            msgbuf(1) = dfnorm!!cDt!!
            msgbuf(2) = fnorm!!cDt!!
!           ------------------------------------------------------------
	if ( typ .eq. 2) then
           call GlobalSumInGroup(GroupID, msgbuf,2)!!cDt!!
!! Sum over the sums of |dR|^2 and |R|^2 in each proc & returns resultant to all procs
!!xDx!! If the potential is modelled as a vector then based on the discretized space, each index of the vector would correspond to a unique point
!!xDx!! In the whole space. Since we are calculating dR and R locally around each atom, we need to calculate the actual norm by summing over the 
!!xDx!! potential form the other atoms. Thus the Global sum in Group. The sum over all the atoms in each proc is reflected in the calc being performed
!!xDx!! in the idq loop.
	endif
!        ------------------------------------------------------------
            dfnorm = sqrt( msgbuf(1) )!!cDt!!
            fnorm  = sqrt( msgbuf(2) )!!cDt!!
!           ============================================================
!           set: vector(2) := alpha*dR/|dR| + (vector(1) - vold)/|dR|
!                vold      := vector(1) 
!                vector(1) := dR/|dR|
!           ============================================================
!
            fac2 = one/dfnorm!!cDt!!
            fac1 = alpha*fac2
            do k = 1,vlen
               pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k)) 
               pvold(k)   = pvect_old(k)
               pvect_old(k) = fac2*pdf(k)
            enddo  
            call broy_sav_r( pu, pvt, pvect_old, pvect_new, iter_broy-1, & 
                                NumRPMIter_broy, vlen )!!cDt!!
!           ============================================================
!!xDx!! nn -> #iters-1, equiv to m in DD Johnson.
            do j = 1,nn - 1          ! off diagonal elements of a(i,j)
               do i = j+1,nn
                  aij = zero
                     do k = 1,vlen
                        aij = aij + pvt(k,j)*pvt(k,i)!!cDt!!
                     enddo
!                 ------------------------------------------------------
!                 call GlobalSumInGroup(GroupID,aij)
!                 ------------------------------------------------------
                  a_r(i,j) = aij
                  a_r(j,i) = aij
               enddo
            enddo
!
            do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
               aij = zero
               cmj = zero
                  do k = 1,vlen!!cDt!!x2
                     cmj = cmj + pvt(k,i)*pf(k) 
                     aij = aij + pvt(k,i)*pvt(k,i)
                  enddo
               a_r(i,i) = aij
               cm_r(i) = cmj
            enddo

	if ( typ .eq. 2) then
!!xDx!! the a_r(i,j) is a_{ij}/(w_iw_j) of Eqn. (13a). |||y for cm_r.
!           ------------------------------------------------------------
            call GlobalSumInGroup(GroupID,a_r,NumRPMIter_broy,NumRPMIter_broy)!!cDt!!
!           ------------------------------------------------------------
            call GlobalSumInGroup(GroupID,cm_r,NumRPMIter_broy)!!cDt!!
!           ------------------------------------------------------------
	endif
!           ============================================================
!           shift down weights in stack
!
!           (TCS, bug fixed 8/5/97: replace  RPMIter by  RPMIter-1 -> see broy_sav)
!           ============================================================
               if (  iter_broy-1 > NumRPMIter_broy ) then
                  do i = 1,NumRPMIter_broy-1
                     pw(i) = pw(i+1)
                  enddo
               endif
!
               rms = 0
               wtmp = zero
               if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
               if ( wtmp < ONE ) wtmp = ONE
               if (  iter_broy > NumRPMIter_broy ) then
                  pw(NumRPMIter_broy) = wtmp
               else
                  pw(lastIter) = wtmp       !w(lastm1)=wtmp
               endif
!              =========================================================
!              now calculate the b-matrix:
!              b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1 !!xDx!! Eqn. (13a). Assuming this is just the standard calc of \beta.
!                 ======================================================
!                 use lapack !!nDn!! everything is local now.
!                 ======================================================
                  do i = 1,nn
                     do j = 1,nn
                        b_r(j,i) = a_r(j,i)*pw(j)*pw(i)!!cDt!!
                     enddo
                     b_r(i,i) = w0**2 + a_r(i,i)*pw(i)*pw(i)!!cDt!!
                  enddo
                  if ( .not.allocated(ipiv) ) then
                     allocate( ipiv(NumBroydenIter) )
                  endif
!                 ------------------------------------------------------
                  call dgetrf( nn, nn, b_r, NumRPMIter_broy,ipiv, info ) !!cDt!!
!                 ------------------------------------------------------
                  call dgetri( nn, b_r, NumRPMIter_broy, ipiv, d_r, nn, info )!!cDt!!
!                 ------------------------------------------------------
!                 write(6,*) ' optimum lwork', d(1,1)
!
               do k = 1,vlen
                  pvect_new(k) = pvold(k) + alpha*pf(k)!!cDt!!
               enddo
               do i = 1,nn
                  gmi = zero
                  do j = 1,nn!!cDt!!
                     gmi = gmi + cm_r(j)*b_r(j,i)*pw(j) 
                  enddo
                  do k = 1,vlen
                     pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)!!cDt!!
                  enddo
               enddo
            

   endif   

enddo ITER

!!Deallocation
deallocate(pvold)
deallocate(pvect_old)
deallocate(pvect_new)
deallocate(pdf)
deallocate(pf)
deallocate(pu)
deallocate(pw)
deallocate(pvt)
deallocate( a_r )
deallocate( b_r )
deallocate( d_r )
deallocate( cm_r )

!!Testing for 4.!!
!!Testing for 4.!!	print *,'Output for type=',typ
!!Testing for 4.!!
!!Testing for 4.!!	   do ii=0,NumPEs-1
!!Testing for 4.!!	    if (MyPE .eq. ii) then
!!Testing for 4.!!	     do k=1,MaxBasis
!!Testing for 4.!!	      print *, 'pn,k,MyPE',pn(k),k,MyPE
!!Testing for 4.!!	      print *, 'po,k,MyPE',pn(k),k,MyPE
!!Testing for 4.!!	     enddo
!!Testing for 4.!!	    endif
!!Testing for 4.!!	   enddo
!!Testing for 4.!!
!!Testing for 4.!!	call syncAllPEs()


deallocate(pn)
deallocate(po)

enddo TYPITER

	   do ii=0,NumPEs-1
	    if (MyPE .eq. ii) then
	     do k=1,MaxBasis
	      print *, 'pn,k,MyPE',pn(k),k,MyPE
	      print *, 'po,k,MyPE',pn(k),k,MyPE
	     enddo
	    endif
	   enddo



!!3.!!	  do k=1,MaxBasis                         
!!3.!!	   do ii=0,NumPEs-1                       
!!3.!!	    if (MyPE .eq. ii) then                
!!3.!!	      print *, 'pn,k,MyPE',pn(k),k,MyPE   
!!3.!!	    endif                                 
!!3.!!	   enddo                                  
!!3.!!	  enddo                                   


call endGroupComm()
call endMPP()

      end program broyden



   subroutine broy_sav_r( fins, fots, vector_old, vector_new, itscf,  & 
                          istore, ivsiz )
!  ==================================================================
!!xDx!! broy_sav_r( pu, pvt, pvect_old, pvect_new, iter-1, NumBroydenIter, vlen)
!!xDx!! where the subroutine assigns value to fins and fots.
!!xDx!! Assigns the matrix u=[u^(1),u^(2),u^(3),...,u^(m)] where m=NumRPMIter_broy
!!xDx!! and u^(i) is as in (13b.iii) in DD Johnson. 
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : zero, one, two, three, ten2m2, ten2m3,ten2m8, &
                               ten2m9, czero, cone, PI4, THIRD, sqrtm1

   implicit none
!
   integer(kind=IntKind), intent(in) :: itscf
   integer(kind=IntKind), intent(in) :: ivsiz
   integer(kind=IntKind), intent(in) :: istore
!
   real(kind=RealKind), intent(in) :: vector_old(1:ivsiz)
   real(kind=RealKind), intent(in) :: vector_new(1:ivsiz)
   real(kind=RealKind), intent(inout) :: fins(1:ivsiz,1:10) !!1:vlen and 1:MaxRPMIter
   real(kind=RealKind), intent(inout) :: fots(1:ivsiz,1:10)
!
   integer(kind=IntKind) :: i, j

!
!  ==================================================================
!  write(6,'('' IN BROY_SAV: istore,itscf '',2i5)') istore,itscf
!  ==================================================================
   if ( itscf <= istore ) then !!xDx!! if ( iter-1 <= NumRPMIter_broy )
!     ===============================================================
!     Load the first istore iterations in increasing iteration count
!     ===============================================================
      do i = 1,ivsiz 
         fins(i,itscf) = vector_new(i) 
      enddo
!!xDx!! pu(1:vlen,iter-1)=pvect_new(1:vlen)
!
      do i = 1,ivsiz
         fots(i,itscf) = vector_old(i)
      enddo
!
!!xDx!! pvt(1:vlen,iter-1)=pvect_old(1:vlen). This is now stored as (dR/|dR|)
   else
!     ===============================================================
!     Re-load so that the ordering is in increasing iteration count
!     ===============================================================
      do j = 1,istore - 1
!        write(6,'('' IN BROY_SAV: j,j+1 '',2i5)') j,j+1
         do i = 1,ivsiz
            fins(i,j) = fins(i,j+1)
         enddo
!!xDx!! pu(1:vlen,j)=pu(1:vlen,j+1) shifting each colmn one left.
!
         do i = 1,ivsiz
            fots(i,j) = fots(i,j+1)
         enddo
!!xDx!! pvt(1:vlen,j)=pvt(1:vlen,j+1) shifting each colmn one left.
!
      enddo
!     ===============================================================
!     Load current charge densities in the last storage location
!     ===============================================================
      do i = 1,ivsiz
         fins(i,istore) = vector_new(i)
      enddo
!!xDx!! loading the last column as the latest entry.
!
      do i = 1,ivsiz
         fots(i,istore) = vector_old(i)
      enddo
!!xDx!! loading the last column as the latest entry.
!
   endif
!
   end subroutine broy_sav_r



