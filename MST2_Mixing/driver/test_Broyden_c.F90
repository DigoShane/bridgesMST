!!xDx!! This is to test out Broyden Mixing Scheme for complex valued functions.
!!xDx!!========================================================================
!!xDx!!         The program no. is controlled by the variable 'typ'
!!xDx!!========================================================================
!!xDx!!1. We just write a simple bryden Mixing process in order to see whether
!!xDx!!   the magnitude of the output was just as high as the magnitude of the
!!xDx!!   input. We are just interested in the output after one iteration.
!!xDx!!   We want to run this for different vector sizes, thus the input
!!xDx!!   are
!!xDx!!   1. vlen=10, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha= 0.1, beta=1-alpha
!!xDx!!   2. vlen=2, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha= 0.1, beta=1-alpha
!!xDx!!   N.B. :- example of typ would be 11, 12,etc.
!!xDx!!2. In order to be sure its working fine, we test it against a Function
!!xDx!!   Besure to chekc out the mathematica file, Broyden_Test.nb in my laptop.
!!xDx!!   1. {y+5==x}&&{x^2-8==y}. For this, the variables are:-
!!xDx!!      vlen=2, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha=0.1, beta=1-alpha
!!xDx!!   2. {y*x^2+5==x}&&{y^3+x^2-8==y}. For this, the variables are:-
!!xDx!!      vlen=2, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha=0.1, beta=1-alpha
!!xDx!!   3. {Iy*x^2+5==x}&&{y^3+x^2-8==y}  or {(y*x^2+5)^0.02==x}&&{(y^3+x^2-8)^0.2==y}.
!!xDx!!      For both cases, the variables are:-
!!xDx!!      vlen=2, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha=0.1, beta=1-alpha
!!xDx!!   4. {x^2-2==x}. For this, the variables are:-
!!xDx!!      vlen=1, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha=0.1, beta=1-alpha
!!xDx!!   5. {x^2-2==x}&&{y^2-2==y}. For this, the variables are:-
!!xDx!!      vlen=2, Basis=10, MaxRPMIter=10, NumRPMIter_broy=10, NumBroydenIter=10, alpha=0.1, beta=1-alpha
!!xDx!!      First a real vlaued i/p and then a complex valued i/p.

      program broyden
   use KindParamModule, only : IntKind, RealKind, CmplxKind, QuadCmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : zero, one, two, three, ten2m2, ten2m3,ten2m8, &
                               ten2m9, ten2m5, czero, cone, PI4, THIRD, sqrtm1

    use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MyPE, GlobalSum, &
                         getCommunicator
   use GroupCommModule, only : initGroupComm, endGroupComm
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
   use GroupCommModule, only : syncAllPEsInGroup
   use GroupCommModule, only : createProcGrid, createBoxGroupFromGrid, &
                               getProcGridID, GlobalSumInGroup

   implicit none

  !!xDx!! Looping/Choice variables
   integer (kind=IntKind) :: ii, jj, j, k, i, typ=25
   integer(kind=IntKind) :: idt, id, idq, iter_broy
  !!xDx!! Parameters
  integer (kind=IntKind) :: Basis, vlen, MaxRPMIter
  real(kind=RealKind) :: tol
  !!xDx!! Iter variables
   integer (kind=IntKind) :: lastIter, nn, lastm1, NumRPMIter_broy, NumBroydenIter
  !!xDx!! Mixing Parameters 
   real(kind=RealKind) :: alpha, beta
  !!xDx!! Used as pointer to original variable and pointer to unstable subspace variable.
   complex(kind=CmplxKind), allocatable :: pvect_old(:), pvect_new(:)
   complex(kind=CmplxKind), allocatable :: pdf(:), pf(:)
   complex(kind=CmplxKind), allocatable :: pvold(:)
   complex(kind=CmplxKind), allocatable :: pu(:,:), pw(:), pvt(:,:)
   complex(kind=CmplxKind), allocatable :: pn(:), po(:)
   !!xDx!! Variables related to Broyden part
   real (kind=RealKind) :: w0, wtmp, dfnorm, fnorm, fac2, fac1
   complex(kind=CmplxKind) :: aij, gmi, cmj
 !!xDx!! Message passing
   complex (kind=CmplxKind) :: msgbuf(1:2)
   complex (kind=CmplxKind), allocatable :: msgbuf1(:) 
   complex (kind=CmplxKind), allocatable :: msgbuf2(:,:) 
   complex (kind=CmplxKind), pointer :: pad(:), pbd(:)
   integer (kind=IntKind) :: INFO
   integer (kind=IntKind) :: dim(1), box(1), color, key, my_rank, group_size
   integer (kind=IntKind) :: itmp, ictxt1, ictxt2, grid, group
 !
   integer (kind=IntKind), allocatable :: ipiv(:)
!
   complex(kind=CmplxKind), allocatable:: a_c(:,:), b_c(:,:), cm_c(:)
   complex(kind=CmplxKind), allocatable:: work(:)
   integer (kind=IntKind) :: lwork
   complex (kind=CmplxKind) :: tmp
!
  integer(kind=IntKind) :: GroupID

  real (kind=RealKind) :: rms
!
  character (len=50) :: PEStr



!!xDx!! Testing assigning & printing variables
  complex(kind=CmplxKind) :: test
!!xDx!! Testing assigning & printing variables

!! MPI initialization
   call initMPP()

   call initGroupComm()

!!xDx!!Allocatin variables
if (typ == 11) then
vlen=10; Basis=10; MaxRPMIter=10; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
elseif (typ == 12) then                       
vlen=2; Basis=10; MaxRPMIter=10; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
elseif (typ == 21) then
vlen=2; Basis=10; MaxRPMIter=40; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
tol = ten2m9
elseif (typ == 22) then
vlen=2; Basis=10; MaxRPMIter=40; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
tol = ten2m9
elseif (typ == 23) then
vlen=2; Basis=10; MaxRPMIter=40; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
tol = ten2m5
elseif (typ == 24) then
vlen=1; Basis=10; MaxRPMIter=40; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
tol = ten2m9
elseif (typ == 25) then
vlen=2; Basis=10; MaxRPMIter=40; NumRPMIter_broy=10; NumBroydenIter=10; alpha= 0.1D0; beta=1-alpha
tol = ten2m9
else
print *, 'Please eneter a vaild typ'
endif
!!xDx!!Allocation arrays
allocate(pn(MaxRPMIter))
allocate(po(MaxRPMIter))
allocate(pvold(vlen))
allocate(pvect_old(vlen))
allocate(pvect_new(vlen))
allocate(pdf(vlen))
allocate(pf(vlen))
allocate(pu(vlen,MaxRPMIter))
allocate(pw(MaxRPMIter))
allocate(pvt(vlen,MaxRPMIter))
      allocate( a_c( NumBroydenIter,NumBroydenIter ) )
      allocate( cm_c( NumBroydenIter ) )
      a_c = CZERO; cm_c = CZERO

   write(PEStr, '(I2)') MyPE+10

!!xDx!! Allocating input -----------------------------------------------------------------------
   pvect_old=CZERO
   pvect_new=CZERO

  if ( typ .gt. 10 .and. typ .lt. 20) then

    pn(1) =cmplx(-13447766.7262306D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(2) =cmplx(-13288002.1753515D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(3) =cmplx(-13130135.5771736D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(4) =cmplx(-12974144.3846144D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(5) =cmplx(-12820006.3184435D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(6) =cmplx(-12667699.3641010D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(7) =cmplx(-12517201.7685530D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(8) =cmplx(-12368492.0371847D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(9) =cmplx(-12221548.9307306D+000,0.000000000000000D+000, kind=CmplxKind)
    pn(10)=cmplx(-12076351.4622412D+000,0.000000000000000D+000, kind=CmplxKind)
    po(1) =cmplx(-13447766.7308236D+000,1.966346749934326D-015, kind=CmplxKind)
    po(2) =cmplx(-13288002.1799445D+000,1.966350770803877D-015, kind=CmplxKind)
    po(3) =cmplx(-13130135.5817666D+000,1.957557378979791D-015, kind=CmplxKind)
    po(4) =cmplx(-12974144.3892073D+000,1.965821466841742D-015, kind=CmplxKind)
    po(5) =cmplx(-12820006.3230365D+000,1.970609412568363D-015, kind=CmplxKind)
    po(6) =cmplx(-12667699.3686940D+000,1.974348450617597D-015, kind=CmplxKind)
    po(7) =cmplx(-12517201.7731460D+000,1.964812199147513D-015, kind=CmplxKind)
    po(8) =cmplx(-12368492.0417776D+000,1.966671533894002D-015, kind=CmplxKind)
    po(9) =cmplx(-12221548.9353236D+000,1.957293682164603D-015, kind=CmplxKind)
    po(10)=cmplx(-12076351.4668341D+000,1.965604261488643D-015, kind=CmplxKind)

    do i=1,vlen
      pvect_new(i)=pn(i)
      pvect_old(i)=po(i)
    enddo
  elseif ( typ .eq. 21) then
    pvect_new(1)=1D0 !!x
    pvect_new(2)=-1D0 !!y
  elseif ( typ .eq. 22) then
   ! pvect_new(1) = cmplx( 1.3D0,-0.6D0) !!x
   ! pvect_new(2) = cmplx(-0.8D0,-1.55D0)!!y
   ! pvect_new(1) = cmplx( 2.7D0, 0.0D0) !!x
   ! pvect_new(2) = cmplx(-0.2D0, 0.0D0)!!y
    pvect_new(1) = dcmplx(-1.4889355994121622D0,-1.1436970024236717D0) !!x
    pvect_new(2) = dcmplx(-0.788132877825017D0 , 1.6949445795623577D0)!!y
  elseif ( typ .eq. 23) then
   ! pvect_new(1) = cmplx(-1.4889355994121622D0,-1.1436970024236717D0) !!x
   ! pvect_new(2) = cmplx(-0.788132877825017D0 , 1.6949445795623577D0)!!y
   ! pvect_new(1) = cmplx( 3.077790D0,0D0) !!x
   ! pvect_new(2) = cmplx(-0.852735D0, 0D0)!!y
    pvect_new(1) = dcmplx( 1.0370D0,-0.0027D0) !!x
    pvect_new(2) = dcmplx( 1.2975D0,-0.773364D0)!!y
  elseif ( typ .eq. 24) then
    pvect_new(1)=1D0 !!x
  elseif ( typ .eq. 25) then
   ! pvect_new(1)=1D0 !!x
   ! pvect_new(2)=0D0 !!y
    pvect_new(1) = dcmplx(-1.0370000000D0,-0.0027000000D0) !!x
    pvect_new(2) = dcmplx( 1.2975000000D0,-0.773364000000D0)!!y
    print *, 'pvect_new(1)',pvect_new(1)
    print *, 'pvect_new(2)',pvect_new(2)
    test= dcmplx( 1.0370D0, -0.00270D0)
    !print *, 'test',test
    !print *, 'alpha',alpha
    !print *, 'beta',beta
  endif

!!xDx!!------------------------------------------------------------------------------------------------ 
!!iDi!! Debugging
             do k=0,NumPEs-1 
              if ( k .eq. MyPE ) then 
               open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
               write (MyPE+10,1120) '========================'
               write (MyPE+10,1340) 'MyPE','pvect_guess(:),iter',MyPE,0
               write (MyPE+10,1120) '========================'
               do i =1,vlen
                  write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
               enddo
               write (MyPE+10,1120) '========================'
               write (MyPE+10,1340) 'MyPE','pvect_old(:),iter',MyPE,0
               write (MyPE+10,1120) '========================'
               do i =1,vlen
                  write (MyPE+10,1241) 'i,MyPE','pvect_old',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
               enddo
              endif
             enddo
1120 format(X,A)
1340 format(2(X,A),(X,I4),(X,I4))
1241 format(2(X,A),2(X,I4),(X,F20.10,SP,F16.10,A))
!!iDi!! Debugging
!!xDx!!------------------------------------------------------------------------------------------------ 


ITER:   do iter_broy=1,MaxRPMIter!!1a.!!1,MaxRPMIter, is for the rest.

!!cDc!!        print *, 'iter_broy=',iter_broy

!!xDx!! setting pvect_new=F[pvect_old]
   if( typ .eq. 21 ) then
     pvect_old = pvect_new
     pvect_new(1)=pvect_old(2)+5
     pvect_new(2)=pvect_old(1)**2-8
   elseif( typ .eq. 22 ) then
     pvect_old = pvect_new
     pvect_new(1)=pvect_old(2)*pvect_old(1)**2+5
     pvect_new(2)=pvect_old(2)**3+pvect_old(1)**2-8
   elseif( typ .eq. 23 ) then
     pvect_old = pvect_new
     pvect_new(1)=pvect_old(2)*pvect_old(1)**2+5
     pvect_new(1)=pvect_new(1)**0.02
     pvect_new(2)=pvect_old(2)**3+pvect_old(1)**2-8
     pvect_new(2)=pvect_new(2)**0.2
   elseif( typ .eq. 24 ) then
     pvect_old = pvect_new
     pvect_new(1)=pvect_old(1)**2-2
   elseif( typ .eq. 25 ) then
     pvect_old = pvect_new
     pvect_new(1)=pvect_old(1)**2-2
     pvect_new(2)=pvect_old(2)**2-2
     !print *, 'pvect_new(1)',pvect_new(1)
     !print *, 'pvect_new(2)',pvect_new(2)
   endif

!!cDc!!,!!xDx!!----------------------------------------------------------------------------------------------
!!cDc!!,!!iDi!! Debugging
!!cDc!!,             do k=0,NumPEs-1 
!!cDc!!,              if ( k .eq. MyPE ) then 
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','F(v),iter',MyPE,0
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','v,iter',MyPE,0
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_old',i,MyPE,real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,              endif
!!cDc!!,             enddo
!!cDc!!,!!iDi!! Debugging
!!cDc!!,!!xDx!!----------------------------------------------------------------------------------------------

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
               !print *, 'pf(k),k',pf(k),k
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k) !!cDt!!
               !print *, 'pvold(k),k',pvold(k),k
            enddo
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k) !!cDt!!
               !print *, 'alpha*pf(k),k',alpha*pf(k),k
               !print *, 'pvect_new(k),k',pvect_new(k),k
            enddo
!

!!cDc!!,!!xDx!!-----------------------------------------------------------------------------------------------
!!cDc!!,!!iDi!! Debugging
!!cDc!!,             do k=0,NumPEs-1 
!!cDc!!,              if ( k .eq. MyPE ) then 
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pf,iter',MyPE,0
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pf(i)),aimag(pf(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pvect_old,iter',MyPE,0
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_old',i,MyPE,real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pvect_new,iter',MyPE,0
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,              endif
!!cDc!!,             enddo
!!cDc!!,!!iDi!! Debugging
!!cDc!!,!!xDx!!-----------------------------------------------------------------------------------------------


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
                  print *, '1 pf,k,iter',pf(k),k,iter_broy
                  print *, '1 pdf,k,iter',pdf(k),k,iter_broy
                  print *, '1 pvect_new,k,iter',pvect_new(k),k,iter_broy
                  print *, '1 pvect_old,k,iter',pvect_old(k),k,iter_broy
               enddo
               do k = 1,vlen  !!xDx!! These are sums over all the atoms in the proc. 
                  dfnorm = dfnorm + real(pdf(k)*conjg(pdf(k)),kind=RealKind)
                  fnorm  = fnorm  + real(pf(k)*conjg(pf(k)),kind=RealKind) 
               enddo
!
            dfnorm = sqrt( dfnorm )
            print *, '2 dfnorm',dfnorm
            fnorm  = sqrt( fnorm )
            print *, '2 fnorm',fnorm
!           ============================================================

            if ( int(typ/10,kind=IntKind) .eq. 2 .and. fnorm .lt. tol) then
               write (MyPE+10,1120) '========================'
               write (MyPE+10,1120) '!!xDx!!!xDx!!!xDx!!!xDx!!!xDx!!'
               write (MyPE+10,1120) '========================'
               write (MyPE+10,1120) 'Solution Converged with'
               write (MyPE+10,1341) 'MyPE','fnorm,iter',MyPE,fnorm,iter_broy
               write (MyPE+10,1341) 'MyPE','dfnorm,iter',MyPE,dfnorm,iter_broy
               write (MyPE+10,1120) '========================'
               write (MyPE+10,1120) '!!xDx!!!xDx!!!xDx!!!xDx!!!xDx!!'
               write (MyPE+10,1120) '========================'
               EXIT
               1341 format(2(X,A),(X,I4),(X,F20.10),(X,I4))
            endif

!!cDc!!,!!xDx!!-----------------------------------------------------------------------------------------------
!!cDc!!,!!iDi!! Debugging
!!cDc!!,             do k=0,NumPEs-1 
!!cDc!!,              if ( k .eq. MyPE ) then 
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1341) 'MyPE','dfnorm,iter',MyPE,dfnorm,iter_broy
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1341) 'MyPE','fnorm,iter',MyPE, fnorm,iter_broy
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,              endif
!!cDc!!,             enddo
!!cDc!!,!!iDi!! Debugging
!!cDc!!,!!xDx!!-----------------------------------------------------------------------------------------------

!           set: vector(2) := alpha*dR/|dR| + (vector(1) - vold)/|dR|
!                vold      := vector(1) 
!                vector(1) := dR/|dR|
!           ============================================================
!
            fac2 = one/dfnorm
            fac1 = alpha*fac2
            do k = 1,vlen
               pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k)) 
               pvold(k)   = pvect_old(k)
               pvect_old(k) = fac2*pdf(k)
               if (iter_broy .eq. 2) then
                  print *, '3 pvect_new,pvold,pvect_old,k',pvect_new(k),pvold(k),pvect_old(k),k
               endif
            enddo  

!!,!!xDx!!------------------------------------------------------------------------------------------------ 
!!,!!iDi!! Debugging
!!,             do k=0,NumPEs-1 
!!,              if ( k .eq. MyPE ) then 
!!,               write (MyPE+10,1120) '                        '
!!,               write (MyPE+10,1120) '========================'
!!,               write (MyPE+10,1340) 'MyPE','u^n(:),iter',MyPE,iter_broy
!!,               write (MyPE+10,1120) '========================'
!!,               do i =1,vlen
!!,                  write (MyPE+10,1241) 'i,MyPE','u^n',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!,               enddo
!!,               write (MyPE+10,1120) '                        '
!!,               write (MyPE+10,1120) '========================'
!!,               write (MyPE+10,1340) 'MyPE','DF(:),iter',MyPE,iter_broy
!!,               write (MyPE+10,1120) '========================'
!!,               do i =1,vlen
!!,                  write (MyPE+10,1241) 'i,MyPE','DF',i,MyPE,real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!,               enddo
!!,               write (MyPE+10,1120) '                        '
!!,              endif
!!,             enddo
!!,!!iDi!! Debugging
!!,!!xDx!!------------------------------------------------------------------------------------------------ 

!           ============================================================
            call broy_sav_c( pu, pvt, pvect_old, pvect_new, iter_broy-1, & 
                                NumRPMIter_broy, vlen )!!cDt!!
!           ============================================================
!!xDx!! nn -> #iters-1, equiv to m in DD Johnson.
            do j = 1,nn - 1          ! off diagonal elements of a(i,j)
               do i = j+1,nn
                  aij = czero
                     do k = 1,vlen
                       aij = aij + conjg(pvt(k,j))*pvt(k,i)
                     enddo
                  a_c(i,j) = aij
                  a_c(j,i) = conjg(aij)
               enddo
            enddo
!
            do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
               aij = czero
               cmj = czero
                  do k = 1,vlen!!cDt!!x2
                     cmj = cmj + conjg(pvt(k,i))*pf(k) 
                     aij = aij + conjg(pvt(k,i))*pvt(k,i)
                  enddo
               a_c(i,i) = aij
               cm_c(i) = cmj
            enddo

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
               rms = dfnorm !!cDc!! ZERO
               wtmp = ZERO
               if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
               if ( wtmp < ONE ) wtmp = ONE
               !!cDc!!,!!xDx!! Trying
               !!cDc!!,    wtmp = 43.47
               !!cDc!!,!!xDx!! Trying
               if (  iter_broy > NumRPMIter_broy ) then
                  pw(NumRPMIter_broy) = wtmp
               else
                  pw(lastIter) = wtmp       !w(lastm1)=wtmp
               endif
                 !! do i = 1,NumRPMIter_broy
                 !!    print *,'pw,i',pw(i),i
                 !! enddo
                 !! print *,'rms',rms
!              =========================================================
!              now calculate the b-matrix:
!              b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1 !!xDx!! Eqn. (13a)
!              ======================================================
               allocate( b_c( nn,nn ) )
               b_c = CZERO
               do i = 1,nn
                  do j = 1,nn
                     b_c(j,i) = a_c(j,i)*pw(j)*pw(i)!!cDt!!
                  enddo
                  b_c(i,i) = w0**2 + a_c(i,i)*pw(i)*pw(i)!!cDt!!
               enddo
               if ( .not.allocated(ipiv) ) then
                  allocate( ipiv(nn) )
               endif
!              ------------------------------------------------------
               call zgetrf( nn, nn, b_c, nn, ipiv, info ) !!cDt!!
!              ------------------------------------------------------
               call zgetri( nn, b_c, nn, ipiv, tmp, -1, info )
               lwork=int(real(tmp,kind=RealKind))
               allocate( work(1:lwork) )
               call zgetri( nn, b_c, nn, ipiv, work, lwork, info) 
               deallocate( work )
!

!!xDx!!------------------------------------------------------------------------------------------------ 
!!iDi!! Debugging
             do k=0,NumPEs-1 
              if ( k .eq. MyPE ) then 
               write (MyPE+10,1120) '                        '
               write (MyPE+10,1120) '========================'
               write (MyPE+10,1340) '4 MyPE','b_c,iter',MyPE,iter_broy
               write (MyPE+10,1120) '========================'
               do i =1,nn
                  write (MyPE+10,1230) ( b_c(i,j), j=1,nn )
               enddo
               write (MyPE+10,1120) '                        '
              endif
             enddo
1230 format(14f16.6)
!!iDi!! Debugging
!!xDx!!------------------------------------------------------------------------------------------------ 

               do k = 1,vlen
                  pvect_new(k) = pvold(k) + alpha*pf(k)!!cDt!!
                  if (iter_broy .eq. 2) then
                     print *, '5 pvect_new,k',pvect_new(k),k
                  endif
               enddo
               do i = 1,nn
                  gmi = zero
                  do j = 1,nn
                     gmi = gmi + cm_c(j)*b_c(j,i)*pw(j) 
                  enddo
                     if (iter_broy .eq. 2) then
                        print *, '5 gmi',gmi
                     endif
                  do k = 1,vlen
                     pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)!!cDt!!
                     if (iter_broy .eq. 2) then
                        print *, '5 pvect_new,k',pvect_new(k),k
                     endif
                  enddo
               enddo
            
               deallocate( b_c )

!!cDc!!,!!xDx!!------------------------------------------------------------------------------------------------ 
!!cDc!!,!!iDi!! Debugging
!!cDc!!,             do k=0,NumPEs-1 
!!cDc!!,              if ( k .eq. MyPE ) then 
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pvect_new(:),iter',MyPE,iter_broy
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pvect_old(:),iter',MyPE,iter_broy
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_old',i,MyPE,real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               write (MyPE+10,1120) '                        '
!!cDc!!,              endif
!!cDc!!,             enddo
!!cDc!!,!!iDi!! Debugging
!!cDc!!,!!xDx!!------------------------------------------------------------------------------------------------ 

   endif   

enddo ITER

!!cDc!!,!!xDx!!------------------------------------------------------------------------------------------------ 
!!cDc!!,!!iDi!! Debugging
!!cDc!!,             do k=0,NumPEs-1 
!!cDc!!,              if ( k .eq. MyPE ) then 
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pvect_new(:),iter',MyPE,iter_broy+1
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               write (MyPE+10,1340) 'MyPE','pvect_old(:),iter',MyPE,iter_broy+1
!!cDc!!,               write (MyPE+10,1120) '========================'
!!cDc!!,               do i =1,vlen
!!cDc!!,                  write (MyPE+10,1241) 'i,MyPE','pvect_old',i,MyPE,real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!cDc!!,               enddo
!!cDc!!,               close ( MyPE+10)
!!cDc!!,              endif
!!cDc!!,             enddo
!!cDc!!,!!iDi!! Debugging
!!cDc!!,!!xDx!!------------------------------------------------------------------------------------------------ 

!!Deallocation
deallocate(pvold)
deallocate(pn)
deallocate(po)
deallocate(pvect_old)
deallocate(pvect_new)
deallocate(pdf)
deallocate(pf)
deallocate(pu)
deallocate(pw)
deallocate(pvt)
deallocate( a_c )
deallocate( cm_c )

call endGroupComm()
call endMPP()

      end program broyden

   subroutine broy_sav_c( fins, fots, vector_old, vector_new, itscf,  & 
                          istore, ivsiz )
!  ==================================================================
!!xDx!! broy_sav_c( pu, pvt, pvect_old, pvect_new, iter-1, NumBroydenIter, vlen)
!!xDx!! where the subroutine assigns value to fins and fots.
!!xDx!! Assigns the matrix u=[u^(1),u^(2),u^(3),...,u^(m)] where m=NumRPMIter_broy
!!xDx!! and u^(i) is as in (13b.iii) in DD Johnson. 
   use KindParamModule, only : IntKind, RealKind, CmplxKind, QuadCmplxKind
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use MathParamModule, only : zero, one, two, three, ten2m2, ten2m3,ten2m8, &
                               ten2m9, czero, cone, PI4, THIRD, sqrtm1

   implicit none
!
   integer(kind=IntKind), intent(in) :: itscf
   integer(kind=IntKind), intent(in) :: ivsiz
   integer(kind=IntKind), intent(in) :: istore
!
   complex(kind=CmplxKind), intent(in) :: vector_old(1:ivsiz)
   complex(kind=CmplxKind), intent(in) :: vector_new(1:ivsiz)
   complex(kind=CmplxKind), intent(inout) :: fins(1:ivsiz,1:10) !!1:vlen and 1:MaxRPMIter
   complex(kind=CmplxKind), intent(inout) :: fots(1:ivsiz,1:10)
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
   end subroutine broy_sav_c

