!!xDx!! This is to test out RPM Mixing Scheme for complex valued functions.
!!xDx!! We were able to get RPM Mixing to work in matlab, we now want to try
!!xDx!! out RPm Mixing on fortran.
!!xDx!!========================================================================
!!xDx!!         The program no. is controlled by the variable 'typ'
!!xDx!!========================================================================
!!xDx!!1. In order to be sure its working fine, we test it against a Function
!!xDx!!   We want to solve equations of the form F[u]=u
!!xDx!!   Besure to chekc out the mathematica file, Broyden_Test.nb in my laptop.
!!xDx!!   1. F[u,v,w,s]= {u^2-2,v^2-2,w^2-2,s^2-2}  
!!xDx!!       Doesnt convergenge using Fixed point  
!!xDx!!   2. F[u,v,w,s]= {v^2-8,u+5,sqrt(w+2),sqrt(s+2)}
!!xDx!!       last two have fixed points first 2 don't  
!!xDx!!   2. F[u,v,w,s]= {uv^2+5w,v^3+u^2-8w,uw-v,sqrt(s+2)}
!!xDx!!       last has fixed points first 3 don't.  
!!xDx!!   3. F[u,v,w,s,t,x,y,z]= {uv^2+5w,v^3+u^2-8w,uw-v,sqrt(s+2), sqrt[t + 4], Sqrt[y], Sqrt[1 - x^2], Sqrt[z + 3v]}
!!xDx!!       last 6 have fixed points first 2 don't. Technically the 3rd should also not have a fixed point but I guess
!!xDx!!       it has for the initial conditions chosen in MEqn1 in RPM_Test.nb 
!!xDx!!   4. F[u,v,w,s,t,x,y,z]= {I*uv^2+5w,v^3+u^2-8w,uw-v,sqrt(s+2), sqrt[t + 4], Sqrt[y], Sqrt[1 - x^2], Sqrt[z + 3v]}
!!xDx!!       last 6 have fixed points first 2 don't. Technically the 3rd should also not have a fixed point but I guess
!!xDx!!       it has for the initial conditions chosen in MEqn2 in RPM_Test.nb 
program RPM
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

#ifdef USE_SCALAPACK

   interface 
      function BLACS_PNUM(ic,j,k) result(i)
         integer :: ic, j, k
         integer :: i
      end function BLACS_PNUM
   end interface
!
   interface 
      function NUMROC(ic,j,k,l,m) result(i)
         integer :: ic, j, k, l, m
         integer :: i
      end function NUMROC
   end interface
!
   interface 
      function indxg2p(ic,j,k,l,m) result(i)
         integer, intent(in) :: ic, j, k, l, m
         integer :: i
      end function indxg2p
   end interface
!
   interface 
      function indxl2g(ic,j,k,l,m) result(i)
         integer, intent(in) :: ic, j, k, l, m
         integer :: i
      end function indxl2g
   end interface
!
   interface 
      subroutine infog1l(ig,j,k,l,m,n,o)
         integer, intent(in) :: ig, j, k, l, m
         integer, intent(out):: n, o 
      end subroutine infog1l
   end interface
!
   interface 
      subroutine infog2l(igr,igc,desc,nprow,npcol,r,c,lrind,lcind,rsrc,csrc)
         integer, intent(in) :: igr, igc, nprow, npcol, r, c
         integer, intent(in) :: desc(9)
         integer, intent(out):: lrind, lcind, rsrc, csrc
      end subroutine infog2l
   end interface

#endif


  !!xDx!! Looping/Choice variables
   integer (kind=IntKind) :: ii, jj, j, k, i, typ=14   
   integer(kind=IntKind) :: idt, id, idq, vlenb4, vlen, iter, iteration   
   integer(kind=IntKind) :: RPMIter=0 !! Used to increment Basis Size    
   integer(kind=IntKind) :: iter_broy=0 !! iter for Broy part of RPM    
  !!xDx!! Parameters    
  integer (kind=IntKind) :: Basis, MaxRPMIter, TotalRPMIter   
  integer (kind=IntKind) :: MaxBasis, NumAtoms, NumTotalMix   
  real(kind=RealKind) :: tol    
  !!xDx!! Iter variables    
   integer (kind=IntKind) :: lastIter, nn, lastm1, NumRPMIter_broy, NumBroydenIter    
  !!xDx!! Mixing Parameters     
   real(kind=RealKind) :: alpha, beta    
  !!xDx!! Used as pointer to original variable and pointer to unstable subspace variable.    
   complex(kind=CmplxKind), allocatable :: pvect_old(:), pvect_new(:)    
   complex(kind=CmplxKind), allocatable :: pdf(:), pf(:)     
   complex(kind=CmplxKind), allocatable :: pvold(:), pv_st_old(:), A(:,:)     
   complex(kind=CmplxKind), allocatable :: pu(:,:), pw(:), pvt(:,:)    
   complex(kind=CmplxKind), allocatable :: pn(:), po(:)     
   real(kind=RealKind) :: normp   
!!xDx!! Projection variables    
   complex(kind=CmplxKind), allocatable :: G(:,:), R(:,:), G_g_c(:,:)      
   complex (kind=CmplxKind) :: Ri11, Ri22
!!xDx!! Subspace variables    
   complex(kind=CmplxKind), allocatable :: q_old(:), q_new(:)      
   complex(kind=CmplxKind), allocatable :: p_old(:), p_new(:)    
   complex(kind=CmplxKind), allocatable :: z_old(:), z_new(:)     
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
   integer (kind=Intkind) :: new_comm, old_comm, old_group, new_group
!                   
   complex(kind=CmplxKind), allocatable:: a_c(:,:), b_c(:,:), cm_c(:)        
!
  integer(kind=IntKind) :: GroupID     
      
  real (kind=RealKind) :: rms      
!          
  character (len=50) :: PEStr     
  character (len=19) fmt !!xDx!! 19 is no. of char in (F7.4,"+",F7.4,"i")

#ifdef USE_SCALAPACK
!  ===================================================================
!  *****      ScaLAPACK parameters for communication
!  ===================================================================
   integer (kind=IntKind), parameter :: DLEN_ = 9
   integer (kind=IntKind) :: ICTXT
   integer (kind=IntKind) :: NPROW
   integer (kind=IntKind) :: NPCOL
   integer (kind=IntKind) :: MYROW
   integer (kind=IntKind) :: MYCOL

!!rDr!!   integer (kind=IntKind) :: DESC_Z( DLEN_ )
   integer (kind=IntKind) :: DESC_G( DLEN_ )
   integer (kind=IntKind) :: DESC_R( DLEN_ )

!!xDx!! ScaLAPACK parameters for local Matrix definition
   integer (kind=IntKind) :: M_Z, N_Z, MB_Z, NB_Z !!nDn!! Check in detail, especially M.
   integer (kind=IntKind) :: M_G, N_G, MB_G, NB_G !!nDn!! Check in detail, especially M.
   integer (kind=IntKind) :: M_R, N_R, MB_R, NB_R !!nDn!! Check in detail, especially M.
   integer (kind=IntKind) :: CSRC_Z=0, RSRC_Z=0, IZ=1, JZ=1
   integer (kind=IntKind) :: CSRC_G=0, RSRC_G=0, IG=1, JG=1
   integer (kind=IntKind) :: CSRC_R=0, RSRC_R=0, IR=1, JR=1
   integer (kind=IntKind) :: ROCSRC, COCSRC, LindR, LindC
#endif
!!xDx!! PArameters for ScaLAPACK & LAPACK routines
   integer (kind=IntKind) :: LWork
   real (kind=RealKind) :: tmp 
   complex (kind=cmplxKind), allocatable :: TAU(:), work(:)
   integer(kind=IntKind), allocatable :: ipiv(:)
!!xDx!! Infor Parameters 
   integer (kind=IntKind) :: INFO_G
   integer (kind=IntKind) :: INFO_R
!  ===================================================================
!!xDx!! Defining the required parallel variable.
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup
!          
           
!!xDx!! Testing assigning & printing variables        
  complex(kind=CmplxKind) :: test       
!!xDx!! Testing assigning & printing variables      
        
!! MPI initialization          
   call initMPP()        
          
   call initGroupComm()       
!      
!! The Equation to solve F(u)=u        
if (typ == 11) then  
   vlenb4=4; Basis=0; MaxRPMIter=6; NumBroydenIter=10; alpha= 0.1; beta=1-alpha; 
   tol = 1D-09; TotalRPMIter=40; 
elseif (typ == 12) then
   vlenb4=2; Basis=0; MaxRPMIter=8; NumBroydenIter=10; alpha= 1D-01; beta=1-alpha;
   tol = 1D-09; TotalRPMIter=100; 
elseif (typ == 13) then
   vlenb4=2; Basis=0; MaxRPMIter=8; NumBroydenIter=10; alpha= 1D-01; beta=1-alpha;
   tol = 1D-09; TotalRPMIter=100;
elseif (typ == 13) then
   vlenb4=2; Basis=0; MaxRPMIter=8; NumBroydenIter=10; alpha= 1D-01; beta=1-alpha;
   tol = 1D-09; TotalRPMIter=100;
elseif (typ == 14) then
   vlenb4=2; Basis=0; MaxRPMIter=8; NumBroydenIter=10; alpha= 1D-01; beta=1-alpha;
   tol = 1D-09; TotalRPMIter=100;
else
print *, 'Please eneter a vaild typ'
endif

   NumAtoms=NumPEs
   NumTotalMix=1 
   vlen=NumTotalMix*vlenb4
   MaxBasis=4!NumAtoms*vlen
   NumRPMIter_broy=MaxRPMIter-2

!!xDx!!Allocation arrays
allocate(pn(NumAtoms*vlenb4))
allocate(po(NumAtoms*vlenb4))
allocate(pvold(MaxBasis))
allocate(pvect_old(vlenb4))
allocate(pvect_new(vlenb4))
allocate(pv_st_old(vlenb4))
allocate(pdf(MaxBasis))
allocate(pf(MaxBasis))
allocate(pu(MaxBasis,NumRPMIter_broy))
allocate(pw(MaxRPMIter))
allocate(pvt(MaxBasis,NumRPMIter_broy))
allocate(A(NumTotalMix*vlenb4,2))

   write(PEStr, '(I2)') MyPE+10

!!xDx!! Allocating input -----------------------------------------------------------------------
   po=CZERO
   pn=CZERO

if  (typ .eq. 11) then
    pn(1) = cmplx(-1.0370D0,-0.0027D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(2) = cmplx( 1.2975D0,-0.773364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(3) = cmplx(-1.0370D0,-0.0027D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(4) = cmplx( 1.2975D0,-0.773364D0)!complex(1.2975000143051100,-0.7733640074729920)
elseif  (typ .eq. 12) then    
    pn(1) = cmplx(-4.30278D0, 0.0007D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(2) = cmplx(-3.30278D0, 0.073364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(3) = cmplx( 3.0070D0, 0.0027D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(4) = cmplx( 0.9975D0, 0.073364D0)!complex(1.2975000143051100,-0.7733640074729920)
elseif  (typ .eq. 13) then    
    pn(1) = cmplx(-1.36278D0, 0.1807D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(2) = cmplx( 0.33278D0,-1.973364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(3) = cmplx(-0.1970D0, 0.5927D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(4) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(5) = cmplx( 2.36278D0, 0.0807D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(6) = cmplx( 0.60278D0,-0.973364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(7) = cmplx( 0.5970D0, 0.0927D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(8) = cmplx( 2.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
elseif  (typ .eq. 14) then    
    pn(1) = cmplx(-1.36278D0, 0.1807D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(2) = cmplx( 0.33278D0,-1.973364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(3) = cmplx(-0.1970D0, 0.5927D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(4) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(5) = cmplx( 2.36278D0, 0.0807D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(6) = cmplx( 0.60278D0,-0.973364D0)!complex(1.2975000143051100,-0.7733640074729920)
    pn(7) = cmplx( 0.5970D0, 0.0927D0)!complex(-1.037099957466130,-2.700000070035458e-003)
    pn(8) = cmplx( 2.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
endif

!!xDx!! Declaring Global variables
M_Z =NumAtoms*vlen  !!?D?!! What the dickens happened to the L values.
N_Z =MaxBasis
MB_Z=NumTotalMix*vlen !! Assuming NumTotalMix is the Local no. of atoms
NB_Z=MaxBasis
M_G =N_Z
N_G =M_Z
MB_G=NB_Z
NB_G=MB_Z
M_R =2 
N_R =NumAtoms*vlen !!?D?!! What the dickens happened to the L values.
MB_R=2
NB_R=NumTotalMix*vlen !! Assuming NumTotalMix is the Local no. of atoms
!q_old; q_new; p_old; p_new; z_old; z_new

#ifdef USE_SCALAPACK
   call initGroupComm()
   dim(1) = NumPEs
   call createProcGrid(1,'1-d Grid',dim)
   grid = getProcGridID('1-d Grid')
   box(1)=NumPEs
   call createBoxGroupFromGrid(grid,box,'Full')
   GroupID = getGroupID('Full')
   new_comm = getGroupCommunicator(GroupID)
   old_comm = getCommunicator()
   call MPI_COMM_GROUP(old_comm,old_group,info)
   call MPI_COMM_GROUP(new_comm,new_group,info)
   call syncAllPEs()
!
   ICTXT = new_comm
!
   call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEs )
   call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!  -------------------------------------------------------------------
   call DESCINIT( DESC_G, M_G, N_G, MB_G, NB_G, RSRC_G, CSRC_G, ICTXT, max( 1, NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW) ) , INFO_G )
   call DESCINIT( DESC_R, M_R, N_R, MB_R, NB_R, RSRC_R, CSRC_R, ICTXT, max( 1, NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW) ) , INFO_R )
!!xDx!! Allocating Matrices related to projection
   allocate ( G( NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW), NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL) ) )
   allocate ( R( NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW), NUMROC(N_R,NB_R,MYCOL,CSRC_R,NPCOL) ) )
#endif

ITERA:   do RPMIter = 1,TotalRPMIter!!1a.!!1,MaxRPMIter, is for the rest.

   !defining F in u=F(u) 
   if( typ .eq. 11 ) then
     po = pn;
     pn(1)=po(1)**2-2
     pn(2)=po(2)**2-2
     pn(3)=po(3)**2-2
     pn(4)=po(4)**2-2
   elseif( typ .eq. 12 ) then 
     po = pn
     pn(1)=po(2)**2-8
     pn(2)=po(1)+5
     pn(3)=sqrt(po(3)+2)
     pn(4)=sqrt(po(4)+2)
   elseif( typ .eq. 13 ) then 
     po = pn
     pn(1)=po(2)*(po(1)**2)+5*po(3)
     pn(2)=(po(2)**3)+(po(1)**2)-8*po(3)
     pn(3)=po(3)*po(1)-po(2)
     pn(4)=sqrt(po(4)+2)
     pn(5)=sqrt(po(5)+4)
     pn(6)=sqrt(po(7))
     pn(7)=sqrt(1-po(6)**2)
     pn(8)=sqrt(po(8)+3)
   elseif( typ .eq. 14 ) then 
     po = pn
     pn(1)=cmplx(0,1)*po(2)*(po(1)**2)+5*po(3)
     pn(2)=(po(2)**3)+(po(1)**2)-8*po(3)
     pn(3)=po(3)*po(1)-po(2)
     pn(4)=sqrt(po(4)+2)
     pn(5)=sqrt(po(5)+4)
     pn(6)=sqrt(po(7))
     pn(7)=sqrt(1-po(6)**2)
     pn(8)=sqrt(po(8)+3)
   endif     
   ! distributing pvect to the different PEs so that we can iterate.
      pvect_new(1:vlen) = pn(MyPE*vlen+1:(MyPE+1)*vlen)
      pvect_old(1:vlen) = po(MyPE*vlen+1:(MyPE+1)*vlen)

   if (MyPE .eq. 0) then
      write (MyPE+10,1122) '========================'
      write (MyPE+10,1342) 'MyPE','pn(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
      write (MyPE+10,1122) '========================'
      do i =1,NumAtoms*vlenb4
         write (MyPE+10,1241) 'i,MyPE','pn',i,MyPE,real(pn(i)),aimag(pn(i)),'*i'
      enddo
      write (MyPE+10,1122) '========================'
      write (MyPE+10,1342) 'MyPE','po(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
      write (MyPE+10,1122) '========================'
      do i =1,NumAtoms*vlenb4
         write (MyPE+10,1241) 'i,MyPE','po',i,MyPE,real(po(i)),aimag(po(i)),'*i'
      enddo
1122 format(X,A)
1342 format(2(X,A),3(X,I4))
1241 format(2(X,A),2(X,I4),(X,F20.10,SP,F16.10,A))
   endif


!!xDx!! Allocating local variables                                                           
   allocate ( q_old( NumTotalMix*vlen ) )
   allocate ( q_new( NumTotalMix*vlen ) )
   allocate ( p_old( NumTotalMix*vlen ) )
   allocate ( p_new( NumTotalMix*vlen ) )
   allocate ( z_old( MaxBasis ) )
   allocate ( z_new( MaxBasis ) )
   if ( RPMIter .eq. 1) then
      allocate(G_g_c(MaxBasis,NumAtoms*vlen))
      G_g_c=CZERO
   endif
  q_old = CZERO
  q_new = CZERO
  p_old = CZERO
  p_new = CZERO
  z_old = CZERO
  z_new = CZERO
  R = CZERO
  G = CZERO

#ifdef USE_SCALAPACK
!!xDx!! initializing G from G_g_c
   do ii=1,M_G
   do jj=1,N_G
   call infog2l ( ii, jj, DESC_G, NPROW, NPCOL, MYROW, MYCOL, LindR, LindC, ROCSRC, COCSRC) ! This gives us loc coord for global G
   if ( MYROW .eq. ROCSRC .and. MYCOL .eq. COCSRC) then
   G(LindR,LindC) = G_g_c(ii,jj)
   end if
   enddo
   enddo
#endif

   p_old = pvect_old
   p_new = pvect_new

   iter = RPMIter

   if (RPMIter .ge. MaxRPMIter)then 
       iter_broy= mod( RPMIter, MaxRPMIter )
   endif

!!xDx!! Defining corresponding corresponding subspace variables which will be passed -------------------------------------------------------------------------
            allocate(msgbuf1(1:MaxBasis))
            msgbuf1=matmul(G,p_old)!! z<-Gu
            call GlobalSumInGroup(GroupID, msgbuf1, MaxBasis)
            z_old=msgbuf1
            deallocate(msgbuf1)
            allocate(msgbuf1(1:MaxBasis))
            msgbuf1=matmul(G,p_new)!! z<-Gu
            call GlobalSumInGroup(GroupID, msgbuf1, MaxBasis)
            z_new=msgbuf1                                            
            deallocate(msgbuf1)
            p_old=matmul(transpose(conjg(G)),z_old)!! p<-Gz
            p_new=matmul(transpose(conjg(G)),z_new)!! p<-Gz

      !------ Defining the stable subspace projection
            q_old = pvect_old-p_old  ! q<- u-p
            q_new = pvect_new-p_new  ! q<- u-p
!!cDc!!,
!!cDc!!,!!cDc!!
!!cDc!!,	do k=0,NumPEs-1
!!cDc!!,	   if ( MyPE .eq. k) then
!!cDc!!,	      do ii=1, NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW)
!!cDc!!,	      do jj=1, NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL)
!!cDc!!,	         print *, 'k, G,ii,jj',k, G(ii,jj),ii,jj
!!cDc!!,	      enddo
!!cDc!!,	      enddo
!!cDc!!,	   endif
!!cDc!!,	enddo
!!cDc!!,!!cDc!!
!!cDc!!,
      !------ Defining the unstable subspace projection as pvect
            !pvect_old=CZERO                                               
            !pvect_new=CZERO                                               
            !pvect_old(1:MaxBasis)=z_old( 1:MaxBasis ) !!iDi!!
            !pvect_new(1:MaxBasis)=z_new( 1:MaxBasis ) !!iDi!!
!!cDc!!     
!!cDc!! do k=0,NumPEs-1 
!!cDc!!    if (MyPE .eq. k) then
!!cDc!!       write (MyPE+10,1122) '========================'
!!cDc!!       write (MyPE+10,1342) '2 MyPE','q_new(:),RPMIter',MyPE,RPMIter
!!cDc!!       write (MyPE+10,1122) '========================'
!!cDc!!       do i =1,vlen
!!cDc!!          write (MyPE+10,1241) 'i,MyPE','q_new',i,MyPE,real(q_new(i)),aimag(q_new(i)),'*i'
!!cDc!!       enddo
!!cDc!!       write (MyPE+10,1122) '========================'
!!cDc!!       write (MyPE+10,1342) '2 MyPE','q_old(:),RPMIter',MyPE,RPMIter
!!cDc!!       write (MyPE+10,1122) '========================'
!!cDc!!       do i =1,vlen
!!cDc!!          write (MyPE+10,1241) 'i,MyPE','q_old',i,MyPE,real(q_old(i)),aimag(q_old(i)),'*i'
!!cDc!!       enddo
!!cDc!!    endif
!!cDc!! enddo
!!cDc!! 
!! Performing Simple Mixing ---------------------------------------------------------------------
            if ( alpha == 1 ) then
            q_new = alpha*q_new
            else
            q_new = alpha*q_new+(1-alpha)*q_old;
            endif

!!cDc!!,!!cDc!!
!!cDc!!,	print *, 'Hello, RPMITER, iter_broy, Basis', RPMIter, iter_broy, Basis
!!cDc!!,!!cDc!!

!! Performing Broyden Mixing -------------------------------------------------------------------------

   if ( iter_broy .eq. 1) then

!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
            pvold = CZERO
            pu    = CZERO
            pvt   = CZERO
            pdf   = CZERO

            pf = z_new - z_old 
            pvold = z_old 
            z_new = z_old + alpha * pf 
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
!!cDc!!
!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (*,1122) '========================'
!!cDc!!      write (*,1042) 'nn', 'iter_broy','lastIter','NumRPMIter_broy',  nn, iter_broy,  lastIter, NumRPMIter_broy
!!cDc!!      write (*,1122) '========================'
!!cDc!!1042 format(4(X,A),4(X,I4))
!!cDc!!   endif
!!cDc!!
!

   if (MyPE .eq. 0) then
      write (MyPE+10,1122) '========================'
      write (MyPE+10,1342) 'MyPE','z_new(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
      write (MyPE+10,1122) '========================'
      do i =1,MaxBasis
         write (MyPE+10,1241) 'i,MyPE','z_new',i,MyPE,real(z_new(i)),aimag(z_new(i)),'*i'
      enddo
      write (MyPE+10,1122) '========================'
      write (MyPE+10,1342) 'MyPE','z_old(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
      write (MyPE+10,1122) '========================'
      do i =1,MaxBasis
         write (MyPE+10,1241) 'i,MyPE','z_old',i,MyPE,real(z_old(i)),aimag(z_old(i)),'*i'
      enddo
   endif

         w0= ten2m2
            dfnorm = zero
            fnorm = zero
            if ( MaxBasis==1 ) then
               z_new(1) = z_old(1) + alpha * pf(1) !!xDx!! ??? why ???
            endif
            do k = 1,MaxBasis !!cDt!!x2
               pdf(k) = z_new(k) - z_old(k) - pf(k)!! dR^(iter)=R^(iter)-R^(iter-1)
               pf(k)  = z_new(k) - z_old(k) !!R^(iter)
            enddo
            do k = 1,MaxBasis  !!xDx!! These are sums over all the atoms in the proc. 
               dfnorm = dfnorm + real(pdf(k)*conjg(pdf(k)),kind=RealKind)
               fnorm  = fnorm  + real(pf(k)*conjg(pf(k)),kind=RealKind) 
            enddo
!
            dfnorm = sqrt( dfnorm )
            fnorm  = sqrt( fnorm )
!!cDc!!
!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1345) 'MyPE','pdf(:),RPMIter, dfnorm', MyPE, RPMIter, dfnorm
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','pdf',i,MyPE,real(pdf(i)),aimag(pdf(i)),'*i'
!!cDc!!      enddo
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1345) 'MyPE','pf(:),RPMIter, fnorm', MyPE, RPMIter, fnorm
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','pf',i,MyPE,real(pf(i)),aimag(pf(i)),'*i'
!!cDc!!      enddo
!!cDc!!   1345 format(2(X,A),2(X,I4),(X,F16.8))
!!cDc!!   endif
!!cDc!!
!
            fac2 = one/dfnorm
            fac1 = alpha*fac2
            !!cDc!! vlen  = MaxBasis  !!??!! RPMMix(id)%vlen
            !!cDc!! if ( vlen==1 ) then
            !!cDc!!    cycle
            !!cDc!! endif
            do k = 1,MaxBasis
               z_new(k) = fac1*pdf(k) + fac2*(z_old(k) - pvold(k)) 
               pvold(k)   = z_old(k)
               z_old(k) = fac2*pdf(k)
            enddo  
!!cDc!!
!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'pu: MyPE','z_new(:), RPMIter, iter_broy', MyPE, RPMIter, iter_broy
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','z_new',i,MyPE,real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!      enddo
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'pvt: MyPE','z_old(:), RPMIter, iter_broy', MyPE, RPMIter, iter_broy
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','z_old',i,MyPE,real(z_old(i)),aimag(z_old(i)),'*i'
!!cDc!!      enddo
!!cDc!!   endif
!!cDc!!

!           ============================================================
            call broy_sav_c( pu, pvt, z_old, z_new, iter_broy-1, & 
                                NumRPMIter_broy, MaxBasis )!!cDt!!

!           ============================================================
            !!xDx!! allocating the Matrices aM_c and bM_c=aM_c^-1
            allocate (a_c(1:nn,1:nn))
            allocate (cm_c(1:nn))
!!cDc!!
!!cDc!!	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!	if ( MyPE .eq. 0) then
!!cDc!!	write (mype+10,1122) '                        '
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	write (mype+10,1342) 'mype','pu(:,:),rpmiter, iter_broy', mype, rpmiter, iter_broy
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	do ii=1,maxbasis
!!cDc!!	   do jj=1,NumRPMIter_broy
!!cDc!!              fmt(8:8) = merge('+',' ',imag(pu(ii,jj)).gt.0)
!!cDc!!              write(mype+10,fmt, advance="no") pu(ii,jj)
!!cDc!!	   enddo
!!cDc!!	   write(mype+10, fmt="(a)") " "
!!cDc!!	enddo
!!cDc!!	write (MyPE+10,1122) '                        '
!!cDc!!	!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	write (mype+10,1342) 'mype','pvt(:,:),rpmiter, iter_broy', mype, rpmiter, iter_broy
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	do ii=1,maxbasis
!!cDc!!	   do jj=1,NumRPMIter_broy
!!cDc!!              fmt(8:8) = merge('+',' ',imag(pvt(ii,jj)).gt.0)
!!cDc!!              write(mype+10,fmt, advance="no") pvt(ii,jj)
!!cDc!!	   enddo
!!cDc!!	   write(mype+10, fmt="(a)") " "
!!cDc!!	enddo
!!cDc!!	write (mype+10,1122) '                        '
!!cDc!!	endif
!!cDc!!

!!xDx!! nn -> #iters-1, equiv to m in DD Johnson.
            do j = 1,nn - 1          ! off diagonal elements of a(i,j)
               do i = j+1,nn
                  aij = czero
                     do k = 1,MaxBasis
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
                  do k = 1,MaxBasis!!cDt!!x2
                     cmj = cmj + conjg(pvt(k,i))*pf(k) 
                     aij = aij + conjg(pvt(k,i))*pvt(k,i)
                  enddo
               a_c(i,i) = aij
               cm_c(i) = cmj
            enddo
!!cDc!!
!!cDc!!	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!	if ( MyPE .eq. 0) then
!!cDc!!	write (mype+10,1122) '                        '
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	write (mype+10,1342) 'mype','a_c(:,:),rpmiter, iter_broy', mype, rpmiter, iter_broy
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	do ii=1,nn
!!cDc!!	   do jj=1,nn
!!cDc!!              fmt(8:8) = merge('+',' ',imag(a_c(ii,jj)).gt.0)
!!cDc!!              write(mype+10,fmt, advance="no") a_c(ii,jj)
!!cDc!!	   enddo
!!cDc!!	   write(mype+10, fmt="(a)") " "
!!cDc!!	enddo
!!cDc!!	write (MyPE+10,1122) '                        '
!!cDc!!	!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!!cDc!!	write (MyPE+10,1122) '========================'
!!cDc!!	write (MyPE+10,1342) '2 MyPE','cm_c(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
!!cDc!!	write (MyPE+10,1122) '========================'
!!cDc!!	do i =1,nn
!!cDc!!	   write (MyPE+10,1241) 'i,MyPE','cm_c',i,MyPE,real(cm_c(i)),aimag(cm_c(i)),'*i'
!!cDc!!	enddo
!!cDc!!	endif
!!cDc!!

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
               if (  iter_broy > NumRPMIter_broy ) then
                  pw(NumRPMIter_broy) = wtmp
               else
                  pw(lastIter) = wtmp       !w(lastm1)=wtmp
               endif
!!cDc!!
!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1002) '2 pw',' dfnorm','MyPE','NumRPMIter_broy', 'iter_broy',dfnorm, MyPE, NumRPMIter_broy, iter_broy
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!	do i =1,nn
!!cDc!!	   write (MyPE+10,1241) 'i,MyPE','pw',i,MyPE,real(pw(i)),aimag(pw(i)),'*i'
!!cDc!!	enddo
!!cDc!!1002 format(5(X,A),(X,F10.8),3(X,I4))
!!cDc!!   endif
!!cDc!!
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
!!cDc!!
!!cDc!!	fmt = '(F9.1,"+",F9.1,"i")'
!!cDc!!	if ( MyPE .eq. 0) then
!!cDc!!	write (mype+10,1122) '                        '
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	write (mype+10,1342) ' b4 inv mype','b_c(:,:),rpmiter, iter_broy', mype, rpmiter, iter_broy
!!cDc!!	write (mype+10,1122) '========================'
!!cDc!!	do ii=1,nn
!!cDc!!	   do jj=1,nn
!!cDc!!              fmt(8:8) = merge('+',' ',imag(b_c(ii,jj)).gt.0)
!!cDc!!              write(mype+10,fmt, advance="no") b_c(ii,jj)
!!cDc!!	   enddo
!!cDc!!	   write(mype+10, fmt="(a)") " "
!!cDc!!	enddo
!!cDc!!	write (MyPE+10,1122) '                        '
!!cDc!!	endif
!!cDc!!
!
               if ( .not.allocated(ipiv) ) then
                  allocate( ipiv(nn) )
               endif
!              ------------------------------------------------------
               call zgetrf( nn, nn, b_c, nn,ipiv, info ) !!cDt!!
!              ------------------------------------------------------
               call zgetri( nn, b_c, nn, ipiv, tmp, -1, info )!!cDt!!
!              ------------------------------------------------------
               LWork = int(real(tmp,kind=RealKind))
               allocate( Work(1:LWork) )
!              ------------------------------------------------------
               call zgetri( nn, b_c, nn, ipiv, Work, LWork, info )!!cDt!!
!              ------------------------------------------------------
               deallocate (Work)
!              ------------------------------------------------------
               deallocate (ipiv)
!              ------------------------------------------------------
!

               do k = 1,MaxBasis
                  z_new(k) = pvold(k) + alpha*pf(k)!!cDt!!
               enddo

               do i = 1,nn
                  gmi = zero
                  do j = 1,nn
                     gmi = gmi + cm_c(j)*b_c(j,i)*pw(j) 
                  enddo
                  do k = 1,MaxBasis
                     z_new(k) = z_new(k) - gmi*pu(k,i)*pw(i)!!cDt!!
                  enddo
               enddo

               deallocate( b_c )
               deallocate( a_c )
               deallocate( cm_c )

   endif   
!!xDx!! End of Broyden Mixing -------------------------------------------------------------------------------------------------------------------------------

         !!xDx!! projecting local back to global variables
            pvect_new(1:vlen)= matmul( Transpose(conjg( G(1:MaxBasis,1:vlen) ) ), z_new(1:MaxBasis) ) 
            pvect_old(1:vlen)= matmul( Transpose(conjg( G(1:MaxBasis,1:vlen) ) ), z_old(1:MaxBasis) ) 
         !!xDx!! Combining the resultant to return to main ------------------------------------------------------------------------------------------------------
               pvect_new=q_new+pvect_new!! u<- q+p

!
!!cDc!!
!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'after MyPE','pvect_new(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!      enddo
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'after MyPE','pvect_old(:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','pvect_old',i,MyPE,real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!cDc!!      enddo
!!cDc!!   endif
!!cDc!!

!!xDx!! Increasing size of basis --------------------------------------------------------------------------------------
            if ( RPMIter .eq. 1) then
               pv_st_old(:) = q_old(:)
            elseif ( RPMIter .eq. 2) then
               A(:,1)=q_old(:)-pv_st_old(:) 
               pv_st_old(:)=q_old(:)
            else
               A(:,2)=A(:,1)
               A(:,1)=q_old(:)-pv_st_old(:) 
               pv_st_old(:)=q_old(:)
            endif

            if ( mod( RPMIter, MaxRPMIter ) .eq. 0 ) then !! Inc basis size.
               R = transpose(conjg(A))
            endif
!!cDc!!,
!!cDc!!,!!cDc!!
!!cDc!!,	do k=0,NumPEs-1
!!cDc!!,	   if ( MyPE .eq. k) then
!!cDc!!,	   do jj=1,NumTotalMix*vlen
!!cDc!!,	      do ii=1,2
!!cDc!!,	         print *, 'k, A(ii,jj),ii,jj',k, A(ii,jj), ii, jj
!!cDc!!,	      enddo
!!cDc!!,	   enddo
!!cDc!!,	   endif
!!cDc!!,	enddo
!!cDc!!,!!cDc!!


	fmt = '(F9.6,"+",F9.6,"i")'
	do k=0,NumPEs-1
	if ( MyPE .eq. k .and. mod( RPMIter, MaxRPMIter ) .eq. 0) then
	write (MyPE+10,1122) '                        '
	write (MyPE+10,1122) '========================'
	write (MyPE+10,1342) ' MyPE','R(:,:),RPMIter, iter_broy', MyPE, RPMIter, iter_broy
	write (MyPE+10,1122) '========================'
	do ii=1,NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW)
	   do jj=1,NUMROC(N_R,NB_R,MYCOL,CSRC_R,NPCOL)
              fmt(8:8) = MERGE('+',' ',imag(R(ii,jj)).gt.0)
              write(MyPE+10,fmt, advance="no") R(ii,jj)
	   enddo
	   write(MyPE+10, fmt="(a)") " "
	enddo
	write (MyPE+10,1122) '                        '
	endif
	enddo


!                                                   
!!xDx!! Performing the calc to increase the Basis size                                                                                                          
      if ( mod( RPMIter, MaxRPMIter ) .eq. 0 ) then !!nDn!! iter_broy= 0 & RPMIter> MaxRPMIter

      !!xDx!! HouseHolder LQ---------------------------------------------------------------- 
#ifdef USE_SCALAPACK
            !------- Allocating Tau -------------------------------------
                  allocate( TAU( 1:NUMROC(IR+min(M_R,N_R)-1, MB_R , MYROW, 0, NPROW) ) )

            !----------------------- Work Place Qwery ----------------------------
                  call PZGELQF( M_R , N_R, R, IR, JR, DESC_R, TAU, tmp, -1, INFO_R )
            !-------------------------------------------------------------------
            LWork = int(tmp)
            allocate( WORK(1:LWork) )
           
            !  ------------------------- Actual QR  ------------------------------
                  call PZGELQF( M_R, N_R, R, IR, JR, DESC_R, TAU, WORK, LWork, INFO_R )
            !  -------------------------------------------------------------------
           
            deallocate (WORK)

       !!xDx!! Storing R_{11} and R_{22}
            Ri11=CZERO
            Ri22=CZERO
       !!nDn!! Ri11 and Ri22 could have non-zero values on different procs.
            call infog2l ( 1, 1, DESC_R, NPROW, NPCOL, MYROW, MYCOL, LindR, LindC, ROCSRC, COCSRC )
            if ( MYROW .eq. ROCSRC .and. MYCOL .eq. COCSRC ) then
            Ri11=R(LindR,LindC)
            endif
            call infog2l ( 2, 2, DESC_R, NPROW, NPCOL, MYROW, MYCOL, LindR, LindC, ROCSRC, COCSRC )
            if ( MYROW .eq. ROCSRC .and. MYCOL .eq. COCSRC ) then
            Ri22=R(LindR,LindC)
            endif
       !!xDx!! Distributing the values
            msgbuf(1) = Ri11 
            msgbuf(2) = Ri22
            call GlobalSumInGroup(GroupID, msgbuf,2)
            Ri11 = msgbuf(1)
            Ri22 = msgbuf(2)
            !  ------------------------- Work Place Qwery  -----------------------
                  call PZUNGLQ( M_R, N_R, M_R, R, IR, JR, DESC_R, TAU, tmp, -1, INFO_R )
            !  -------------------------------------------------------------------
            LWork =int(tmp)
            allocate( WORK(1:LWork) )
            !  ------------------------- Actual Q  -------------------------------
                  call PZUNGLQ( M_R, N_R, M_R, R, IR, JR, DESC_R, TAU, WORK, LWork, INFO_R )
            !  -------------------------------------------------------------------
               deallocate (WORK)
               deallocate (TAU)
!!nDn!! write the LAPACK version.
#endif
	fmt = '(F9.4,"+",F9.4,"i")'
	if ( MyPE .eq. 0 ) then
	write (MyPE+10,1122) '                        '
	fmt(8:8) = MERGE('+',' ',imag(Ri11).gt.0)
	write (MyPE+10,1122) '========================'
	write (MyPE+10,1122) 'MyPE','Ri11'
	write (MyPE+10,1122) '========================'
	write (MyPE+10,fmt) Ri11
	write (MyPE+10,1122) '                        '
	fmt(8:8) = MERGE('+',' ',imag(Ri22).gt.0)
	write (MyPE+10,1122) '========================'
	write (MyPE+10,1122) 'MyPE','Ri22'
	write (MyPE+10,1122) '========================'
	write (MyPE+10,fmt) Ri22
	write (MyPE+10,1122) '                        '
	endif

!
!!xDx!! Based on the Magnitude of the eigen value, we will keep only 1 basis.--------------------------------------------------------------------------------------
         if ( abs(Ri11) .gt. abs(Ri22)*10**3 ) then
            if ( Basis .lt. MaxBasis ) then
              Basis=Basis+1
              G(Basis,:)=R(1,:) !! R is of size vlen, Z_g is of size NumAtoms*vlen, need to allocate properly.
            else
              do jj=2,MaxBasis !!nDn!! Shifting every column left by 1.
                G(jj-1,:)=G(jj,:)
              enddo
              G(MaxBasis,:)=R(1,:)
            endif
         elseif ( abs(Ri22) .gt. abs(Ri11)*10**3 ) then
            if ( Basis .lt. MaxBasis ) then
              Basis=Basis+1
              G(Basis,:)=R(2,:)
            else
              do jj=2,MaxBasis !!nDn!! Shifting every column left by 1.
                G(jj-1,:)=G(jj,:)
              enddo
              G(MaxBasis,:)=R(2,:)
            endif
         else
            if ( Basis .lt. MaxBasis-1 ) then
              Basis=Basis+2
              G(Basis-1,:) = R(1,:) 
              G(Basis  ,:) = R(2,:) 
            else
              do jj=3,MaxBasis !!nDn!! Shifting every column left by 1.
                G(jj-2,:)=G(jj,:)
              enddo
              G(MaxBasis-1,:)=R(1,:)
              G(MaxBasis,:)  =R(2,:)
            endif
         endif

!!cDc!!
!!cDc!!!	fmt = '(F10.5,"+",F10.5,"i")'
!!cDc!!	fmt = '(F7.4,"+",F7.4,"i")'
!!cDc!!	do k=0,NumPEs-1
!!cDc!!	if ( MyPE .eq. k) then
!!cDc!!	write (MyPE+10,1122) '                        '
!!cDc!!	write (MyPE+10,1122) '========================'
!!cDc!!	write (MyPE+10,1342) ' b4 reorth MyPE','G(:,:),RPMIter, Basis', MyPE, RPMIter, Basis
!!cDc!!	write (MyPE+10,1122) '========================'
!!cDc!!	do ii=1,NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW)
!!cDc!!	   do jj=1,NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL)
!!cDc!!!              fmt(11:11) = MERGE('+',' ',imag(G(ii,jj)).gt.0)
!!cDc!!              fmt(8:8) = MERGE('+',' ',imag(G(ii,jj)).gt.0)
!!cDc!!              write(MyPE+10,fmt, advance="no") G(ii,jj)
!!cDc!!	   enddo
!!cDc!!	   write(MyPE+10, fmt="(a)") " "
!!cDc!!	enddo
!!cDc!!	write (MyPE+10,1122) '                        '
!!cDc!!	endif
!!cDc!!	enddo
!!cDc!!


!!xDx!! Re-orthogonalize G based on new value of basis------------------------------------------------------------------------------------------------------------
#ifdef USE_SCALAPACK
            if ( Basis .gt. 1) then
            !!xDx!! LQ Decomposition
               !------- Allocating Tau -------------------------------------
                     allocate( TAU( 1:NUMROC(IG+min(M_G,N_G)-1, MB_G, MYROW,0,NPROW) ) )
               !----------------------- Work Place Qwery ----------------------------
                     call PZGELQF( M_G , N_G, G, IG, JG, DESC_G, TAU, tmp, -1, INFO_G )
               !-------------------------------------------------------------------
                  LWork = int(tmp)
                  allocate( WORK(1:LWork) )
               !  ------------------------- Actual LQ  ------------------------------
                     call PZGELQF( M_G , N_G, G, IG, JG, DESC_G, TAU, WORK, LWork, INFO_G )
               !  -------------------------------------------------------------------
                  deallocate (WORK)
               !  ------------------------- Work Place Qwery  -----------------------
                     call PZUNGLQ( M_G , N_G, M_G, G, IG, JG, DESC_G, TAU, tmp, -1, INFO_G )
               !  -------------------------------------------------------------------
                  LWork =int(tmp)
                  allocate( WORK(1:LWork) )
               !  ------------------------- Actual Q  -------------------------------
                     call PZUNGLQ( M_G , N_G, M_G, G, IG, JG, DESC_G, TAU, WORK, LWork, INFO_G )
               !  -------------------------------------------------------------------
               deallocate (WORK)
               deallocate (TAU)
               
                 if ( Basis .lt. MaxBasis ) then
                    G(Basis+1:MaxBasis,:) = CZERO
                 endif

            endif !!Basis .gt. 1

            !!xDx!! Re-Allocating G_g_c from G
                G_g_c=CZERO
                do ii=1,NUMROC(M_G,MB_G,MYROW,0,NPROW)
                do jj=1,NUMROC(N_G,NB_G,MYCOL,0,NPCOL)
                G_g_c( indxl2g( ii, MB_G, MYROW, 0, NPROW ) , indxl2g( jj, NB_G, MYCOL, 0, NPCOL ) ) = G(ii,jj)
                enddo
                enddo
            !!xDx!! Global sum to transfer the whole matrix
                allocate(msgbuf2(1:M_G,1:N_G))
                msgbuf2 = G_g_c
                call GlobalSumInGroup(GroupID, msgbuf2, M_G, N_G)
                G_g_c = msgbuf2
                deallocate(msgbuf2)
#endif



	fmt = '(F9.6,"+",F9.6,"i")'
!	fmt = '(F7.4,"+",F7.4,"i")'
	if ( MyPE .eq. 0) then
	write (MyPE+10,1122) '                        '
	write (MyPE+10,1122) '========================'
	write (MyPE+10,1342) 'MyPE','G_g_c(:,:),RPMIter, Basis', MyPE, RPMIter, Basis
	write (MyPE+10,1122) '========================'
	do ii=1,MaxBasis
	   do jj=1,NumAtoms*vlen
!              fmt(11:11) = MERGE('+',' ',imag(G_g_c(ii,jj)).gt.0)
              fmt(8:8) = MERGE('+',' ',imag(G_g_c(ii,jj)).gt.0)
              write(MyPE+10,fmt, advance="no") G_g_c(ii,jj)
	   enddo
	   write(MyPE+10, fmt="(a)") " "
	enddo
	write (MyPE+10,1122) '                        '
	endif



      endif !!nDn!! iter_broy= 0

!!xDx!! Allocating local variables                                                           
        deallocate (q_old)
        deallocate (q_new)
        deallocate (p_old)
        deallocate (p_new)
        deallocate (z_old)
        deallocate (z_new)
!

   
        ! Allocating initial guess to the different PEs
        pn = CZERO
        pn(MyPE*vlen+1:(MyPE+1)*vlen) = pvect_new(1:vlen)
        allocate(msgbuf1(1:NumAtoms*vlenb4))
        msgbuf1 = pn
        call GlobalSumInGroup(GroupID, msgbuf1, NumAtoms*vlenb4)
        pn = msgbuf1
        deallocate(msgbuf1)

!           ============================================================
            normp=ZERO
            do k = 1,NumAtoms*vlenb4 !!xDx!! These are sums over all the atoms in the proc. 
               normp  = normp  + real((pn(k)-po(k))*conjg(pn(k)-po(k)),kind=RealKind) 
            enddo
!!cDc!!
!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,2341) 'MyPE','normp, RPMITER',MyPE,normp, RPMIter
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,2341) 'MyPE','tol, RPMITER',MyPE,tol, RPMIter
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'MyPE','int(typ/10,kind=IntKind)',MyPE,int(typ/10,kind=IntKind)
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'MyPE','typ/10, RPMIter',MyPE,typ/10, RPMIter
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) 'MyPE','typ',MyPE,typ
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      2341 format(2(X,A),(X,I4),(X,F20.16),(X,I4))
!!cDc!!   endif
!!cDc!!

            if ( (int(typ/10,kind=IntKind) .eq. 1) .and. (normp .lt. tol) ) then
               if ( MyPE .eq. 0 ) then
                  write (*,1121) '========================'
                  write (*,1121) '!!xDx!!!xDx!!!xDx!!!xDx!!!xDx!!'
                  write (*,1121) '========================'
                  write (*,1121) 'Solution Converged with'
                  write (*,1341) 'MyPE','normp,iter',MyPE,normp,RPMIter
                  write (*,1342) 'MyPE','Basis',MyPE,Basis
                  write (*,1122) '========================'
                  write (*,1122) '========================'
                  write (*,1342) 'MyPE','pn(:),RPMIter',MyPE,RPMIter
                  write (*,1122) '========================'
                  do k = 1,NumAtoms*vlenb4 
                     write (*,5241) 'k','pn',k,real(pn(k)),aimag(pn(k)),'*i'
                  enddo
                  write (*,1121) '========================'
                  write (*,1121) '!!xDx!!!xDx!!!xDx!!!xDx!!!xDx!!'
                  write (*,1121) '========================'
               endif
               call syncAllPEs()
               EXIT ITERA
               1121 format(X,A)
               1341 format(2(X,A),(X,I4),(X,F20.10),(X,I4))
               5241 format(2(X,A),(X,I4),(X,F20.10,SP,F16.10,A))
            endif

!           ============================================================


!!cDc!!   if (MyPE .eq. 0) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) '3 MyPE','pn(:),RPMIter',MyPE,RPMIter
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,NumAtoms*vlenb4
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','pn',i,MyPE,real(pn(i)),aimag(pn(i)),'*i'
!!cDc!!      enddo
!!cDc!!   endif

!!cDc!!    
!!cDc!!do k=0,NumPEs-1 
!!cDc!!   if (MyPE .eq. k) then
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) '2 MyPE','pvect_new(:),RPMIter',MyPE,RPMIter
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,vlen
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','pvect_new',i,MyPE,real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!      enddo
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      write (MyPE+10,1342) '2 MyPE','q_old(:),RPMIter',MyPE,RPMIter
!!cDc!!      write (MyPE+10,1122) '========================'
!!cDc!!      do i =1,vlen
!!cDc!!         write (MyPE+10,1241) 'i,MyPE','q_old',i,MyPE,real(q_old(i)),aimag(q_old(i)),'*i'
!!cDc!!      enddo
!!cDc!!   endif
!!cDc!!enddo
!!cDc!!


enddo ITERA

!!Deallocation
deallocate(pvold)
deallocate(pn)
deallocate(po)
deallocate(pvect_old)
deallocate(pvect_new)
deallocate(pv_st_old)
deallocate(pdf)
deallocate(pf)
deallocate(pu)
deallocate(pw)
deallocate(pvt)
deallocate(G)
deallocate(R)
deallocate(A)

call BLACS_GRIDEXIT( ICTXT )
call endGroupComm()
call endMPP()


   stop 'Ok'
end program RPM

subroutine broy_sav_c( fins, fots, vector_old, vector_new, itscf,  & 
                          istore, ivsiz )
!  ==================================================================
!!xDx!! broy_sav_c( pu, pvt, z_old, z_new, iter-1, NumBroydenIter, MaxBasis)
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
   complex(kind=CmplxKind), intent(inout) :: fins(1:ivsiz,1:istore) !!1:vlen and 1:MaxRPMIter
   complex(kind=CmplxKind), intent(inout) :: fots(1:ivsiz,1:istore)
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

