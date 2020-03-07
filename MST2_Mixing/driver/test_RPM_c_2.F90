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

      module MixingModule
      !
         use KindParamModule, only : IntKind, RealKind, CmplxKind
         use ErrorHandlerModule, only : ErrorHandler, WarningHandler
         use MathParamModule, only : zero, one, two, three, ten2m2, ten2m3,ten2m8, &
                                     ten2m9, czero, cone, PI4, THIRD, sqrtm1
      !!xDx!! Communicator Modules
         use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MyPE, GlobalSum, &
                               getCommunicator
         use GroupCommModule, only : initGroupComm, endGroupComm
         use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
         use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
         use GroupCommModule, only : syncAllPEsInGroup
         use GroupCommModule, only : createProcGrid, createBoxGroupFromGrid, &
                                     getProcGridID, GlobalSumInGroup
      !!xDx!!
         use PublicTypeDefinitionsModule, only : MixListRealStruct, &
                                                 MixListCmplxStruct
!      !!xDx!! Mainly need NumAtoms
!         use SystemModule, only : getNumAtoms, getNumAlloyElements, getAlloyElementContent
      !
      public :: initMixing,         &
                mixValues,          &
                calRPM_c,           &
                setRPMSpace_c,      &
                setRPMMixing,       &
                getNumAtoms,        &
                setNumAtoms,        &
                endMixing
      !
      !!xDx!! Some real RPM Parameters made Public
         integer(kind=IntKind) :: MaxBasis, Basis=0, NumAtoms 
      !!xDx!! Global G matrix. = Transpose of Global Z matrix and global storage for Projection Matrix.
         real(kind=RealKind), allocatable  :: G_g(:,:)
         complex(kind=CmplxKind), allocatable  :: G_g_c(:,:)
!
         integer(kind=IntKind) :: GroupID     
      !
!      private
      !
      !  General parameters
      !
         logical :: Initialized = .false.
         logical :: InitializedWkSpace = .false.
      !
         integer(kind=IntKind) :: iter_count = 0
         integer(kind=IntKind) :: NumTotalMix
         integer(kind=IntKind) :: NumTypes
         integer(kind=IntKind), allocatable :: NumQuantities(:)
      !!xDx!!
         integer(kind=IntKind) :: lastIter = 0
      !!xDx!!
         integer(kind=IntKind) :: MaxVlen
         integer(kind=IntKind) :: MixingMethod =  3 ! -1 - No method have been set yet
                                                    !  0 - Simple
                                                    !  1 - Broyden 
                                                    !  2 - D.G. Anderson
                                                    !  3 - Recursive Projection Method (RPM)
      !
      !!xDx!! Real RPM Mixing Parameters
         integer(kind=IntKind) :: MaxNumRPMIter_broy !!#History used for Broyden Mixing, MaxNumRPMIter_broy<=MaxRPMIter-2!!
         integer(kind=IntKind) :: MaxRPMIter !! Used to increment Basis Size
         integer(kind=IntKind) :: MaxSizeRPMVect
         integer(kind=IntKind) :: RPMInvMethod = 1
      !
      integer(kind=IntKind) :: NumRPMIter_broy !! Used for Broyden Mixing
      integer(kind=IntKind) :: SizeRPMVect
      !
      type RPMStruct
      logical :: Initialized
      !
      integer(kind=IntKind) :: vlen
      !
      real(kind=RealKind) :: alpha
      !
      real(kind=RealKind) :: rms
      !
      real(kind=RealKind), pointer :: mesh(:)
      !
      real(kind=RealKind), pointer :: vector_old_r(:)
      real(kind=RealKind), pointer :: vector_new_r(:)
      !
      complex(kind=CmplxKind), pointer :: vector_old_c(:)
      complex(kind=CmplxKind), pointer :: vector_new_c(:)
      !
      real(kind=RealKind), pointer :: u_r(:,:)
      real(kind=RealKind), pointer :: vt_r(:,:)
      real(kind=RealKind), pointer :: f_r(:)
      real(kind=RealKind), pointer :: A_r(:,:)
      real(kind=RealKind), pointer :: df_r(:)
      real(kind=RealKind), pointer :: w_r(:)
      real(kind=RealKind), pointer :: vold_r(:) 
      real(kind=RealKind), pointer :: v_st_old_r(:) !!xDx!! vold_r for Stable subspace 
      !
      complex(kind=CmplxKind), pointer :: u_c(:,:)
      complex(kind=CmplxKind), pointer :: vt_c(:,:)
      complex(kind=CmplxKind), pointer :: f_c(:)
      complex(kind=CmplxKind), pointer :: A_c(:,:)
      complex(kind=CmplxKind), pointer :: df_c(:)
      complex(kind=CmplxKind), pointer :: w_c(:)
      complex(kind=CmplxKind), pointer :: vold_c(:)
      complex(kind=CmplxKind), pointer :: v_st_old_c(:) !!xDx!! vold_c for Stable subspace 
      end type RPMStruct
      !
      type (RPMStruct), allocatable :: RPMMix(:)
      !!xDx!! Not sure how to obtain NumLocalAtoms
      !!xDx!!

      !!xDx!! Defining the parameters for Parallel Linear Algebra
!      #ifdef USE_SCALAPACK
      !!xDx!! scaLAPACK Parameters
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
!      #endif
      !!xDx!!
      contains


      subroutine initMixing( nt, nq, CmplxArrayList )
      !  ===================================================================
      use GroupCommModule, only : getGroupID
      !
      implicit none
      !
      integer(kind=IntKind), intent(in) :: nt, nq(nt)
      integer(kind=IntKind) :: idt, idq
      !
      type (MixListCmplxStruct), target :: CmplxArrayList
      type (MixListCmplxStruct), pointer :: p_CAL
      !
      if ( nt<1 ) then
      call ErrorHandler('initMixing','Number of Types is less than 1')
      endif
      !
      NumTypes = nt
      allocate( NumQuantities(NumTypes) )
      NumQuantities = nq
      GroupID = getGroupID('Full')
      !
      NumTotalMix = 0
      p_CAL => CmplxArrayList
      nullify( p_CAL%next )
      !
      do idt = 1,NumTypes
      if ( NumQuantities(idt)<1 ) then
       call ErrorHandler('initMixing', &
            'Number of Quntities is less than 1 sor type',idt)
      endif
      do idq = 1,NumQuantities(idt)
       NumTotalMix = NumTotalMix +1
       if (idq==NumQuantities(idt)) cycle
       allocate(p_CAL%next)
       p_CAL => p_CAL%next
       nullify( p_CAL%next )
      enddo
      enddo
      !   nullify( p_CAL%next )
      !
      Initialized = .true.
      !
      end subroutine initMixing

      subroutine endMixing()
      !  ===================================================================
      !
      implicit none
      !
      if ( Initialized ) then
      deallocate( NumQuantities )
!      call delMixing()
      if ( MixingMethod == 3 ) then
       deallocate( RPMMix )
      endif
      Initialized = .false.
      endif
      !
      MixingMethod = -1
      !
      end subroutine endMixing

      subroutine mixValues( list_q )
      !  ===================================================================
      !
      implicit none
      !
      type(MixListCmplxStruct) :: list_q
      !
      iter_count =iter_count+1
      !
      if ( MixingMethod == 3 ) then
      call calRPM_c(list_q)
      endif
      !
      end subroutine mixValues

      subroutine  setRPMSpace_c( id, size )
!  ==================================================================
!
      implicit  none
!
      integer(kind=IntKind), intent(in) :: id, size
!
      allocate( RPMMix(id)%u_c(size,NumRPMIter_broy) )
      allocate( RPMMix(id)%vt_c(size,NumRPMIter_broy) )
      allocate( RPMMix(id)%f_c(size) )
      allocate( RPMMix(id)%A_c(2,size) )
      allocate( RPMMix(id)%df_c(size) )
      allocate( RPMMix(id)%w_c(NumRPMIter_broy) )
      allocate( RPMMix(id)%vold_c(size) )
      allocate( RPMMix(id)%v_st_old_c(size) )
      RPMMix(id)%u_c = CZERO
      RPMMix(id)%vt_c = CZERO
      RPMMix(id)%f_c = CZERO
      RPMMix(id)%A_c = CZERO
      RPMMix(id)%df_c = CZERO
      RPMMix(id)%w_c = CZERO
      RPMMix(id)%vold_c = CZERO
      RPMMix(id)%v_st_old_c = CZERO
      RPMMix(id)%vlen = size
!
      nullify( RPMMix(id)%u_r )
      nullify( RPMMix(id)%vt_r )
      nullify( RPMMix(id)%f_r )
      nullify( RPMMix(id)%A_r )
      nullify( RPMMix(id)%df_r )
      nullify( RPMMix(id)%w_r )
      nullify( RPMMix(id)%vold_r )
      nullify( RPMMix(id)%v_st_old_r )
!
      end subroutine setRPMSpace_c

      function getNumAtoms() result(NmAtms)
      !  ===================================================================
      
      implicit none
      !
      integer(kind=IntKind) ::  NmAtms
      !
      NmAtms=NumAtoms
      !
      end function getNumAtoms

      subroutine setNumAtoms( NmAtms )
      !  ===================================================================
      
      implicit none
      !
      integer(kind=IntKind), intent(in) ::  NmAtms
      !
      NumAtoms=NmAtms
      !
      end subroutine setNumAtoms

      subroutine setRPMMixing( idt, idq, amix, MaxBasis, MaxRPMIter,    &
                               MaxNumRPMIter_broy)
      !  ===================================================================
      use RadialGridModule, only : getNumRmesh
      use ScfDataModule, only    : inputRPMParam
      
      implicit none
      !
      integer(kind=IntKind), intent(in) :: idt, idq, MaxBasis
      integer(kind=IntKind), intent(in) :: MaxRPMIter 
      integer(kind=IntKind), intent(in) :: MaxNumRPMIter_broy
      !
      real(kind=RealKind), intent(in) :: amix
      !
      integer(kind=IntKind) :: id
      !
      if ( .not.Initialized ) then
      call ErrorHandler('setRPMMix',                              &
        'The MixingModule needs to be initialized first')
      endif
      !
      if ( idt<0 .or. idt>NumTypes ) then
      call ErrorHandler('setRPMMixing',"Invalid type index",idq)
      endif
      if ( idq<0 .or. idq>NumQuantities(idt) ) then
      call ErrorHandler('setRPMMixing',"Invalid quantity index",idq)
      endif
      !
      if ( .not.allocated(RPMMix) ) then
      allocate( RPMMix(NumTotalMix) )
      do id = 1, NumTotalMix
       RPMMix(id)%Initialized = .false.
      enddo
      endif  
!      call inidq!putRPMParam(MaxBasis, MaxRPMIter, MaxNumRPMIter_broy)
!      NumAtoms=getLocalNumAtoms()
      !allocate(G_g(MaxBasis,NumAtoms*1031))
      !!xDx!! Assuming same Rmesh used for all atoms.
      !!???!!   allocate(G_g(MaxBasis,NumAtoms*(1+getNumRmesh(1))))
      !!???!!   G_g=ZERO
      !
      NumRPMIter_broy = MaxNumRPMIter_broy
      !
      id = idq!getMixID(idt,idq)
      if ( .not.RPMMix(id)%Initialized ) then
      RPMMix(id)%alpha = amix
      endif
      !
      if ( MaxNumRPMIter_broy .ge. MaxRPMIter-1) then
      call WarningHandler('setRPMMixing', &
         'Max #terms used for history must be less than #RPM steps-1')
      endif 
      !
      if (RPMInvMethod<0 .or. RPMInvMethod>1) then
      call WarningHandler('initMixing', &
         'Invalid inversion method switching to default(invert1)')
      RPMInvMethod = 0
      endif 
      !
      if ( MixingMethod /= 3 ) then
      MixingMethod = 3
      endif
      !
      end subroutine setRPMMixing
      !  ===================================================================
      !
      !  *******************************************************************
      !
      !  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!
      subroutine delRPMMixing( )
      !  ===================================================================
      implicit none
      !
      integer(kind=IntKind) :: id
      !
      if ( allocated(RPMMix) ) then
      do id = 1,NumTotalMix
       if ( .not.RPMMix(id)%Initialized ) cycle
       RPMMix(id)%Initialized = .false.
       if ( associated(RPMMix(id)%u_r) ) deallocate( RPMMix(id)%u_r )
       if ( associated(RPMMix(id)%vt_r) ) deallocate( RPMMix(id)%vt_r )
       if ( associated(RPMMix(id)%f_r) ) deallocate( RPMMix(id)%f_r )
       if ( associated(RPMMix(id)%A_r) ) deallocate( RPMMix(id)%A_r )
       if ( associated(RPMMix(id)%df_r) ) deallocate( RPMMix(id)%df_r )
       if ( associated(RPMMix(id)%w_r) ) deallocate( RPMMix(id)%w_r )
       if ( associated(RPMMix(id)%vold_r) ) deallocate( RPMMix(id)%vold_r )
       if ( associated(RPMMix(id)%v_st_old_r) ) deallocate( RPMMix(id)%v_st_old_r )
       if ( associated(RPMMix(id)%u_c) ) deallocate( RPMMix(id)%u_c )
       if ( associated(RPMMix(id)%vt_c) ) deallocate( RPMMix(id)%vt_c )
       if ( associated(RPMMix(id)%f_c) ) deallocate( RPMMix(id)%f_c )
       if ( associated(RPMMix(id)%A_c) ) deallocate( RPMMix(id)%A_c )
       if ( associated(RPMMix(id)%df_c) ) deallocate( RPMMix(id)%df_c )
       if ( associated(RPMMix(id)%w_c) ) deallocate( RPMMix(id)%w_c )
       if ( associated(RPMMix(id)%vold_c) ) deallocate( RPMMix(id)%vold_c )
       if ( associated(RPMMix(id)%v_st_old_c) ) deallocate( RPMMix(id)%v_st_old_c )
      enddo
       deallocate(G_g)
       deallocate(G_g_c)
      endif
      !
      end subroutine delRPMMixing

      subroutine calRPM_c(list_q) !!xDx!! argument is the potnetial.
     !  ==================================================================
     !!nDn!! Nota Bene:- For the Broyden part we did everything wrt the
     !!nDn!! unstable subspace projections, z_old and z_new. Thus we got
     !!nDn!! rid of the idt and idq loops. Not entirely sure, but I think
     !!nDn!! that it shouldn't matter. 
     !!nDn!! Nota Bene:- We have done the boryden mixing independently in
     !!nDn!! each proc as size of MaxBasis is chosen to be small. Would 
     !!nDn!! have to rewrite the global sum in group parts. Which isn't much.
     !!xDx!!
        use printDataModule, only : printDataReal, printDataInteger,                &
                                    printDataPointer,printDataRealArray,            &
                                    printDataComplex,printDataComplexArray           
     !!xDx!!
     
     !!xDx!! Communicators for ScaLAPACK
      use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MyPE, GlobalSum, &
                            getCommunicator
      use GroupCommModule, only : initGroupComm, endGroupComm
      use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
      use GroupCommModule, only : getNumClustersOfGroup
      use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
      use GroupCommModule, only : syncAllPEsInGroup, bcastMessageInGroup
      use GroupCommModule, only : createProcGrid, createBoxGroupFromGrid, &
                                  getProcGridID, GlobalSumInGroup
!   
      use PublicTypeDefinitionsModule, only : MixListCmplxStruct
!
      use KindParamModule, only : IntKind, RealKind, CmplxKind
   
      implicit none
   
!   #ifdef USE_SCALAPACK
   
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
  
!   #endif
   
     !
      type ( MixListCmplxStruct ), target :: list_q
      integer (kind=IntKind) :: NumAtoms
     !!xDx!! Looping variables
      integer (kind=IntKind) :: ii, jj, j, k, i
      integer(kind=IntKind) :: idt, id, idq, vlen, iter, iteration
      integer(kind=IntKind) :: RPMIter=0 !! Used to increment Basis Size
      integer(kind=IntKind) :: iter_broy=0 !! iter for Broy part of RPM
     !!xDx!! Mixing Parameters 
      real(kind=RealKind) :: alpha, beta
     !!xDx!! Used as pointer to original variable and pointer to unstable subspace variable.
      complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:), pdf(:), pf(:), A(:,:) !A is for LQ Decomposititon
      complex(kind=CmplxKind), pointer :: pvold(:), pv_st_old(:),pu(:,:), pw(:), pvt(:,:) !pvold is used as a pointer for storing the previous iteration for Q and P subspace iteration.
     !!xDx!! Subspace variables
      complex(kind=CmplxKind), allocatable :: q_old(:,:), q_new(:,:) 
      complex(kind=CmplxKind), allocatable :: p_old(:), p_new(:)
      complex(kind=CmplxKind), allocatable :: z_old(:), z_new(:) 
     !!xDx!! MAtrices related to basis, Projection,
      complex(kind=CmplxKind), allocatable ::  G(:,:), R(:,:)
     !!xDx!! Pointers
      type ( MixListCmplxStruct ), pointer :: plq_tmp, qloc_tmp, zloc_tmp, xiloc_tmp 
     !!xDx!! Variables related to Broyden part
      integer (kind=IntKind) :: lastm1,nn
      real (kind=RealKind) :: w0, wtmp, dfnorm, fnorm, fac2, fac1
      real (kind=RealKind) :: rms !!xDx!!Rms of pot. allocated in setupMixRealArrayList
      complex(kind=CmplxKind) :: aij, gmi, cmj
      complex(kind=CmplxKind), allocatable :: aM_c(:,:), bM_c(:,:), cMM_C(:)
     !!xDx!! Message passing
      complex(kind=CmplxKind) :: Ri11, Ri22
      complex(kind=CmplxKind) :: msgbuf(1:2)
      complex(kind=CmplxKind), allocatable :: msgbuf1(:) 
      complex(kind=CmplxKind), allocatable :: msgbuf2(:,:) 
      complex(kind=CmplxKind), pointer :: pad(:), pbd(:)
      integer (kind=IntKind) :: INFO
     !!xDx!! Misc/Debugging/Trash variables
      real (kind=RealKind) :: Mux
      character (len=19) fmt !!xDx!! 19 is no. of char in (F7.4,"+",F7.4,"i")
      character (len=21) fmt1 
   
!   #ifdef USE_SCALAPACK
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
!   #endif
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
!  General parameters
!
      logical :: Initialized = .false.
      logical :: InitializedWkSpace = .false.

      iter = RPMIter
      plq_tmp => list_q
   
!??      if ( .not.InitializedWkSpace ) then
!??         call setWorkingSpace_c(MixingMethod)
         InitializedWkSpace = .true.
!??      endif
     !!xDx!! we are assigning all the local pointers to the global variable here
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = idq!getMixID(idt,idq)
            vlen = plq_tmp%size
!            RPMMix(id)%rms = plq_tmp%rms
            RPMMix(id)%vector_old_c => plq_tmp%vector_old(1:vlen) !! n_{in}
            RPMMix(id)%vector_new_c => plq_tmp%vector_new(1:vlen) !! n_{out}
            if ( .not.RPMMix(id)%Initialized ) then
               call setRPMSpace_c(id,vlen)
            endif
            plq_tmp => plq_tmp%next
         enddo                     
      enddo                        
   
	print *, 'Hey1'
!   #ifdef USE_SCALAPACK
     !  ===================================================================
     !  Initialize ScaLAPACK and set up matrix distribution
     !  ===================================================================
     
     !!xDx!! Setting some variables
     NumAtoms = getNumAtoms()

     !!xDx!! initializing some variables
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
     !!nDn!! The values have been chosen by back-calculation so each proc has what it calculates.
     !
      ICTXT = getGroupCommunicator(GroupID)
      NumPEsInGroup = getNumPEsInGroup(GroupID)
     !
	print *, 'M_Z,N_Z,MB_Z,NB_Z,M_G,N_G,MB_G,NB_G,M_R,N_R,MB_R,NB_R',      &
                  M_Z,N_Z,MB_Z,NB_Z,M_G,N_G,MB_G,NB_G,M_R,N_R,MB_R,NB_R
      call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEsInGroup )
	print *, 'NumAtoms', NumAtoms
      call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)
	print *, 'Hey4'
     !  -------------------------------------------------------------------
      call DESCINIT( DESC_G, M_G, N_G, MB_G, NB_G, RSRC_G, CSRC_G, ICTXT &
           , max( 1, NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW) ) , INFO_G )
	print *, 'Hey5'
      call DESCINIT( DESC_R, M_R, N_R, MB_R, NB_R, RSRC_R, CSRC_R, ICTXT, max( 1, NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW) ) , INFO_R )
     !!xDx!! Allocating Matrices related to projection
      allocate ( G( NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW), NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL) ) )
      allocate ( R( NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW), NUMROC(N_R,NB_R,MYCOL,CSRC_R,NPCOL) ) )
     !!xDx!! Allocating local variables                                                           
	print *, 'calRPM 3'
      allocate ( q_old( NumTotalMix, vlen ) )
      allocate ( q_new( NumTotalMix, vlen ) )
      allocate ( p_old( NumTotalMix*vlen ) )
      allocate ( p_new( NumTotalMix*vlen ) )
      allocate ( z_old( MaxBasis ) )
      allocate ( z_new( MaxBasis ) )
      if ( iter .eq. 1) then
         allocate(G_g_c(MaxBasis,NumAtoms*vlen))
         G_g_c=CZERO
      endif
	print *, 'calRPM 4'
     q_old = CZERO
     q_new = CZERO
     p_old = CZERO
     p_new = CZERO
     z_old = CZERO
     z_new = CZERO
     R = CZERO
     G = CZERO
    
	print *, 'calRPM 5'
     !!xDx!! initializing G from G_g_c
      do ii=1,M_G
      do jj=1,N_G
      call infog2l ( ii, jj, DESC_G, NPROW, NPCOL, MYROW, MYCOL, LindR, LindC, ROCSRC, COCSRC) ! This gives us loc coord for global G
      if ( MYROW .eq. ROCSRC .and. MYCOL .eq. COCSRC) then
      G(LindR,LindC) = G_g_c(ii,jj)
      end if
      enddo
      enddo
!   #endif
   
	print *, 'calRPM 6'
         do idt = 1,NumTypes
            do idq = 1,NumQuantities(idt)
               id = idq!getMixID(idt,idq)
               alpha = RPMMix(id)%alpha
               vlen =RPMMix(id)%vlen
               A => RPMMix(id)%A_c(1:vlen,1:2)!!qDq!!
               pf => RPMMix(id)%f_c(1:vlen)
               pdf => RPMMix(id)%df_c(1:vlen)
               pvect_old => RPMMix(id)%vector_old_c(1:vlen)
               pvect_new => RPMMix(id)%vector_new_c(1:vlen)
               pvold    => RPMMix(id)%vold_c(1:vlen) !! RPMMix(id)%vold_c stores old for G update.
               do k = 1,vlen
                  p_old((id-1)*vlen+k) = pvect_old(k)
                  p_new((id-1)*vlen+k) = pvect_new(k)
               enddo
            enddo
         enddo
   
!!cDc!!     !!iDi!! Debugging
!!cDc!!     !!            if ( RPMIter .eq. 2*MaxRPMIter+2 )then 
!!cDc!!                do k=0,NumPEs-1 
!!cDc!!                 if ( k .eq. MyPE ) then 
!!cDc!!                  call printDataInteger('pvect_new(:),iter',iter,MYPE,1)
!!cDc!!                  call printDataComplexArray('pvect_new(:)',pvect_new(1:vlen),vlen,MyPE, 2) 
!!cDc!!                  call printDataInteger('pvect_old(:),iter',iter,MYPE,2)
!!cDc!!                  call printDataComplexArray('pvect_old(:)',pvect_old(1:vlen),vlen,MyPE, 2) 
!!cDc!!                 endif 
!!cDc!!                enddo 
!!cDc!!     !!            endif 
!!cDc!!     !!iDi!! Debugging
   
   
     !!xDx!! Incrementing iter counters----------------------------------------------------------------------------------------------------------------------------
               RPMIter=RPMIter+1 !! Incrementing RPMIter
               if (RPMIter .ge. MaxRPMIter) then!D!Only after the first MaxRPMIter we generate an unstable subspace
               iter_broy= mod( RPMIter, MaxRPMIter )!D!After every MaxRPMIter steps. we get a new subspace
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
     !--------------------------------------------------
   
         do idt = 1,NumTypes
            do idq = 1,NumQuantities(idt)
               id = idq!getMixID(idt,idq)
               alpha = RPMMix(id)%alpha
               vlen =RPMMix(id)%vlen
               A => RPMMix(id)%A_c(1:vlen,1:2)!!qDq!!
               pf => RPMMix(id)%f_c(1:vlen)
               pdf => RPMMix(id)%df_c(1:vlen)
               pvect_old => RPMMix(id)%vector_old_c(1:vlen)
               pvect_new => RPMMix(id)%vector_new_c(1:vlen)
               pvold    => RPMMix(id)%vold_c(1:vlen) !! RPMMix(id)%vold_c stores old for G update???
         !xDx!------ Defining the stable subspace projection
               q_old(id,1:vlen) =pvect_old(1:vlen)-p_old((id-1)*vlen+1:id*vlen)!! q<- u-p
               q_new(id,1:vlen) =pvect_new(1:vlen)-p_new((id-1)*vlen+1:id*vlen)!! q<- u-p
   
     !!xDx!! Performing Simple Mixing ----------------------------------------------------------------------------------------------------------------------------
               if ( alpha .ge. ONE ) then
               q_new(id,1:vlen)=alpha*q_new(id,1:vlen)
               else
               q_new(id,1:vlen)=alpha*q_new(id,1:vlen)+(1-alpha)*q_old(id,1:vlen)
               endif
            enddo
         enddo
   
     !!xDx!! Performing Broyden Mixing ---------------------------------------------------------------------------------------------------------------------------
   
     !!riDri!! SNGLPROCBROY:    if ( MYROW .eq. 0 .and. MYCOL .eq. 0 ) then !!iDi!!
      if ( iter_broy .eq. 1) then
   
     !     ===============================================================
     !     first iteration: preform linear mixing, load f and vold, set
     !                      different pointers and variables
     !     ===============================================================
         lastIter = 0             ! initialize pointers
         lastm1 = lastIter - 1
     !
               id = 1
               alpha = RPMMix(id)%alpha
               pf => RPMMix(id)%f_c(1:MaxBasis)
               pvold    => RPMMix(id)%vold_c(1:MaxBasis) !! BroydenMix(id)%vold_c used to store pvold for P-subspace
               !!xDx!! everytime we restart the broyden mixing, we allocate everything to zero.
                  pu    => RPMMix(id)%u_c(1:MaxBasis,1:NumRPMIter_broy)!!xDx!! u as in eqn (13b.iii)
                  pvt   => RPMMix(id)%vt_c(1:MaxBasis,1:NumRPMIter_broy)!!*D*!!
                  pdf   => RPMMix(id)%df_c(1:MaxBasis)
                  pvold = CZERO
                  pu    = CZERO
                  pvt   = CZERO
                  pdf   = CZERO
     !
                  pf = z_new - z_old
                  pvold = z_old
                  z_new = z_old + alpha*pf
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
   
!!cDc!!     !!iDi!! Debugging
!!cDc!!     !!            if ( RPMIter .eq. 2*MaxRPMIter+2 )then 
!!cDc!!                do k=0,NumPEs-1 
!!cDc!!                 if ( k .eq. MyPE ) then 
!!cDc!!                  call printDataInteger('z_new(:),iter',iter,MYPE,1)
!!cDc!!                  call printDataComplexArray('z_new(:)',z_new(1:vlen),MaxBasis,MyPE, 2) 
!!cDc!!                  call printDataInteger('z_old(:),iter',iter,MYPE,2)
!!cDc!!                  call printDataComplexArray('z_old(:)',z_old(1:vlen),MaxBasis,MyPE, 2) 
!!cDc!!                 endif 
!!cDc!!                enddo 
!!cDc!!     !!            endif 
!!cDc!!     !!iDi!! Debugging
   
            w0= ten2m2
               dfnorm = zero
               fnorm = zero
               id = 1
               alpha = RPMMix(id)%alpha
     	!!xDx!!
     	print *, 'alpha', RPMMix(1)%alpha
     	!!xDx!!
               pf  => RPMMix(id)%f_c(1:MaxBasis) !!This is F^(iter-1):=F[n^(iter-1)]=n_{out}[n^(iter-1)]-n^(iter-1).
            !!nDn!! based on what follows, we want to use pvect_{new,old} as temporary placeholders for p_{new,old}, their projected counterparts.
               pdf => RPMMix(id)%df_c(1:MaxBasis)
               if ( MaxBasis==1 ) then
                  z_new(1) = z_old(1) + alpha * pf(1) !!xDx!! ??? why ???
               endif
               do k = 1,MaxBasis !!cDt!!x2
                  pdf(k) = z_new(k) - z_old(k) - pf(k)!! dR^(itr)=R^(itr)-R^(itr-1)
                  pf(k)  = z_new(k) - z_old(k) !!R^(itr)
               enddo
     !        
               do k = 1,MaxBasis !!xDx!! sums over all the atoms in the proc
               dfnorm = dfnorm + real(pdf(k)*conjg(pdf(k)),kind=RealKind)
               fnorm  = fnorm  + real(pf(k)*conjg(pf(k)),kind=RealKind)
               enddo
     !
     !           ------------------------------------------------------------
               dfnorm = sqrt( dfnorm )!!cDt!!
               fnorm  = sqrt( fnorm  )!!cDt!!
     !           ============================================================
     !
               fac2 = one/dfnorm!!cDt!!
                  id = 1
                  alpha = RPMMix(id)%alpha
                  fac1 = alpha*fac2
                  pu=> RPMMix(id)%u_c(1:MaxBasis,1:NumRPMIter_broy)!!xDx!!u in eqn (13b.iii)
                  pvt=> RPMMix(id)%vt_c(1:MaxBasis,1:NumRPMIter_broy)!!*D*!!
                  pdf=> RPMMix(id)%df_c(1:MaxBasis)
                  pvold=> RPMMix(id)%vold_c(1:MaxBasis) !! BroydenMix(id)%vold_c used for BRoyden method
                  do k = 1,MaxBasis!!cDt!!x3
                     z_new(k) = fac1*pdf(k) + fac2*(z_old(k) - pvold(k)) 
                     !! aI(dR/|dR|)+(p^(2)-p^(1))/|dR|, Eqn (13b)
                     pvold(k)   = z_old(k)
                     z_old(k) = fac2*pdf(k)  !!=(dR/|dR|) storeds the part of the unit vector of dF.
                  enddo  
     !              ---------------------------------------------------------
                  call broy_sav_c( pu, pvt, z_old, z_new, iter_broy-1, & 
                                   NumRPMIter_broy, MaxBasis )!!cDt!!
     !           ============================================================
               !!xDx!! allocating the Matrices aM_c and bM_c=aM_c^-1
               allocate (aM_c(1:nn,1:nn))
               allocate (cMM_c(1:nn))
               allocate (bM_c(1:nn,1:nn))
     !
               !!xDx!! nn -> #iters-1, equiv to m in DD Johnson.
               do j = 1,nn - 1          ! off diagonal elements of a(i,j)
                  do i = j+1,nn
                     aij = CZERO
                        id = 1
                        pvt => RPMMix(id)%vt_c(1:MaxBasis,1:NumRPMIter_broy)
                        do k = 1,MaxBasis
                           aij = aij + conjg(pvt(k,j))*pvt(k,i)
                        enddo
                     aM_c(i,j) = aij
                     aM_c(j,i) = conjg(aij)
                  enddo
               enddo
     !
               do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
                  aij = CZERO
                  cmj = CZERO
                     id = 1
                     pvt => RPMMix(id)%vt_c(1:MaxBasis,1:NumRPMIter_broy)
                     pf  => RPMMix(id)%f_c(1:MaxBasis)
                     do k = 1,MaxBasis!!cDt!!x2
                        cmj = cmj + conjg(pvt(k,i))*pf(k) !!xDx!! c^m_k/w_k in Eqn (15b). Evaluated for fixed m, ie, for current iteration.
                        aij = aij + conjg(pvt(k,i))*pvt(k,i)
                     enddo
                  aM_c(i,i) = aij
                  cMM_c(i) = cmj
               enddo
     !
   
     !           ============================================================
     !
                  id = 1
                  alpha = RPMMix(id)%alpha
                  pvold => RPMMix(id)%vold_c(1:MaxBasis)
                  pw => RPMMix(id)%w_c(1:NumRPMIter_broy) !!xDx!! Weights,
                  pu => RPMMix(id)%u_c(1:MaxBasis,1:NumRPMIter_broy)
                  pf => RPMMix(id)%f_c(1:MaxBasis)
     !
                  if (  iter_broy-1 > NumRPMIter_broy ) then
                     do i = 1,NumRPMIter_broy-1
                        pw(i) = pw(i+1)
                     enddo
                  endif
     !
                  rms = dfnorm !!cDc!!RPMMix(id)%rms
                  wtmp = zero
                  if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
                  if ( wtmp < ONE ) wtmp = ONE
                  if (  iter_broy > NumRPMIter_broy ) then
                     pw(NumRPMIter_broy) = wtmp
                  else
                     pw(lastIter) = wtmp       !w(lastm1)=wtmp
                  endif
   
   	if ( MyPE .eq. 0 ) then
   	write (MyPE+30,1122) '                        '
   	write (MyPE+30,1122) '========================'
   	write (MyPE+30,1122) 'MyPE','wtmp'
   	write (MyPE+30,1122) '========================'
   	write (MyPE+30,'(X,F20.10)') wtmp
   	write (MyPE+30,1122) '                        '
   	endif
   
     !              =========================================================
     !              now calculate the b-matrix:
     !              b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1 !!xDx!! Eqn. (13a). Assuming this is just the standard calc of \beta.
     !              =========================================================
     !                 ======================================================
     !                 use lapack !!nDn!! everything is local now.
     !                 ======================================================
                     bM_c =CZERO
                     do i = 1,nn
                        do j = 1,nn
                           bM_c(j,i) = aM_c(j,i)*pw(j)*pw(i)!!cDt!!
                        enddo
                        bM_c(i,i) = w0**2 + aM_c(i,i)*pw(i)*pw(i)!!cDt!!
                     enddo
                     if ( .not.allocated(ipiv) ) then
                        allocate( ipiv(nn) )
                     endif
     !                 ------------------------------------------------------
                     call zgetrf( nn, nn, bM_c, nn,ipiv, info ) !!cDt!!
     !                 ------------------------------------------------------
                     call zgetri( nn, bM_c, nn, ipiv, tmp, -1, info )!!cDt!!
     !                 ------------------------------------------------------
                     LWork = int(real(tmp,kind=RealKind))
                     allocate( Work(1:LWork) )
     !                 ------------------------------------------------------
                     call zgetri( nn, bM_c, nn, ipiv, Work, LWork, info )!!cDt!!
     !                 ------------------------------------------------------
                     deallocate (Work)
     !                 ------------------------------------------------------
                     deallocate (ipiv)
     !                 ------------------------------------------------------
     !
                     do k = 1,MaxBasis
                        z_new(k) = pvold(k) + alpha*pf(k)!!cDt!!
                     enddo
                     do i = 1,nn
                        gmi = CZERO
                        do j = 1,nn!!cDt!!
                           gmi = gmi + cMM_c(j)*bM_c(j,i)*pw(j) !!xDx!! \gamma_{ml} in Eqn. (15b), the cMM_c is wihtout w_j
                        enddo
                        do k = 1,MaxBasis
                           z_new(k) = z_new(k) - gmi*pu(k,i)*pw(i)!!cDt!!
                        enddo
                     enddo
     !
               !!xDx!! Deallocating the MAtrices aM_c and bM_c=aM_c^-1
               deallocate(aM_c)
               deallocate(bM_c)
               deallocate(cMM_c)
   
!!cDc!!     !!iDi!! Debugging
!!cDc!!     !!            if ( RPMIter .eq. 2*MaxRPMIter+2 )then 
!!cDc!!                  do k=0,NumPEs-1 
!!cDc!!                   if ( k .eq. MyPE ) then 
!!cDc!!                    call printDataInteger('z_new(:),iter, aftr',iter,MYPE,1)
!!cDc!!                    call printDataComplexArray('z_new(:)',z_new(1:vlen),MaxBasis,MyPE, 2) 
!!cDc!!                    call printDataInteger('z_old(:),iter, aftr',iter,MYPE,2)
!!cDc!!                    call printDataComplexArray('z_old(:)',z_old(1:vlen),MaxBasis,MyPE, 2) 
!!cDc!!                   endif 
!!cDc!!                  enddo 
!!cDc!!     !!            endif 
!!cDc!!     !!iDi!! Debugging
   
   
      endif   
   
     !!xDx!! End of Broyden Mixing -------------------------------------------------------------------------------------------------------------------------------
         do idt = 1,NumTypes
            do idq = 1,NumQuantities(idt)
               id = idq!getMixID(idt,idq)
               alpha = RPMMix(id)%alpha
               vlen =RPMMix(id)%vlen
               A => RPMMix(id)%A_c(1:vlen,1:2)
               pf => RPMMix(id)%f_c(1:vlen)
               pdf => RPMMix(id)%df_c(1:vlen)
               pvect_old => RPMMix(id)%vector_old_c(1:vlen)
               pvect_new => RPMMix(id)%vector_new_c(1:vlen)
               pvold    => RPMMix(id)%vold_c(1:vlen)
               pv_st_old    => RPMMix(id)%v_st_old_c(1:vlen)
            !!xDx!! projecting local back to global variables
               pvect_new(1:vlen)= matmul( Transpose(conjg( G(1:MaxBasis,(id-1)*vlen+1:id*vlen) ) ), z_new(1:MaxBasis) ) !!iDi!!
               pvect_old(1:vlen)= matmul( Transpose(conjg( G(1:MaxBasis,(id-1)*vlen+1:id*vlen) ) ), z_old(1:MaxBasis) ) !!iDi!!
            !!xDx!! Combining the resultant to return to main -------------------------------------------------------------------------------------------------------
     !            RPMMix(id)%vector_new_c(1:vlen)=q_new(id,1:vlen)+pvect_new!! u<- q+p
               do k=1,vlen
                  pvect_new(k)=q_new(id,k)+pvect_new(k)!! u<- q+p
               enddo
     !
!!cDc!!   		print *, 'id',id
!!cDc!!     !!iDi!! Debugging
!!cDc!!     !!            if ( RPMIter .eq. 2*MaxRPMIter+2 )then 
!!cDc!!                do k=0,NumPEs-1 
!!cDc!!                 if ( k .eq. MyPE ) then 
!!cDc!!                  call printDataInteger('pvect_new(:),iter,aftr',iter,MYPE,1)
!!cDc!!                  call printDataComplexArray('pvect_new(:)',pvect_new(1:vlen),vlen,MyPE, 2) 
!!cDc!!                  call printDataInteger('pvect_old(:),iter,aftr',iter,MYPE,2)
!!cDc!!                  call printDataComplexArray('pvect_old(:)',pvect_old(1:vlen),vlen,MyPE, 2) 
!!cDc!!                 endif 
!!cDc!!                enddo 
!!cDc!!     !!            endif 
!!cDc!!     !!iDi!! Debugging
   
     !!xDx!! Increasing size of basis --------------------------------------------------------------------------------------
               if ( RPMIter .eq. 1) then
                  do k = 1,vlen
                     pv_st_old(k) = q_old(id,k)
                  enddo
               elseif ( RPMIter .eq. 2) then
                  do k = 1,vlen
                     A(k,1)=q_old(id,k)-pv_st_old(k) 
                     pv_st_old(k)=q_old(id,k)
                  enddo
               else
                  A(:,2)=A(:,1)
                  do k = 1,vlen
                     A(k,1)=q_old(id,k)-pv_st_old(k) 
                     pv_st_old(k)=q_old(id,k)
                  enddo
               endif
   
               if ( mod( RPMIter, MaxRPMIter ) .eq. 0 ) then !! Inc basis size.
                  do ii=1,2
                     do jj=1,vlen
                        R( ii, (id-1)*vlen+jj ) = conjg(A(jj,ii))
                     enddo
                  enddo
               endif
   
   	fmt1 = '(F19.4,"+",F19.4,"i")'
   	do k=0,NumPEs-1
   	if ( MyPE .eq. k .and. mod( RPMIter, MaxRPMIter ) .eq. 0) then
   	write (MyPE+40,1122) '                        '
   	write (MyPE+40,1122) '========================'
   	write (MyPE+40,1342) ' MyPE','R^T(:,:),RPMIter, id', MyPE, RPMIter, id
   	write (MyPE+40,1122) '========================'
   	do jj=1,vlen
   	   do ii=1,2
                 fmt1(9:9) = MERGE('+',' ',imag(conjg(A(jj,ii))).gt.0)
                 write(MyPE+40,fmt1, advance="no") conjg(A(jj,ii))
   	   enddo
   	   write(MyPE+40, "(a)") " "
   	enddo
   	write (MyPE+40,1122) '                        '
   	endif
   	enddo
   
   
            enddo
         enddo
     !                                                   
     !!xDx!! Performing the calc to increase the Basis size                                                                                                          
         if ( mod( RPMIter, MaxRPMIter ) .eq. 0 ) then !!nDn!! iter_broy= 0 & RPMIter> MaxRPMIter
   
         !!xDx!! HouseHolder LQ---------------------------------------------------------------- 
!   #ifdef USE_SCALAPACK
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
!   #endif
   
   	fmt1 = '(F19.6,"+",F19.6,"i")'
   	if ( MyPE .eq. 0 ) then
   	write (MyPE+30,1122) '                        '
   	fmt1(9:9) = MERGE('+',' ',imag(Ri11).gt.0)
   	write (MyPE+30,1122) '========================'
   	write (MyPE+30,1122) 'MyPE','Ri11'
   	write (MyPE+30,1122) '========================'
   	write (MyPE+30,fmt1) Ri11
   	write (MyPE+30,1122) '                        '
   	fmt1(9:9) = MERGE('+',' ',imag(Ri22).gt.0)
   	write (MyPE+30,1122) '========================'
   	write (MyPE+30,1122) 'MyPE','Ri22'
   	write (MyPE+30,1122) '========================'
   	write (MyPE+30,fmt1) Ri22
   	write (MyPE+30,1122) '                        '
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
   
     !!xDx!! Re-orthogonalize G based on new value of basis------------------------------------------------------------------------------------------------------------
!   #ifdef USE_SCALAPACK
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
!   #endif
   
   	fmt = '(F9.6,"+",F9.6,"i")'
   	if ( MyPE .eq. 0) then
   	write (MyPE+20,1122) '                        '
   	write (MyPE+20,1122) '========================'
   	write (MyPE+20,1342) ' MyPE','G_g_c^T(:,:),RPMIter, Basis', MyPE, RPMIter, Basis
   	write (MyPE+20,1122) '========================'
   	do jj=1,NumAtoms*vlen
   	   do ii=1,MaxBasis
                 fmt(8:8) = MERGE('+',' ',imag(G_g_c(ii,jj)).gt.0)
                 write(MyPE+20,fmt, advance="no") G_g_c(ii,jj)
   	   enddo
   	   write(MyPE+20, fmt="(a)") " "
   	enddo
   	write (MyPE+20,1122) '                        '
        1122 format(X,A)
        1342 format(2(X,A),3(X,I4))
   	endif
   
   
         endif !!nDn!! iter_broy= 0
   
     !!xDx!! Deallocate variables
           deallocate(G)
           deallocate(R)
     !!xDx!! Allocating local variables                                                           
           deallocate (q_old)
           deallocate (q_new)
           deallocate (p_old)
           deallocate (p_new)
           deallocate (z_old)
           deallocate (z_new)
     !
      nullify(plq_tmp)
     !
      end subroutine calRPM_c


      subroutine broy_sav_c( fins, fots, vector_old, vector_new, itscf,  & 
                             istore, ivsiz )
!  ==================================================================
      implicit none
!
      integer(kind=IntKind), intent(in) :: itscf
      integer(kind=IntKind), intent(in) :: ivsiz
      integer(kind=IntKind), intent(in) :: istore
!
      complex(kind=CmplxKind), target :: vector_old(:)
      complex(kind=CmplxKind), target :: vector_new(:)
      complex(kind=CmplxKind), target :: fins(:,:)
      complex(kind=CmplxKind), target :: fots(:,:)
!
      integer(kind=IntKind) :: i, j
!
!  ==================================================================
!  write(6,'('' IN BROY_SAV: istore,itscf '',2i5)') istore,itscf
!  ==================================================================
      if ( itscf <= istore ) then
!     ===============================================================
!     Load the first istore iterations in increasing iteration count
!     ===============================================================
         do i = 1,ivsiz
            fins(i,itscf) = vector_new(i)
         enddo
!
         do i = 1,ivsiz
            fots(i,itscf) = vector_old(i)
         enddo
!
      else
!     ===============================================================
!     Re-load so that the ordering is in increasing iteration count
!     ===============================================================
         do j = 1,istore - 1
!        write(6,'('' IN BROY_SAV: j,j+1 '',2i5)') j,j+1
            do i = 1,ivsiz
               fins(i,j) = fins(i,j+1)
            enddo
!
            do i = 1,ivsiz
               fots(i,j) = fots(i,j+1)
            enddo
!
         enddo
!     ===============================================================
!     Load current charge densities in the last storage location
!     ===============================================================
         do i = 1,ivsiz
            fins(i,istore) = vector_new(i)
         enddo
!
         do i = 1,ivsiz
            fots(i,istore) = vector_old(i)
         enddo
!
      endif
!
      end subroutine broy_sav_c

      end module MixingModule

!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

      program RPM_2
      use KindParamModule, only : IntKind, RealKind, CmplxKind, QuadCmplxKind
      use ErrorHandlerModule, only : ErrorHandler, WarningHandler
      use MathParamModule, only : zero, one, two, three, ten2m2, ten2m3,ten2m8, &
                                  ten2m9, ten2m5, czero, cone, PI4, THIRD, sqrtm1
!   
      use PublicTypeDefinitionsModule, only : MixListCmplxStruct
!
       use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MyPE, GlobalSum, &
                            getCommunicator
      use GroupCommModule, only : initGroupComm, endGroupComm
      use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
      use GroupCommModule, only : GlobalMaxInGroup, getGroupCommunicator
      use GroupCommModule, only : syncAllPEsInGroup
      use GroupCommModule, only : createProcGrid, createBoxGroupFromGrid, &
                               getProcGridID, GlobalSumInGroup
!   
      use DataServiceCenterModule, only : getDataStorage,                 &
                                       ComplexType, ComplexMark,       &
                                       createDataStorage,              &
                                       isDataStorageExisting,          &
                                       setDataStorageLDA,              &
                                       setDataStorage2Value,           &
                                       deleteDataStorage
!
      use MixingModule, only : initMixing, endMixing, mixValues, setRPMMixing
!
!
      implicit none

      integer (kind=IntKind) :: typ=14, MaxRPMIter, vlen, Basis
      integer (kind=IntKind) :: NumBroydenIter, TotalRPMIter
      integer (kind=IntKind) :: LocalNumAtoms, NumAtoms, NumTotalMix
      integer (kind=IntKind) :: MaxBasis, NumRPMIter_broy
      integer (kind=IntKind) :: RPMIter
!
      integer (kind=IntKind) :: k
!
      real(kind=RealKind) :: normp, alpha, beta, tol   
      complex(kind=CmplxKind), allocatable :: pn(:), po(:)     
!
      character (len=50) :: PEStr     
!
      integer (kind=IntKind) :: id, ia, nr, is
      integer (kind=IntKind) :: data_size, data_size_ns, data_size_save
      integer (kind=IntKind) :: lmax, jl, jl_nonZero, ind_jl, jmax, MaxSpecies
      integer (kind=IntKind) :: NumMix(1)
!
      complex (kind=CmplxKind), pointer :: pStore_old(:,:), pStore_new(:,:)
      type (MixListCmplxStruct), pointer :: p_CAL
!
      integer (kind=IntKind) :: dim(1), box(1), color, key, my_rank, group_size      
      integer (kind=IntKind) :: itmp, ictxt, ictxt1, ictxt2, grid, group     
      integer (kind=Intkind) :: new_comm, old_comm, old_group, new_group
      integer (kind=IntKind) :: info                
!
      integer(kind=IntKind) :: GroupID     
!
!
!  ===================================================================
!  Mixing quantity
!  ===================================================================
      type (MixListCmplxStruct), target :: CmplxArrayList
!
        
!! MPI initialization          
     call initMPP()        
          
     call initGroupComm()       
!#ifdef USE_SCALAPACK
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
!!cDc!!     call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEs )
!!cDc!!     call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!!cDc!!!  -------------------------------------------------------------------
!!cDc!!     call DESCINIT( DESC_G, M_G, N_G, MB_G, NB_G, RSRC_G, CSRC_G, ICTXT, max( 1, NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW) ) , INFO_G )
!!cDc!!     call DESCINIT( DESC_R, M_R, N_R, MB_R, NB_R, RSRC_R, CSRC_R, ICTXT, max( 1, NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW) ) , INFO_R )
!!cDc!!!!xDx!! Allocating Matrices related to projection
!!cDc!!     allocate ( G( NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW), NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL) ) )
!!cDc!!     allocate ( R( NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW), NUMROC(N_R,NB_R,MYCOL,CSRC_R,NPCOL) ) )
!#endif

!! The Equation to solve F(u)=u        
      if (typ == 14) then
         vlen=2; Basis=0; MaxRPMIter=8; NumBroydenIter=10; alpha= 1D-01;
         beta=1-alpha; tol = 1D-09; TotalRPMIter=100; LocalNumAtoms=2
      else
         print *, 'Please eneter a vaild typ'
      endif

      NumAtoms=NumPEs*LocalNumAtoms
      NumTotalMix=LocalNumAtoms !! by defn also =NumMix, i/p to initMixing_c 
      nr=vlen
      MaxBasis=4
      NumRPMIter_broy=MaxRPMIter-2

!!nDn!! Based on how the LSMS code is set up, we have 
!!nDn!! N= #(Atoms/box)x#(RadialMesh)x#(Lvalues)x#(PEs/box)
!!nDn!! for typ=4 (N=8), we want the parameters to be
!!nDn!! #(Atoms/box)=2, #(RadialMesh)=2, #(Lvalues)=1, #(PEs/box)=2

!!xDx!!Allocation arrays
      allocate(pn(NumPEs*LocalNumAtoms*vlen))
      allocate(po(NumPEs*LocalNumAtoms*vlen))


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
   
         NumMix(1)=LocalNumAtoms

         call setNumAtoms(NumAtoms)
 
         call initMixing( 1, NumMix, CmplxArrayList )

!!xDx!! setupMixingScheme
      is = 1
!  do is = 1,n_spin_pola
      do id=1,LocalNumAtoms
            call setRPMMixing(is, id, alpha, MaxBasis, MaxRPMIter,        &
                              NumRPMIter_broy)
      enddo

      ITERAA:   do RPMIter = 1,TotalRPMIter!!1a.!!1,MaxRPMIter, is for the rest.
      
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

     !!xDx!! setupMixCmplxArrayList
            data_size = 0
            MaxSpecies = 0
            do id = 1,LocalNumAtoms
               nr = vlen !!nr=2 for typ=4
               lmax = 0 !! set to 0
               jl_nonZero = 1
               jmax = (lmax+1)*(lmax+2)/2 !! lmax=0 =>jmax=1
               data_size = max(data_size,nr*jl_nonZero) !! shld be 2. so that each PE has 2x1 vector.
               MaxSpecies = 1 !?!max(MaxSpecies, getLocalNumSpecies(id)) !!set to 1
            enddo
            data_size_ns = data_size ! data_size*n_spin_pola+n_spin_pola !! WTSet =data_size
         !
            data_size_save = data_size_ns
         !
            if (.not.isDataStorageExisting('MixingVectorOld')) then
         !     ----------------------------------------------------------------
               call createDataStorage('MixingVectorOld',                        &
                                      LocalNumAtoms*MaxSpecies*data_size_ns     &
                                      ,ComplexType)
               call setDataStorage2Value('MixingVectorOld',CZERO)
            endif
            if (.not.isDataStorageExisting('MixingVectorNew')) then
               call createDataStorage('MixingVectorNew',                        &
                                      LocalNumAtoms*MaxSpecies*data_size_ns,ComplexType)
               call setDataStorage2Value('MixingVectorNew',CZERO)
            endif
         !
            pStore_old => getDataStorage('MixingVectorOld',data_size_ns,        &
                                         LocalNumAtoms,ComplexMark)
            pStore_new => getDataStorage('MixingVectorNew',data_size_ns,        &
                                         LocalNumAtoms,ComplexMark)
         !
            p_CAL => CmplxArrayList
            do id = 1, LocalNumAtoms
               nr = vlen !!=2 (typ=4)
               p_CAL%size = data_size_ns
               p_CAL%vector_old => pStore_old(1:data_size_ns,id)
               p_CAL%vector_new => pStore_new(1:data_size_ns,id)
               p_CAL%vector_old(:) = CZERO
               p_CAL%vector_new(:) = CZERO
               p_CAL%rms = ZERO
               p_CAL%vector_old(1:vlen) = po(MyPE*vlen*LocalNumAtoms+(id-1)*vlen+1:MyPE*vlen*LocalNumAtoms+id*vlen) 
               p_CAL%vector_new(1:vlen) = pn(MyPE*vlen*LocalNumAtoms+(id-1)*vlen+1:MyPE*vlen*LocalNumAtoms+id*vlen) 
      !
               if (associated(p_CAL%next)) then
                  p_CAL => p_CAL%next
               else if ( id/=LocalNumAtoms ) then
                  call ErrorHandler('setupMixCmplxArrayList',          &
                            'CmplxArrayList is not set up properlyi',id)
               endif
            enddo
      !!xDx!!----------------------------------------------------------------------------------------------
  	print *, 'Hello 5' 
      !!xDx!! calling Mixing-------------------------------------------------------------------------------
            call mixValues(CmplxArrayList)
!            call calRPM_c(CmplxArrayList)
	print *, 'Hello 6' 
      !!xDx!! update Mix CmplxValues----------------------------------------------------------------------
         data_size = 0
         MaxSpecies = 0
         do id = 1,LocalNumAtoms
            nr = vlen
            lmax = 0
            jl_nonZero = 1
            data_size = max(data_size,nr*jl_nonZero)
            MaxSpecies = 1 !?!max(MaxSpecies,getLocalNumSpecies(id))
         enddo
         data_size_ns = data_size !data_size*n_spin_pola+n_spin_pola
         if (.not.isDataStorageExisting('MixingVectorOld')) then
            call ErrorHandler('updateMixCmplxValues','MixingVectorOld not  &
            defined')
         endif
         if (.not.isDataStorageExisting('MixingVectorNew')) then
            call ErrorHandler('updateMixCmplxValues','MixingVectorNew not  &
            defined')
         endif
      !
         pStore_old => getDataStorage('MixingVectorOld',data_size_ns, &
                                      LocalNumAtoms,ComplexMark)
         pStore_new => getDataStorage('MixingVectorNew',data_size_ns, &
                                      LocalNumAtoms,ComplexMark)
         p_CAL => CmplxArrayList
         do id = 1, LocalNumAtoms
            nr = vlen
            p_CAL%vector_old => pStore_old(1:data_size_ns,id)
            p_CAL%vector_new => pStore_new(1:data_size_ns,id)
            pn = CZERO
            pn(MyPE*vlen*LocalNumAtoms+(id-1)*vlen+1:MyPE*vlen*            & 
            LocalNumAtoms+id*vlen) = p_CAL%vector_new(1:vlen)
            p_CAL => p_CAL%next
         enddo
      !!xDx!!
	print *, 'Hello 7'
   
!           ============================================================
            normp=ZERO
            do k = 1,NumPEs*LocalNumAtoms*vlen !!xDx!! These are sums over all the atoms in the proc. 
               normp  = normp  + real((pn(k)-po(k))*conjg(pn(k)-po(k)),kind=RealKind) 
            enddo

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
                  do k = 1,NumPEs*LocalNumAtoms*vlen 
                     write (*,5241) 'k','pn',k,real(pn(k)),aimag(pn(k)),'*i'
                  enddo
                  write (*,1121) '========================'
                  write (*,1121) '!!xDx!!!xDx!!!xDx!!!xDx!!!xDx!!'
                  write (*,1121) '========================'
               endif
               call syncAllPEs()
               EXIT ITERAA
               1121 format(X,A)
               1122 format(X,A)
               1341 format(2(X,A),(X,I4),(X,F20.10),(X,I4))
               1342 format(2(X,A),3(X,I4))
               5241 format(2(X,A),(X,I4),(X,F20.10,SP,F16.10,A))
            endif

!           ============================================================

      enddo ITERAA

      !!Deallocation
      deallocate(pn)
      deallocate(po)
      
      call endGroupComm()
      call endMPP()
      
         stop 'Ok'
      end program RPM_2



