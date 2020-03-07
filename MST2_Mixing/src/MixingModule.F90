!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  History: Original code written by D.D. Johnson (see PRB 38, 12807)
!                  note:  there are a few typos in that paper but 
!                  the code is working!
!              Rewritten by W. A. Shelton for LSMS code 6/21/94
!                  this version is easy to read (no goto!!!! more comments ...)
!                  and is setup for MPP machines (not tested)
!              Rewritten by T. C. Schulthess, ORNL, March 97
!                  this version should work for any code (see comments below)
!
!     Bug fixes:   TCS, 8/5/97 see comments below 
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     further comments on how to use this subroutine:
!     (if nothing useful stands here I had no time yet to update these
!     comments, please consult usage in lkkr code version 0.6.3 or later,
!     or call Thomas Schulthess (423) 5768067)
!
!      vector(r,i) -> i=1: old vector (input), scratch (ouput)
!                  -> i=2: new vector (input), mixed vector (output)
!      vlen        -> length of vector
!      alpha       -> linear mixing factor
!      rms         -> RMS difference between old and new vector
!      iter        -> iteration number (if 1, linear mixing, broyden reset)
!      broylen     -> number of iterations that are used for mixing
!                     (<=MaxNumBroydenIter)
!      u, vt, f, df, w and vold are working arrays that need to be saved
!                    between call to this subroutine
!      a, b, d, cm, and ipiv are working arrays that need not be saved
!      MaxNumBroydenIter    -> maximum number of iterations that can be saved
!      MaxSizeBroydenVect       -> maximum length of vectors
!
!      See declaration for exact dimentsions and types
!
!      There are two options for matrix inversions, a Gaussian
!      elimination routine called invert1 and calls to lapack routines
!      with pivoting (see comments "using invert1" and "using lapack").
!      Obviously only one can be used, comment out the other one.
!
!      When using this subroutine in a parallel code in which only 
!      parts of the vectors are known on every node, make sure that 
!      the calls to gldsum (global sum) are correct (LKKR and LSMS 
!      codes have different calls).
!
!      In a serial code, either comment out the calls to GlobalSum or
!      provide a dummy subroutine
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!xDx!! This module receives q_list which has pvect_old and pvect_new
!!xDx!! each of size vlen. pvect_new is the output of the Kohn-Sham map
!!xDx!! pvect_old is the value of the potential that was used as input to
!!xDx!! the KS map. The variable vlen is of size (jmax)x(NumRpts). For MT
!!xDx!! calc jmax is 1 (only spherical component). When there is a symmetry
!!xDx!! then some components are zero due to it. For cubic symm all even values
!!xDx!! of l (except l=2) and values of m divisible by 4 are non-zero. We can
!!xDx!! identify this using an lflag.
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
!!xDx!! Mainly need NumAtoms
   use SystemModule, only : getNumAtoms, getNumAlloyElements, getAlloyElementContent
!
public :: initMixing,         &
          mixValues,          &
          setRPMMixing,       &
          setBroydenMixing,   &
          setSimpleMixing,    &
          setDGAMixing,       &
          setMixingAlpha,     &
          resetMixing,        &
          endMixing
!
   interface initMixing
      module procedure initMixing_r, initMixing_c
   end interface initMixing
!
   interface mixValues
      module procedure mixValues_r, mixValues_c
   end interface mixValues
!
   interface resetMixing
      module procedure resetMixing_r, resetMixing_c
   end interface resetMixing
!
!!xDx!! Some real RPM Parameters made Public
   integer(kind=IntKind) :: MaxBasis, Basis=0, NumAtoms 
!!xDx!! Global G matrix. = Transpose of Global Z matrix and global storage for Projection Matrix.
   real(kind=RealKind), allocatable  :: G_g(:,:)
   complex(kind=CmplxKind), allocatable  :: G_g_c(:,:)
!
private
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
   integer(kind=IntKind) :: MixingMethod = -1 ! -1 - No method have been set yet
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
!  Simple Mixing Parameters
!
type SimpleStruct
logical :: Initialized
!
integer (kind=IntKind) :: vlen
!
real(kind=RealKind) :: alpha
!
real(kind=RealKind), pointer :: vector_old_r(:)
real(kind=RealKind), pointer :: vector_new_r(:)
!
complex(kind=CmplxKind), pointer :: vector_old_c(:)
complex(kind=CmplxKind), pointer :: vector_new_c(:)
end type SimpleStruct
!
type (SimpleStruct), allocatable :: SimpleMix(:)
!
!  Broyden Mixing Parameters
!
integer(kind=IntKind) :: MaxNumBroydenIter = 10
integer(kind=IntKind) :: MaxSizeBroydenVect
integer(kind=IntKind) :: BroydenInvMethod = 1
!
!   integer(kind=IntKind) :: lastIter = 0 !!xDx!! Commented becuase its definitin have been moved above.
!
integer(kind=IntKind) :: NumBroydenIter
integer(kind=IntKind) :: SizeBroydenVect
!
type BroydenStruct
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
real(kind=RealKind), pointer :: df_r(:)
real(kind=RealKind), pointer :: w_r(:)
real(kind=RealKind), pointer :: vold_r(:)
!
complex(kind=CmplxKind), pointer :: u_c(:,:)
complex(kind=CmplxKind), pointer :: vt_c(:,:)
complex(kind=CmplxKind), pointer :: f_c(:)
complex(kind=CmplxKind), pointer :: df_c(:)
complex(kind=CmplxKind), pointer :: w_c(:)
complex(kind=CmplxKind), pointer :: vold_c(:)
end type BroydenStruct
!
type (BroydenStruct), allocatable :: BroydenMix(:)
!
!  D.G.Anderson mixing Parameters
!
integer(kind=IntKind) :: MaxNumDGAIter = 10
integer(kind=IntKind) :: NumDGAIter
integer(kind=IntKind) :: GroupID
!
type DGAStruct
logical :: Initialized
!
integer (kind=IntKind) :: vlen
!
real(kind=RealKind) :: alpha
!
real(kind=RealKind), pointer :: f_old_r(:,:)
real(kind=RealKind), pointer :: f_new_r(:,:)
!
complex(kind=CmplxKind), pointer :: f_old_c(:,:)
complex(kind=CmplxKind), pointer :: f_new_c(:,:)
!
real(kind=RealKind), pointer :: vector_old_r(:)
real(kind=RealKind), pointer :: vector_new_r(:)
!
complex(kind=CmplxKind), pointer :: vector_old_c(:)
complex(kind=CmplxKind), pointer :: vector_new_c(:)
end type DGAStruct
!
type (DGAStruct), allocatable :: DGAMix(:)
!
!  Working Space
!
integer (kind=IntKind), allocatable, target :: ipiv(:)
!
real(kind=RealKind), allocatable, target :: a_r(:,:), b_r(:,:)
real(kind=RealKind), allocatable, target :: d_r(:,:), cm_r(:)
!
complex(kind=CmplxKind), allocatable, target :: a_c(:,:), b_c(:,:)
complex(kind=CmplxKind), allocatable, target :: d_c(:,:), cm_c(:)
!
!!xDx!! Defining the parameters for Parallel Linear Algebra
#ifdef USE_SCALAPACK
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
#endif
!!xDx!!
contains
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initMixing_r( nt, nq, RealArrayList )
!  ===================================================================
use GroupCommModule, only : getGroupID
!
implicit none
!
integer(kind=IntKind), intent(in) :: nt, nq(nt) 
integer(kind=IntKind) :: idt, idq
!
type (MixListRealStruct), target :: RealArrayList
type (MixListRealStruct), pointer :: p_RAL
!
if ( nt<1 ) then
call ErrorHandler('initMixing','Number of Types is less than 1')
endif
!
NumTypes = nt
allocate( NumQuantities(NumTypes) )
NumQuantities = nq
GroupID = getGroupID('Unit Cell') 
!!xDx!! This is paralellized over the no. of atoms in the unit cell and L values. Check RPM Notes
!
NumTotalMix = 0
p_RAL => RealArrayList
nullify( p_RAL%next )
!
do idt = 1,NumTypes
if ( NumQuantities(idt)<1 ) then
 call ErrorHandler('initMixing', &
      'Number of Quntities is less than 1 sor type',idt)
endif
do idq = 1,NumQuantities(idt)
 NumTotalMix = NumTotalMix +1
 if ( idq==NumQuantities(idt)) cycle
 allocate(p_RAL%next)
 p_RAL => p_RAL%next
 nullify( p_RAL%next )
enddo
enddo
!   nullify( p_RAL%next )
!
Initialized = .true.
!
end subroutine initMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initMixing_c( nt, nq, CmplxArrayList )
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
GroupID = getGroupID('Unit Cell')
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
end subroutine initMixing_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine endMixing()
!  ===================================================================
!
implicit none
!
if ( Initialized ) then
deallocate( NumQuantities )
call delWorkingSpace()
call delMixing()
if ( MixingMethod == 0 ) then
 deallocate( SimpleMix )
else if ( MixingMethod == 1 ) then
 deallocate( BroydenMix )
else if ( MixingMethod == 2 ) then
 deallocate( DGAMix )
!!xDx!! modified for test_RPM_c_3 typ=14
else if ( MixingMethod == 3 ) then
 deallocate( RPMMix )
!!xDx!! modified for test_RPM_c_3 typ=14
endif
Initialized = .false.
endif
!
MixingMethod = -1
!
end subroutine endMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mixValues_r( list_q )
!  ===================================================================
!
implicit none
!
type(MixListRealStruct) :: list_q
!
iter_count =iter_count+1
!
if ( MixingMethod == 0 ) then
call calSimpleMixing_r(list_q)
else if ( MixingMethod == 1 ) then
call calBroydenMixing_r(list_q)
else if ( MixingMethod == 3 ) then
call calRPM_r(list_q)
else
call calDGAMixing_r(list_q)
endif
!
end subroutine mixValues_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mixValues_c( list_q )
!  ===================================================================
!
implicit none
!
type(MixListCmplxStruct) :: list_q
!
iter_count =iter_count+1
!
if ( MixingMethod == 0 ) then
call calSimpleMixing_c(list_q)
else if ( MixingMethod == 1 ) then
call calBroydenMixing_c(list_q)
else if ( MixingMethod == 3 ) then
call calRPM_c(list_q)
else 
call calDGAMixing_c(list_q)
endif
!
end subroutine mixValues_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setSimpleMixing( idt, idq, amix )
!  ===================================================================
implicit none
!
integer(kind=IntKind), intent(in) :: idt, idq
!
real(kind=RealKind), intent(in) :: amix
!
integer(kind=IntKind) :: id
!
if ( .not.Initialized ) then
call ErrorHandler('setSimpleMixing',                            &
  'The MixingModule needs to be initialized first')
endif
!
if ( idt<0 .or. idt>NumTypes ) then
call ErrorHandler('setSimpleMixing',"Invalid type index",idt)
endif
if ( idq<0 .or. idq>NumQuantities(idt) ) then
call ErrorHandler('setSimpleMixing',"Invalid quantity index",idq)
endif
!
if ( .not.allocated(SimpleMix) ) then
allocate( SimpleMix(NumTotalMix) )
do id = 1,NumTotalMix
 SimpleMix(id)%Initialized = .false.
enddo 
endif  
!
id = getMixID(idt,idq)
if ( .not.SimpleMix(id)%Initialized ) then
SimpleMix(id)%alpha = amix
endif
!
if ( MixingMethod /=0 ) then
MixingMethod = 0
endif
!
end subroutine setSimpleMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine delSimpleMixing()
!  ===================================================================
implicit none
!
if ( allocated(SimpleMix) ) then
!      deallocate( SimpleMix )
endif
!
end subroutine delSimpleMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calSimpleMixing_r( list_q)
!  ===================================================================

!!xDx!!
use printDataModule, only : printDataReal, printDataInteger,                &
		       printDataPointer,printDataRealArray

implicit none
!
type ( MixListRealStruct ), target :: list_q
!
integer(kind=IntKind) :: j, idq, vlen
integer(4) :: iter
!
real(kind=RealKind) :: alpha, beta
real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:)
!
type ( MixListRealStruct ), pointer :: plq_tmp
!
plq_tmp => list_q
do idq = 1,NumTotalMix
vlen = plq_tmp%size
if ( .not.SimpleMix(idq)%Initialized ) then
 call setSimpleSpace(idq,vlen)
endif
!
SimpleMix(idq)%vector_old_r => plq_tmp%vector_old(1:vlen)
SimpleMix(idq)%vector_new_r => plq_tmp%vector_new(1:vlen)
!
alpha = SimpleMix(idq)%alpha
beta = ONE-alpha
if ( beta < zero ) beta = ZERO
!
vlen = SimpleMix(idq)%vlen
pvect_old => SimpleMix(idq)%vector_old_r(1:vlen)
pvect_new => SimpleMix(idq)%vector_new_r(1:vlen)
do j = 1,vlen
 pvect_new(j) = alpha*pvect_new(j) + beta*pvect_old(j)
enddo
plq_tmp => plq_tmp%next
enddo
!
!!ciDic!!
!!ciDic!!!!iDi!! Debugging
!!ciDic!!             iter = int(iter_count,4)
!!ciDic!!             do j=0,NumPEs-1 
!!ciDic!!              if ( j .eq. MyPE ) then
!!ciDic!!               if (iter_count .eq. 1)then  
!!ciDic!!                  call printDataInteger('pvect_new(:), iter',iter,MYPE,1)
!!ciDic!!               else
!!ciDic!!                  call printDataInteger('pvect_new(:), iter',iter,MYPE,iter_count)
!!ciDic!!               endif
!!ciDic!!                  call printDataRealArray('pvect_new(:)',pvect_new(1:vlen),vlen,MYPE, 2) 
!!ciDic!!               if (iter_count .eq. 1)then  
!!ciDic!!                  call printDataInteger('pvect_old(:), iter',iter,MYPE,2)
!!ciDic!!               else
!!ciDic!!                  call printDataInteger('pvect_old(:), iter',iter,MYPE,iter_count)
!!ciDic!!               endif
!!ciDic!!                  call printDataRealArray('pvect_old(:)',pvect_old(1:vlen),vlen,MYPE, 2) 
!!ciDic!!              endif 
!!ciDic!!             enddo 
!!ciDic!!!!iDi!! Debugging
!!ciDic!!
!
nullify(plq_tmp)
!
end subroutine calSimpleMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calSimpleMixing_c( list_q)
!  ===================================================================
implicit none
!
type ( MixListCmplxStruct ), target :: list_q
!
integer(kind=IntKind) :: j, idq, vlen 
!
real(kind=RealKind) :: alpha, beta
complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:)
!
type ( MixListCmplxStruct ), pointer :: plq_tmp
!
plq_tmp => list_q
do idq = 1,NumTotalMix
vlen = plq_tmp%size
if ( .not.SimpleMix(idq)%Initialized ) then
 call setSimpleSpace(idq,vlen)
endif
!
SimpleMix(idq)%vector_old_c => plq_tmp%vector_old(1:vlen)
SimpleMix(idq)%vector_new_c => plq_tmp%vector_new(1:vlen)
!
alpha = SimpleMix(idq)%alpha
beta = ONE-alpha
if ( beta < zero ) beta = ZERO
!
vlen = SimpleMix(idq)%vlen
pvect_old => SimpleMix(idq)%vector_old_c(1:vlen)
pvect_new => SimpleMix(idq)%vector_new_c(1:vlen)
do j = 1,vlen
 pvect_new(j) = alpha*pvect_new(j) + beta*pvect_old(j)
enddo
plq_tmp => plq_tmp%next
enddo
!
nullify(plq_tmp)
!
end subroutine calSimpleMixing_c

!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setDGAMixing( idt, idq, amix )
!  ==================================================================
implicit none
!
integer(kind=IntKind), intent(in) :: idt, idq
!  integer(kind=IntKind), optional   :: ndgaiter
!
integer(kind=IntKind) :: id
!
real(kind=RealKind), intent(in) :: amix
!
if ( .not.Initialized ) then
call ErrorHandler('setDGAMix',                                 &
  'The MixingModule needs to be initialized first')
endif
!
if ( idt<0 .or. idt>NumTypes ) then
call ErrorHandler('setDGAMixing',"Invalid type index",idt)
endif
if ( idq<0 .or. idq>NumQuantities(idt) ) then
call ErrorHandler('setDGAMixing',"Invalid quantity index",idq)
endif
!
if ( .not.allocated(DGAMix) ) then
allocate( DGAMix(NumTotalMix) )
do id = 1,NumTotalMix
 DGAMix(id)%Initialized = .false.
enddo
endif
!
NumDGAIter = MaxNumDGAIter
!  if ( present(ndgaiter) ) then
!     if ( ndgaiter > MaxNumDGAIter ) then
!        call WarningHandler('setDGAMixing', &
!            'The number of DGA iterations exeedes the maximum limit')
!     endif 
!     NumDGAIter = ndgaiter
!  endif
!
id = getMixID(idt,idq)
if ( .not.DGAMix(id)%Initialized ) then
DGAMix(id)%alpha = amix
endif
!
if ( MixingMethod /=2 ) then
MixingMethod = 2
endif
!
end subroutine setDGAMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine delDGAMixing()
!  ===================================================================
implicit none
!
integer(kind=IntKind) :: iq
!
if ( allocated(DGAMix) ) then
do iq = 1,NumTotalMix
 if ( .not.DGAMix(iq)%Initialized ) cycle
 DGAMix(iq)%Initialized = .false.
 if ( associated(DGAMix(iq)%f_old_r) ) deallocate( DGAMix(iq)%f_old_r )
 if ( associated(DGAMix(iq)%f_new_r) ) deallocate( DGAMix(iq)%f_new_r )
 if ( associated(DGAMix(iq)%f_old_c) ) deallocate( DGAMix(iq)%f_old_c )
 if ( associated(DGAMix(iq)%f_new_c) ) deallocate( DGAMix(iq)%f_new_c )
enddo 
!      deallocate( DGAMix )
endif
!
end subroutine delDGAMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calDGAMixing_r( list_q )
!  ===================================================================
!
implicit none
!
type(MixListRealStruct), target :: list_q
!
integer(kind=IntKind) :: i, ii, jj, n, iter, isz, vlen
integer(kind=IntKind) :: idt, idq, id
!
real(kind=RealKind) :: alpha
real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:) , pb(:)
real(kind=RealKind), pointer :: pf_old(:,:), pf_new(:,:)
real(kind=RealKind), pointer :: pfi_t(:), pfo_t(:)
!
type(MixListRealStruct), pointer :: plq_tmp
!
iter = min( iter_count, NumDGAIter )
n    = mod( iter_count-1, NumDGAIter ) + 1
!
if ( .not.InitializedWkSpace ) then
call setWorkingSpace_r(MixingMethod)
endif
!
plq_tmp => list_q
do idt = 1,NumTypes
a_r = zero
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 vlen = plq_tmp%size
 if ( .not.DGAMix(id)%Initialized ) then
    call setDGASpace_r(id,vlen)
 endif
 DGAMix(id)%vector_old_r => plq_tmp%vector_old(1:vlen)
 DGAMix(id)%vector_new_r => plq_tmp%vector_new(1:vlen)
 plq_tmp => plq_tmp%next
enddo
!
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 pvect_old => DGAMix(id)%vector_old_r(1:vlen)
 pf_old    => DGAMix(id)%f_old_r(1:NumDGAIter,1:vlen)
 pvect_new => DGAMix(id)%vector_new_r(1:vlen)
 pf_new    => DGAMix(id)%f_new_r(1:NumDGAIter,1:vlen)
 vlen = DGAMix(id)%vlen
 do i = 1,vlen
    pf_old(n,i)  = pvect_old(i)
    pf_new(n,i) = pvect_new(i)
 enddo
!
 do jj = 1,iter
    do ii = 1,jj
!              -------------------------------------------------------
       a_r(ii,jj) = a_r(ii,jj) +                               &
		  dgarms_r(pf_old,pf_new,ii,jj,vlen)
!              -------------------------------------------------------
       if ( ii.ne.jj ) then
	  a_r(jj,ii) = a_r(jj,ii) + a_r(ii,jj)
       endif
    enddo
 enddo
!
enddo
!     ================================================================
!     sum the a(i,j) matrix over all other atoms......................
!     ================================================================
!     ----------------------------------------------------------------
call GlobalSumInGroup(GroupID,a_r,NumDGAIter,NumDGAIter)
!     ----------------------------------------------------------------
!
isz = iter+1
pb => b_r(1:isz,1)
do i = 1,iter
 a_r(i,isz) = one
 a_r(isz,i) = one
 pb(i) = zero
enddo
a_r(isz,isz) = zero
pb(isz) = one
!
!     ================================================================
!     call the solver.................................................
!     ================================================================
!     ----------------------------------------------------------------
call dgaleq_r(a_r,b_r,isz,NumDGAIter)
!     ----------------------------------------------------------------
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 vlen = DGAMix(id)%vlen
 alpha = DGAMix(id)%alpha
!        =============================================================
!        perform the mixing...........................................
!        =============================================================
 pf_old    => DGAMix(id)%f_old_r(1:NumDGAIter,1:vlen)
 pf_new    => DGAMix(id)%f_new_r(1:NumDGAIter,1:vlen)
 pvect_new => DGAMix(id)%vector_new_r(1:vlen)
 do i = 1,vlen
    pfi_t => pf_old(1:iter,i)
    pfo_t => pf_new(1:iter,i)
!           ----------------------------------------------------------
    pvect_new(i) = dgasad_r( NumDGAIter, pfi_t, pfo_t, &
			   pb, iter, alpha )
!           ----------------------------------------------------------
 enddo
enddo
enddo
!
end subroutine calDGAMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calDGAMixing_c( list_q )
!  ===================================================================
!
implicit none
!
type(MixListCmplxStruct), target :: list_q
!
integer(kind=IntKind) :: i, ii, jj, n, iter, isz, vlen
integer(kind=IntKind) :: idt, idq, id
!
real(kind=RealKind) :: alpha
complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:) , pb(:)
complex(kind=CmplxKind), pointer :: pf_old(:,:), pf_new(:,:)
complex(kind=CmplxKind), pointer :: pfi_t(:), pfo_t(:)
!
type(MixListCmplxStruct), pointer :: plq_tmp
!
iter = min( iter_count, NumDGAIter )
n    = mod( iter_count-1, NumDGAIter ) + 1
!
if ( .not.InitializedWkSpace ) then
call setWorkingSpace_c(MixingMethod)
endif
!
plq_tmp => list_q
do idt = 1,NumTypes
a_c = czero
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 vlen = plq_tmp%size
 if ( .not.DGAMix(id)%Initialized ) then
    call setDGASpace_c(id,vlen)
 endif
 DGAMix(id)%vector_old_c => plq_tmp%vector_old(1:vlen)
 DGAMix(id)%vector_new_c => plq_tmp%vector_new(1:vlen)
 plq_tmp => plq_tmp%next
enddo
!
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 pvect_old => DGAMix(id)%vector_old_c(1:vlen)
 pf_old    => DGAMix(id)%f_old_c(1:NumDGAIter,1:vlen)
 pvect_new => DGAMix(id)%vector_new_c(1:vlen)
 pf_new    => DGAMix(id)%f_new_c(1:NumDGAIter,1:vlen)
 vlen = DGAMix(id)%vlen
 do i = 1,vlen
    pf_old(n,i) = pvect_old(i)
    pf_new(n,i) = pvect_new(i)
 enddo
!
 do jj = 1,iter
    do ii = 1,jj
!              -------------------------------------------------------
       a_c(ii,jj) = a_c(ii,jj) +                                  &
		  dgarms_c(pf_old,pf_new,ii,jj,vlen)
!              -------------------------------------------------------
       if ( ii.ne.jj ) then
	  a_c(jj,ii) = a_c(jj,ii) + a_c(ii,jj)
       endif
    enddo
 enddo
!
enddo
!     ================================================================
!     sum the a(i,j) matrix over all other atoms......................
!     ================================================================
!     ----------------------------------------------------------------
call GlobalSumInGroup(GroupID,a_c,NumDGAIter,NumDGAIter)
!     ----------------------------------------------------------------
!
isz = iter+1
pb => b_c(1:isz,1)
do i = 1,iter
 a_c(i,isz) = cone
 a_c(isz,i) = cone
 pb(i) = czero
enddo
a_c(isz,isz) = zero
pb(isz) = one
!
!     ================================================================
!     call the solver.................................................
!     ================================================================
!     ----------------------------------------------------------------
call dgaleq_c(a_c,b_c,isz,NumDGAIter)
!     ----------------------------------------------------------------
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 vlen = DGAMix(id)%vlen
 alpha = DGAMix(id)%alpha
!        =============================================================
!        perform the mixing...........................................
!        =============================================================
 pf_old    => DGAMix(id)%f_old_c(1:NumDGAIter,1:vlen)
 pf_new    => DGAMix(id)%f_new_c(1:NumDGAIter,1:vlen)
 pvect_new => DGAMix(id)%vector_new_c(1:vlen)
 do i = 1,vlen
    pfi_t => pf_old(1:iter,i)
    pfo_t => pf_new(1:iter,i)
!           ----------------------------------------------------------
    pvect_new(i) = dgasad_c( NumDGAIter, pfi_t, pfo_t, &
			   pb, iter, alpha )
!           ----------------------------------------------------------
 enddo
enddo
enddo
!
end subroutine calDGAMixing_c
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
subroutine setRPMMixing( idt, idq, amix )
!  ===================================================================
use RadialGridModule, only : getNumRmesh
use ScfDataModule, only    : inputRPMParam

implicit none
!
integer(kind=IntKind), intent(in) :: idt, idq
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
!!xDx!! removed for test_RPM_c_3 typ=14
call inputRPMParam(MaxBasis, MaxRPMIter, MaxNumRPMIter_broy)
NumAtoms=getNumAtoms()
!allocate(G_g(MaxBasis,NumAtoms*1031))
!!xDx!! Assuming same Rmesh used for all atoms.
!!???!!   allocate(G_g(MaxBasis,NumAtoms*(1+getNumRmesh(1))))
!!???!!   G_g=ZERO
!
!!xDx!! removed for test_RPM_c_3 typ=14
NumRPMIter_broy = MaxNumRPMIter_broy
!
id = getMixID(idt,idq)
if ( .not.RPMMix(id)%Initialized ) then
RPMMix(id)%alpha = amix
endif
!
!!xDx!! removed for test_RPM_c_3 typ=14
if ( MaxNumRPMIter_broy .ge. MaxRPMIter-1) then
call WarningHandler('setRPMMixing', &
   'Max #terms used for history must be less than #RPM steps-1')
MaxNumRPMIter_broy = MaxRPMIter-2
endif 
!!xDx!! removed for test_RPM_c_3 typ=14
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
!!xDx!! commented for test_RPM_c_3.F90 typ=14
! deallocate(G_g)
!!xDx!! commented for test_RPM_c_3.F90 typ=14
!!xDx!! TEST
 if ( allocated(G_g_c))then
    deallocate(G_g_c)
 endif
endif
!
end subroutine delRPMMixing
!  ===================================================================
!
!  *******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!! This is an edited version. the actual one is commented below
!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine calRPM_r(list_q) !!xDx!! argument is the potnetial.
!  ==================================================================

!!xDx!!
use printDataModule, only : printDataReal, printDataInteger,                &
		       printDataPointer,printDataRealArray
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

!
type ( MixListRealStruct ), target :: list_q
!!xDx!! Looping variables
integer (kind=IntKind) :: ii, jj, j, k, i
integer(kind=IntKind) :: idt, id, idq, vlen, iter!!xDx!!vlen=1031
integer(kind=IntKind) :: RPMIter=0 !! Used to increment Basis Size
integer(kind=IntKind) :: iter_broy=0 !! iter for Broy part of RPM
!!xDx!! Mixing Parameters 
real(kind=RealKind) :: alpha, beta
!!xDx!! Used as pointer to original variable and pointer to unstable subspace variable.
real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:), pdf(:), pf(:), A(:,:) !A is for LQ Decomposititon
real(kind=RealKind), pointer :: pvold(:), pv_st_old(:),pu(:,:), pw(:), pvt(:,:) !pvold is used as a pointer for storing the previous iteration for Q and P subspace iteration.
!!xDx!! Subspace variables
real(kind=RealKind), allocatable :: q_old(:,:), q_new(:,:) 
real(kind=RealKind), allocatable :: p_old(:), p_new(:)
real(kind=RealKind), allocatable :: z_old(:), z_new(:) 
!!xDx!! MAtrices related to basis, Projection,
real(kind=RealKind), allocatable :: ipiv(:),work(:), P(:,:),G(:,:), R(:,:) !!rDr!!, Z(:,:)
!!xDx!! Pointers
type ( MixListRealStruct ), pointer :: plq_tmp, qloc_tmp, zloc_tmp, xiloc_tmp 
!!xDx!! Variables related to Broyden part
integer (kind=IntKind) :: lastm1,nn
real (kind=RealKind) :: w0, wtmp, dfnorm, fnorm, fac2, fac1
real (kind=RealKind) :: aij, cmj, gmi
real (kind=RealKind) :: rms !!xDx!!Rms of pot. allocated in setupMixRealArrayList
real(kind=RealKind), allocatable :: aM_r(:,:), bM_r(:,:), cMM_R(:)
!!xDx!! Message passing
real (kind=RealKind) :: Ri11=0, Ri22=0
real (kind=RealKind) :: msgbuf(1:2)
real (kind=RealKind), allocatable :: msgbuf1(:) 
real (kind=RealKind), allocatable :: msgbuf2(:,:) 
real(kind=RealKind), pointer :: pad(:), pbd(:)
integer (kind=IntKind) :: INFO

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
integer (kind=IntKind) :: DESC_P( DLEN_ )
integer (kind=IntKind) :: DESC_R( DLEN_ )

!!xDx!! ScaLAPACK parameters for local Matrix definition
integer (kind=IntKind) :: M_Z, N_Z, MB_Z, NB_Z !!nDn!! Check in detail, especially M.
integer (kind=IntKind) :: M_G, N_G, MB_G, NB_G !!nDn!! Check in detail, especially M.
integer (kind=IntKind) :: M_P, N_P, MB_P, NB_P !!nDn!! Check in detail, especially M.
integer (kind=IntKind) :: M_R, N_R, MB_R, NB_R !!nDn!! Check in detail, especially M.
integer (kind=IntKind) :: CSRC_Z=0, RSRC_Z=0, IZ=1, JZ=1
integer (kind=IntKind) :: CSRC_G=0, RSRC_G=0, IG=1, JG=1
integer (kind=IntKind) :: CSRC_P=0, RSRC_P=0, IP=1, JP=1
integer (kind=IntKind) :: CSRC_R=0, RSRC_R=0, IR=1, JR=1
integer (kind=IntKind) :: ROCSRC, COCSRC, LindR, LindC
#else
!!xDx!! LAPACK Parameters
#endif
!!xDx!! PArameters for ScaLAPACK routines
integer (kind=IntKind) :: LWork
real (kind=RealKind) :: tmp 
real (kind=RealKind), allocatable :: TAU(:)
!!xDx!! Infor Parameters 
integer (kind=IntKind) :: INFO_G
integer (kind=IntKind) :: INFO_P
integer (kind=IntKind) :: INFO_R
!  ===================================================================
!!xDx!! Defining the required parallel variable.
integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup
!
iter = iter_count
plq_tmp => list_q

if ( .not.InitializedWkSpace ) then
call setWorkingSpace_r(MixingMethod)
endif
!!xDx!! we are assigning all the local pointers to the global variable here
do idt = 1,NumTypes
do idq = 1,NumQuantities(idt)
 id = getMixID(idt,idq)
 vlen = plq_tmp%size
!!cDc!!	print *, 'vlen, MyPE', vlen, MyPE
 RPMMix(id)%rms = plq_tmp%rms
 RPMMix(id)%vector_old_r => plq_tmp%vector_old(1:vlen) !! n_{in}
!!This is the potential that is used along wiht g_0 to calculate the tau matrix.
 RPMMix(id)%vector_new_r => plq_tmp%vector_new(1:vlen) !! n_{out}
!!after calculating the charge density using the contour integration of the Greens fn.
 if ( .not.RPMMix(id)%Initialized ) then
    call setRPMSpace_r(id,vlen)
 endif
 plq_tmp => plq_tmp%next
enddo                     
enddo                        


#ifdef USE_SCALAPACK
!  ===================================================================
!  Initialize ScaLAPACK and set up matrix distribution
!  ===================================================================
!!xDx!! initializing some variables
M_Z =NumAtoms*vlen  !!?D?!! What the dickens happened to the L values.
N_Z =MaxBasis
MB_Z=NumTotalMix*vlen !! Assuming NumTotalMix is the Local no. of atoms
NB_Z=MaxBasis
M_G =N_Z
N_G =M_Z
MB_G=NB_Z
NB_G=MB_Z
M_P =M_Z
N_P =M_Z
MB_P=MB_Z
NB_P=MB_Z
M_R =2 
N_R =NumAtoms*vlen !!?D?!! What the dickens happened to the L values.
MB_R=2
NB_R=NumTotalMix*vlen !! Assuming NumTotalMix is the Local no. of atoms
!!nDn!! The values have been chosen by back-calculation so each proc has what it calculates.
!
ICTXT = getGroupCommunicator(GroupID)
NumPEsInGroup = getNumPEsInGroup(GroupID)
!
call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEsInGroup )
call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)

!  -------------------------------------------------------------------
call DESCINIT( DESC_G, M_G, N_G, MB_G, NB_G, RSRC_G, CSRC_G, ICTXT, max( 1, NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW) ) , INFO_G )
call DESCINIT( DESC_P, M_P, N_P, MB_P, NB_P, RSRC_P, CSRC_P, ICTXT, max( 1, NUMROC(M_P,MB_P,MYROW,RSRC_P,NPROW) ) , INFO_P )
call DESCINIT( DESC_R, M_R, N_R, MB_R, NB_R, RSRC_R, CSRC_R, ICTXT, max( 1, NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW) ) , INFO_R )
!!xDx!! Allocating Matrices related to projection
allocate ( G( NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW), NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL) ) )
allocate ( P( NUMROC(M_P,MB_P,MYROW,RSRC_P,NPROW), NUMROC(N_P,NB_P,MYCOL,CSRC_P,NPCOL) ) )
allocate ( R( NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW), NUMROC(N_R,NB_R,MYCOL,CSRC_R,NPCOL) ) )
P=0                                                                                       
!!xDx!! Allocating local variables                                                           
allocate ( q_old( NumTotalMix, vlen ) )
allocate ( q_new( NumTotalMix, vlen ) )
allocate ( p_old( NumTotalMix*vlen ) )
allocate ( p_new( NumTotalMix*vlen ) )
allocate ( z_old( NumTotalMix*MaxBasis ) )
allocate ( z_new( NumTotalMix*MaxBasis ) )
if ( iter .eq. 1) then
allocate(G_g(MaxBasis,NumAtoms*vlen))
G_g=ZERO
endif
!!xDx!! initializing G from G_g
do ii=1,M_G
do jj=1,N_G
call infog2l ( ii, jj, DESC_G, NPROW, NPCOL, MYROW, MYCOL, LindR, LindC, ROCSRC, COCSRC) ! This gives us loc coord for global G
if ( MYROW .eq. ROCSRC .and. MYCOL .eq. COCSRC) then
G(LindR,LindC) = G_g(ii,jj)
end if
enddo
enddo
#endif

do idt = 1,NumTypes
 do idq = 1,NumQuantities(idt)
    id = getMixID(idt,idq)
    alpha = RPMMix(id)%alpha
    vlen =RPMMix(id)%vlen
    A => RPMMix(id)%A_r(1:vlen,1:2)!!qDq!!
    pf => RPMMix(id)%f_r(1:vlen)
    pdf => RPMMix(id)%df_r(1:vlen)
    pvect_old => RPMMix(id)%vector_old_r(1:vlen)
    pvect_new => RPMMix(id)%vector_new_r(1:vlen)
    pvold    => RPMMix(id)%vold_r(1:vlen) !! RPMMix(id)%vold_r stores old for G update.
    do k = 1,vlen
       p_old((id-1)*vlen+k) = pvect_old(k)
       p_new((id-1)*vlen+k) = pvect_new(k)
	if (MyPE .eq. 0) then
!!cDc!!	 print *, 'k,pvect_old',k,pvect_old(k)
	endif
    enddo
 enddo
enddo

!!xDx!! Incrementing iter counters----------------------------------------------------------------------------------------------------------------------------
    RPMIter=RPMIter+1 !! Incrementing RPMIter
            if (RPMIter .ge. MaxRPMIter) then!D!Only after the first MaxRPMIter we generate an unstable subspace
            iter_broy= mod( RPMIter, MaxRPMIter )!D!After every MaxRPMIter steps. we get a new subspace
            endif

!!xDx!! Defining corresponding corresponding subspace variables which will be passed -------------------------------------------------------------------------
#ifdef USE_SCALAPACK
            call PDGEMM( 'T', 'N', M_P, N_P, N_Z, 1.0D0, G, IG, JG,     &  !! N_Z=M_G
            DESC_G, G, IG, JG, DESC_G, 0.0D0, P, IP, JP, DESC_P )
!!?D?!! Wouldn't it be better for large calculations to store P instead of G?
#else
            P = matmul(transpose(G),G)!!nDn!! Equivalent Lapack version of the above needs to be written.
#endif



      !!xDx!! Local Multiplication of Pu----------------------------------------------------------------------------------------------------------------------------
            allocate(msgbuf1(1:M_P))
            msgbuf1=matmul(P,p_old)!! p<-Pu
            call GlobalSumInGroup(GroupID, msgbuf1, M_P)
            !!nDn!! ith proc should store the ith block.
            p_old(1:NumTotalMix*vlen)=                                  &
            msgbuf1(MYCOL*NumTotalMix*vlen+1:(MYCOL+1)*NumTotalMix*vlen)
            deallocate(msgbuf1)
            allocate(msgbuf1(1:M_P))
            msgbuf1=matmul(P,p_new)!! p<-Pu
            call GlobalSumInGroup(GroupID, msgbuf1, M_P)
            p_new(1:NumTotalMix*vlen)=                                  &
            msgbuf1(MYCOL*NumTotalMix*vlen+1:(MYCOL+1)*NumTotalMix*vlen)
            deallocate(msgbuf1)
      !!xDx!! Allocating local coordinates for the unstable subspace !!iDi!!
            allocate(msgbuf1(1:MaxBasis))                            !!iDi!!
            msgbuf1=matmul(G,p_old)!! p<-Pu                          !!iDi!!
            call GlobalSumInGroup(GroupID, msgbuf1, MaxBasis)        !!iDi!!
            z_old=msgbuf1                                            !!iDi!!
            deallocate(msgbuf1)                                      !!iDi!!
            allocate(msgbuf1(1:MaxBasis))                            !!iDi!!
            msgbuf1=matmul(G,p_new)!! p<-Pu                          !!iDi!!
            call GlobalSumInGroup(GroupID, msgbuf1, MaxBasis)        !!iDi!!
            z_new=msgbuf1                                            !!iDi!!
            deallocate(msgbuf1)                                      !!iDi!!
            deallocate(P) !! we dont use it later.
!!nDn!! We can actually use MPI_REDUCE with SUM as the operation.

      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
            vlen =RPMMix(id)%vlen
!!cDc!!		print *, 'vlen, MyPE', vlen, MyPE
            A => RPMMix(id)%A_r(1:vlen,1:2)!!qDq!!
            pf => RPMMix(id)%f_r(1:vlen)
            pdf => RPMMix(id)%df_r(1:vlen)
            pvect_old => RPMMix(id)%vector_old_r(1:vlen)
            pvect_new => RPMMix(id)%vector_new_r(1:vlen)
            pvold    => RPMMix(id)%vold_r(1:vlen) !! RPMMix(id)%vold_r stores old for G update???
      !xDx!------ Defining the stable subspace projection
            q_old(id,1:vlen) =pvect_old(1:vlen)-p_old((id-1)*vlen+1:id*vlen)!! q<- u-p
            q_new(id,1:vlen) =pvect_new(1:vlen)-p_new((id-1)*vlen+1:id*vlen)!! q<- u-p
      !xDx!------ Defining the unstable subspace projection as pvect
            pvect_old=ZERO                                               !!iDi!!
            pvect_new=ZERO                                               !!iDi!!
            if ( id .eq. 1 )then                           !!iDi!!
                 pvect_old(1:MaxBasis)=z_old( 1:MaxBasis ) !!iDi!!
                 pvect_new(1:MaxBasis)=z_new( 1:MaxBasis ) !!iDi!!
            endif                                          !!iDi!!
!!xDx!! Performing Simple Mixing ----------------------------------------------------------------------------------------------------------------------------
            if ( alpha .ge. ONE ) then
            q_new(id,1:vlen)=alpha*q_new(id,1:vlen)
            else
            q_new(id,1:vlen)=alpha*q_new(id,1:vlen)+(1-alpha)*q_old(id,1:vlen)
            endif
         enddo
      enddo


!!xDx!! Performing Broyden Mixing ---------------------------------------------------------------------------------------------------------------------------
!!nDn!! pvect_{old,new} contains the unstable subspace.

!!riDri!! SNGLPROCBROY:    if ( MYROW .eq. 0 .and. MYCOL .eq. 0 ) then !!iDi!!
   if ( iter_broy .eq. 1) then

!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
            vlen = MaxBasis !?!RPMMix(id)%vlen
            pf => RPMMix(id)%f_r(1:vlen)
            pvect_old => RPMMix(id)%vector_old_r(1:vlen)
            pvect_new => RPMMix(id)%vector_new_r(1:vlen) !!nDn!! Note that this pvect_. is now the variable projected on the unstable subspace.
            pvold    => RPMMix(id)%vold_r(1:vlen) !! BroydenMix(id)%vold_r used to store pvold for P-subspace
            !!xDx!! everytime the broyden mixing, we allocate everything to zero.
               pu    => RPMMix(id)%u_r(1:vlen,1:NumRPMIter_broy)!!xDx!! u as in eqn (13b.iii)
               pvt   => RPMMix(id)%vt_r(1:vlen,1:NumRPMIter_broy)!!*D*!!
               pdf   => RPMMix(id)%df_r(1:vlen)
               pvold = ZERO
               pu    = ZERO
               pvt   = ZERO
               pdf   = ZERO
!
            do k = 1,vlen
               pf(k) = pvect_new(k) - pvect_old(k) !!cDt!!
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k) !!cDt!!
            enddo
!
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k) !!cDt!!
            enddo
         enddo
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
      if (  iter_broy > NumRPMIter_broy ) then !!xDx!!=MaxNumRPMIter(=4) 
         nn = NumRPMIter_broy
      else
         nn = lastIter
      endif
!
   LoopType: do idt = 1,NumTypes
         w0= ten2m2
            dfnorm = zero
            fnorm = zero
            do idq = 1,NumQuantities(idt)
               id = getMixID(idt,idq)
               vlen = MaxBasis  !!??!! RPMMix(id)%vlen
               alpha = RPMMix(id)%alpha
               pvect_old => RPMMix(id)%vector_old_r(1:vlen)
               pvect_new => RPMMix(id)%vector_new_r(1:vlen)
               pf  => RPMMix(id)%f_r(1:vlen) !!This is F^(iter-1):=F[n^(iter-1)]=n_{out}[n^(iter-1)]-n^(iter-1).
            !!nDn!! based on what follows, we want to use pvect_{new,old} as temporary placeholders for p_{new,old}, their projected counterparts.
               pdf => RPMMix(id)%df_r(1:vlen)
               if ( vlen==1 ) then
                  pvect_new(1) = pvect_old(1) + alpha * pf(1) !!xDx!! ??? why ???
                  cycle !!cDt!!
               endif
               do k = 1,vlen !!cDt!!x2
                  pdf(k) = pvect_new(k) - pvect_old(k) - pf(k)!! dR^(iter)=R^(iter)-R^(iter-1)
                  pf(k)  = pvect_new(k) - pvect_old(k) !!R^(iter)
               enddo
!              
               do k = 1,vlen  !!xDx!! These are sums over all the atoms in the proc. 
                  dfnorm = dfnorm + pdf(k)*pdf(k)!!cDt!!
                  fnorm  = fnorm  + pf(k)*pf(k)!!cDt!!
               enddo
            enddo
! 
            msgbuf(1) = dfnorm!!cDt!!
            msgbuf(2) = fnorm!!cDt!!
!           ------------------------------------------------------------
            dfnorm = sqrt( msgbuf(1) )!!cDt!!
            fnorm  = sqrt( msgbuf(2) )!!cDt!!
!           ============================================================
!
            fac2 = one/dfnorm!!cDt!!
            do idq = 1,NumQuantities(idt)
               id = getMixID(idt,idq)
               alpha = RPMMix(id)%alpha
               vlen  = MaxBasis  !!??!! RPMMix(id)%vlen
               if ( vlen==1 ) then
                  cycle
               endif
               fac1 = alpha*fac2
               pvect_old => RPMMix(id)%vector_old_r(1:vlen)
               pvect_new => RPMMix(id)%vector_new_r(1:vlen)
               pu    => RPMMix(id)%u_r(1:vlen,1:NumRPMIter_broy)!!xDx!! u as in eqn (13b.iii)
               pvt   => RPMMix(id)%vt_r(1:vlen,1:NumRPMIter_broy)!!*D*!!
               pdf   => RPMMix(id)%df_r(1:vlen)
               pvold => RPMMix(id)%vold_r(1:vlen) !! BroydenMix(id)%vold_r used for BRoyden method
               do k = 1,vlen!!cDt!!x3
                  pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k)) !! aI(dR/|dR|)+(p^(2)-p^(1))/|dR|, Eqn (13b)
                  pvold(k)   = pvect_old(k)
                  pvect_old(k) = fac2*pdf(k)  !!=(dR/|dR|) storeds the part of the unit vector of dF.
               enddo  
!              ---------------------------------------------------------
               call broy_sav_r( pu, pvt, pvect_old, pvect_new, iter_broy-1, & 
                                NumRPMIter_broy, vlen )!!cDt!!
            enddo
!           ============================================================
            !!xDx!! allocating the Matrices aM_r and bM_r=aM_r^-1
            allocate (aM_r(1:nn,1:nn))
            allocate (cMM_r(1:nn))
            allocate (bM_r(1:nn,1:nn))
!
            !!xDx!! nn -> #iters-1, equiv to m in DD Johnson.
            do j = 1,nn - 1          ! off diagonal elements of a(i,j)
               do i = j+1,nn
                  aij = zero
                  do idq = 1,NumQuantities(idt)
                     id = getMixID(idt,idq)
                     vlen = MaxBasis  !!??!! RPMMix(id)%vlen
                     if ( vlen==1 ) then
                        cycle
                     endif
                     pvt => RPMMix(id)%vt_r(1:vlen,1:NumRPMIter_broy)
                     do k = 1,vlen
                        aij = aij + pvt(k,j)*pvt(k,i)!!cDt!!
                     enddo
                  enddo
                  aM_r(i,j) = aij
                  aM_r(j,i) = aij
               enddo
            enddo
!
            do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
               aij = zero
               cmj = zero
               do idq = 1, NumQuantities(idt)
                  id = getMixID(idt,idq)
                  vlen = MaxBasis  !!??!! RPMMix(id)%vlen
                  if ( vlen==1 ) then
                     cycle
                  endif
                  pvt => RPMMix(id)%vt_r(1:vlen,1:NumRPMIter_broy)
                  pf  => RPMMix(id)%f_r(1:vlen)
                  do k = 1,vlen!!cDt!!x2
                     cmj = cmj + pvt(k,i)*pf(k) !!xDx!! c^m_k/w_k in Eqn (15b). Eavaluated for fixed m, ie, for the current iteration.
                     aij = aij + pvt(k,i)*pvt(k,i)
                  enddo
               enddo
               aM_r(i,i) = aij
               cMM_r(i) = cmj
            enddo
!
!           ============================================================
            do idq = 1,NumQuantities(idt)
!
               id = getMixID(idt,idq)
               vlen = MaxBasis  !!??!! RPMMix(id)%vlen
               if ( vlen==1 ) then
                  cycle
               endif
               alpha = RPMMix(id)%alpha
               pvect_old => RPMMix(id)%vector_old_r(1:vlen)
               pvect_new => RPMMix(id)%vector_new_r(1:vlen)
!
               pvold => RPMMix(id)%vold_r(1:vlen)
               pw => RPMMix(id)%w_r(1:NumRPMIter_broy) !!xDx!! Weights,
               pu => RPMMix(id)%u_r(1:vlen,1:NumRPMIter_broy)
               pf => RPMMix(id)%f_r(1:vlen)
!
               if (  iter_broy-1 > NumRPMIter_broy ) then
                  do i = 1,NumRPMIter_broy-1
                     pw(i) = pw(i+1)
                  enddo
               endif
!
               rms = RPMMix(id)%rms
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
!              =========================================================
!                 ======================================================
!                 use lapack !!nDn!! everything is local now.
!                 ======================================================

                  do i = 1,nn
                     do j = 1,nn
                        bM_r(j,i) = aM_r(j,i)*pw(j)*pw(i)!!cDt!!
                     enddo
                     bM_r(i,i) = w0**2 + aM_r(i,i)*pw(i)*pw(i)!!cDt!!
                  enddo
                  allocate( ipiv(1:nn) )
!                 ------------------------------------------------------
                  call dgetrf( nn, nn, bM_r, nn,ipiv, info ) !!cDt!!
!                 ------------------------------------------------------
                  call dgetri( nn, bM_r, nn, ipiv, tmp, -1, info )!!cDt!!
!                 ------------------------------------------------------
                  LWork = int(real(tmp,kind=RealKind))
                  allocate( Work(1:LWork) )
!                 ------------------------------------------------------
                  call dgetri( nn, bM_r, nn, ipiv, Work, LWork, info )!!cDt!!
!                 ------------------------------------------------------
                  deallocate (Work)
!                 ------------------------------------------------------
                  deallocate (ipiv)
!                 ------------------------------------------------------
!
               do k = 1,vlen
                  pvect_new(k) = pvold(k) + alpha*pf(k)!!cDt!!
               enddo
               do i = 1,nn
                  gmi = zero
                  do j = 1,nn!!cDt!!
                     gmi = gmi + cMM_r(j)*bM_r(j,i)*pw(j) !!xDx!! \gamma_{ml} in Eqn. (15b), the cMM_r is wihtout w_j
                  enddo
                  do k = 1,vlen
                     pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)!!cDt!!
                  enddo
               enddo
            !!xDx!! temporarily storing the variables in z_new and z_old.
               if (idq .eq. 1) then                                                       !!iDi!!
                 z_new(1:MaxBasis) = pvect_new(1:MaxBasis)                                  !!iDi!!
                 z_old(1:MaxBasis) = pvect_old(1:MaxBasis)                                  !!iDi!!
               endif                                                                      !!iDi!!
            enddo
            !!xDx!! Deallocating the MAtrices aM_r and bM_r=aM_r^-1
            deallocate(aM_r)
            deallocate(bM_r)
            deallocate(cMM_r)

      enddo LoopType

   endif   

!!xDx!! End of Broyden Mixing -------------------------------------------------------------------------------------------------------------------------------
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
            vlen =RPMMix(id)%vlen
            A => RPMMix(id)%A_r(1:vlen,1:2)
            pf => RPMMix(id)%f_r(1:vlen)
            pdf => RPMMix(id)%df_r(1:vlen)
            pvect_old => RPMMix(id)%vector_old_r(1:vlen)
            pvect_new => RPMMix(id)%vector_new_r(1:vlen)
            pvold    => RPMMix(id)%vold_r(1:vlen)
            pv_st_old    => RPMMix(id)%v_st_old_r(1:vlen)
         !!xDx!! projecting local back to global variables
            pvect_new(1:vlen)= matmul( Transpose( G(1:MaxBasis,(id-1)*vlen+1:id*vlen) ), z_new(1:MaxBasis) ) !!iDi!!
            pvect_old(1:vlen)= matmul( Transpose( G(1:MaxBasis,(id-1)*vlen+1:id*vlen) ), z_old(1:MaxBasis) ) !!iDi!!
         !!xDx!! Combining the resultant to return to main -------------------------------------------------------------------------------------------------------
!            RPMMix(id)%vector_new_r(1:vlen)=q_new(id,1:vlen)+pvect_new!! u<- q+p
            do k=1,vlen
               pvect_new(k)=q_new(id,k)+pvect_new(k)!! u<- q+p
            enddo
!
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
               R( ii, (id-1)*vlen+jj ) =A( jj, ii )
               enddo
               enddo
            endif

         enddo
      enddo
!                                                   
!!xDx!! Performing the calc to increase the Basis size                                                                                                          
      if ( mod( RPMIter, MaxRPMIter ) .eq. 0 ) then !!nDn!! iter_broy= 0, RPMIter> MaxRPMIter

      !!xDx!! HouseHolder LQ---------------------------------------------------------------- 
#ifdef USE_SCALAPACK
            !------- Allocating Tau -------------------------------------
                  allocate( TAU( 1:NUMROC(IR+min(M_R,N_R)-1, MB_R , MYROW, 0, NPROW) ) )

            !----------------------- Work Place Qwery ----------------------------
                  call PDGELQF( M_R , N_R, R, IR, JR, DESC_R, TAU, tmp, -1, INFO_R )
            !-------------------------------------------------------------------
            LWork = int(tmp)
            allocate( WORK(1:LWork) )
           
            !  ------------------------- Actual QR  ------------------------------
                  call PDGELQF( M_R, N_R, R, IR, JR, DESC_R, TAU, WORK, LWork, INFO_R )
            !  -------------------------------------------------------------------
           
            deallocate (WORK)

       !!xDx!! Storing R_{11} and R_{22}
            Ri11=0
            Ri22=0
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
                  call PDORGLQ( M_R, N_R, M_R, R, IR, JR, DESC_R, TAU, tmp, -1, INFO_R )
            !  -------------------------------------------------------------------
            LWork =int(tmp)
            allocate( WORK(1:LWork) )
            !  ------------------------- Actual Q  -------------------------------
                  call PDORGLQ( M_R, N_R, M_R, R, IR, JR, DESC_R, TAU, WORK, LWork, INFO_R )
            !  -------------------------------------------------------------------
               deallocate (WORK)
               deallocate (TAU)
!!nDn!! write the LAPACK version.
#endif
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
#ifdef USE_SCALAPACK
            !------- Allocating Tau -------------------------------------
                  allocate( TAU( 1:NUMROC(IG+min(M_G,N_G)-1, MB_G, MYROW,0,NPROW) ) )
            !----------------------- Work Place Qwery ----------------------------
                  call PDGELQF( M_G , N_G, G, IG, JG, DESC_G, TAU, tmp, -1, INFO_G )
            !-------------------------------------------------------------------
               LWork = int(tmp)
               allocate( WORK(1:LWork) )
            !  ------------------------- Actual QR  ------------------------------
                  call PDGELQF( M_G , N_G, G, IG, JG, DESC_G, TAU, WORK, LWork, INFO_G )
            !  -------------------------------------------------------------------
               deallocate (WORK)
            !  ------------------------- Work Place Qwery  -----------------------
                  call PDORGLQ( M_G , N_G, M_G, G, IG, JG, DESC_G, TAU, tmp, -1, INFO_G )
            !  -------------------------------------------------------------------
               LWork =int(tmp)
               allocate( WORK(1:LWork) )
            !  ------------------------- Actual Q  -------------------------------
                  call PDORGLQ( M_G , N_G, M_G, G, IG, JG, DESC_G, TAU, WORK, LWork, INFO_G )
            !  -------------------------------------------------------------------
               deallocate (WORK)
               deallocate (TAU)

            !!xDx!! Re-Allocating G_g from G
                G_g=0
                do ii=1,NUMROC(M_G,MB_G,MYROW,0,NPROW)
                do jj=1,NUMROC(N_G,NB_G,MYCOL,0,NPCOL)
                G_g( indxl2g( ii, MB_G, MYROW, 0, NPROW ) , indxl2g( jj, NB_G, MYCOL, 0, NPCOL ) ) = G(ii,jj)
                enddo
                enddo
            !!xDx!! Global sum to transfer the whole matrix
                allocate(msgbuf2(1:M_G,1:N_G))
                msgbuf2 = G_g
                call GlobalSumInGroup(GroupID, msgbuf2, M_G, N_G)
                G_g = msgbuf2
                deallocate(msgbuf2)
#endif

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
end subroutine calRPM_r
!  ===================================================================
!
!  *******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calRPM_c(list_q) !!xDx!! argument is the potnetial.
!  ==================================================================

!!xDx!! removed for test_RPM_c_3 typ=14
!!xDx!!
!   use printDataModule, only : printDataReal, printDataInteger,                &
!                               printDataPointer,printDataRealArray,            &
!                               printDataComplex,printDataComplexArray           
!!xDx!!
!!xDx!! removed for test_RPM_c_3 typ=14

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

!
   type ( MixListCmplxStruct ), target :: list_q
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

!!xDx!! removed for test_RPM_c_3 typ=14
!NumAtoms=4
!MaxBasis=6
!MaxRPMIter=6
!NumRPMIter_broy=MaxRPMIter-2
!!xDx!! removed for test_RPM_c_3 typ=14

   iter = iter_count
   plq_tmp => list_q

   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_c(MixingMethod)
   endif
!!xDx!! we are assigning all the local pointers to the global variable here
   do idt = 1,NumTypes
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         RPMMix(id)%rms = plq_tmp%rms
         RPMMix(id)%vector_old_c => plq_tmp%vector_old(1:vlen) !! n_{in}
         RPMMix(id)%vector_new_c => plq_tmp%vector_new(1:vlen) !! n_{out}
         if ( .not.RPMMix(id)%Initialized ) then
            call setRPMSpace_c(id,vlen,MaxBasis)
         endif
         plq_tmp => plq_tmp%next
      enddo                     
   enddo                        

#ifdef USE_SCALAPACK
!  ===================================================================
!  Initialize ScaLAPACK and set up matrix distribution
!  ===================================================================
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
   call BLACS_GRIDINIT( ICTXT, 'R', 1, NumPEsInGroup )
   call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)
!  -------------------------------------------------------------------
   call DESCINIT( DESC_G, M_G, N_G, MB_G, NB_G, RSRC_G, CSRC_G, ICTXT, max( 1, NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW) ) , INFO_G )
   call DESCINIT( DESC_R, M_R, N_R, MB_R, NB_R, RSRC_R, CSRC_R, ICTXT, max( 1, NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW) ) , INFO_R )
!!xDx!! Allocating Matrices related to projection
   allocate ( G( NUMROC(M_G,MB_G,MYROW,RSRC_G,NPROW), NUMROC(N_G,NB_G,MYCOL,CSRC_G,NPCOL) ) )
   allocate ( R( NUMROC(M_R,MB_R,MYROW,RSRC_R,NPROW), NUMROC(N_R,NB_R,MYCOL,CSRC_R,NPCOL) ) )
!!xDx!! Allocating local variables                                                           
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
  q_old = CZERO
  q_new = CZERO
  p_old = CZERO
  p_new = CZERO
  z_old = CZERO
  z_new = CZERO
  R = CZERO
  G = CZERO
 
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

      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
            vlen =RPMMix(id)%vlen
!!xDx!! commented for test_RPM_c_3 typ=14
!            A => RPMMix(id)%A_c(1:vlen,1:2)!!qDq!!
!            pf => RPMMix(id)%f_c(1:vlen)
!            pdf => RPMMix(id)%df_c(1:vlen)
!            pvold    => RPMMix(id)%vold_c(1:vlen) !! RPMMix(id)%vold_c stores old for G update.
!!xDx!! commented for test_RPM_c_3 typ=14
            pvect_old => RPMMix(id)%vector_old_c(1:vlen)
            pvect_new => RPMMix(id)%vector_new_c(1:vlen)
            do k = 1,vlen
               p_old((id-1)*vlen+k) = pvect_old(k)
               p_new((id-1)*vlen+k) = pvect_new(k)
            enddo

!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (MyPE+10,1122) '========================'
!!cDc!!,       write (MyPE+10,1342) 'pvect_new(:)','MyPE,RPMIter, id', MyPE, RPMIter+1, id
!!cDc!!,       write (MyPE+10,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+10,1242) real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+10,1122) '========================'
!!cDc!!,       write (MyPE+10,1342) 'pvect_old(:)','MyPE,RPMIter, id', MyPE, RPMIter+1, id
!!cDc!!,       write (MyPE+10,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+10,1242) real(pvect_old(i)),aimag(pvect_old(i)),'*i'
!!cDc!!,       enddo
!!cDc!!, 1122 format(X,A)
!!cDc!!, 1342 format(2(X,A),3(X,I4))
!!cDc!!, 1241 format(2(X,A),2(X,I4),(X,F20.10,SP,F16.10,A))
!!cDc!!,  1242 format((X,F20.10,SP,F16.10,A))
!!cDc!!,    endif
!!cDc!!, enddo
!!cDc!!, !!cDc!!

         enddo
      enddo

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
!!cDc!!, !!cDc!!
!!cDc!!,    if (MyPE .eq. 0) then
!!cDc!!,       write (MyPE+12,1122) '========================'
!!cDc!!,       write (MyPE+12,1342) 'b4 z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+12,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+12,1241) 'i,MyPE','z_new',i,MyPE,real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, !!cDc!!

      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
            vlen =RPMMix(id)%vlen
!!xDx!! commented for test_RPM_c_3 typ=14
!            A => RPMMix(id)%A_c(1:vlen,1:2)!!qDq!!
!            pf => RPMMix(id)%f_c(1:vlen)
!            pdf => RPMMix(id)%df_c(1:vlen)
!            pvold    => RPMMix(id)%vold_c(1:vlen) !! RPMMix(id)%vold_c stores old for G update??
!!xDx!! commented for test_RPM_c_3 typ=14
            pvect_old => RPMMix(id)%vector_old_c(1:vlen)
            pvect_new => RPMMix(id)%vector_new_c(1:vlen)
      !xDx!------ Defining the stable subspace projection
            q_old(id,1:vlen) =pvect_old(1:vlen)-p_old((id-1)*vlen+1:id*vlen)!! q<- u-p
            q_new(id,1:vlen) =pvect_new(1:vlen)-p_new((id-1)*vlen+1:id*vlen)!! q<- u-p

!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (MyPE+14,1122) '========================'
!!cDc!!,       write (MyPE+14,1342) 'after q_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+14,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+14,1241) 'i,MyPE','q_new',i,MyPE,real(q_new(id,i)),aimag(q_new(id,i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+14,1122) '========================'
!!cDc!!,       write (MyPE+14,1342) 'after q_old(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+14,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+14,1241) 'i,MyPE','q_old',i,MyPE,real(q_old(id,i)),aimag(q_old(id,i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, enddo
!!cDc!!, !!cDc!!

!!cDc!!      !xDx!------ Defining the unstable subspace projection as pvect
!!cDc!!            pvect_old=CZERO                                               !!iDi!!
!!cDc!!            pvect_new=CZERO                                               !!iDi!!
!!cDc!!            if ( id .eq. 1 )then                           !!iDi!!
!!cDc!!                 pvect_old(1:MaxBasis)=z_old( 1:MaxBasis ) !!iDi!!
!!cDc!!                 pvect_new(1:MaxBasis)=z_new( 1:MaxBasis ) !!iDi!!
!!cDc!!            endif                                          !!iDi!!
      !xDx! Performing Simple Mixing ----------------------------------------------------------------
            if ( alpha .ge. ONE ) then
            q_new(id,1:vlen)=alpha*q_new(id,1:vlen)
            else
            q_new(id,1:vlen)=alpha*q_new(id,1:vlen)+(1-alpha)*q_old(id,1:vlen)
            endif
         enddo
      enddo
!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!, do id=1,2
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (MyPE+16,1122) '========================'
!!cDc!!,       write (MyPE+16,1342) 'after q_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+16,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+16,1241) 'i,MyPE','q_new',i,MyPE,real(q_new(id,i)),aimag(q_new(id,i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, enddo
!!cDc!!, enddo
!!cDc!!, !!cDc!!

!!xDx!! Performing Broyden Mixing ---------------------------------------------------------------------------------------------------------------------------
!!nDn!! pvect_{old,new} contains the unstable subspace.

!!riDri!! SNGLPROCBROY:    if ( MYROW .eq. 0 .and. MYCOL .eq. 0 ) then !!iDi!!
   if ( iter_broy .eq. 1) then

!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
!!cDc!!      do idt = 1,NumTypes
!!cDc!!         do idq = 1,NumQuantities(idt)
            id = 1 !!cDc!!getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
!!cDc!!            vlen = MaxBasis !?!RPMMix(id)%vlen
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
!!cDc!!, !!cDc!!
!!cDc!!,    if (MyPE .eq. 0) then
!!cDc!!,       write (MyPE+20,1122) '========================'
!!cDc!!,       write (MyPE+20,1342) 'b4 z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+20,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+20,1241) 'i,MyPE','z_new',i,MyPE,real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+20,1122) '========================'
!!cDc!!,       write (MyPE+20,1342) 'b4 pvold(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+20,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+20,1241) 'i,MyPE','pvold',i,MyPE,real(pvold(i)),aimag(pvold(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, !!cDc!!
            pf = z_new - z_old
            pvold = z_old
            z_new = z_old + alpha*pf
!!cDc!!, !!cDc!!
!!cDc!!,    if (MyPE .eq. 0) then
!!cDc!!,       write (MyPE+21,1122) '========================'
!!cDc!!,       write (MyPE+21,1342) 'after pf(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+21,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+21,1242) real(pf(i)),aimag(pf(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+21,1122) '========================'
!!cDc!!,       write (MyPE+21,1342) 'after z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+21,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+21,1242) real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+21,1122) '========================'
!!cDc!!,       write (MyPE+21,1342) 'after z_old(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+21,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+21,1242) real(z_old(i)),aimag(z_old(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, !!cDc!!

!
!!cDc!!         enddo
!!cDc!!      enddo
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
               id = 1 
               alpha = RPMMix(id)%alpha
               pf  => RPMMix(id)%f_c(1:MaxBasis) 
               pdf => RPMMix(id)%df_c(1:MaxBasis)
               if ( MaxBasis==1 ) then
                  z_new(1) = z_old(1) + alpha * pf(1) 
               endif
               do k = 1,MaxBasis 
                  pdf(k) = z_new(k)-z_old(k)-pf(k)
                  pf(k)  = z_new(k)-z_old(k) 
               enddo
!          

               do k = 1,MaxBasis!!xDx!!These r sums over all the atoms in the proc
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
               pu    => RPMMix(id)%u_c(1:MaxBasis,1:NumRPMIter_broy)
               pvt   => RPMMix(id)%vt_c(1:MaxBasis,1:NumRPMIter_broy)
               pdf   => RPMMix(id)%df_c(1:MaxBasis)
               pvold => RPMMix(id)%vold_c(1:MaxBasis) 
!!cDc!!, !!cDc!!
!!cDc!!,    if (MyPE .eq. 0) then
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       write (MyPE+30,1342) 'b4 pf(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+30,1242) real(pf(i)),aimag(pf(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       write (MyPE+30,1342) 'b4 z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+30,1242) real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       write (MyPE+30,1342) 'b4 z_old(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+30,1242) real(z_old(i)),aimag(z_old(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       write (MyPE+30,1342) 'b4 pvold(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+30,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+30,1242) real(pvold(i)),aimag(pvold(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, !!cDc!!
               do k = 1,MaxBasis
                  z_new(k) = fac1*pdf(k) + fac2*(z_old(k) - pvold(k)) 
                  pvold(k)   = z_old(k)
                  z_old(k) = fac2*pdf(k)  
               enddo 
!!cDc!!, !!cDc!!
!!cDc!!,    if (MyPE .eq. 0) then
!!cDc!!,       write (MyPE+31,1122) '========================'
!!cDc!!,       write (MyPE+31,1342) 'after pvold(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+31,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+31,1242) real(pvold(i)),aimag(pvold(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+31,1122) '========================'
!!cDc!!,       write (MyPE+31,1342) 'after z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+31,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+31,1242) real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,       write (MyPE+31,1122) '========================'
!!cDc!!,       write (MyPE+31,1342) 'after z_old(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+31,1122) '========================'
!!cDc!!,       do i =1,MaxBasis
!!cDc!!,          write (MyPE+31,1242) real(z_old(i)),aimag(z_old(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, !!cDc!!

!              ---------------------------------------------------------
               call broy_sav_c( pu, pvt, z_old, z_new, iter_broy-1, & 
                                NumRPMIter_broy, MaxBasis )!!cDt!!
!!cDc!!, !!cDc!!
!!cDc!!, 	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!, 	if ( MyPE .eq. 0) then
!!cDc!!, 	write (MyPE+32,1122) '                        '
!!cDc!!, 	write (MyPE+32,1122) '========================'
!!cDc!!, 	write (MyPE+32,1342) 'MyPE','pu(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+32,1122) '========================'
!!cDc!!, 	do ii=1,MaxBasis
!!cDc!!, 	   do jj=1,NumRPMIter_broy
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(pu(ii,jj)).gt.0)
!!cDc!!,               write(MyPE+32,fmt, advance="no") pu(ii,jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+32, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+32,1122) '========================'
!!cDc!!, 	write (MyPE+32,1342) 'MyPE','pvt(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+32,1122) '========================'
!!cDc!!, 	do ii=1,MaxBasis
!!cDc!!, 	   do jj=1,NumRPMIter_broy
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(pvt(ii,jj)).gt.0)
!!cDc!!,               write(MyPE+32,fmt, advance="no") pvt(ii,jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+32, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+32,1122) '                        '
!!cDc!!, 	endif
!!cDc!!, !!cDc!!


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
                     cmj = cmj + conjg(pvt(k,i))*pf(k) 
                     aij = aij + conjg(pvt(k,i))*pvt(k,i)
                  enddo
               aM_c(i,i) = aij
               cMM_c(i) = cmj
            enddo
!
!!cDc!!, !!cDc!!
!!cDc!!, 	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!, 	if ( MyPE .eq. 0) then
!!cDc!!, 	write (MyPE+33,1122) '                        '
!!cDc!!, 	write (MyPE+33,1122) '========================'
!!cDc!!, 	write (MyPE+33,1342) 'MyPE','aM_c(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+33,1122) '========================'
!!cDc!!, 	do ii=1,nn
!!cDc!!, 	   do jj=1,nn
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(aM_c(ii,jj)).gt.0)
!!cDc!!,               write(MyPE+33,fmt, advance="no") aM_c(ii,jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+33, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+33,1122) '========================'
!!cDc!!, 	write (MyPE+33,1342) 'MyPE','cMM_c(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+33,1122) '========================'
!!cDc!!, 	do ii=1,nn
!!cDc!!,            fmt(8:8) = MERGE('+',' ',imag(cMM_c(ii)).gt.0)
!!cDc!!,            write(MyPE+33,fmt, advance="no") cMM_c(ii)
!!cDc!!, 	   write(MyPE+33, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+33,1122) '                        '
!!cDc!!, 	endif
!!cDc!!, !!cDc!!
!           ============================================================
!
               id = 1 
               alpha = RPMMix(id)%alpha
!
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
               rms = dfnorm !RPMMix(id)%rms
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
!              b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1 
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
!!cDc!!, !!cDc!!
!!cDc!!, 	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!, 	if ( MyPE .eq. 0) then
!!cDc!!, 	write (MyPE+34,1122) '                        '
!!cDc!!, 	write (MyPE+34,1122) '========================'
!!cDc!!, 	write (MyPE+34,1342) 'MyPE','by inv bM_c(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+34,1122) '========================'
!!cDc!!, 	do ii=1,nn
!!cDc!!, 	   do jj=1,nn
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(bM_c(ii,jj)).gt.0)
!!cDc!!,               write(MyPE+34,fmt, advance="no") bM_c(ii,jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+34, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+34,1122) '========================'
!!cDc!!, 	write (MyPE+34,1342) 'MyPE','pw(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+34,1122) '========================'
!!cDc!!, 	do ii=1,nn
!!cDc!!,            fmt(8:8) = MERGE('+',' ',imag(pw(ii)).gt.0)
!!cDc!!,            write(MyPE+34,fmt, advance="no") pw(ii)
!!cDc!!, 	   write(MyPE+34, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+34,1122) '                        '
!!cDc!!, 	endif
!!cDc!!, !!cDc!!
                  allocate( ipiv(1:nn) )
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
!!cDc!!, !!cDc!!
!!cDc!!, 	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!, 	if ( MyPE .eq. 0) then
!!cDc!!, 	write (MyPE+35,1122) '                        '
!!cDc!!, 	write (MyPE+35,1122) '========================'
!!cDc!!, 	write (MyPE+35,1342) 'MyPE','bM_c(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+35,1122) '========================'
!!cDc!!, 	do ii=1,nn
!!cDc!!, 	   do jj=1,nn
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(bM_c(ii,jj)).gt.0)
!!cDc!!,               write(MyPE+35,fmt, advance="no") bM_c(ii,jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+35, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+35,1122) '                        '
!!cDc!!, 	endif
!!cDc!!, !!cDc!!
!
               do k = 1,MaxBasis
                  z_new(k) = pvold(k) + alpha*pf(k)!!cDt!!
               enddo
!!cDc!!, !!cDc!!
!!cDc!!, if (MyPE .eq. 0) then
!!cDc!!,    write (MyPE+36,1122) '========================'
!!cDc!!,    write (MyPE+36,1342) 'b4 z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,    write (MyPE+36,1122) '========================'
!!cDc!!,    do i =1,vlen
!!cDc!!,       write (MyPE+36,1242) real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,    enddo
!!cDc!!, endif
!!cDc!!, !!cDc!!
               do i = 1,nn
                  gmi = CZERO
                  do j = 1,nn!!cDt!!
                     gmi = gmi + cMM_c(j)*bM_c(j,i)*pw(j) 
                  enddo
                  do k = 1,MaxBasis
                     z_new(k) = z_new(k) - gmi*pu(k,i)*pw(i)!!cDt!!
                  enddo
               enddo
!!cDc!!, !!cDc!!
!!cDc!!, if (MyPE .eq. 0) then
!!cDc!!,    write (MyPE+37,1122) '========================'
!!cDc!!,    write (MyPE+37,1342) 'aftr z_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,    write (MyPE+37,1122) '========================'
!!cDc!!,    do i =1,vlen
!!cDc!!,       write (MyPE+37,1242) real(z_new(i)),aimag(z_new(i)),'*i'
!!cDc!!,    enddo
!!cDc!!, endif
!!cDc!!, !!cDc!!
            !!xDx!! Deallocating the MAtrices aM_c and bM_c=aM_c^-1
            deallocate(aM_c)
            deallocate(bM_c)
            deallocate(cMM_c)

   endif   

!!xDx!! End of Broyden Mixing -------------------------------------------------------------------------------------------------------------------------------
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = RPMMix(id)%alpha
            vlen =RPMMix(id)%vlen
            A => RPMMix(id)%A_c(1:vlen,1:2)
!!xDx!! commented for testRPM_c_3 typ=14
!            pvold    => RPMMix(id)%vold_c(1:vlen)
!            pf => RPMMix(id)%f_c(1:vlen)
!            pdf => RPMMix(id)%df_c(1:vlen)
!!xDx!! commented for testRPM_c_3 typ=14
            pvect_old => RPMMix(id)%vector_old_c(1:vlen)
            pvect_new => RPMMix(id)%vector_new_c(1:vlen)
            pv_st_old    => RPMMix(id)%v_st_old_c(1:vlen)
         !!xDx!! projecting local back to global variables
            pvect_new(1:vlen)= matmul( Transpose(conjg( G(1:MaxBasis,(id-1)*vlen+1:id*vlen) ) ), z_new(1:MaxBasis) ) !!iDi!!
            pvect_old(1:vlen)= matmul( Transpose(conjg( G(1:MaxBasis,(id-1)*vlen+1:id*vlen) ) ), z_old(1:MaxBasis) ) !!iDi!!
!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (MyPE+40,1122) '========================'
!!cDc!!,       write (MyPE+40,1342) 'p_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+40,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+40,1242) real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, enddo
!!cDc!!, !!cDc!!
         !!xDx!! Combining the resultant to return to main -------------------------------------------------------------------------------------------------------
!            RPMMix(id)%vector_new_c(1:vlen)=q_new(id,1:vlen)+pvect_new!! u<- q+p
            do k=1,vlen
               pvect_new(k)=q_new(id,k)+pvect_new(k)!! u<- q+p
            enddo
!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (MyPE+42,1122) '========================'
!!cDc!!,       write (MyPE+42,1342) 'pvect_new(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,       write (MyPE+42,1122) '========================'
!!cDc!!,       do i =1,vlen
!!cDc!!,          write (MyPE+42,1242) real(pvect_new(i)),aimag(pvect_new(i)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, enddo
!!cDc!!, !!cDc!!
!
!!xDx!! Increasing size of basis --------------------------------------------------------------------------------------
!!xDx!! modified for test_RPM_c_3 typ=14
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
               R( ii, (id-1)*vlen+jj ) = conjg(A( jj, ii ))
               enddo
               enddo
!!cDc!!, !!cDc!!
!!cDc!!, 	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!, 	do k=0,NumPEs-1
!!cDc!!, 	if ( MyPE .eq. k) then
!!cDc!!, 	write (MyPE+50,1122) '                        '
!!cDc!!, 	write (MyPE+50,1122) '========================'
!!cDc!!, 	write (MyPE+50,1342) 'MyPE','R(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+50,1122) '========================'
!!cDc!!, 	do ii=1,2
!!cDc!!, 	   do jj=1,vlen
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(R(ii,(id-1)*vlen+jj)).gt.0)
!!cDc!!,               write(MyPE+50,fmt, advance="no") R(ii,(id-1)*vlen+jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+50, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+50,1122) '                        '
!!cDc!!, 	endif
!!cDc!!, 	enddo
!!cDc!!, !!cDc!!

            endif
!!xDx!! modified for test_RPM_c_3 typ=14

!!cDc!!, !!cDc!!
!!cDc!!,  do k=0,NumPEs-1
!!cDc!!,     if (MyPE .eq. k) then
!!cDc!!,        write (MyPE+52,1122) '                        '
!!cDc!!,        write (MyPE+52,1122) '========================'
!!cDc!!,        write (MyPE+52,1342) 'q_old(:)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,        write (MyPE+52,1122) '========================'
!!cDc!!,        do i =1,vlen
!!cDc!!,           write (MyPE+52,1242) real(q_old(id,i)),aimag(q_old(id,i)),'*i'
!!cDc!!,        enddo
!!cDc!!,        write (MyPE+52,1122) '========================'
!!cDc!!,        write (MyPE+52,1342) 'A(:,1)','MyPE,RPMIter, id', MyPE, RPMIter, id
!!cDc!!,        write (MyPE+52,1122) '========================'
!!cDc!!,        do i =1,vlen
!!cDc!!,           write (MyPE+52,1242) real(A(i,1)),aimag(A(i,1)),'*i'
!!cDc!!,        enddo
!!cDc!!, 
!!cDc!!,        write (MyPE+52,1122) '                        '
!!cDc!!,     endif
!!cDc!!,  enddo
!!cDc!!, !!cDc!!

         enddo
      enddo
!                                                   
!!xDx!! Performing the calc to increase the Basis size                                                                                                          
      if ( mod( RPMIter, MaxRPMIter ) .eq. 0 ) then !!nDn!! iter_broy= 0 & RPMIter> MaxRPMIter

      !!xDx!! HouseHolder LQ---------------------------------------------------------------- 
#ifdef USE_SCALAPACK
            !------- Allocating Tau -------------------------------------
                  allocate( TAU( 1:NUMROC(IR+min(M_R,N_R)-1, MB_R , MYROW, 0, NPROW) ) )

!!xDx!!Test
        print *, 'MyPE', MyPE
        print *, 'M_R', M_R
        print *, 'N_R', N_R
        print *, 'INFO_R', INFO_R
!!xDx!!Test

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
!!cDc!!	print *, 'RPMIter,Ri11,Ri22',RPMIter,Ri11,Ri22
!!nDn!! write the LAPACK version.
#endif

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

!!cDc!!, !!cDc!!
!!cDc!!, 	fmt = '(F9.6,"+",F9.6,"i")'
!!cDc!!, 	do k=0,NumPEs-1
!!cDc!!, 	if ( MyPE .eq. k) then
!!cDc!!, 	write (MyPE+60,1122) '                        '
!!cDc!!, 	write (MyPE+60,1122) '========================'
!!cDc!!, 	write (MyPE+60,1342) 'MyPE','G(:,:),RPMIter, id', MyPE, RPMIter, id
!!cDc!!, 	write (MyPE+60,1122) '========================'
!!cDc!!, 	do ii=1,Basis
!!cDc!!, 	   do jj=1,NumQuantities(1)*vlen
!!cDc!!,               fmt(8:8) = MERGE('+',' ',imag(G(ii,jj)).gt.0)
!!cDc!!,               write(MyPE+60,fmt, advance="no") G(ii,jj)
!!cDc!!, 	   enddo
!!cDc!!, 	   write(MyPE+60, fmt="(a)") " "
!!cDc!!, 	enddo
!!cDc!!, 	write (MyPE+60,1122) '                        '
!!cDc!!, 	endif
!!cDc!!, 	enddo
!!cDc!!, !!cDc!!
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
               !  ------------------------- Actual QR  ------------------------------
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

!!xDx!! edited for test_RPM_c_3 typ=14
                 if ( Basis .lt. MaxBasis ) then
                    G(Basis+1:MaxBasis,:) = CZERO
                 endif
            endif !!Basis .gt. 1
!!xDx!! edited for test_RPM_c_3 typ=14

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
!!xDx!! removed for test_RPM_c_3 typ=14
!   nullify(plq_tmp)
!!xDx!! removed for test_RPM_c_3 typ=14
!
end subroutine calRPM_c

!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!xDx!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine setBroydenMixing( idt, idq, amix )
!  ===================================================================
   implicit none
!
   integer(kind=IntKind), intent(in) :: idt, idq
!  integer(kind=IntKind), optional   :: nbroyiter
!  integer(kind=IntKind), optional   :: inv_method
!
   real(kind=RealKind), intent(in) :: amix
!
   integer(kind=IntKind) :: id
!
   if ( .not.Initialized ) then
      call ErrorHandler('setBroydenMix',                              &
          'The MixingModule needs to be initialized first')
   endif
!
   if ( idt<0 .or. idt>NumTypes ) then
      call ErrorHandler('setBroydenMixing',"Invalid type index",idq)
   endif
   if ( idq<0 .or. idq>NumQuantities(idt) ) then
      call ErrorHandler('setBroydenMixing',"Invalid quantity index",idq)
   endif
!
   if ( .not.allocated(BroydenMix) ) then
      allocate( BroydenMix(NumTotalMix) )
      do id = 1, NumTotalMix
         BroydenMix(id)%Initialized = .false.
      enddo
   endif  
!
   NumBroydenIter = MaxNumBroydenIter
!  if ( present(nbroyiter) ) then
!     NumBroydenIter = nbroyiter
!     if ( nbroyiter > MaxNumBroydenIter ) then
!        call WarningHandler('setBroydenMixing', &
!            'The number of Broyden iterations exeedes the maximum limit')
!        NumBroydenIter = MaxNumBroydenIter
!     endif
!  endif
!
   id = getMixID(idt,idq)
   if ( .not.BroydenMix(id)%Initialized ) then
      BroydenMix(id)%alpha = amix
   endif
!
!  if (present(inv_method)) BroydenInvMethod = inv_method
!
   if (BroydenInvMethod<0 .or. BroydenInvMethod>1) then
      call WarningHandler('initMixing', &
           'Invalid inversion method switching to default(invert1)')
      BroydenInvMethod = 0
   endif 
!
   if ( MixingMethod /= 1 ) then
      MixingMethod = 1
   endif
!
   end subroutine setBroydenMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine delBroydenMixing( )
!  ===================================================================
   implicit none
!
   integer(kind=IntKind) :: id
!
   if ( allocated(BroydenMix) ) then
      do id = 1,NumTotalMix
         if ( .not.BroydenMix(id)%Initialized ) cycle
         BroydenMix(id)%Initialized = .false.
         if ( associated(BroydenMix(id)%u_r) ) deallocate( BroydenMix(id)%u_r )
         if ( associated(BroydenMix(id)%vt_r) ) deallocate( BroydenMix(id)%vt_r )
         if ( associated(BroydenMix(id)%f_r) ) deallocate( BroydenMix(id)%f_r )
         if ( associated(BroydenMix(id)%df_r) ) deallocate( BroydenMix(id)%df_r )
         if ( associated(BroydenMix(id)%w_r) ) deallocate( BroydenMix(id)%w_r )
         if ( associated(BroydenMix(id)%vold_r) ) deallocate( BroydenMix(id)%vold_r )
         if ( associated(BroydenMix(id)%u_c) ) deallocate( BroydenMix(id)%u_c )
         if ( associated(BroydenMix(id)%vt_c) ) deallocate( BroydenMix(id)%vt_c )
         if ( associated(BroydenMix(id)%f_c) ) deallocate( BroydenMix(id)%f_c )
         if ( associated(BroydenMix(id)%df_c) ) deallocate( BroydenMix(id)%df_c )
         if ( associated(BroydenMix(id)%w_c) ) deallocate( BroydenMix(id)%w_c )
         if ( associated(BroydenMix(id)%vold_c) ) deallocate( BroydenMix(id)%vold_c )
      enddo
!      deallocate( BroydenMix )
   endif
!
   end subroutine delBroydenMixing
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBroydenMixing_r( list_q )
!  ===================================================================
!!xDx!! Expalnation of how the ocde works
!!xDx!! 1. n is the input potential at each iteration. 
!!xDx!! This is denoted as pvect_new when setting up the potential and pvect_old after the Greens function has been calculated.
!!xDx!! 2. This is used to calculate the charge density by calculating the greens function and integrating over energy
!!xDx!! 3. This charge is then used as the input for a poisson solver to get the new potential (n_{out}[n]). Denoted pvect_new
!!xDx!! 4. The mixing scheme is used to mix pvect_new and pvect_old, this is stored in pvect_new and sent as input for the next iteration.

!!xDx!! Notes for Digo:-
!!xDx!! 1. Wrt DDjohnson, n_out[n] is pvect_new. 
!!xDx!! 2. Wrt DDjohnson, n is pvect_old. 
!!xDx!! 3. pf:= pointer to F, where F[n]:=n_out[n]-n as defined in DD Johson.
!!xDx!! 4. For each iter, pvect_old is the pvect_new as returned by the previous iteration.
!!xDx!! 

!!xDx!! removed for test_RPM_c_3 typ=14
!!xDx!!
!   use printDataModule, only : printDataReal, printDataInteger,                &
!                               printDataPointer,printDataRealArray
!!xDx!! removed for test_RPM_c_3 typ=14

!
   implicit none
!
   type(MixListRealStruct), target :: list_q
!
   integer(kind=IntKind) :: i, j, k, info, iter, idq, idt, id
   integer(kind=IntKind) :: lastm1, nn, vlen
!
   real(kind=RealKind) :: fac1, fac2, fnorm, dfnorm, w0
   real(kind=RealKind) :: aij, gmi, cmj, wtmp
   real(kind=RealKind) :: msgbuf(2)
   real(kind=RealKind) :: rms
   real(kind=RealKind) :: alpha
!
   real(kind=RealKind), pointer :: pf(:), pdf(:), pu(:,:), pw(:), pvt(:,:)
   real(kind=RealKind), pointer :: pvect_old(:), pvect_new(:), pvold(:)
   real(kind=RealKind), pointer :: pad(:), pbd(:)
!
   type(MixListRealStruct), pointer :: plq_tmp
!
   iter = iter_count
   plq_tmp => list_q
!
   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_r(MixingMethod)
   endif
!
   do idt = 1,NumTypes
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         BroydenMix(id)%rms = plq_tmp%rms
         BroydenMix(id)%vector_old_r => plq_tmp%vector_old(1:vlen)
         BroydenMix(id)%vector_new_r => plq_tmp%vector_new(1:vlen)
         if ( .not.BroydenMix(id)%Initialized ) then
            call setBroydenSpace_r(id,vlen)
         endif
         plq_tmp => plq_tmp%next
      enddo
   enddo
!
   if ( iter_count == 1) then
!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!

      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen =BroydenMix(id)%vlen
            pf => BroydenMix(id)%f_r(1:vlen)
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
!
            pvold    => BroydenMix(id)%vold_r(1:vlen)
            do k = 1,vlen
               pf(k) = pvect_new(k) - pvect_old(k)
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k)
            enddo
!
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k)
            enddo
         enddo
      enddo
!
   else
!     ===============================================================
!     iter_count > 1: this is where the non-linear mixing is done
!     ===============================================================
      lastIter = lastIter + 1      ! update pointers
      lastm1   = lastIter - 1
!
!     ===============================================================
!     set current lenght of broyden cycle
!     ===============================================================
!
      if ( iter > NumBroydenIter ) then 
         nn = NumBroydenIter
      else
         nn = lastIter
      endif
!
      LoopType: do idt = 1,NumTypes
!        ============================================================
!        set weighting factor for the zeroth iteration
!        ============================================================
         w0= ten2m2
!
!        ============================================================
!        find: f[i] := vector(2)[i] - vector(1)[i]
!              df   := f[i] - f[i-1]
!        ============================================================
!!xDx!! calculated new F and dF.  
         dfnorm = zero
         fnorm = zero
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
            pf  => BroydenMix(id)%f_r(1:vlen) !!This is F^(iter-1):=F[n^(iter-1)]=n_{out}[n^(iter-1)]-n^(iter-1).
            pdf => BroydenMix(id)%df_r(1:vlen)
            if ( vlen==1 ) then
               pvect_new(1) = pvect_old(1) + alpha * pf(1) !!xDx!! ??? why ???
               cycle
            endif
            do k = 1,vlen
               pdf(k) = pvect_new(k) - pvect_old(k) - pf(k)!! dF^(iter)=F^(iter)-F^(iter-1)
               pf(k)  = pvect_new(k) - pvect_old(k) !!F^(iter)
            enddo

!           =========================================================
!           find: fnorm  := |f|
!                 dfnorm := |df|
!           =========================================================
            do k = 1,vlen  !!xDx!! These are sums over all the atoms in the proc. 
               dfnorm = dfnorm + pdf(k)*pdf(k) 
               fnorm  = fnorm  + pf(k)*pf(k)
            enddo
         enddo
!
         msgbuf(1) = dfnorm
         msgbuf(2) = fnorm
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID, msgbuf,2) !! Sum over the sums of ||dF||^2 and ||F||^2 in each proc & returns resultant to all procs
!!xDx!! If the potential is modelled as a vector then based on the discretized space, each index of the vector would correspond to a unique point
!!xDx!! In the whole space. Since we are calculating dF and F locally around each atom, we need to calculate the actual norm by summing over the 
!!xDx!! potential form the other atoms. Thus the Global sum in Group. The sum over all the atoms in each proc is reflected in the calc being performed
!!xDx!! in the idq loop.
!        ------------------------------------------------------------
         dfnorm = sqrt( msgbuf(1) )  
         fnorm  = sqrt( msgbuf(2) )  
!        ============================================================
!        set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
!             vold      := vector(1) 
!             vector(1) := df/|df|
!        ============================================================
!
         fac2 = one/dfnorm
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen  = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            fac1 = alpha*fac2
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
            pu    => BroydenMix(id)%u_r(1:vlen,1:NumBroydenIter)
            pvt   => BroydenMix(id)%vt_r(1:vlen,1:NumBroydenIter)
            pdf   => BroydenMix(id)%df_r(1:vlen)
            pvold => BroydenMix(id)%vold_r(1:vlen) !!n^(1), n^(iter-1), temp storage of pvect_old

            do k = 1,vlen
               pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k)) !! aI(dF/|dF|)+(n^(2)-n^(1))/|dF|, Eqn (13b)
               pvold(k)   = pvect_old(k)
               pvect_old(k) = fac2*pdf(k)  !! storeds the part of the unit vector of dF.
            enddo !! pvect_new now changes to |u^(n)> 
!           =========================================================
!           store vector(1) and vector(2) in the stacks u and vt
!           v(i) = v(i+1)
!           u(i) = u(i+1)
!           u(nn) = vector(1)
!           v(nn) = vector(2)
!           =========================================================
!           ---------------------------------------------------------
            call broy_sav_r( pu, pvt, pvect_old, pvect_new, iter-1, & 
                             NumBroydenIter, vlen ) 
!!xDx!! if iter-1<=NumBroydenIter, pu(:,iter-1)=pvect_new(:) calculated @ iter and
!!xDx!! pvt(:,iter-1)=pvect_old(:) calculated @ iter.
!!xDx!! else iter-1>NumBroydenIter and each column of pu and pvt contain pvect_new
!!xDx!! and pvect_old for all the previous iterations with dedicated colums. Shift
!!xDx!! each column to the left and store the latest pvect_. in the the final coln.
!           ---------------------------------------------------------
         enddo
!        ============================================================
!        calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:
!           a(i,j<i) = v(j)*v(l)
!           a(i,i) = v(i)*v(i)
!           a(j,i) = a(i,j)
!           cm(i) = sum_(l=1,i) v(l)*f
!        ============================================================
!!xDx!! nn -> #iters-1, equiv to m in DD Johnson.
         do j = 1,nn - 1          ! off diagonal elements of a(i,j)
            do i = j+1,nn
               aij = zero
               do idq = 1,NumQuantities(idt)
                  id = getMixID(idt,idq)
                  vlen = BroydenMix(id)%vlen
                  if ( vlen==1 ) then
                     cycle
                  endif
                  pvt => BroydenMix(id)%vt_r(1:vlen,1:NumBroydenIter)
                  do k = 1,vlen
                     aij = aij + pvt(k,j)*pvt(k,i)
                  enddo
               enddo
!              ------------------------------------------------------
!              call GlobalSumInGroup(GroupID,aij)
!              ------------------------------------------------------
               a_r(i,j) = aij
               a_r(j,i) = aij
            enddo
         enddo
!
         do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
            aij = zero
            cmj = zero
            do idq = 1, NumQuantities(idt)
               id = getMixID(idt,idq)
               vlen = BroydenMix(id)%vlen
               if ( vlen==1 ) then
                  cycle
               endif
               pvt => BroydenMix(id)%vt_r(1:vlen,1:NumBroydenIter)
               pf  => BroydenMix(id)%f_r(1:vlen)
               do k = 1,vlen
                  cmj = cmj + pvt(k,i)*pf(k) !!xDx!! c^m_k/w_k in Eqn (15b). Eavaluated for fixed m, ie, for the current iteration.
                  aij = aij + pvt(k,i)*pvt(k,i)
               enddo
            enddo
!           msgbuf(1) = aij
!           msgbuf(2) = cmj
!           ---------------------------------------------------------
!           call GlobalSumInGroup(GroupID,msgbuf,2)
!           ---------------------------------------------------------
!           aij = msgbuf(1)
!           cmj = msgbuf(2)
            a_r(i,i) = aij
            cm_r(i) = cmj
         enddo
!!xDx!! the a_r(i,j) is a_{ij}/(w_iw_j) of Eqn. (13a). |||y for cm_r.
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,a_r,NumBroydenIter,NumBroydenIter)
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,cm_r,NumBroydenIter)
!        ------------------------------------------------------------
!        ============================================================
!        shift down weights in stack
!
!        (TCS, bug fixed 8/5/97: replace iter by iter-1 -> see broy_sav)
!        ============================================================
         do idq = 1,NumQuantities(idt)
!
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_r(1:vlen)
!
            pvold => BroydenMix(id)%vold_r(1:vlen)
            pw => BroydenMix(id)%w_r(1:NumBroydenIter) !!?D?!!
            pu => BroydenMix(id)%u_r(1:vlen,1:NumBroydenIter)
            pf => BroydenMix(id)%f_r(1:vlen)
!
            if ( iter-1 > NumBroydenIter ) then
               do i = 1,NumBroydenIter-1
                  pw(i) = pw(i+1)
               enddo
            endif
!
            rms = BroydenMix(id)%rms
            wtmp = zero
            if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
            if ( wtmp < ONE ) wtmp = ONE
            if ( iter > NumBroydenIter ) then
               pw(NumBroydenIter) = wtmp
            else
               pw(lastIter) = wtmp       !w(lastm1)=wtmp
            endif
!           =========================================================
!           now calculate the b-matrix:
!           b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1 !!xDx!! Eqn. (13a). Assuming this is just the standard calc of \beta.
!           =========================================================
            if ( BroydenInvMethod==0 .and. vlen>3*NumBroydenIter ) then
!              ======================================================
!              use invert1
!              ======================================================
               do i = 1,nn
                  do j = 1,nn
                     d_r(j,i) = a_r(j,i)*pw(j)*pw(i)
                     b_r(j,i) = zero
                  enddo
                  b_r(i,i) = one
                  d_r(i,i) = w0**2 + a_r(i,i)*pw(i)*pw(i)
               enddo
!
!              this is very unlikely
!
               if ( 3*NumBroydenIter > vlen ) then !!?D?!! this will be false based on #1604.
                  call ErrorHandler( "callBroyden", & 
                      "need larger dimension for MaxSizeBroydenVect")
               endif
               pvect_old => BroydenMix(id)%vector_old_r(1:vlen)
               pad =>                                                 &
               BroydenMix(id)%vector_old_r(NumBroydenIter+1:2*NumBroydenIter)
               pbd =>                                                 &
               BroydenMix(id)%vector_old_r(2*NumBroydenIter+1:3*NumBroydenIter)
!              ------------------------------------------------------
               call invert1_r(d_r,b_r,nn,pvect_old,pad,pbd,NumBroydenIter)
!              ------------------------------------------------------
            else
!              ======================================================
!              use lapack
!              ======================================================
               do i = 1,nn
                  do j = 1,nn
                     b_r(j,i) = a_r(j,i)*pw(j)*pw(i)
                  enddo
                  b_r(i,i) = w0**2 + a_r(i,i)*pw(i)*pw(i)
               enddo
               if ( .not.allocated(ipiv) ) then
                  allocate( ipiv(NumBroydenIter) )
               endif
!              ------------------------------------------------------
               call dgetrf( nn, nn, b_r, NumBroydenIter,ipiv, info )
!              ------------------------------------------------------
               call dgetri( nn, b_r, NumBroydenIter, ipiv, d_r, nn, info )
!              ------------------------------------------------------
!              write(6,*) ' optimum lwork', d(1,1)
            endif
!
            do k = 1,vlen
               pvect_new(k) = pvold(k) + alpha*pf(k)
            enddo
            do i = 1,nn
               gmi = zero
               do j = 1,nn
                  gmi = gmi + cm_r(j)*b_r(j,i)*pw(j) !!xDx!! \gamma_{ml} in Eqn. (15b), the cm_r is wihtout w_j
               enddo
               do k = 1,vlen
                  pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)
               enddo
            enddo
         enddo
      enddo LoopType
   endif
!
   end subroutine calBroydenMixing_r
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calBroydenMixing_c( list_q )
!  ===================================================================

!!xDx!! removed for test_RPM_c_3 typ=14
!!xDx!!
!   use printDataModule, only : printDataReal, printDataInteger,                &
!                               printDataPointer,printDataRealArray,            &
!                               printDataComplex,printDataComplexArray           
!!xDx!!
!!xDx!! removed for test_RPM_c_3 typ=14
!
   implicit none
!
   type(MixListCmplxStruct), target :: list_q
!
   integer(kind=IntKind) :: i, j, k, info, iter, idq, idt, id
   integer(kind=IntKind) :: lastm1, nn, vlen
!
   real(kind=RealKind) :: fac1, fac2, fnorm, dfnorm, w0
   real(kind=RealKind) :: wtmp
   real(kind=RealKind) :: msgbuf(2)
   real(kind=RealKind) :: rms
   real(kind=RealKind) :: alpha
!
   complex(kind=CmplxKind) :: aij, gmi, cmj
   complex(kind=CmplxKind), pointer :: pf(:), pdf(:), pu(:,:), pw(:), pvt(:,:)
   complex(kind=CmplxKind), pointer :: pvect_old(:), pvect_new(:), pvold(:)
   complex(kind=CmplxKind), pointer :: pad(:), pbd(:)
!
   type(MixListCmplxStruct), pointer :: plq_tmp
!
   iter = iter_count
   plq_tmp => list_q
!
   if ( .not.InitializedWkSpace ) then
      call setWorkingSpace_c(MixingMethod)
   endif
!
   do idt = 1,NumTypes
      do idq = 1,NumQuantities(idt)
         id = getMixID(idt,idq)
         vlen = plq_tmp%size
         BroydenMix(id)%rms = plq_tmp%rms
         BroydenMix(id)%vector_old_c => plq_tmp%vector_old(1:vlen)
         BroydenMix(id)%vector_new_c => plq_tmp%vector_new(1:vlen)
         if ( .not.BroydenMix(id)%Initialized ) then
            call setBroydenSpace_c(id,vlen)
         endif
         plq_tmp => plq_tmp%next
      enddo
   enddo
!
   if ( iter_count == 1) then
!     ===============================================================
!     first iteration: preform linear mixing, load f and vold, set
!                      different pointers and variables
!     ===============================================================
      lastIter = 0             ! initialize pointers
      lastm1 = lastIter - 1
!
      do idt = 1,NumTypes
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen = BroydenMix(id)%vlen
            pf => BroydenMix(id)%f_c(1:vlen)
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
!
            pvold    => BroydenMix(id)%vold_c(1:vlen)
            do k = 1,vlen
               pf(k) = pvect_new(k) - pvect_old(k)
            enddo
            do k = 1,vlen
               pvold(k) = pvect_old(k)
            enddo
!
            do k = 1,vlen          ! this is the linear mixing
               pvect_new(k) = pvect_old(k) + alpha * pf(k)
            enddo

         enddo
      enddo
!
   else
!     ===============================================================
!     iter_count > 1: this is where the non-linear mixing is done
!     ===============================================================
      lastIter = lastIter + 1      ! update pointers
      lastm1   = lastIter - 1
!
!     ===============================================================
!     set current lenght of broyden cycle
!     ===============================================================
      if ( iter > NumBroydenIter ) then
         nn = NumBroydenIter
      else
         nn = lastIter
      endif
!
      LoopType: do idt = 1,NumTypes
!        ============================================================
!        set weighting factor for the zeroth iteration
!        ============================================================
         w0= ten2m2
!
!        ============================================================
!        find: f[i] := vector(2)[i] - vector(1)[i]
!              df   := f[i] - f[i-1]
!        ============================================================
         dfnorm = zero
         fnorm = zero
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
            pf  => BroydenMix(id)%f_c(1:vlen)
            pdf => BroydenMix(id)%df_c(1:vlen)
            if ( vlen==1 ) then
               pvect_new(1) = pvect_old(1) + alpha * pf(1)
               cycle
            endif
            do k = 1,vlen
               pdf(k) = pvect_new(k) - pvect_old(k) - pf(k)
               pf(k)  = pvect_new(k) - pvect_old(k)
            enddo
!           =========================================================
!           find: fnorm  := |f|
!                 dfnorm := |df|
!           =========================================================
            do k = 1,vlen
               dfnorm = dfnorm + real(pdf(k)*conjg(pdf(k)),kind=RealKind)
               fnorm  = fnorm  + real(pf(k)*conjg(pf(k)),kind=RealKind)
            enddo
         enddo
!
         msgbuf(1) = dfnorm
         msgbuf(2) = fnorm
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID, msgbuf,2)
!        ------------------------------------------------------------
         dfnorm = sqrt( msgbuf(1) )
         fnorm  = sqrt( msgbuf(2) )
!        ============================================================
!        set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|
!             vold      := vector(1) 
!             vector(1) := df/|df|
!        ============================================================
!
         fac2 = one/dfnorm
         do idq = 1,NumQuantities(idt)
            id = getMixID(idt,idq)
            alpha = BroydenMix(id)%alpha
            vlen  = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            fac1 = alpha*fac2
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
            pu    => BroydenMix(id)%u_c(1:vlen,1:NumBroydenIter)
            pvt   => BroydenMix(id)%vt_c(1:vlen,1:NumBroydenIter)
            pdf   => BroydenMix(id)%df_c(1:vlen)
            pvold => BroydenMix(id)%vold_c(1:vlen)
!
            do k = 1,vlen
               pvect_new(k) = fac1*pdf(k) + fac2*(pvect_old(k) - pvold(k))
               pvold(k)   = pvect_old(k)
               pvect_old(k) = fac2*pdf(k)
            enddo
!           =========================================================
!           store vector(1) and vector(2) in the stacks u and vt
!           v(i) = v(i+1)
!           u(i) = u(i+1)
!           u(nn) = vector(1)
!           v(nn) = vector(2)
!           =========================================================
!           ---------------------------------------------------------
            call broy_sav_c( pu, pvt, pvect_old, pvect_new, iter-1, & 
                             NumBroydenIter, vlen )
!           ---------------------------------------------------------
         enddo
!        ============================================================
!        calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:
!           a(i,j<i) = v(j)*v(l)
!           a(i,i) = v(i)*v(i)
!           a(j,i) = a(i,j)
!           cm(i) = sum_(l=1,i) v(l)*f
!        ============================================================
         do j = 1,nn - 1          ! off diagonal elements of a(i,j)
            do i = j+1,nn
               aij = czero
               do idq = 1,NumQuantities(idt)
                  id = getMixID(idt,idq)
                  vlen = BroydenMix(id)%vlen
                  if ( vlen==1 ) then
                     cycle
                  endif
                  pvt => BroydenMix(id)%vt_c(1:vlen,1:NumBroydenIter)
                  do k = 1,vlen
                     aij = aij + conjg(pvt(k,j))*pvt(k,i)
                  enddo
               enddo
!              ------------------------------------------------------
!              call GlobalSumInGroup(GroupID,aij)
!              ------------------------------------------------------
               a_c(i,j) = aij
               a_c(j,i) = conjg(aij)
            enddo
         enddo
!
         do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
            aij = czero
            cmj = czero
            do idq = 1, NumQuantities(idt)
               id = getMixID(idt,idq)
               vlen = BroydenMix(id)%vlen
               if ( vlen==1 ) then
                  cycle
               endif
               pvt => BroydenMix(id)%vt_c(1:vlen,1:NumBroydenIter)
               pf  => BroydenMix(id)%f_c(1:vlen)
               do k = 1,vlen
                  cmj = cmj + conjg(pvt(k,i))*pf(k)
                  aij = aij + conjg(pvt(k,i))*pvt(k,i)
               enddo
            enddo
!           msgbuf(1) = aij
!           msgbuf(2) = cmj
!           ---------------------------------------------------------
!           call GlobalSumInGroup(GroupID,msgbuf,2)
!           ---------------------------------------------------------
!           aij = msgbuf(1)
!           cmj = msgbuf(2)
            a_c(i,i) = aij
            cm_c(i) = cmj
         enddo
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,a_c,NumBroydenIter,NumBroydenIter)
!        ------------------------------------------------------------
         call GlobalSumInGroup(GroupID,cm_c,NumBroydenIter)
!        ------------------------------------------------------------
!        ============================================================
!        shift down weights in stack
!
!        (TCS, bug fixed 8/5/97: replace iter by iter-1 -> see broy_sav)
!        ============================================================
         do idq = 1,NumQuantities(idt)
!
            id = getMixID(idt,idq)
            vlen = BroydenMix(id)%vlen
            if ( vlen==1 ) then
               cycle
            endif
            alpha = BroydenMix(id)%alpha
            pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
            pvect_new => BroydenMix(id)%vector_new_c(1:vlen)
!
            pvold => BroydenMix(id)%vold_c(1:vlen)
            pw => BroydenMix(id)%w_c(1:NumBroydenIter)
            pu => BroydenMix(id)%u_c(1:vlen,1:NumBroydenIter)
            pf => BroydenMix(id)%f_c(1:vlen)
!
            if ( iter-1 > NumBroydenIter ) then
               do i = 1,NumBroydenIter-1
                  pw(i) = pw(i+1)
               enddo
            endif
!

            rms = BroydenMix(id)%rms
            wtmp = zero
            if ( rms > ten2m9 ) wtmp = TWO*sqrt(ten2m2/rms)
            if ( wtmp < ONE ) wtmp = ONE
            if ( iter > NumBroydenIter ) then
               pw(NumBroydenIter) = wtmp
            else
               pw(lastIter) = wtmp       !w(lastm1)=wtmp
            endif
!           =========================================================
!           now calculate the b-matrix:
!           b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
!           =========================================================
            if ( BroydenInvMethod==0 .and. vlen>3*NumBroydenIter ) then
!              ======================================================
!              use invert1
!              ======================================================
               do i = 1,nn
                  do j = 1,nn
                     d_c(j,i) = a_c(j,i)*pw(j)*pw(i)
                     b_c(j,i) = czero
                  enddo
                  b_c(i,i) = cone
                  d_c(i,i) = w0**2 + a_c(i,i)*pw(i)*pw(i)
               enddo
!
!              this is very unlikely
!
               if ( 3*NumBroydenIter > vlen ) then 
                  call ErrorHandler( "callBroyden", & 
                      "need larger dimension for MaxSizeBroydenVect")
               endif
               pvect_old => BroydenMix(id)%vector_old_c(1:vlen)
               pad =>                                                 &
               BroydenMix(id)%vector_old_c(NumBroydenIter+1:2*NumBroydenIter)
               pbd =>                                                 &
               BroydenMix(id)%vector_old_c(2*NumBroydenIter+1:3*NumBroydenIter)
!              ------------------------------------------------------
               call invert1_c(d_c,b_c,nn,pvect_old,pad,pbd,NumBroydenIter)
!              ------------------------------------------------------
            else
!              ======================================================
!              use lapack
!              ======================================================
               b_c = CZERO
               do i = 1,nn
                  do j = 1,nn
                     b_c(j,i) = a_c(j,i)*pw(j)*pw(i)
                  enddo
                  b_c(i,i) = w0**2 + a_c(i,i)*pw(i)*pw(i)
               enddo
               if ( .not.allocated(ipiv) ) then
                  allocate( ipiv(NumBroydenIter) )
                  ipiv = 0
               endif
!              ------------------------------------------------------
               call zgetrf( nn, nn, b_c, NumBroydenIter,ipiv, info )
!              ------------------------------------------------------
               call zgetri( nn, b_c, NumBroydenIter, ipiv, d_c, nn, info )
!              ------------------------------------------------------
!              write(6,*) ' optimum lwork', d(1,1)
            endif
!
            do k = 1,vlen
               pvect_new(k) = pvold(k) + alpha*pf(k)
            enddo
            do i = 1,nn
               gmi = czero
               do j = nn, 1, -1
                  gmi = gmi + cm_c(j)*b_c(j,i)*pw(j)
               enddo
               do k = 1,vlen
                  pvect_new(k) = pvect_new(k) - gmi*pu(k,i)*pw(i)
               enddo
            enddo
         enddo
      enddo LoopType
   endif
!
   end subroutine calBroydenMixing_c

!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine broy_sav_r( fins, fots, vector_old, vector_new, itscf,  & 
                          istore, ivsiz )
!  ==================================================================
!!xDx!! broy_sav_r( pu, pvt, pvect_old, pvect_new, iter-1, NumBroydenIter, vlen)
!!xDx!! where the subroutine assigns value to fins and fots.
!!xDx!! Assigns the matrix u=[u^(1),u^(2),u^(3),...,u^(m)] where m=NumRPMIter_broy
!!xDx!! and u^(i) is as in (13b.iii) in DD Johnson. 
   implicit none
!
   integer(kind=IntKind), intent(in) :: itscf
   integer(kind=IntKind), intent(in) :: ivsiz
   integer(kind=IntKind), intent(in) :: istore
!
   real(kind=RealKind), target :: vector_old(:)
   real(kind=RealKind), target :: vector_new(:)
   real(kind=RealKind), target :: fins(:,:)
   real(kind=RealKind), target :: fots(:,:)
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
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
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
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine invert1_r( pa, pb, m, ptd, pad, pbd, mm ) 
!  ==================================================================
!  ==================================================================
!  temporary for testting, should be lapack routine in future!
!  ==================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: m, mm
!
   real(kind=RealKind), target :: pa(:,:)
   real(kind=RealKind), target :: pb(:,:)
   real(kind=RealKind), target :: ptd(:)
   real(kind=RealKind), target :: pad(:)
   real(kind=RealKind), target :: pbd(:)
!
   integer(kind=IntKind) :: n, i, j, k
!
   real(kind=RealKind) :: atmp
!
!  ==================================================================
!  parameter (mm=5)
!  ==================================================================
!
!  ==================================================================
!  subroutine to preform gaussian elimination
!             no zeros along the diagonal
!  ==================================================================
!
   n = m
   if ( n > mm ) then
      call ErrorHandler('invert1',' invert: matrix a too large')
   endif
!
   do i = 1,n
      atmp = pa(i,i)
      if ( abs(atmp) < ten2m8 ) then
         call ErrorHandler("invert1","matrix has zero diagonal(row)",i)
      endif
   enddo
!
   if ( n /= 1 ) then
!
      do i = 1,n
!
         do j = 1,n
            ptd(j) = pa(j,i)/pa(i,i)
         enddo
!
         ptd(i) = zero
!
         do k = 1,n
            pbd(k) = pb(i,k)
            pad(k) = pa(i,k)
         enddo
!
         do k = 1,n
            do j = 1,n
               pb(j,k) = pb(j,k)-(ptd(j)*pbd(k))
               pa(j,k) = pa(j,k)-(ptd(j)*pad(k))
            enddo
         enddo
!
      enddo
!
      do i = 1,n
         do j = 1,n
            pb(j,i) = pb(j,i)/pa(j,j)
         enddo
      enddo
!
   else
      pb(1,1) = one/pa(1,1)
   endif
!
   end subroutine invert1_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine invert1_c( pa, pb, m, ptd, pad, pbd, mm ) 
!  ==================================================================
!  ==================================================================
!  temporary for testting, should be lapack routine in future!
!  ==================================================================
!
   implicit none
!
   integer(kind=IntKind), intent(in) :: m, mm
!
   complex(kind=CmplxKind), target :: pa(:,:)
   complex(kind=CmplxKind), target :: pb(:,:)
   complex(kind=CmplxKind), target :: ptd(:)
   complex(kind=CmplxKind), target :: pad(:)
   complex(kind=CmplxKind), target :: pbd(:)
!
   integer(kind=IntKind) :: n, i, j, k
!
   complex(kind=CmplxKind) :: atmp
!
!  ==================================================================
!  parameter (mm=5)
!  ==================================================================
!
!  ==================================================================
!  subroutine to preform gaussian elimination
!             no zeros along the diagonal
!  ==================================================================
!
   n = m
   if ( n > mm ) then
      call ErrorHandler('invert1',' invert: matrix a too large')
   endif
!
   do i = 1,n
      atmp = pa(i,i)
      if ( abs(atmp) < ten2m8 ) then
         call ErrorHandler("invert1","matrix has zero diagonal(row)",i)
      endif
   enddo
!
   if ( n /= 1 ) then
!
      do i = 1,n
!
         do j = 1,n
            ptd(j) = pa(j,i)/pa(i,i)
         enddo
!
         ptd(i) = czero
!
         do k = 1,n
            pbd(k) = conjg(pb(i,k))
            pad(k) = conjg(pa(i,k))
         enddo
!
         do k = 1,n
            do j = 1,n
               pb(j,k) = pb(j,k) - ptd(j)*pbd(k)
               pa(j,k) = pa(j,k) - ptd(j)*pad(k)
            enddo
         enddo
!
      enddo
!
      do i = 1,n
         do j = 1,n
            pb(j,i) = pb(j,i)/pa(j,j)
         enddo
      enddo
!
   else
      pb(1,1) = one/pa(1,1)
   endif
!
   end subroutine invert1_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine dgaleq_r( a, y, n, ipits)
!  ==================================================================
!
   implicit   none
!
   integer(kind=IntKind), intent(in) :: n, ipits
!
   real(kind=RealKind), intent(inout) :: a(ipits+1,ipits+1)
   real(kind=RealKind), intent(inout) :: y(ipits+1)
!
   integer(kind=IntKind) :: i, j, k, ij
!
   real(kind=RealKind) :: f1,f2
!
!  ==================================================================
!  ******************************************************************
!  Important note: the following codes may will result a problem:
!         a(i-1,i-1) becomes zero!!!
!  It needs to be fixed.
!  ******************************************************************
!  ==================================================================
   do i = 2,n
      f1 = -one/a(i-1,i-1)
      do j = i,n
         f2 = f1*a(j,i-1)
         do k = 1,n
            a(j,k) = a(j,k)+f2*a(i-1,k)
         enddo
         y(j) = y(j)+f2*y(i-1)
      enddo
   enddo
   y(n) = y(n)/a(n,n)
   do ij = 1,n-1
      i = n-ij
      do j = 1,ij
         y(i) = y(i)-y(i+j)*a(i,i+j)
      enddo
      y(i) = y(i)/a(i,i)
   enddo
!
   end subroutine dgaleq_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine dgaleq_c( a, y, n, ipits)
!  ==================================================================
!
   implicit   none
!
   integer(kind=IntKind), intent(in) :: n, ipits
!
   complex(kind=CmplxKind), intent(inout) :: a(ipits+1,ipits+1)
   complex(kind=CmplxKind), intent(inout) :: y(ipits+1)
!
   integer(kind=IntKind) :: i, j, k, ij
!
   complex(kind=CmplxKind) :: f1,f2
!
!  ==================================================================
!  ******************************************************************
!  Important note: the following codes may will result a problem:
!         a(i-1,i-1) becomes zero!!!
!  It needs to be fixed.
!  ******************************************************************
!  ==================================================================
   do i = 2,n
      f1 = -cone/a(i-1,i-1)
      do j = i,n
         f2 = f1*a(j,i-1)
         do k = 1,n
            a(j,k) = a(j,k)+f2*conjg(a(i-1,k))
         enddo
         y(j) = y(j)+f2*y(i-1)
      enddo
   enddo
   y(n) = y(n)/a(n,n)
   do ij = 1,n-1
      i = n-ij
      do j = 1,ij
         y(i) = y(i)-y(i+j)*conjg(a(i,i+j))
      enddo
      y(i) = y(i)/a(i,i)
   enddo
!
   end subroutine dgaleq_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   function  dgarms_r(fins,fots,i,j,npts)         result(dga_rms)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: i, j, npts
!
   real(kind=RealKind), target :: fins(:,:), fots(:,:)
!
   integer(kind=IntKind) :: k
!
   real(kind=RealKind) :: dga_rms
!
   dga_rms = zero
   do k = 1,npts
      dga_rms = dga_rms + (fins(i,k)-fots(i,k))*(fins(j,k)-fots(j,k))
   enddo
!
   end function dgarms_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   function  dgarms_c(fins,fots,i,j,npts)         result(dga_rms)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: i, j, npts
!
   complex(kind=CmplxKind), target :: fins(:,:), fots(:,:)
!
   integer(kind=IntKind) :: k
!
   complex(kind=CmplxKind) :: dga_rms
!
   dga_rms = czero
   do k = 1,npts
      dga_rms = dga_rms + (fins(i,k)-fots(i,k))*conjg(fins(j,k)-fots(j,k))
   enddo
!
   end function dgarms_c
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  dgasad_r(ipits,cin,cot,b,niter,amix)       result(dga_sad)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: ipits, niter
!
   real(kind=RealKind), intent(in) :: amix
   real(kind=RealKind), target :: cin(:), cot(:), b(:)
!
   integer(kind=IntKind) :: i, iter
!
   real(kind=RealKind) :: optin, optot
!
   real(kind=RealKind) :: dga_sad
!
   iter = min(niter,ipits)
   optin = zero
   optot = zero
   do i = 1,iter
      optin = optin + b(i)*cin(i)
      optot = optot + b(i)*cot(i)
   enddo
   dga_sad = optin + amix*( optot - optin )
!
   end function dgasad_r
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  dgasad_c(ipits,cin,cot,b,niter,amix)       result(dga_sad)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: ipits, niter
!
   real(kind=RealKind), intent(in) :: amix
   complex(kind=CmplxKind), target :: cin(:), cot(:), b(:)
!
   integer(kind=IntKind) :: i, iter
!
   complex(kind=CmplxKind) :: optin, optot
!
   complex(kind=CmplxKind) :: dga_sad
!
   iter = min(niter,ipits)
   optin = czero
   optot = czero
   do i = 1,iter
      optin = optin + b(i)*cin(i)
      optot = optot + b(i)*cot(i)
   enddo
   dga_sad = optin + amix*( optot - optin )
!
   end function dgasad_c
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function  getMixID( idt, idq )                          result(id)
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: idt, idq
!
   integer(kind=IntKind) :: i
!
   integer(kind=IntKind) :: id
!
   id = 0
   do i = 1,idt-1
      id = id + NumQuantities(idt)
   enddo
   id = id + idq
!
   end function getMixID
!  ==================================================================
!
!  ******************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine  setMixingAlpha( idt, idq , amix )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: idq, idt
!
   real(kind=RealKind), intent(in) :: amix 
!
   integer(kind=IntKind) :: id
!
   id = getMixID(idt,idq)
   if ( MixingMethod == 0) then
      SimpleMix(id)%alpha = amix
   else if ( MixingMethod == 1 ) then
      BroydenMix(id)%alpha = amix
   else if ( MixingMethod == 2 ) then
      DGAMix(id)%alpha = amix
   endif
!
   end subroutine setMixingAlpha
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine setSimpleSpace( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   SimpleMix(id)%vlen = size
!
   SimpleMix(id)%Initialized = .true.
!
   end subroutine setSimpleSpace
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setDGASpace_r( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( DGAMix(id)%f_old_r(NumDGAIter,size) )
   allocate( DGAMix(id)%f_new_r(NumDGAIter,size) )
   DGAMix(id)%vlen = size
!
   nullify( DGAMix(id)%f_old_c )
   nullify( DGAMix(id)%f_new_c )
!
   DGAMix(id)%Initialized = .true.
!
   end subroutine setDGASpace_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setDGASpace_c( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( DGAMix(id)%f_old_c(NumDGAIter,size) )
   allocate( DGAMix(id)%f_new_c(NumDGAIter,size) )
   DGAMix(id)%vlen = size
!
   nullify( DGAMix(id)%f_old_r )
   nullify( DGAMix(id)%f_new_r )
!
   DGAMix(id)%Initialized = .true.
!
   end subroutine setDGASpace_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!
   subroutine  setRPMSpace_r( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( RPMMix(id)%u_r(size,NumRPMIter_broy) )
   allocate( RPMMix(id)%vt_r(size,NumRPMIter_broy) )
   allocate( RPMMix(id)%f_r(size) )
   allocate( RPMMix(id)%A_r(size,2) ) !!xDx!! for 2D Krylov subspace.
   allocate( RPMMix(id)%df_r(size) )
   allocate( RPMMix(id)%w_r(NumRPMIter_broy) )
   allocate( RPMMix(id)%vold_r(size) )
   allocate( RPMMix(id)%v_st_old_r(size) )
   RPMMix(id)%u_r = ZERO
   RPMMix(id)%vt_r = ZERO
   RPMMix(id)%f_r = ZERO
   RPMMix(id)%A_r = ZERO
   RPMMix(id)%df_r = ZERO
   RPMMix(id)%w_r = ZERO
   RPMMix(id)%vold_r = ZERO
   RPMMix(id)%v_st_old_r = ZERO
   RPMMix(id)%vlen = size
!
   nullify( RPMMix(id)%u_c )
   nullify( RPMMix(id)%vt_c )
   nullify( RPMMix(id)%f_c )
   nullify( RPMMix(id)%A_c )
   nullify( RPMMix(id)%df_c )
   nullify( RPMMix(id)%w_c )
   nullify( RPMMix(id)%vold_c )
   nullify( RPMMix(id)%v_st_old_c )
!
   RPMMix(id)%Initialized = .true.
!
   end subroutine setRPMSpace_r

!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
!!xDx!! modified because of error testRPM_c_3 typ=14
   subroutine  setRPMSpace_c( id, size1, size2 )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size1, size2
!
   allocate( RPMMix(id)%u_c(size2,NumRPMIter_broy) )
   allocate( RPMMix(id)%vt_c(size2,NumRPMIter_broy) )
   allocate( RPMMix(id)%f_c(size2) )
!!xDx!! modified because of error testRPM_c_3 typ=14
   allocate( RPMMix(id)%A_c(size1,2) ) !!A was 2xvlen at first
!!xDx!! modified because of error testRPM_c_3 typ=14
   allocate( RPMMix(id)%df_c(size2) )
   allocate( RPMMix(id)%w_c(NumRPMIter_broy) )
   allocate( RPMMix(id)%vold_c(size2) )
   allocate( RPMMix(id)%v_st_old_c(size1) )
   RPMMix(id)%u_c = CZERO
   RPMMix(id)%vt_c = CZERO
   RPMMix(id)%f_c = CZERO
   RPMMix(id)%A_c = CZERO
   RPMMix(id)%df_c = CZERO
   RPMMix(id)%w_c = CZERO
   RPMMix(id)%vold_c = CZERO
   RPMMix(id)%v_st_old_c = CZERO
   RPMMix(id)%vlen = size1
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
!!xDx!! modified for test_RPM_c_3 typ=14
   RPMMix(id)%Initialized = .true.
!!xDx!! modified for test_RPM_c_3 typ=14
!
   end subroutine setRPMSpace_c
!!xDx!! modified because of error testRPM_c_3 typ=14
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================

!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!xDx!!
   subroutine  setBroydenSpace_r( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( BroydenMix(id)%u_r(size,NumBroydenIter) )
   allocate( BroydenMix(id)%vt_r(size,NumBroydenIter) )
   allocate( BroydenMix(id)%f_r(size) )
   allocate( BroydenMix(id)%df_r(size) )
   allocate( BroydenMix(id)%w_r(NumBroydenIter) )
   allocate( BroydenMix(id)%vold_r(size) )
   BroydenMix(id)%u_r = ZERO
   BroydenMix(id)%vt_r = ZERO
   BroydenMix(id)%f_r = ZERO
   BroydenMix(id)%df_r = ZERO
   BroydenMix(id)%w_r = ZERO
   BroydenMix(id)%vold_r = ZERO
   BroydenMix(id)%vlen = size
!
   nullify( BroydenMix(id)%u_c )
   nullify( BroydenMix(id)%vt_c )
   nullify( BroydenMix(id)%f_c )
   nullify( BroydenMix(id)%df_c )
   nullify( BroydenMix(id)%w_c )
   nullify( BroydenMix(id)%vold_c )
!
   BroydenMix(id)%Initialized = .true.
!
   end subroutine setBroydenSpace_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setBroydenSpace_c( id, size )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: id, size
!
   allocate( BroydenMix(id)%u_c(size,NumBroydenIter) )
   allocate( BroydenMix(id)%vt_c(size,NumBroydenIter) )
   allocate( BroydenMix(id)%f_c(size) )
   allocate( BroydenMix(id)%df_c(size) )
   allocate( BroydenMix(id)%w_c(NumBroydenIter) )
   allocate( BroydenMix(id)%vold_c(size) )
   BroydenMix(id)%u_c = CZERO
   BroydenMix(id)%vt_c = CZERO
   BroydenMix(id)%f_c = CZERO
   BroydenMix(id)%df_c = CZERO
   BroydenMix(id)%w_c = CZERO
   BroydenMix(id)%vold_c = CZERO
   BroydenMix(id)%vlen = size
!
   nullify( BroydenMix(id)%u_r )
   nullify( BroydenMix(id)%vt_r )
   nullify( BroydenMix(id)%f_r )
   nullify( BroydenMix(id)%df_r )
   nullify( BroydenMix(id)%w_r )
   nullify( BroydenMix(id)%vold_r )
!
   BroydenMix(id)%Initialized = .true.
!
   end subroutine setBroydenSpace_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setWorkingSpace_r( MixId )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: MixId 
!
   if (MixId == 1 .or. MixId ==3) then
      allocate( a_r( NumBroydenIter,NumBroydenIter ) )
      allocate( b_r( NumBroydenIter,NumBroydenIter ) )
      allocate( d_r( NumBroydenIter,NumBroydenIter ) )
      allocate( cm_r( NumBroydenIter ) )
      a_r = ZERO; b_r = ZERO; d_r = ZERO; cm_r = ZERO
!
      if ( BroydenInvMethod == 1 ) then
         allocate( ipiv(NumBroydenIter) )
      endif
   else if (MixId == 2) then
      allocate( a_r(NumDGAIter+1,NumDGAIter+1), b_r(NumDGAIter+1,NumDGAIter+1) )
      a_r = ZERO; b_r = ZERO
   endif
!
   InitializedWkSpace = .true.
!
   end subroutine setWorkingSpace_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  setWorkingSpace_c( MixId )
!  ==================================================================
!
   implicit  none
!
   integer(kind=IntKind), intent(in) :: MixId
!
   if (MixId == 1) then
      allocate( a_c( NumBroydenIter,NumBroydenIter ) )
      allocate( b_c( NumBroydenIter,NumBroydenIter ) )
      allocate( d_c( NumBroydenIter,NumBroydenIter ) )
      allocate( cm_c( NumBroydenIter ) )
      a_c = CZERO; b_c = CZERO; d_c = CZERO; cm_c = CZERO
!
      if ( BroydenInvMethod == 1 ) then
         allocate( ipiv(NumBroydenIter) )
      endif
   else if (MixId == 2) then
      allocate( a_c(NumDGAIter+1,NumDGAIter+1), b_c(NumDGAIter+1,NumDGAIter+1) )
      a_c = CZERO; b_c = CZERO
   endif
!
   InitializedWkSpace = .true.
!
   end subroutine setWorkingSpace_c
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  delWorkingSpace()
!  ==================================================================
!
   implicit  none
!
   if ( MixingMethod == 1 ) then
      if ( allocated(a_r) ) deallocate( a_r, b_r, d_r, cm_r )
      if ( allocated(a_c) ) deallocate( a_c, b_c, d_c, cm_c )
!
      if ( BroydenInvMethod == 1 ) then
         if ( allocated(ipiv) ) deallocate( ipiv )
      endif
   else if (MixingMethod == 2) then
      if ( allocated(a_r) ) deallocate( a_r, b_r )
      if ( allocated(a_c) ) deallocate( a_c, b_c )
   endif
!
   InitializedWkSpace = .false.
!
   end subroutine delWorkingSpace
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  delMixing()
!  ==================================================================
!
   implicit  none
!
   call delWorkingSpace()
!
   if ( MixingMethod == 0 ) then
      return
   else if ( MixingMethod == 1 ) then
      call delBroydenMixing()
!!xDx!! modified for test_RPM_c_3 typ=14
   else if ( MixingMethod == 3 ) then
      call delRPMMixing()
!!xDx!! modified for test_RPM_c_3 typ=14
   else
      call delDGAMixing()
   endif
!
   end subroutine delMixing
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  resetMixing_r(pMixingList)
!  ==================================================================
!
   implicit none
!
   type (MixListRealStruct), target:: pMixingList
!
   type (MixListRealStruct), pointer:: pML
!
   pML => pMixingList
   do while ( associated(pML) )
      pML%size = 0
      pML%rms = ZERO
      nullify(pML%mesh)
      nullify(pML%vector_old)
      nullify(pML%vector_new)
      pML => pML%next
   end do
!
   call delWorkingSpace()
   call delMixing()
!
   end subroutine resetMixing_r
!  ==================================================================
!
!  ******************************************************************
!
!  ==================================================================
   subroutine  resetMixing_c(pMixingList)
!  ==================================================================
!
   implicit none
!
   type (MixListCmplxStruct), target:: pMixingList
!
   type (MixListCmplxStruct), pointer:: pML
!
   pML => pMixingList
   do while ( associated(pML) )
      pML%size = 0
      pML%rms = ZERO
      nullify(pML%mesh)
      nullify(pML%vector_old)
      nullify(pML%vector_new)
      pML => pML%next
   end do
!
   call delWorkingSpace()
   call delMixing()
   iter_count = 0
!
   end subroutine resetMixing_c
!  ==================================================================
end module MixingModule
