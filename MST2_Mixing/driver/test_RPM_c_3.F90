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
!!xDx!!================================================================================================
!!xDx!!------------------------------Important Note----------------------------------------------------
!!xDx!!================================================================================================
!!xDx!!1.  In order to run this piece of code, we needed to make changes to the MixingModule. This
!!xDx!!    Mixing has been saved under MixingModule_RPM.F90
!!xDx!!    Need to add MixingModule_RPM.F90 to the makefile of src, else the program will not run 
!!xDx!!    because it wont' find an .o file for MixingModule_RPM. 
!!xDx!!    (Most prob, MixingModule_RPM.o)
!!xDx!!
!!xDx!!2.  As a temporary fix to get the RPM Module to work, removed test_RPM_c_3 from the makefile
!!xDx!!    in driver. In order to run it again, here is the necessary steps.
!!xDx!!    a. Under TEST_CASES = ....
!!xDx!!       $(ODIR)/test_RPM_c_3 \
!!xDx!!    b. Under TEST_CASES_O = ....
!!xDx!!       $(ODIR)/test_RPM_c_3.o \
!!xDx!!    c. uncomment the occurance of '$(ODIR)/test_RPM_c_3...', ie, the main make.
             


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
      use MixingModule, only : initMixing, endMixing, mixValues,       &
                               setRPMMixing, setBroydenMixing
!
!
      implicit none

      integer (kind=IntKind) :: typ=14, MaxRPMIter, vvlen, Basis
      integer (kind=IntKind) :: NumBroydenIter, TotalRPMIter
      integer (kind=IntKind) :: LocalNumAtoms, NumAtoms, NumTotalMix
      integer (kind=IntKind) :: MaxBasis, NumRPMIter_broy
      integer (kind=IntKind) :: RPMIIter
!
      integer (kind=IntKind) :: k, ii
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
 !!xDx!! Message passing                     
      integer(kind=IntKind) :: GroupID     
      complex (kind=CmplxKind), allocatable :: msgbuf1(:)          
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
     call createBoxGroupFromGrid(grid,box,'Unit Cell')
     GroupID = getGroupID('Unit Cell')
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
         vvlen=4; Basis=0; MaxRPMIter=6; NumBroydenIter=10; alpha= 1D-01;
         beta=1-alpha; tol = 1D-09; TotalRPMIter=500; LocalNumAtoms=2
      else
         print *, 'Please eneter a vaild typ'
      endif

      NumAtoms=NumPEs*LocalNumAtoms
      NumTotalMix=LocalNumAtoms !! by defn also =NumMix, i/p to initMixing_c 
      nr=vvlen
      MaxBasis=6
      NumRPMIter_broy=MaxRPMIter-2

!!nDn!! Based on how the LSMS code is set up, we have 
!!nDn!! N= #(Atoms/box)x#(RadialMesh)x#(Lvalues)x#(PEs/box)
!!nDn!! for typ=4 (N=8), we want the parameters to be
!!nDn!! #(Atoms/box)=2, #(RadialMesh)=2, #(Lvalues)=1, #(PEs/box)=2

!!xDx!!Allocation arrays
      allocate(pn(NumPEs*LocalNumAtoms*vvlen))
      allocate(po(NumPEs*LocalNumAtoms*vvlen))


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
          pn(9) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(10) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(11) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(12) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(13) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(14) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(15) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
          pn(16) = cmplx( 0.9975D0, 0.003364D0)!complex(1.2975000143051100,-0.7733640074729920)
      endif
   
         NumMix(1)=LocalNumAtoms

         call initMixing( 1, NumMix, CmplxArrayList )

!!xDx!! setupMixingScheme
      is = 1
!  do is = 1,n_spin_pola
      do id=1,LocalNumAtoms
            call setRPMMixing(is, id, alpha )
!            call setBroydenMixing(is, id, alpha )
      enddo

      ITERAA:   do RPMIIter = 1,TotalRPMIter!!1a.!!1,MaxRPMIter, is for the rest.
     
	!!cDc!!
		print *, 'iter',RPMIIter
	!!cDc!!
 
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
           pn(9)=sqrt(po(9)+2)
           pn(10)=sqrt(po(10)+2)
           pn(11)=sqrt(po(11)+2)
           pn(12)=sqrt(po(12)+2)
           pn(13)=sqrt(po(13)+2)
           pn(14)=sqrt(po(14)+2)
           pn(15)=sqrt(po(15)+2)
           pn(16)=sqrt(po(16)+2)
         endif     

!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (10,1122) '========================'
!!cDc!!,       write (10,1342) 'main b4 pn(:)','MyPE,RPMIIter, id', MyPE, RPMIIter, id
!!cDc!!,       write (10,1122) '========================'
!!cDc!!,       do ii =1,NumPEs*vvlen*LocalNumAtoms
!!cDc!!,          write (10,1241) 'ii,MyPE','pn',ii,MyPE,real(pn(ii)),aimag(pn(ii)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, 1241 format(2(X,A),2(X,I4),(X,F20.10,SP,F16.10,A))
!!cDc!!, enddo
!!cDc!!, !!cDc!!


     !!xDx!! setupMixCmplxArrayList
            data_size = 0
            MaxSpecies = 0
            do id = 1,LocalNumAtoms
               nr = vvlen !!nr=2 for typ=4
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
               nr = vvlen !!=2 (typ=4)
               p_CAL%size = data_size_ns
               p_CAL%vector_old => pStore_old(1:data_size_ns,id)
               p_CAL%vector_new => pStore_new(1:data_size_ns,id)
               p_CAL%vector_old(:) = CZERO
               p_CAL%vector_new(:) = CZERO
               p_CAL%rms = ZERO
               p_CAL%vector_old(1:vvlen) = po(MyPE*vvlen*LocalNumAtoms+(id-1)*vvlen+1:MyPE*vvlen*LocalNumAtoms+id*vvlen) 
               p_CAL%vector_new(1:vvlen) = pn(MyPE*vvlen*LocalNumAtoms+(id-1)*vvlen+1:MyPE*vvlen*LocalNumAtoms+id*vvlen) 
      !
               if (associated(p_CAL%next)) then
                  p_CAL => p_CAL%next
               else if ( id/=LocalNumAtoms ) then
                  call ErrorHandler('setupMixCmplxArrayList',          &
                            'CmplxArrayList is not set up properlyi',id)
               endif
            enddo
      !!xDx!!----------------------------------------------------------------------------------------------
      !!xDx!! calling Mixing-------------------------------------------------------------------------------
            call mixValues(CmplxArrayList)
      !!xDx!! update Mix CmplxValues----------------------------------------------------------------------
         data_size = 0
         MaxSpecies = 0
         do id = 1,LocalNumAtoms
            nr = vvlen
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
         pn = CZERO
         do id = 1, LocalNumAtoms
            nr = vvlen
            p_CAL%vector_old => pStore_old(1:data_size_ns,id)
            p_CAL%vector_new => pStore_new(1:data_size_ns,id)
            pn(MyPE*vvlen*LocalNumAtoms+(id-1)*vvlen+1:MyPE*vvlen*            & 
            LocalNumAtoms+id*vvlen) = p_CAL%vector_new(1:vvlen)
            p_CAL => p_CAL%next
         enddo
      !!xDx!!
  
            ! Allocating initial guess to the different PEs
            allocate(msgbuf1(1:NumPEs*LocalNumAtoms*vvlen))
            msgbuf1 = pn
            call GlobalSumInGroup(GroupID, msgbuf1, NumPEs*LocalNumAtoms*vvlen)
            pn = msgbuf1
            deallocate(msgbuf1)

!!cDc!!, !!cDc!!
!!cDc!!, do k=0,NumPEs-1
!!cDc!!,    if (MyPE .eq. k) then
!!cDc!!,       write (10,1122) '========================'
!!cDc!!,       write (10,1342) 'main aftr pn(:)','MyPE,RPMIIter, id', MyPE, RPMIIter, id
!!cDc!!,       write (10,1122) '========================'
!!cDc!!,       do ii =1,NumPEs*vvlen*LocalNumAtoms
!!cDc!!,          write (10,1241) 'ii,MyPE','pn',ii,MyPE,real(pn(ii)),aimag(pn(ii)),'*i'
!!cDc!!,       enddo
!!cDc!!,    endif
!!cDc!!, enddo
!!cDc!!, !!cDc!!

!           ============================================================
            normp=ZERO
            do k = 1,NumPEs*LocalNumAtoms*vvlen !!xDx!! These are sums over all the atoms in the proc. 
               normp  = normp  + real((pn(k)-po(k))*conjg(pn(k)-po(k)),kind=RealKind) 
            enddo

            if ( (int(typ/10,kind=IntKind) .eq. 1) .and. (normp .lt. tol) ) then
               if ( MyPE .eq. 0 ) then
                  write (*,1121) '========================'
                  write (*,1121) '!!xDx!!!xDx!!!xDx!!!xDx!!!xDx!!'
                  write (*,1121) '========================'
                  write (*,1121) 'Solution Converged with'
                  write (*,1341) 'MyPE','normp,iter',MyPE,normp,RPMIIter
                  write (*,1342) 'MyPE','Basis',MyPE,Basis
                  write (*,1122) '========================'
                  write (*,1122) '========================'
                  write (*,1342) 'MyPE','pn(:),RPMIIter',MyPE,RPMIIter
                  write (*,1122) '========================'
                  do k = 1,NumPEs*LocalNumAtoms*vvlen 
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

      !!xDx!! De-allocating the Mixing Parameters.
      call endMixing()

      !!Deallocation
      deallocate(pn)
      deallocate(po)
      
      call endGroupComm()
      call endMPP()
      
         stop 'Ok'
      end program RPM_2



