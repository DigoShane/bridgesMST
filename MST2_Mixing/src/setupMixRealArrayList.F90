!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupMixRealArrayList( LocalNumAtoms, n_spin_pola,       &
                                     RealArrayList, r_rms, p_rms )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MathParamModule, only : ZERO, HALF
!
   use MPPModule, only : GlobalSum
!
   use ScfDataModule, only : isPotentialMixing, isChargeMixing
!
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
   use DataServiceCenterModule, only : getDataStorage, RealType, RealMark, &
                                       createDataStorage,               &
                                       setDataStorage2Value,            &
                                       isDataStorageExisting,           &
                                       deleteDataStorage
!
   use PotentialModule, only : getOldSphPotr => getSphPotr
   use PotentialModule, only : getOldVdif => getVdif
!
   use AtomModule, only : getLocalNumSpecies
!
   use PotentialGenerationModule, only : getNewSphPotr => getSphPotr
   use PotentialGenerationModule, only : getNewVdif => getVdif
!
   use ChargeDensityModule, only : getSphChargeDensity, &
                                   getSphMomentDensity
!
   use PublicTypeDefinitionsModule, only : MixListRealStruct
!
   use MixingModule, only : resetMixing
!
   use WriteMatrixModule, only : writeMatrix
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms,n_spin_pola
   real (kind=RealKind), intent(in) :: r_rms(n_spin_pola,LocalNumAtoms)
   real (kind=RealKind), intent(in) :: p_rms(n_spin_pola,LocalNumAtoms)
!
   type (MixListRealStruct), target :: RealArrayList
!
   integer (kind=IntKind), save :: data_size_save = 0
!
   integer (kind=IntKind) :: id, ia, nr, is, data_size, MaxSpecies, reset_flag
!
   real (kind=RealKind) :: factor
   real (kind=RealKind), pointer :: pStore_old(:,:,:), pStore_new(:,:,:)
   real (kind=RealKind), pointer :: ptmp1(:), ptmp2(:)
!
   type (MixListRealStruct), pointer :: p_RAL
!
!!xDx!! initializing variables for Data storage
   data_size = 0; MaxSpecies = 0
   do id = 1,LocalNumAtoms
      nr = getNumRmesh(id)+1
      data_size = max(data_size,nr)
      MaxSpecies = max(MaxSpecies,getLocalNumSpecies(id))
   enddo
   data_size = data_size*n_spin_pola
!
!!xDx!! flag check to ensure proper data storage
   reset_flag = 0
   if ( data_size_save/=0 .and. data_size_save/=data_size) then
      reset_flag = 1
   endif
!  -------------------------------------------------------------------
   call GlobalSum(reset_flag)
!  -------------------------------------------------------------------
   if (reset_flag > 0) then
      call resetMixing(RealArrayList)
      call deleteDataStorage("MixingVectorOld")
      call deleteDataStorage("MixingVectorNew")
   endif
   data_size_save = data_size
!
!!xDx!! Checking existence of space for mixing vectors.
   if (.not.isDataStorageExisting("MixingVectorOld")) then
      call createDataStorage('MixingVectorOld',LocalNumAtoms*MaxSpecies*data_size,RealType)
      call setDataStorage2Value('MixingVectorOld',ZERO)
   endif
   if (.not.isDataStorageExisting("MixingVectorNew")) then
      call createDataStorage('MixingVectorNew',LocalNumAtoms*MaxSpecies*data_size,RealType)
      call setDataStorage2Value('MixingVectorNew',ZERO)
   endif
!!!,,!!,!!cDc!! Creating a DataStorage For the RPM variable qloc 
!!!,,!!,   if (.not.isDataStorageExisting("MixingVectorOld")) then
!!!,,!!,      call createDataStorage('MixingVectorOld',LocalNumAtoms*MaxSpecies*data_size,RealType)
!!!,,!!,      call setDataStorage2Value('MixingVectorOld',ZERO)
!!!,,!!,   endif
!!!,,!!,   if (.not.isDataStorageExisting("MixingVectorNew")) then
!!!,,!!,      call createDataStorage('MixingVectorNew',LocalNumAtoms*MaxSpecies*data_size,RealType)
!!!,,!!,      call setDataStorage2Value('MixingVectorNew',ZERO)
!!!,,!!,   endif
!!!,,!!,!!cDc!!
!
!!xDx!! Point to appropriate data.
   pStore_old => getDataStorage("MixingVectorOld",data_size, &
                                MaxSpecies,LocalNumAtoms,RealMark)
   pStore_new => getDataStorage("MixingVectorNew",data_size, &
                                MaxSpecies,LocalNumAtoms,RealMark)
!!xDx!! initializing the targets
   p_RAL => RealArrayList
   do id = 1, LocalNumAtoms
      nr = getNumRmesh(id)+1
      p_RAL%mesh => getRmesh(id)
      p_RAL%size = nr*n_spin_pola
      do ia = 1, getLocalNumSpecies(id)
         p_RAL%vector_old => pStore_old(1:nr*n_spin_pola,ia,id)
         p_RAL%vector_new => pStore_new(1:nr*n_spin_pola,ia,id)
         p_RAL%vector_old(:) = ZERO
         p_RAL%vector_new(:) = ZERO
         p_RAL%rms = ZERO
         do is = 1, n_spin_pola
            if (isPotentialMixing()) then
               factor = real(2-is,kind=RealKind)
               p_RAL%vector_old(nr*(is-1)+1:nr*is-1) = getOldSphPotr(id,ia,is)
               p_RAL%vector_old(nr*is:nr*is) = factor*getOldVdif()
               p_RAL%vector_new(nr*(is-1)+1:nr*is-1) = getNewSphPotr(id,ia,is)
               p_RAL%vector_new(nr*is:nr*is) = factor*getNewVdif()
!              if (id == 1) then
!                 write(6,*) 'New Vr: spin =',is
!                 call writeMatrix('New Vr',p_RAL%vector_new(nr*(is-1)+1:nr*is),nr)
!              endif
               if (n_spin_pola==1) then
                 p_RAL%rms = p_rms(is,id)
               else
                 p_RAL%rms = (p_rms(1,id)+p_rms(2,id))*half
               endif
            else
               factor = real(3-is*2,kind=RealKind)
               ptmp1 => getSphChargeDensity('TotalOld', id, ia)
               if (n_spin_pola==2) then
                  ptmp2 => getSphMomentDensity('TotalOld', id, ia)
               endif
               if (n_spin_pola==2) then
                  p_RAL%vector_old(nr*(is-1)+1:nr*is) = (ptmp1+factor*ptmp2)*half
               else
                  p_RAL%vector_old(1:nr) = ptmp1(1:nr)
               endif
               ptmp1 => getSphChargeDensity('TotalNew', id, ia)
               if (n_spin_pola==2) then
                  ptmp2 => getSphMomentDensity('TotalNew', id, ia)
               endif
               if (n_spin_pola==2) then
                  p_RAL%vector_new(nr*(is-1)+1:nr*is) = (ptmp1+factor*ptmp2)*half
               else
                  p_RAL%vector_new(1:nr) = ptmp1(1:nr)
               endif
               if (n_spin_pola==1) then
                 p_RAL%rms = r_rms(is,id)
               else
                 p_RAL%rms = r_rms(1,id)*half
               endif
            endif
         enddo
         if (associated(p_RAL%next)) then
            p_RAL => p_RAL%next
         else if ( id /= LocalNumAtoms ) then
            call ErrorHandler('setupMixCmplxArrayList',                  &
                              'CmplxArrayList is not set up properly')
         endif
      enddo
   enddo
!
   end subroutine setupMixRealArrayList
