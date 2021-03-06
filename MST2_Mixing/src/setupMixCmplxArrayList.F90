!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setupMixCmplxArrayList( LocalNumAtoms, n_spin_pola,      &
                                      CmplxArrayList, r_rms, p_rms )
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use MathParamModule, only : ZERO, HALF, CZERO
!
   use MPPModule, only : GlobalSum
!
   use ScfDataModule, only : isPotentialMixing, isChargeMixing
!
   use RadialGridModule, only : getNumRmesh, getRmesh, getGrid
!
   use DataServiceCenterModule, only : getDataStorage,                 &
                                       ComplexType, ComplexMark,       &
                                       createDataStorage,              &
                                       isDataStorageExisting,          &
                                       setDataStorageLDA,              &
                                       setDataStorage2Value,           &
                                       deleteDataStorage
!
   use PotentialModule, only : getOldPotential => getPotential
   use PotentialModule, only : getOldVdif => getVdif
!
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
!
   use ChargeDensityModule, only : isChargeComponentZero, getRhoLmax
   use ChargeDensityModule, only : getChargeComponentFlag, getChargeDensity
   use ChargeDensityModule, only : getMomentDensity
   use ChargeDensityModule, only : getVPCharge, getVPMomSize
!
   use PotentialGenerationModule, only : getNewPotential => getPotential
   use PotentialGenerationModule, only : getNewVdif => getVdif
   use PotentialGenerationModule, only : isPotComponentZero, getPotLmax
   use PotentialGenerationModule, only : getPotComponentFlag
!
   use PublicTypeDefinitionsModule, only : MixListCmplxStruct
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
   type (MixListCmplxStruct), target :: CmplxArrayList
!
   integer (kind=IntKind), save :: data_size_save = 0
!
   integer (kind=IntKind) :: id, ia, nr, is, data_size, data_size_ns, reset_flag
   integer (kind=IntKind) :: lmax, jl, jl_nonZero, ind_jl, jmax, MaxSpecies
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind) :: factor
   real (kind=RealKind), pointer :: rptmp(:)
!
   complex (kind=CmplxKind), pointer :: ptmp1(:), ptmp2(:)
   complex (kind=CmplxKind), pointer :: pStore_old(:,:,:), pStore_new(:,:,:)
   type (MixListCmplxStruct), pointer :: p_CAL
!
   data_size = 0
   MaxSpecies = 0
   do id = 1,LocalNumAtoms
      nr = getNumRmesh(id)
      if ( isPotentialMixing() ) then
         lmax = getPotLmax(id)
         flag_jl => getPotComponentFlag(id)
      else
         lmax = getRhoLmax(id)
         flag_jl => getChargeComponentFlag(id)
      endif
      jl_nonZero = 0
      jmax = (lmax+1)*(lmax+2)/2
      do jl = 1, jmax
         if ( flag_jl(jl) /= 0 ) then
            jl_nonZero = jl_nonZero + 1
         endif
      enddo
!      jl_nonZero = jmax
      data_size = max(data_size,nr*jl_nonZero)
      MaxSpecies = max(MaxSpecies, getLocalNumSpecies(id))
   enddo
   data_size_ns = data_size*n_spin_pola+n_spin_pola
!
   reset_flag = 0
   if ( data_size_save/=0 .and. data_size_save/=data_size_ns) then
!     write(6,*) "WARNING: Mixing have been reseted !!!", data_size_save, data_size_ns
      reset_flag = 1
   endif
!  -------------------------------------------------------------------
   call GlobalSum(reset_flag)
!  -------------------------------------------------------------------
   if (reset_flag > 0) then
!     ----------------------------------------------------------------
      call resetMixing(CmplxArrayList)
      call deleteDataStorage("MixingVectorOld")
      call deleteDataStorage("MixingVectorNew")
!     ----------------------------------------------------------------
   endif
   data_size_save = data_size_ns
!
   if (.not.isDataStorageExisting('MixingVectorOld')) then
!     ----------------------------------------------------------------
      call createDataStorage('MixingVectorOld',                        &
                             LocalNumAtoms*MaxSpecies*data_size_ns,ComplexType)
      call setDataStorage2Value('MixingVectorOld',CZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
!        -------------------------------------------------------------
         call setDataStorageLDA('MixingVectorOld',data_size_ns)
!        -------------------------------------------------------------
      enddo
   endif
   if (.not.isDataStorageExisting('MixingVectorNew')) then
      call createDataStorage('MixingVectorNew',                        &
                             LocalNumAtoms*MaxSpecies*data_size_ns,ComplexType)
      call setDataStorage2Value('MixingVectorNew',CZERO)
!     ----------------------------------------------------------------
      do id = 1,LocalNumAtoms
!        -------------------------------------------------------------
         call setDataStorageLDA('MixingVectorNew',data_size_ns)
!        -------------------------------------------------------------
      enddo
   endif
!
   pStore_old => getDataStorage('MixingVectorOld',data_size_ns,        &
                                MaxSpecies,LocalNumAtoms,ComplexMark)
   pStore_new => getDataStorage('MixingVectorNew',data_size_ns,        &
                                MaxSpecies,LocalNumAtoms,ComplexMark)
!
   p_CAL => CmplxArrayList
   do id = 1, LocalNumAtoms
      nr = getNumRmesh(id)
      p_CAL%size = data_size_ns
      p_CAL%mesh => getRmesh(id)
      do ia = 1, getLocalNumSpecies(id)
         p_CAL%vector_old => pStore_old(1:data_size_ns,ia,id)
         p_CAL%vector_new => pStore_new(1:data_size_ns,ia,id)
         p_CAL%vector_old(:) = CZERO
         p_CAL%vector_new(:) = CZERO
         p_CAL%rms = ZERO
         flag_jl => getPotComponentFlag(id)
         do is = 1, n_spin_pola
            ind_jl = data_size*(is-1)+is-1
            if (isPotentialMixing()) then
               lmax = getPotLmax(id)
               jmax = (lmax+1)*(lmax+2)/2
               factor = real(2-is,kind=RealKind)
               do jl = 1, jmax
                  if ( flag_jl(jl) /= 0 ) then
                     ptmp1 => getOldPotential(id,ia,is,jl)
                     p_CAL%vector_old(ind_jl+1:ind_jl+nr) = ptmp1(1:nr)
                     ptmp1 => getNewPotential("Total",id,ia,is,jl)
                     p_CAL%vector_new(ind_jl+1:ind_jl+nr) = ptmp1(1:nr)
                     ind_jl = ind_jl+ nr
                  endif
               enddo
               rptmp => getOldVdif()
               p_CAL%vector_old(ind_jl+1) = factor*rptmp(1)
               rptmp => getNewVdif()
               p_CAL%vector_new(ind_jl+1) = factor*rptmp(1)
               if (n_spin_pola==1) then
                  p_CAL%rms = p_rms(is,id)
               else
                  p_CAL%rms = (p_rms(1,id)+p_rms(2,id))*half
               endif
            else
               lmax = getRhoLmax(id)
               jmax = (lmax+1)*(lmax+2)/2
               factor = real(3-is*2,kind=RealKind)
               flag_jl => getChargeComponentFlag(id)
               do jl = 1, jmax
                  if ( flag_jl(jl) /=0 ) then
                     ptmp1 => getChargeDensity('TotalOld', id, ia, jl)
                     if (n_spin_pola==2) then
                        ptmp2 => getMomentDensity('TotalOld', id, ia, jl)
                        p_CAL%vector_old(ind_jl+1:ind_jl+nr) = (ptmp1(1:nr)+ &
                                                    factor*ptmp2(1:nr))*half
                     else
                        p_CAL%vector_old(ind_jl+1:ind_jl+nr) = ptmp1(1:nr)
                     endif
                     ptmp1 => getChargeDensity('TotalNew', id, ia, jl)
                     if (n_spin_pola==2) then
                        ptmp2 => getMomentDensity('TotalNew', id, ia, jl)
                        p_CAL%vector_new(ind_jl+1:ind_jl+nr) = (ptmp1(1:nr)+ &
                                                    factor*ptmp2(1:nr))*half
                     else
                        p_CAL%vector_new(ind_jl+1:ind_jl+nr) = ptmp1(1:nr)
                     endif
                     ind_jl = ind_jl+ nr
                  endif
               enddo
               if (n_spin_pola==2) then
                  p_CAL%vector_old(ind_jl+1) = half*(getVPCharge('Old',id,ia)+ &
                                               factor*getVPMomSize('Old',id,ia))
                  p_CAL%vector_new(ind_jl+1) = half*(getVPCharge('New',id,ia)+ &
                                               factor*getVPMomSize('New',id,ia))
                  p_CAL%rms = r_rms(1,id)*half
               else
                  p_CAL%vector_old(ind_jl+1) = getVPCharge('Old',id,ia)
                  p_CAL%vector_new(ind_jl+1) = getVPCharge('New',id,ia)
                  p_CAL%rms = r_rms(is,id)
               endif
            endif
         enddo
         if (associated(p_CAL%next)) then
            p_CAL => p_CAL%next
         else if ( id/=LocalNumAtoms ) then
            call ErrorHandler('setupMixCmplxArrayList',                  &
                              'CmplxArrayList is not set up properly', id)
         endif
      enddo
   enddo
!   nullify(p_CAL, ptmp1, ptmp2, rptmp, flag_jl, pStore_old, pStore_new)
!
   end subroutine setupMixCmplxArrayList
