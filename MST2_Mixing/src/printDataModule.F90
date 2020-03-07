!!xDx!! This module is written by me so that I can get the procs to print
!!xDx!! form only 1 PE. This way th outputs don't get muddled up.
!!xDx!! An important point to note is that this has to be included
!!xDx!! any point above the module which is calling this in the "makefile"

!!xDx!! Things to Edit later
!!xDx!! 1. Add interface to spearate Data (integer, real, Complex)

module printDataModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use GroupCommModule, only : getMyClusterIndex, getMyPEinGroup
   use MPPModule, only : initMPP, syncAllPEs, endMPP, NumPEs, MYPE

public :: printDataReal,    &
          printDataComplex, &
          printDataRealArray,    &
          printDataComplexArray,    &
          printDataInteger, &
          printDataPointer

interface printDataRealArray
   module procedure printDataRealArray1, printDataRealArray2
end interface printDataRealArray 

interface printDataComplexArray
   module procedure printDataComplexArray1, printDataComplexArray2
end interface printDataComplexArray 


contains

	subroutine printDataInteger(Text,Data2pnt,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
	integer(4), intent(in) :: Data2pnt
!	integer (kind=IntKind), intent(in) :: Data2pnt
        integer (kind=IntKind), intent(in) :: PE, iter
	character (len=50) :: PEStr

   write(PEStr, '(I2)') MyPE+10

 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1120) '========================'
 elseif ( iter .gt. 1) then
  write (MyPE+10,1120) '========================'
 else
  print*, 'Error from printDataInteger'
 endif

   if (PE .eq. MYPE) then
        write (MyPE+10,1340) 'MyPE',Text,MyPE,Data2pnt
   endif


1190 format(X,A)
1120 format(X,A)
1340 format(2(X,A),(X,I4),(X,I4))


	end subroutine printDataInteger
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!


	subroutine printDataReal(Text,Data2pnt,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
	real (kind=RealKind), intent(in) :: Data2pnt
        integer (kind=IntKind), intent(in) :: PE, iter
	character (len=50) :: PEStr

   write(PEStr, '(I2)') MyPE+10

 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1190) '========================'
 elseif ( iter .gt. 1) then
!  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1190) '========================'
 else
  print*, 'Error from printDataRealArray1'
 endif

   if (PE .eq. MYPE) then
        write (MyPE+10,1340) 'MyPE',Text,MyPE,Data2pnt
   endif


1190 format(X,A)
1340 format(2(X,A),(X,I4),(X,F16.10))


	end subroutine printDataReal
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

	subroutine printDataComplex(Text,Data2pnt,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
	complex (kind=CmplxKind), intent(in) :: Data2pnt
        integer (kind=IntKind), intent(in) :: PE, iter
	character (len=50) :: PEStr

   write(PEStr, '(I2)') MyPE+10

 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1191) '========================'
 elseif ( iter .gt. 1) then
!  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1191) '========================'
 else
  print*, 'Error from printDataRealArray1'
 endif

   if (PE .eq. MYPE) then
        write (MyPE+10,1341) 'MyPE',Text,MyPE,real(Data2pnt),aimag(Data2pnt),'*i'
   endif


1191 format(X,A)
1341 format(2(X,A),(X,I4),(X,F16.10,SP,F16.10,A))

	end subroutine printDataComplex
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

	subroutine printDataRealArray1(Text,Data2pnt,n,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
        integer (kind=IntKind), intent(in) :: n
	real (kind=RealKind), intent(in) :: Data2pnt(1:n)
        integer (kind=IntKind), intent(in) :: PE, iter
        integer (kind=IntKind) :: i
	character (len=50) :: PEStr

   write(PEStr, '(I2)') MyPE+10

 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
 elseif ( iter .gt. 1) then
!  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1190) '========================'
 else
  print*, 'Error from printDataRealArray1'
 endif

   if (PE .eq. MYPE) then
    do i =1,n
        write (MyPE+10,1240) 'i,MyPE',Text,i,MyPE,Data2pnt(i)
    enddo
   endif


1190 format(X,A)
1240 format(2(X,A),2(X,I4),(X,F16.10))
	end subroutine printDataRealArray1
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

	subroutine printDataComplexArray1(Text,Data2pnt,n,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
        integer (kind=IntKind), intent(in) :: n
	complex (kind=CmplxKind), intent(in) :: Data2pnt(1:n)
        integer (kind=IntKind), intent(in) :: PE, iter
        integer (kind=IntKind) :: i
	character (len=50) :: PEStr

   write(PEStr, '(I2)') MyPE+10

 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
	print *, 'Hello form Complx Array1, iter=1'
 elseif ( iter .gt. 1) then
!  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1191) '========================'
 else
  print*, 'Error from printDataComplexArray1'
 endif

   if (PE .eq. MYPE) then
    do i =1,n
        write (MyPE+10,1241) 'i,MyPE',Text,i,MyPE,real(Data2pnt(i)),aimag(Data2pnt(i)),'*i'
    enddo
   endif


1191 format(X,A)
1241 format(2(X,A),2(X,I4),(X,F20.10,SP,F16.10,A))
	end subroutine printDataComplexArray1
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

	subroutine printDataRealArray2(Text,Data2pnt,n1,n2,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
        integer (kind=IntKind), intent(in) :: n1, n2
	real (kind=RealKind), intent(in) :: Data2pnt(1:n1,1:n2)
        integer (kind=IntKind), intent(in) :: PE
        integer (kind=IntKind), intent(in) :: iter
        integer (kind=IntKind) :: i,j
	character (len=50) :: PEStr

	write(PEStr,'(I2)') MyPE+10


 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
 elseif ( iter .gt. 1) then
!  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1190) '========================'
 else
  print*, 'Error from printDataRealArray2'
 endif

   if (PE .eq. MYPE) then
      do i =1,n1
         write (MyPE+10,1230) ( Data2pnt(i,j), j=1,n2 )
      enddo
   endif

!!cDc!!   if (PE .eq. MYPE) then
!!cDc!!      do i =1,n1
!!cDc!!       do j =1,n2
!!cDc!!	write (MyPE+100,1220) 'i,j,MyPE,',Text,i,j,MyPE,Data2pnt(i,j)
!!cDc!!!	write (MyPE+100,1230) ( Data2pnt(i,j), j=1,n2 )
!!cDc!!       enddo
!!cDc!!      enddo
!!cDc!!   endif

1190 format(X,A)
1210 format(X,A,I1,X,A)
1220 format(2(X,A),3(X,I),f16.10)
1230 format(14f16.6)

	end subroutine printDataRealArray2
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

	subroutine printDataComplexArray2(Text,Data2pnt,n1,n2,PE,iter)

	implicit none

	character (len=*), intent(in), optional :: Text
        integer (kind=IntKind), intent(in) :: n1, n2
	complex (kind=CmplxKind), intent(in) :: Data2pnt(1:n1,1:n2)
        integer (kind=IntKind), intent(in) :: PE
        integer (kind=IntKind), intent(in) :: iter
        integer (kind=IntKind) :: i,j
	character (len=50) :: PEStr
	character (len=19) fmt

	write(PEStr,'(I2)') MyPE+10

 if ( iter .eq. 1) then
  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
 elseif ( iter .gt. 1) then
!  open ( MyPE+10,file='MyFilefor'//TRIM(PEStr)//'.out' )
  write (MyPE+10,1191) '========================'
 else
  print*, 'Error from printDataComplexArray2'
 endif

   fmt = '(F7.2,"+",F7.2,"i")'
   if (PE .eq. MYPE) then
      do i =1,n1
        do j =1,n2
           fmt(8:8) = MERGE('+',' ',imag(Data2pnt(i,j)).gt.0)
           write(MyPE+10,fmt, advance="no") Data2pnt(i,j)
!!?!!           write (MyPE+10,1231) ( real(Data2pnt(i,j)), aimag(Data2pnt(i,j)),'*i', j=1,n2 )
        enddo
        write(*, fmt="(a)") " " ! This should be--> write(MyPE+10, fmt="(a)") " "
      enddo
   endif

1191 format(X,A)
1210 format(X,A,I1,X,A)
1220 format(2(X,A),3(X,I),f16.10)
1231 format(14f16.6,SP,14f10.6,A)
	end subroutine printDataComplexArray2
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!


!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!
	subroutine printDataPointer(Text,Data2pnt,gid,print_pe)

	implicit none

	character (len=*), intent(in) :: Text
	real (kind=RealKind), pointer, intent(in) :: Data2pnt(:)   
        integer (kind=IntKind), intent(in) :: gid, print_pe


   if (getMyPEinGroup(gid) == print_pe) then
	print *, Text
	print *, Data2pnt(:)
   endif


	end subroutine printDataPointer
!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!!!xDx!!

end module printDataModule
!!xDx!!
