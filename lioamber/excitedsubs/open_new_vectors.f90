subroutine open_new_vectors(Ra,Rb,Esub,EallA,EallB,Ta,Tb,Eold,&
                            Ndima,Ndimb,Ndim,Sdim,Nstat,M,NCOa,&
                            NCOb,New,conv)
   implicit none

   integer, intent(in) :: Ndima, Ndimb, Ndim, Sdim, Nstat, M, &
                          NCOa, NCOb
   real*8, intent(in) :: Ra(Ndim,Nstat), Rb(Ndim,Nstat)
   real*8, intent(in) :: Esub(Nstat), EallA(M), EallB(M)
   real*8, intent(inout) :: Ta(Ndim,Sdim+Nstat), Tb(Ndim,Sdim+Nstat)
   real*8, intent(inout) :: Eold(Nstat)
   logical, intent(out) :: conv
   integer, intent(out) :: New

   integer :: i, iv, occ, virt, Nvirt, ind, NCOc
   integer :: dim_tot, j, k, append ! puede que galla q borrar esto
   real*8 :: temp, ERROR_A, ERROR_B, ERROR, MAX_ERROR, MAX_ENE, &
             tolv, tole, diffE, prodA, prodB, prodT, norm_A, norm_B
   real*8, dimension(:,:), allocatable :: QvecA, QvecB, Ttot
   integer, dimension(:), allocatable :: valid_id

   allocate(QvecA(Ndim,Nstat),QvecB(Ndim,Nstat))
   allocate(valid_id(Nstat))

   QvecA = 0.0D0; QvecB = 0.0D0
   MAX_ERROR = 0.0D0
   MAX_ENE = 0.0D0
   tolv = 1.0D-5
   tole = 1.0D-5
   New = 0

   do iv=1,Nstat
     diffE = abs(Eold(iv) - Esub(iv))
     if(diffE > MAX_ENE) MAX_ENE = diffE
     call norma(Ra(1:Ndima,iv),Ndima,ERROR_A)
     call norma(Rb(1:Ndimb,iv),Ndimb,ERROR_B)
     
     ERROR = ERROR_A + ERROR_B
     if(ERROR > MAX_ERROR) MAX_ERROR = ERROR

     ! CALCULATE NEW TRIALS VECTORS
     if(ERROR > tolv .or. diffE > tole) then
        ! ALPHA
        New = New + 1
        valid_id(New) = iv
        Nvirt = M - NCOa
        NCOc = NCOa + 1
        do i=1,Ndima
          ind = i - 1
          occ = NCOa - (ind/Nvirt)
          virt = mod(ind,Nvirt) + NCOc
          temp = EallA(virt) - EallA(occ)
          temp = 1.0D0 / (Esub(iv) - temp)
          QvecA(i,iv) = temp * Ra(i,iv)
        enddo
        ! BETA
        Nvirt = M - NCOb
        NCOc = NCOb + 1
        do i=1,Ndimb
          ind = i - 1
          occ = NCOb - (ind/Nvirt)
          virt = mod(ind,Nvirt) + NCOc
          temp = EallB(virt) - EallB(occ)
          temp = 1.0D0 / (Esub(iv) - temp)
          QvecB(i,iv) = temp * Rb(i,iv)
        enddo
     else
        write(*,"(1X,A,I2,1X,A)") "Vector:",iv,"Converged"
     endif
   enddo ! ENDDO STATES

   write(*,8070) MAX_ERROR, tolv, MAX_ENE, tole

!  CHECK CONVERGENCE
   if(Sdim + New >= Ndima+Ndimb) then
      conv = .true.
      return
   endif

   if(MAX_ERROR < tolv .and. MAX_ENE < tole) then
     conv = .true.
   else ! APPEND NEW VECTORS
     conv = .false.
     append = 0

     ! ORTHONORMALIZATION
     do iv=1,New
       do j=1,Sdim+append
         call prod(Ta(1:Ndima,j),QvecA(1:Ndima,valid_id(iv)),norm_A,Ndima)
         call prod(Tb(1:Ndimb,j),QvecB(1:Ndimb,valid_id(iv)),norm_B,Ndimb)
         QvecA(1:Ndima,valid_id(iv))=QvecA(1:Ndima,valid_id(iv)) - (norm_A+norm_B)*Ta(1:Ndima,j)
         QvecB(1:Ndimb,valid_id(iv))=QvecB(1:Ndimb,valid_id(iv)) - (norm_A+norm_B)*Tb(1:Ndimb,j)
       enddo
       call norma(QvecA(1:Ndima,valid_id(iv)),Ndima,ERROR_A)
       call norma(QvecB(1:Ndimb,valid_id(iv)),Ndimb,ERROR_B)
       ERROR = dsqrt(ERROR_A + ERROR_B)
       ERROR = 1.0D0 / ERROR
       Ta(1:Ndima,Sdim+iv) = ERROR * QvecA(1:Ndima,valid_id(iv))
       Tb(1:Ndimb,Sdim+iv) = ERROR * QvecB(1:Ndimb,valid_id(iv))
       append = append + 1
     enddo
   endif

   Eold = Esub
   deallocate(QvecA,QvecB,valid_id)

8070   FORMAT(1X,"VectorsError (crit) = ", F15.7," (",ES9.2,")" &
              " - EnergyError (crit) = ", F15.7," (",ES9.2,")" )
end subroutine open_new_vectors
