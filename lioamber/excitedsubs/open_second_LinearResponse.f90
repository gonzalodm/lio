subroutine open_second_LinearResponse(Calpha0,Cbeta0,Ealpha0,Ebeta0, &
                              M,NCOa,NCOb)
use lrdata, only: lambda_LR, state_LR, Xflr, Eflr, Xslr, Eslr, &
                  nstates, Tdip_flr, Tdip_slr
use garcha_mod, only: rhoalpha, rhobeta
use liosubs_dens, only: builds_densmat, messup_densmat
! C_0  = Molecular Orbitals Unperturbed
! E_0  = Molecular Energies Unperturbed
   implicit none

   integer, intent(in) :: M, NCOa, NCOb
   real*8,  intent(in) :: Calpha0(M,M), Ealpha0(M)
   real*8,  intent(in) :: Cbeta0(M,M), Ebeta0(M)

   character(len=4) :: ii_char, jj_char
   integer :: ii, jj, NCOc, ele, N_slr, Nvirta, Nvirtb, Ndima, Ndimb, Ndim
   real*8, allocatable, dimension(:,:,:) :: Dens1
   real*8, allocatable, dimension(:,:) :: Calpha1, Cbeta1, Tdip
   real*8, allocatable, dimension(:) :: temp, Ene1, Osc

   print*, "Using lambda:", lambda_LR
   ! Deallocate Matrix form with 0
   call open_basis_deinitLR()
  
   Nvirta = M - NCOa
   Nvirtb = M - NCOb
   allocate(temp(M))
   ! Form Perturbed Molecular Orbitals
   ! $$$$ OCCUPIED ALPHA
   allocate(Calpha1(M,M))
   NCOc = NCOa + 1
   do ii=1,NCOa
      temp = 0.0d0
      do jj=1,Nvirta
         ele = (ii - 1) * Nvirta + jj
         temp = temp + Xflr(ele,state_LR,1) * Calpha0(:,NCOa+jj)
      enddo
      Calpha1(:,NCOc-ii) = Calpha0(:,NCOc-ii) + lambda_LR * temp
   enddo
   ! $$$$ VIRTUAL ALPHA
   do ii=1,Nvirta
      temp = 0.0d0
      do jj=1,NCOa
         ele = (jj-1) * Nvirta + ii
         temp = temp + Xflr(ele,state_LR,1) * Calpha0(:,NCOc-jj)
      enddo
      Calpha1(:,NCOa+ii) = Calpha0(:,NCOa+ii) - lambda_LR * temp
    enddo

   ! $$$$ OCCUPIED BETA
   allocate(Cbeta1(M,M))
   NCOc = NCOb + 1
   do ii=1,NCOb
      temp = 0.0d0
      do jj=1,Nvirtb
         ele = (ii - 1) * Nvirtb + jj
         temp = temp + Xflr(ele,state_LR,2) * Cbeta0(:,NCOb+jj)
      enddo
      Cbeta1(:,NCOc-ii) = Cbeta0(:,NCOc-ii) + lambda_LR * temp
   enddo
   ! $$$$ VIRTUAL BETA
   do ii=1,Nvirtb
      temp = 0.0d0
      do jj=1,NCOb
         ele = (jj-1) * Nvirtb + ii
         temp = temp + Xflr(ele,state_LR,2) * Cbeta0(:,NCOc-jj)
      enddo
      Cbeta1(:,NCOb+ii) = Cbeta0(:,NCOb+ii) - lambda_LR * temp
   enddo
   deallocate(temp)

   ! BASIS INITIALIZATION
   call open_basis_initLR(Calpha1,Cbeta1,NCOa,NCOb,M)

   ! FORM PERTURBED DENSITY
   allocate(Dens1(M,M,2)); Dens1 = 0.0d0 ! 1=alpga, 2=beta
   call builds_densmat(M, NCOa, 1.0d0, Calpha1, Dens1(:,:,1))
   call messup_densmat(Dens1(:,:,1))
   call sprepack('L',M,rhoalpha,Dens1(:,:,1))
   call builds_densmat(M, NCOb, 1.0d0, Cbeta1, Dens1(:,:,2))
   call messup_densmat(Dens1(:,:,2))
   call sprepack('L',M,rhobeta,Dens1(:,:,2))

   ! PERFORM SECOND LINEAR RESPONSE
   call open_linear_response(Calpha1,Cbeta1,Ealpha0,Ebeta0,&
                              M,NCOa,NCOb)

   Ndima = NCOa * Nvirta
   Ndimb = NCOb * Nvirtb
   Ndim  = max(Ndima,Ndimb)

   ! OBTAIN CORRECT ORDERING OF STATES
   ! this overwrite Tdip_slr with the corrected transition dipole perturbed
   call wich_state(Xflr,Xslr,Calpha1,Cbeta1,Ndim,M,NCOa,NCOb)

   ! CALCULATED TRANSITION DIPOLE OF INTERESTED STATE
   N_slr = nstates - state_LR
   allocate(Tdip(N_slr,3),Ene1(N_slr))
   do ii=1,N_slr
      ele = state_LR + ii
      Tdip(ii,:) = ( Tdip_slr(ele,:) - Tdip_flr(ele,:) ) / lambda_LR
      !print*, "F", ele, Tdip_flr(ele,1), Tdip_flr(ele,2), Tdip_flr(ele,3)
      !print*, "S", ele, Tdip_slr(ele,1), Tdip_slr(ele,2), Tdip_slr(ele,3)
      Ene1(ii)   = Eflr(ele) - Eflr(state_LR)
   enddo

   ! CALCULATED OSCILLATOR STRENGHT
   allocate(Osc(N_slr)); Osc = 0.0d0
   call ObtainOsc(Tdip,Ene1,Osc,N_slr)
   write (jj_char, '(i4)') state_LR
   print*,""
   write(*,"(2X,A,9X,A,5X,A,13X,A)") "STATES","ENERGY[Ha]","LAMBDA[nm]","F. Osc."
   do ii=1,N_slr
      write (ii_char, '(i4)') state_LR+ii
      write(*,"(1X,A,A,A,5X,F8.4,6X,F12.6,2X,F20.10)") adjustl(jj_char),"-> ",adjustl(ii_char),&
                                                       Ene1(ii),45.56335D0/Ene1(ii),&
                                                       Osc(ii)
   enddo
   print*,""
   deallocate(Calpha1,Cbeta1,Tdip,Osc,Ene1,Dens1)

end subroutine open_second_LinearResponse

subroutine wich_state(X0,X1,Ca1,Cb1,Ndim,M,NCOa,NCOb)
use lrdata, only: nstates

   implicit none

   integer, intent(in)    :: Ndim, M, NCOa, NCOb
   real*8,  intent(in)    :: X0(Ndim,nstates,2), Ca1(M,M), Cb1(M,M)
   real*8,  intent(inout) :: X1(Ndim,nstates,2)

   character(len=4) :: ii_char, jj_char
   integer :: ii, jj, kk
   real*8  :: temp, Ctol, scale_factor
   logical :: asigned
   real*8,  dimension(:), allocatable :: indexes, Ene_fake, OsSt_fake
   real*8,  dimension(:,:), allocatable :: X1temp

   Ctol = 0.001d0
   write(*,"(1X,A,F6.3,1X,A)") "Using Ctol= ",Ctol,"to projection States"

   allocate(X1temp(Ndim,nstates),indexes(nstates))
   ! ALPHA
   do ii=1,nstates ! Perturbed
      asigned = .false.
      do jj=1,nstates ! Unperturbed
         temp = 0.0d0
         do kk=1,Ndim
            temp = temp + X1(kk,ii,1) * X0(kk,jj,1)
         enddo
         if((dabs(temp)>(0.5d0-Ctol)) .and. (dabs(temp)<(0.5d0+Ctol))) then
            asigned = .true.
            scale_factor = 1.0d0
            if(temp < 0.0d0) scale_factor = -1.0d0
            X1temp(:,ii) = scale_factor * X1(:,jj,1)
            if ( ii /= jj ) then
               write (ii_char, '(i4)') ii
               write (jj_char, '(i4)') jj
               write (*,"(1X,A,A,A,A)") "Switch state A ", adjustl(ii_char), "to ", adjustl(jj_char)
            endif
         else
            indexes(jj) = dabs(temp)
         endif
         if (asigned) exit ! exit of jj loop
      enddo
      if(.not. asigned) then
         write (ii_char, '(i4)') ii
         write (jj_char, '(i4)') maxloc(indexes)
         write(*,"(1X,A,A,1X,A)") "WARNING: The State ",adjustl(ii_char),"couldn't asigned"
         write(*,"(1X,A,A,A,A,1X,A,F8.4)") "MAX: |X1_a(",adjustl(ii_char),")*X0_a(",adjustl(jj_char),")^T|=",maxval(indexes)
         print*, ""
         X1temp(:,ii) = X1(:,ii,1)
      endif
   enddo
   X1(:,:,1) = X1temp(:,:); X1temp = 0.0d0

   ! BETA
   do ii=1,nstates ! Perturbed
      asigned = .false.
      do jj=1,nstates ! Unperturbed
         temp = 0.0d0
         do kk=1,Ndim
            temp = temp + X1(kk,ii,2) * X0(kk,jj,2)
         enddo
         if((dabs(temp)>(0.5d0-Ctol)) .and. (dabs(temp)<(0.5d0+Ctol))) then
            asigned = .true.
            scale_factor = 1.0d0
            if(temp < 0.0d0) scale_factor = -1.0d0
            X1temp(:,ii) = scale_factor * X1(:,jj,2)
            if ( ii /= jj ) then
               write (ii_char, '(i4)') ii
               write (jj_char, '(i4)') jj
               write (*,"(1X,A,A,A,A)") "Switch state B ", adjustl(ii_char), "to ", adjustl(jj_char)
            endif
         else
            indexes(jj) = dabs(temp)
         endif
         if (asigned) exit ! exit of jj loop
      enddo
      if(.not. asigned) then
         write (ii_char, '(i4)') ii
         write (jj_char, '(i4)') maxloc(indexes)
         write(*,"(1X,A,A,1X,A)") "WARNING: The State ",adjustl(ii_char),"couldn't asigned"
         write(*,"(1X,A,A,A,A,1X,A,F8.4)") "MAX: |X1_b(",adjustl(ii_char),")*X0_b(",adjustl(jj_char),")^T|=",maxval(indexes)
         print*, ""
         X1temp(:,ii) = X1(:,ii,2)
      endif
   enddo
   X1(:,:,2) = X1temp(:,:); X1temp = 0.0d0
   deallocate(X1temp,indexes)

   ! OBTAIN CORRECT TRANSITION DIPOLE OF PERTURBED VECTORS
   allocate(Ene_fake(nstates),OsSt_fake(nstates))
   Ene_fake = 0.0d0 ! NOT reference
   call open_OscStr(X1(:,:,1),X1(:,:,2),Ene_fake,Ca1,Cb1,OsSt_fake, &
                    M,NCOa,NCOb,Ndim,nstates) ! OVERWRITE DIPOLES
   deallocate(Ene_fake,OsSt_fake)

end subroutine wich_state

