!     Code for Anisotropies in the Microwave Background
!     by Antony Lewis (http://cosmologist.info/) and Anthony Challinor
!     See readme.html for documentation. This is a sample driver routine that reads
!     in one set of parameters and produdes the corresponding output. 

    subroutine camb_wrap(pparams,ell,cltt,clte,clee,clbb)
        use IniFile
        use CAMB
        use LambdaGeneral
        use Lensing
       !krig use RECFAST 
	use Reionization   !krig    
#ifdef NAGF95
        use F90_UNIX
#endif
        implicit none
	!int sfdel(3000)
	!real sfdcl(3000)
      
        Type(CAMBparams) P
        
        character(LEN=Ini_max_string_len) numstr, VectorFileName, &
            InputFile, ScalarFileName, TensorFileName, TotalFileName, LensedFileName
        integer i
        character(LEN=Ini_max_string_len) TransferFileNames(max_transfer_redshifts), &
               MatterPowerFileNames(max_transfer_redshifts), outroot
        real(dl) output_factor, Age
	real(dl), intent(in) :: pparams(6)
	real(dl), intent(out) :: cltt(3000),ell(3000)
	real(dl), intent(out) :: clte(3000),clee(3000),clbb(3000)

#ifdef WRITE_FITS
       character(LEN=Ini_max_string_len) FITSfilename
#endif

#ifndef NAGF95
#ifndef __INTEL_COMPILER_BUILD_DATE
        integer iargc
        !external iargc
#endif        
#endif
        logical bad

        InputFile = ''
   

       ! if (iargc() /= 0)  call getarg(1,InputFile)
       ! if (InputFile == '') stop 'No parameter input file'

       ! call Ini_Open(InputFile, 1, bad, .false.)
       ! if (bad) stop 'Error opening parameter file'

       ! Ini_fail_on_not_found = .false.
    
       ! outroot = Ini_Read_String('output_root')
       ! if (outroot /= '') outroot = trim(outroot) // '_'
        
        call CAMB_SetDefParams(P)

        P%WantScalars = .TRUE.
        P%WantVectors = .FALSE.
        P%WantTensors = .FALSE.
        
        P%OutputNormalization=outNone
        !if (Ini_Read_Logical('COBE_normalize',.false.))  P%OutputNormalization=outCOBE
        output_factor = 7.4311e12

        P%WantCls= P%WantScalars .or. P%WantTensors .or. P%WantVectors

        P%WantTransfer=.TRUE.
        
        P%NonLinear = 0
   
        P%DoLensing = .false.
        if (P%WantCls) then
          if (P%WantScalars  .or. P%WantVectors) then
           P%Max_l = 2000
           P%Max_eta_k = 4000
           if (P%WantScalars) then
             P%DoLensing = .false.
             if (P%DoLensing) lensing_method = 1
           end if
           if (P%WantVectors) then
            if (P%WantScalars .or. P%WantTensors) stop 'Must generate vector modes on their own'
            i = Ini_Read_Int('vector_mode')
            if (i==0) then 
              vec_sig0 = 1
              Magnetic = 0
            else if (i==1) then
              Magnetic = -1
              vec_sig0 = 0
            else
              stop 'vector_mode must be 0 (regular) or 1 (magnetic)'
            end if 
           end if
          end if

          if (P%WantTensors) then
           P%Max_l_tensor = 1500
           P%Max_eta_k_tensor =  3000
          end if
        endif

                
!  Read initial parameters.
       
       w_lam = -1.0
       cs2_lam = 1.0

       P%h0     = pparams(3)
 
       !if (Ini_Read_Logical('use_physical',.false.)) then 

        P%omegab = pparams(1)/(P%H0/100)**2
        P%omegac = pparams(2)/(P%H0/100)**2
        P%omegan = 0.0/(P%H0/100)**2
        !P%omegav = 1- Ini_Read_Double('omk') - P%omegab-P%omegac - P%omegan
	P%omegav = 1- 0.0 - P%omegab-P%omegac - P%omegan
        
	!write(*,*)'ombh2 ',pparams(1),' omch2 ',pparams(2),' h ',P%H0
  
      ! else
       
       ! P%omegab = Ini_Read_Double('omega_baryon')
       ! P%omegac = Ini_Read_Double('omega_cdm')
       ! P%omegav = Ini_Read_Double('omega_lambda')
       ! P%omegan = Ini_Read_Double('omega_neutrino')

       !end if

       P%tcmb   = 2.726
       P%yhe    = 0.24
       P%Num_Nu_massless  = 3.04
       P%Num_Nu_massive   = 0.0
   
       P%nu_mass_splittings = .true.
       P%Nu_mass_eigenstates = 1
       if (P%Nu_mass_eigenstates > max_nu) stop 'too many mass eigenstates'
       numstr = '0'
       if (numstr=='') then
         P%Nu_mass_degeneracies(1)= P%Num_nu_massive
       else
        read(numstr,*) P%Nu_mass_degeneracies(1:P%Nu_mass_eigenstates)
       end if
       numstr = '1'
       if (numstr=='') then
        P%Nu_mass_fractions(1)=1  
        if (P%Nu_mass_eigenstates >1) stop 'must give nu_mass_fractions for the eigenstates'
       else
        read(numstr,*) P%Nu_mass_fractions(1:P%Nu_mass_eigenstates)
       end if

       if (P%NonLinear==NonLinear_lens .and. P%DoLensing) then
          if (P%WantTransfer) &
             write (*,*) 'overriding transfer settings to get non-linear lensing'
          P%WantTransfer  = .true.
          call Transfer_SetForNonlinearLensing(P%Transfer)
          P%Transfer%high_precision=  .false.
       
       else if (P%WantTransfer)  then
        P%Transfer%high_precision=  .false.
        P%transfer%kmax          =  10
        P%transfer%k_per_logint  =  0
        P%transfer%num_redshifts =  1
       ! if (P%transfer%num_redshifts > max_transfer_redshifts) stop 'Too many redshifts'
       ! do i=1, P%transfer%num_redshifts
	
        !     write (numstr,*) i 
        !     numstr=adjustl(numstr)
        !     P%transfer%redshifts(i)  = Ini_Read_Double('transfer_redshift('//trim(numstr)//')',0._dl)
        !     TransferFileNames(i)     = Ini_Read_String('transfer_filename('//trim(numstr)//')')
        !     if (TransferFileNames(i)/= '') &
        !           TransferFileNames(i) = trim(outroot)//TransferFileNames(i)
        !     MatterPowerFilenames(i)  = Ini_Read_String('transfer_matterpower('//trim(numstr)//')')
        !     if (MatterPowerFilenames(i) /= '') &
        !         MatterPowerFilenames(i)=trim(outroot)//MatterPowerFilenames(i)
        
	!end do
        !P%transfer%kmax=P%transfer%kmax*(P%h0/100._dl)
                
       !else
       !  P%transfer%high_precision = .false.
       !endif
  
        !krig P%Reionization = .true.
        !krig P%use_optical_depth = P%Reionization .and. .true.
	P%Reion%use_optical_depth = .true.
  
       ! if ( P%use_optical_depth) then
        
	      P%Reion%optical_depth = pparams(4)
	      P%Reion%delta_redshift=0.5 !krig
        
	 !  else if (P%Reionization) then
         !     P%Reion%redshift = Ini_Read_Double('re_redshift')
         !     P%Reion%fraction = Ini_Read_Double('re_ionization_frac')
       ! end if 

           Ini_fail_on_not_found = .false. 
           
       ! RECFAST_fudge = 1.14

           i = 1
           if (i==2) then
           ! use_Dubrovich = .true.
           else if (i/=1) then
             stop 'Unknown recombination'
           end if

           P%InitPower%nn = 1
           if (P%InitPower%nn>nnmax) stop 'Too many initial power spectra - increase nnmax in InitialPower'
           P%InitPower%rat(:) = 1
           do i=1, P%InitPower%nn
              write (numstr,*) i 
              numstr=adjustl(numstr)
              P%InitPower%an(i) = pparams(5)!&
                   !Ini_Read_Double('scalar_spectral_index('//trim(numstr)//')')

              P%InitPower%n_run(i) = 0!&
                   !Ini_Read_Double('scalar_nrun('//trim(numstr)//')',0._dl)
    
            !  if (P%WantTensors) then
            !     P%InitPower%ant(i) = Ini_Read_Double('tensor_spectral_index('//trim(numstr)//')')
            !     P%InitPower%rat(i) = Ini_Read_Double('initial_ratio('//trim(numstr)//')')
            !  end if              

              P%InitPower%ScalarPowerAmp(i) = pparams(6)
              !Always need this as may want to set tensor amplitude even if scalars not computed
           end do
      
            if (P%WantScalars .or. P%WantTransfer) then
            P%Scalar_initial_condition = 1
            if (P%Scalar_initial_condition == initial_vector) then
                P%InitialConditionVector=0
              numstr = '-1 0 0 0 0'
              read (numstr,*) P%InitialConditionVector(1:initial_iso_neutrino_vel)
            end if
        end if

        
      ! if (P%WantScalars) then
      !    ScalarFileName = trim(outroot)//Ini_Read_String('scalar_output_file')
      !    LensedFileName =  trim(outroot) //Ini_Read_String('lensed_output_file')
      !  end if
      !  if (P%WantTensors) then
      !    TensorFileName =  trim(outroot) //Ini_Read_String('tensor_output_file')
      !   if (P%WantScalars) TotalFileName =  trim(outroot) //Ini_Read_String('total_output_file')
      !  end if
      !  if (P%WantVectors) then
      !    VectorFileName =  trim(outroot) //Ini_Read_String('vector_output_file')
      !  end if
         
#ifdef WRITE_FITS
       ! if (P%WantCls) then
       ! FITSfilename =  trim(outroot) //Ini_Read_String('FITS_filename',.true.)
       ! if (FITSfilename /='') then
       ! inquire(file=FITSfilename, exist=bad)
       ! if (bad) then
       !  open(unit=18,file=FITSfilename,status='old')
       !  close(18,status='delete')
       ! end if
       !end if
       ! end if
#endif        
       

       Ini_fail_on_not_found = .false. 

!optional parameters controlling the computation

       P%AccuratePolarization = .true.
       P%AccurateReionization = .false.
       P%AccurateBB = .false.
        
       !Mess here to fix typo with backwards compatibility
       if (Ini_Read_String('do_late_rad_trunction') /= '') then
         DoLateRadTruncation = .true.
        ! if (Ini_Read_String('do_late_rad_truncation')/='') stop 'check do_late_rad_xxxx'
       else
        DoLateRadTruncation = .true.
       end if
       DoTensorNeutrinos = .false.
       FeedbackLevel = 0
       
       P%MassiveNuMethod  = 3

       ThreadNum      = 0 !Ini_Read_Int('number_of_threads',0)
       AccuracyBoost  = 1 !Ini_Read_Double('accuracy_boost',1.d0)
       lAccuracyBoost = 1 !Ini_Read_Double('l_accuracy_boost',1.d0)
       lSampleBoost   = 1 !Ini_Read_Double('l_sample_boost',1.d0)

       !if (outroot /= '') then
       !  call Ini_SaveReadValues(trim(outroot) //'params.ini',1)
       !end if

       !call Ini_Close

       if (CAMB_ValidateParams(P)) then
       !write(*,*)'it is valid'
#ifdef RUNIDLE
       call SetIdle
#endif 

       if (FeedbackLevel > 0) then
         Age = CAMB_GetAge(P) 
         write (*,'("Age of universe/GYr  = ",f7.3)') Age  
       end if 

	!write(*,*)'about to call get results in inidriver '
       call CAMB_GetResults(P)
    
        if (P%WantTransfer .and. .not. (P%NonLinear==NonLinear_lens .and. P%DoLensing)) then
         !call Transfer_SaveToFiles(MT,TransferFileNames)
         !call Transfer_SaveMatterPower(MT,MatterPowerFileNames)
         if ((P%OutputNormalization /= outCOBE) .or. .not. P%WantCls)  call Transfer_output_sig8(MT)
        end if

        
  
         if (P%OutputNormalization == outCOBE) then

            if (P%WantTransfer) call Transfer_output_Sig8AndNorm(MT)
           
          end if
         
	! open(unit=55,file='sfdoutput.sav',form='formatted',status='unknown')
	! do i=lmin,CP%Max_l
	! write(55,'(1I5,3E15.5)')i ,output_factor*Cl_scalar(i,1,C_Temp:C_Cross)
	! end do
	! close(55)
      
        do i=lmin,CP%Max_l
	  ell(i)=i
	  cltt(i)=output_factor*Cl_scalar(i,1,C_Temp)
	  clee(i)=output_factor*Cl_scalar(i,1,C_E)
	  clte(i)=output_factor*Cl_scalar(i,1,C_Cross)
	  clbb(i)=0.0
	  !write(*,*)'writing ',i,output_factor*Cl_scalar(i,1,C_Temp)
	end do
      
     
      else
       write(*,*)'your parameters make no sense, but I do not care'
       !write(*,*)lmin,CP%Max_l
       do i=1,3000
         !write(*,*)'i ',i
	  ell(i)=i
	  cltt(i)=-100.0
	  clee(i)=-100.0
	  clte(i)=-100.0
	  clbb(i)=-100.0
	  !write(*,*)'writing ',i,output_factor*Cl_scalar(i,1,C_Temp)
	end do
	!write(*,*)'end do'
      end if
     ! write(*,*)'end if'
  end if
  !write(*,*)'leaving camb routine'
        end subroutine camb_wrap


