  module module_mosaic_aerchem_intr


  implicit none


  contains


  !-----------------------------------------------------------------------
  subroutine aerchemistry( it_mosaic, dtchem_in,              & !intent-ins
     pr_atm, rh, te, cair_mol_m3, cair_mol_cc,                &
     jaerosolstate, jaerosolstate_bgn, jhyst_leg,             & !intent-inouts
     rbox, dp_dry_a, dp_wet_a, sigmag_a,                      & 
     gas_avg, gas_netprod_otrproc,                            & 
     mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a,  &
     aH2O_a, gam_ratio, iter_mesa_out                         ) !intent-outs

  use module_data_mosaic_kind, only: r8
  use module_data_mosaic_main, only: &
       m_partmc_mosaic, ntot_max, ntot_used
  use module_data_mosaic_aero, only : &
       mosaic_vars_aa_type, &
       dens_aer_mac, &
       ibc_a, ioc_a, ilim2_a, ioin_a, &
       mw_aer_mac, mw_comp_a, msectional, msize_framework, &
       naer, nbin_a, nbin_a_max, ngas_aerchtot, ngas_volatile, &
       nmax_astem, nmax_mesa, nsalt, &
       use_cam5mam_soa_params, use_cam5mam_accom_coefs
  use module_mosaic_box_aerchem, only: mosaic_box_aerchemistry


  !Subroutine arguments
  integer,  intent(in)  :: it_mosaic
  real(r8), intent(in)  :: dtchem_in, pr_atm, rh, te, cair_mol_m3, cair_mol_cc

  integer, intent(inout),  dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
  integer, intent(inout),  dimension(nbin_a_max) :: jhyst_leg

  real(r8), intent(inout), dimension(ntot_used)     :: rbox
  real(r8), intent(inout), dimension(nbin_a_max)    :: dp_dry_a, dp_wet_a, sigmag_a
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas_avg  ! average gas conc. over dtchem time step (nmol/m3)
  real(r8), intent(inout), dimension(ngas_aerchtot) :: gas_netprod_otrproc
            ! gas_netprod_otrproc = gas net production rate from other processes
            !    such as gas-phase chemistry and emissions (nmol/m3/s)
            ! NOTE - currently in the mosaic box model, gas_netprod_otrproc is set to zero for all
            !        species, so mosaic_aerchemistry does not apply production and condensation together
  real(r8), intent(inout), dimension(nbin_a_max)    :: mass_dry_a_bgn, mass_dry_a
  real(r8), intent(inout), dimension(nbin_a_max)    :: dens_dry_a_bgn, dens_dry_a

  integer,  intent(out),   dimension(nbin_a_max)    :: iter_mesa_out
  real(r8), intent(out),   dimension(nbin_a_max)    :: aH2O_a, gam_ratio

  !Local variables
  character(len=250) :: infile, tmp_str

  logical :: debug_mosaic = .false.

  integer :: ierr, ibin, igas, iaer, istate, iaer_in, istate_in, ibin_in
  integer :: mcall_load_mosaic_parameters, mcall_print_aer_in
  integer :: n
  integer :: unitn

  real(r8) :: dtchem, RH_pc, aH2O, P_atm, T_K, aer_tmp
  real(r8), dimension(naer,3,nbin_a_max) :: aer
  real(r8), dimension(ngas_aerchtot)     :: gas
  real(r8), dimension(nbin_a_max)        :: num_a, water_a, water_a_hyst
  real(r8), dimension(naer)              :: kappa_nonelectro
  real(r8)                               :: uptkrate_h2so4  ! rate of h2so4 uptake by aerosols (1/s)

  type (mosaic_vars_aa_type) :: mosaic_vars_aa


  dtchem      = dtchem_in
  RH_pc       = RH                                    ! RH(%)
  aH2O        = 0.01_r8*RH_pc                             ! aH2O (aerosol water activity)
  P_atm       = pr_atm                                ! P(atm)
  T_K         = te                                      ! T(K)

  kappa_nonelectro(:) = 0.0_r8
  kappa_nonelectro(ibc_a  ) = 0.0001  ! previously kappa_poa = 0.0001
  kappa_nonelectro(ioc_a  ) = 0.0001  ! previously kappa_bc  = 0.0001
  kappa_nonelectro(ilim2_a) = 0.1     ! previously kappa_soa = 0.1
  kappa_nonelectro(ioin_a ) = 0.06    ! previously kappa_oin = 0.06


  ! for box model
  !    on first time step, call load_mosaic_parameters and call print_aer twice
  !    after first time step, just call print_aer once
  ! for cam5 or wrfchem, this routine will loop over multiple grid boxes
  !    call load_mosaic_parameters for first grid box (and all time steps)
  !    call print_aer never
  if (it_mosaic <= 1) then
     mcall_load_mosaic_parameters = 1
     mcall_print_aer_in = 2
  else
     mcall_load_mosaic_parameters = 0
     mcall_print_aer_in = 1
  end if


  ! map variables from rbox (and other) arrays to mosaic aerchem working arrays
  call map_mosaic_species_aerchem_box( 0, jaerosolstate,  &          
       rbox, aer, gas, jhyst_leg, num_a, Dp_dry_a,        &
       sigmag_a, water_a, water_a_hyst, cair_mol_m3       )

  
  !BSINGH - Following block is introduced to reproduce errors Mosaic
  !         model encounters in other models (CAM,WRF etc.). This block repopulate
  !         all the information which is going into the mosaic box (intent-ins and
  !         intent-inouts). It is a binary read to preserve the accuracy.
  if(debug_mosaic) then

     ! set these control variables
     use_cam5mam_soa_params  = 1
     use_cam5mam_accom_coefs = 1

     !Read a binary file which has all the inputs to the mosaic box
     !and stop the model
     
     unitn = 101
     infile = 'mosaic_error_48.bin'
     open( unitn, file=trim(infile), status='old', form='unformatted', CONVERT = 'BIG_ENDIAN' )
     
     read(unitn)aH2O
     read(unitn)T_K
     read(unitn)P_atm
     read(unitn)RH_pc
     read(unitn)dtchem
     

     do ibin = 1, nbin_a_max
        read(unitn)num_a(ibin),water_a(ibin),Dp_dry_a(ibin),        &
             sigmag_a(ibin),dp_wet_a(ibin),jhyst_leg(ibin),          &
             jaerosolstate(ibin)
     end do
     
     
     do igas = 1, ngas_aerchtot
        if (igas <= ngas_volatile) then
           read(unitn) gas(igas), gas_avg(igas), gas_netprod_otrproc(igas)
        else
           gas(igas) = 0.0 ; gas_avg(igas) = 0.0 ; gas_netprod_otrproc(igas) = 0.0
        end if
     enddo
     
     do ibin = 1, nbin_a_max
        do istate = 1, 3
           do iaer = 1 , naer
              read(unitn)iaer_in,istate_in,ibin_in, aer_tmp
              aer(iaer_in,istate_in,ibin_in) = aer_tmp                    
           end do
        end do
     end do
     close(unitn)
  endif
  !BSINGH -----xxx ENDS reading file for debugging mosaic xxxx----
  


  ! calculate gas-aerosol exchange over timestep dtchem
  ! (during this calculation there is no transfer of particles between bins)

 ! aH2O = 0.999 ! min(0.99, aH2O)	! RAZ 2/14/2014
 ! RH_pc = 99.9


  allocate( mosaic_vars_aa%iter_mesa(nbin_a_max), stat=ierr )
  if (ierr /= 0) then
     print *, '*** subr aerchemistry - allocate error for mosaic_vars_aa%iter_mesa'
     stop
  end if
  mosaic_vars_aa%it_host = 0
  mosaic_vars_aa%it_mosaic = it_mosaic
  mosaic_vars_aa%hostgridinfo(:) = 0
  mosaic_vars_aa%f_mos_fail = -1
  mosaic_vars_aa%isteps_astem = 0
  mosaic_vars_aa%isteps_astem_max = 0
  mosaic_vars_aa%jastem_call = 0
  mosaic_vars_aa%jastem_fail = -1
  mosaic_vars_aa%jmesa_call = 0
  mosaic_vars_aa%jmesa_fail = 0
  mosaic_vars_aa%niter_mesa_max = 0
  mosaic_vars_aa%nmax_astem = nmax_astem
  mosaic_vars_aa%nmax_mesa = nmax_mesa
  mosaic_vars_aa%fix_astem_negative = 0
  !BSINGH - flag_itr_kel becomes true when kelvin iteration in mdofule_mosaic_ext.F90 are greater then 100
  mosaic_vars_aa%flag_itr_kel = .false.
  !BSINGH - zero_water_flag becomes .true. if water is zero in liquid phase
  mosaic_vars_aa%zero_water_flag = .false.
  mosaic_vars_aa%cumul_steps_astem = 0.0_r8
  mosaic_vars_aa%niter_mesa = 0.0_r8
  mosaic_vars_aa%xnerr_astem_negative(:,:) = 0.0_r8
  mosaic_vars_aa%iter_mesa(1:nbin_a_max) = 0


  call mosaic_box_aerchemistry(              aH2O,               T_K,            &!Intent-ins
       P_atm,                  RH_pc,        dtchem,                             &
       mcall_load_mosaic_parameters,         mcall_print_aer_in, sigmag_a,       &
       kappa_nonelectro,                                                         &
       jaerosolstate,          aer,                                              &!Intent-inouts
       num_a,                  water_a,      gas,                                &
       gas_avg,                gas_netprod_otrproc,              Dp_dry_a,       &
       dp_wet_a,               jhyst_leg,                                        &
       mosaic_vars_aa,                                                           &
       mass_dry_a_bgn,         mass_dry_a,                                       &!Intent-outs
       dens_dry_a_bgn,         dens_dry_a,   water_a_hyst,       aH2O_a,         &
       uptkrate_h2so4,         gam_ratio,    jaerosolstate_bgn                   )


! if (jASTEM_fail > 0 .or. zero_water_flag .or. f_mos_fail > 0 ) then
  if (mosaic_vars_aa%jastem_fail > 0 .or. mosaic_vars_aa%zero_water_flag .or. mosaic_vars_aa%f_mos_fail > 0) then
     !Write error message and stop the model.
     write(tmp_str,*) 'Error in Mosaic,jASTEM_fail= ', mosaic_vars_aa%jASTEM_fail, &
        ' zero_water_flag: ', mosaic_vars_aa%zero_water_flag, &
        '  f_mos_fail:', mosaic_vars_aa%f_mos_fail
     print*, trim(adjustl(tmp_str))
     print*, 'Fortran Stop in aerosol/aerchemistry.f90'
     stop
  endif

  iter_mesa_out(1:nbin_a_max) = mosaic_vars_aa%iter_mesa(1:nbin_a_max)
 
  if ( sum(mosaic_vars_aa%xnerr_astem_negative(:,:)) > 0.0_r8 ) then
     print '(a,i10)', 'it_mosaic, xnerr_astem_negative =', it_mosaic
     do n = 1, 4
        print '(1p,5e10.2)', mosaic_vars_aa%xnerr_astem_negative(:,n)
     end do
     print '(a)', 'Fortran Stop in aerosol/aerchemistry.f90'
     stop
  end if

  deallocate( mosaic_vars_aa%iter_mesa, stat=ierr )
  if (ierr /= 0) then
     print *, '*** subr aerchemistry - deallocate error for mosaic_vars_aa%iter_mesa'
     stop
  end if


  ! map variables to rbox (and other) arrays from mosaic aerchem working arrays
  call map_mosaic_species_aerchem_box( 1, jaerosolstate, &          
       rbox, aer, gas, jhyst_leg, num_a, Dp_dry_a,       &
       sigmag_a, water_a, water_a_hyst, cair_mol_m3      )


  return
  end subroutine aerchemistry




  !***********************************************************************
  ! maps gas and aerosol information between
  !    rbox, jhyst_leg, ... AND
  !    aerchemistry working arrays (gas, aer, num_a, water_a, ...)
  !
  ! author: Rahul A. Zaveri
  ! update: nov 2001
  !-----------------------------------------------------------------------
  subroutine map_mosaic_species_aerchem_box( imap, jaerosolstate,  &
       rbox, aer, gas, jhyst_leg, num_a, Dp_dry_a,                 &
       sigmag_a, water_a, water_a_hyst, cair_mol_m3                )

    use module_data_mosaic_main, only: &
         r8, naer_tot, ngas_max, m_partmc_mosaic,       & !Parameters
         ntot_used, avogad, mw_air, piover6

    use module_data_mosaic_boxmod, only: &
         kh2so4, khno3, khcl, knh3, kmsa, karo1, karo2, kalk1, kole1, kapi1,      & !TBD
         kapi2, klim1, klim2, knum_a, kdpdry_a, ksigmag_a, kjhyst_a, kwater_a       !TBD

    use module_data_mosaic_aero, only: nbin_a_max,naer,ngas_aerchtot,jtotal,       & !Parameters
         nbin_a,                                                                   & !Input
         ih2so4_g,ihno3_g,ihcl_g,inh3_g,imsa_g,iaro1_g,iaro2_g,ialk1_g,iole1_g,    &
         iapi1_g,iapi2_g,ilim1_g,ilim2_g,                                          & !TBD
         jhyst_lo, jhyst_up, jhyst_undefined,                                      &
         mhyst_method, mhyst_uporlo_waterhyst,                                     &
         mw_aer_mac,                                                               &
         all_solid, all_liquid, mixed, no_aerosol


    ! subr arguments
    integer, intent(in) :: imap
    integer, intent(inout), dimension(nbin_a_max) :: jaerosolstate, jhyst_leg

    real(r8), intent(in) :: cair_mol_m3
    real(r8), intent(inout), dimension(ntot_used) :: rbox
    real(r8), intent(inout), dimension(nbin_a_max) :: num_a, Dp_dry_a, sigmag_a, water_a, water_a_hyst
    real(r8), intent(inout), dimension(ngas_aerchtot) :: gas
    real(r8), intent(inout), dimension(naer,3,nbin_a_max) :: aer

    ! local variables
    integer ibin, iaer, noffset
    real(r8) :: conv_aer, conv_aerinv
    real(r8) :: conv_gas, conv_gasinv
    real(r8) :: conv_num, conv_numinv
    real(r8) :: conv_wat, conv_watinv
    real(r8) :: tmpa



    if ((imap < 0) .or. (imap > 1)) then
       write(*,*) &
       '*** map_mosaic_species_BOX fatal error - bad imap =', imap
       stop
    end if


    ! define conversion factors
    ! BOX
    ! gases -- rbox = umol/mol,   gas = rbox*conv_gas = nmol/m^3
    conv_gas = 1.e3_r8*cair_mol_m3
    conv_gasinv = 1.0_r8/conv_gas

    ! aerosol mass -- rbox = ug/kg,   aer = rbox*conv_aer/mw_aer = nmol/m^3
    conv_aer = mw_air*cair_mol_m3
    conv_aerinv = 1.0_r8/conv_aer

    ! aerosol water -- rbox = ug/kg,   water_a = rbox*conv_wat = kg/m^3
    conv_wat = 1.e-12_r8*mw_air*cair_mol_m3
    conv_watinv = 1.0_r8/conv_wat

    ! aerosol number -- rbox = #/kg,   num_a = rbox*conv_num = #/cm^3
    conv_num = 1.e-9_r8*mw_air*cair_mol_m3
    conv_numinv = 1.0_r8/conv_num


    if (imap == 0) then    
       ! map from host code arrays (rbox in this case) to 
       ! mosaic aerchem working arrays (gas, aer, num_a, etc)
       gas(:) = 0.0_r8
       aer(:,:,:) = 0.0_r8

       ! gases -- rbox = mol/mol,   gas = nmol/m^3
       gas(ih2so4_g) = rbox(kh2so4)*conv_gas
       gas(ihno3_g)  = rbox(khno3)*conv_gas
       gas(ihcl_g)   = rbox(khcl)*conv_gas
       gas(inh3_g)   = rbox(knh3)*conv_gas
       gas(imsa_g)   = rbox(kmsa)*conv_gas
       gas(iaro1_g)  = rbox(karo1)*conv_gas
       gas(iaro2_g)  = rbox(karo2)*conv_gas
       gas(ialk1_g)  = rbox(kalk1)*conv_gas
       gas(iole1_g)  = rbox(kole1)*conv_gas
       gas(iapi1_g)  = rbox(kapi1)*conv_gas
       gas(iapi2_g)  = rbox(kapi2)*conv_gas
       gas(ilim1_g)  = rbox(klim1)*conv_gas
       gas(ilim2_g)  = rbox(klim2)*conv_gas
       
       !print*,'BALLI:in-map:', gas(inh3_g),cnn(knh3),conv1,knh3
       ! aerosol
       do ibin = 1, nbin_a

          noffset = ngas_max + naer_tot*(ibin - 1)
          num_a(ibin)      = rbox(noffset + knum_a)*conv_num    ! aerosol number -- rbox = #/kg,   num_a = #/cm^3
          water_a(ibin)    = rbox(noffset + kwater_a)*conv_wat  ! aerosol water -- rbox = ug/kg,   water_a = kg/m^3

          if (mhyst_method == mhyst_uporlo_waterhyst) then
             ! in this case, rbox holds water_a_hyst
             water_a_hyst(ibin) = rbox(noffset + kjhyst_a)*conv_wat ! rbox = ug/kg,   water_a_hyst = kg/m^3
             ! value of jhyst_leg should not matter, so set it to undefined
             jhyst_leg(ibin) = jhyst_undefined
          else
             ! in this case, use the incoming jhyst_leg value (unchanged)
             ! also, input value of water_a_hyst should not be important, so set it to zero
             water_a_hyst(ibin) = 0.0_r8
          end if

          do iaer = 1, naer
             ! aerosol mass components -- rbox = ug/kg,   aer = nmol/m^3
             !    (for oin, bc, oc, molecular weight = 1.0 so moles = grams)
             aer(iaer,jtotal,ibin) = rbox(noffset+kwater_a+iaer)*conv_aer/mw_aer_mac(iaer)
          enddo

       enddo

    else if (imap == 1) then
       ! map from mosaic aerchem working arrays (gas, aer, num_a, etc)
       ! back to host code arrays (rbox in this case)

       rbox(kh2so4)  = gas(ih2so4_g)*conv_gasinv
       rbox(khno3)   = gas(ihno3_g)*conv_gasinv
       rbox(khcl)    = gas(ihcl_g)*conv_gasinv
       rbox(knh3)    = gas(inh3_g)*conv_gasinv
       rbox(kmsa)    = gas(imsa_g)*conv_gasinv
       rbox(karo1)   = gas(iaro1_g)*conv_gasinv
       rbox(karo2)   = gas(iaro2_g)*conv_gasinv
       rbox(kalk1)   = gas(ialk1_g)*conv_gasinv
       rbox(kole1)   = gas(iole1_g)*conv_gasinv
       rbox(kapi1)   = gas(iapi1_g)*conv_gasinv
       rbox(kapi2)   = gas(iapi2_g)*conv_gasinv
       rbox(klim1)   = gas(ilim1_g)*conv_gasinv
       rbox(klim2)   = gas(ilim2_g)*conv_gasinv

       ! aerosol
       do ibin = 1, nbin_a

          noffset = ngas_max + naer_tot*(ibin - 1)
          rbox(noffset + knum_a)    = num_a(ibin)*conv_numinv
          rbox(noffset + kwater_a)  = water_a(ibin)*conv_watinv

          if (mhyst_method == mhyst_uporlo_waterhyst) then
             ! in this case, rbox holds water_a_hyst
             if ( jaerosolstate(ibin) == all_solid  .or. &
                  jaerosolstate(ibin) == all_liquid .or. &
                  jaerosolstate(ibin) == mixed      ) then
                rbox(noffset + kjhyst_a) = water_a_hyst(ibin)*conv_watinv
             else
                rbox(noffset + kjhyst_a) = 0.0_r8
             end if
             ! value of jhyst_leg should not matter, so leave it unchanged
          else
             ! when mhyst_method /= mhyst_uporlo_waterhyst, do nothing,
             ! leave both jhyst_leg and rbox unchanged
             if ( jaerosolstate(ibin) == all_solid  .or. &
                  jaerosolstate(ibin) == all_liquid .or. &
                  jaerosolstate(ibin) == mixed      ) then
                continue
!               jhyst_leg(ibin) = jhyst_leg(ibin)
             else
                jhyst_leg(ibin) = jhyst_undefined
             end if
          end if

          do iaer = 1, naer
             rbox(noffset+kwater_a+iaer) = aer(iaer,jtotal,ibin)*conv_aerinv*mw_aer_mac(iaer)
          enddo

       enddo

    endif

    return
  end subroutine map_mosaic_species_aerchem_box


  end module module_mosaic_aerchem_intr
