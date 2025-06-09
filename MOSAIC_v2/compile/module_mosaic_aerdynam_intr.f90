  module module_mosaic_aerdynam_intr


  implicit none


  contains


  !-----------------------------------------------------------------------
  subroutine aerosoldynamics( it, t_out, t_in,    & !intent-ins
     pr_atm, rh, te, cair_mol_m3,                 &
     cnn, jaerosolstate, dp_wet_a,                & !intent-inouts
     aH2O_a, gam_ratio, iter_MESA                 ) !intent-outs

  use module_data_mosaic_kind, only: r8
  use module_data_mosaic_main, only: &
       m_partmc_mosaic, ntot_max, ntot_used
  use module_data_mosaic_aero, only : &
       dens_aer_mac, mw_aer_mac, mw_comp_a, &
       msectional, msize_framework, &
       naer, nbin_a, nbin_a_max, ngas_aerchtot, nsalt

  use module_mosaic_aerchem_intr, only: aerchemistry
  use module_mosaic_sect_intr, only: sectional_interface_1


  !Subroutine arguments
  integer,  intent(in)  :: it
  integer,  intent(inout), dimension(nbin_a_max) :: jaerosolstate
  integer,  intent(out),   dimension(nbin_a_max) :: iter_MESA

  real(r8), intent(in) :: t_out, t_in
  real(r8), intent(in) :: pr_atm, rh, te, cair_mol_m3

  real(r8), intent(inout), dimension(ntot_max)   :: cnn
  real(r8), intent(inout), dimension(nbin_a_max) :: dp_wet_a
  real(r8), intent(out),   dimension(nbin_a_max) :: aH2O_a, gam_ratio

  !Local variables
  integer :: it_mosaic
  integer, dimension(nbin_a_max) :: jaerosolstate_bgn
  integer, dimension(nbin_a_max) :: jhyst_leg, jhyst_leg_sv1

  real(r8) :: cair_mol_cc
  real(r8) :: dtchem
  real(r8) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr
  real(r8), dimension(ngas_aerchtot)     :: gas_avg  ! average gas conc. over dtchem time step (nmol/m3)
  real(r8), dimension(ngas_aerchtot)     :: gas_netprod_otrproc
            ! gas_netprod_otrproc = gas net production rate from other processes
            !    such as gas-phase chemistry and emissions (nmol/m3/s)
            ! NOTE - currently in the mosaic box model, gas_netprod_otrproc is set to zero for all
            !        species, so mosaic_aerchemistry does not apply production and condensation together
  real(r8), dimension(nbin_a_max)        :: dens_dry_a_bgn, dens_dry_a
  real(r8), dimension(nbin_a_max)        :: mass_dry_a_bgn, mass_dry_a
  real(r8), dimension(nbin_a_max)        :: dp_dry_a, dp_dry_a_sv1
  real(r8), dimension(nbin_a_max)        :: sigmag_a, sigmag_a_sv1

  real(r8), allocatable, dimension(:)    :: rbox, rbox_sv1

  dtchem = t_out - t_in

  if ( (m_partmc_mosaic <= 0) .and. &
       (msize_framework == msectional) ) then
     allocate( rbox_sv1(1:ntot_used) )
  end if
  allocate( rbox(1:ntot_used) )


  it_mosaic   = it
  cair_mol_cc = cair_mol_m3*1.e-6_r8            ! air conc in mol/cc


  ! map variables from cnn array to rbox (and other) arrays
  call aerodynam_map_mosaic_species( 0, cair_mol_cc, cair_mol_m3, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       cnn, rbox, jaerosolstate, jhyst_leg, dp_dry_a, sigmag_a )


  if ( (m_partmc_mosaic <= 0) .and. &
       (msize_framework == msectional) ) then
     rbox_sv1(1:ntot_used) = rbox(1:ntot_used)
     dp_dry_a_sv1(1:nbin_a) = dp_dry_a(1:nbin_a)
     sigmag_a_sv1(1:nbin_a) = sigmag_a(1:nbin_a)
     jhyst_leg_sv1(1:nbin_a) = jhyst_leg(1:nbin_a)
  end if


  gas_avg = 0.0_r8
  gas_netprod_otrproc = 0.0_r8
  mass_dry_a_bgn = 0.0_r8
  mass_dry_a = 0.0_r8
  dens_dry_a_bgn = 0.0_r8
  dens_dry_a = 0.0_r8


  call aerchemistry( it_mosaic, dtchem,                       & !intent-ins
     pr_atm, rh, te, cair_mol_m3, cair_mol_cc,                &
     jaerosolstate, jaerosolstate_bgn, jhyst_leg,             & !intent-inouts
     rbox, dp_dry_a, dp_wet_a, sigmag_a,                      & 
     gas_avg, gas_netprod_otrproc,                            & 
     mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a,  &
     aH2O_a, gam_ratio, iter_MESA                             ) !intent-outs


  ! for sectional framework, calculate
  !    transfer of particles betweens due to growth/shrinkage
  !    new particle nucleation (optional)
  !    particle coagulation (optional)
  if ( (m_partmc_mosaic <= 0) .and. &
       (msize_framework == msectional) ) then

      call sectional_interface_1( dtchem, &
          rbox_sv1, rbox, &
          it_mosaic, jaerosolstate_bgn, jaerosolstate, &
          jhyst_leg, jhyst_leg_sv1, &
          dp_dry_a, dp_dry_a_sv1, &
          sigmag_a, sigmag_a_sv1, &
          mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
          fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
          pr_atm, rh, te, cair_mol_m3, cair_mol_cc )

      deallocate( rbox_sv1 )

  !else if (msize_framework == mmodal) then
  !   do similar calculations for modal framework, but this not yet implemented
  end if


  ! map variables back to cnn array from rbox (and other) arrays
  call aerodynam_map_mosaic_species( 1, cair_mol_cc, cair_mol_m3, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       cnn, rbox, jaerosolstate, jhyst_leg, dp_dry_a, sigmag_a )


  deallocate( rbox )


  return
  end subroutine aerosoldynamics




  !-----------------------------------------------------------------------
  subroutine aerodynam_map_mosaic_species( imap, cair_mol_cc, cair_mol_m3, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       cnn, rbox, jaerosolstate, jhyst_leg, dp_dry_a, sigmag_a )
    !
    !   maps between cnn array and rbox array
    !
    use module_data_mosaic_kind
    use module_data_mosaic_main, only: &
       avogad, &
       mw_air, naer_tot, ngas_max, ntot_max, ntot_used
    use module_data_mosaic_boxmod, only: &
       kdpdry_a, kjhyst_a, knum_a, ksigmag_a, kwater_a
    use module_data_mosaic_aero, only: &
       jhyst_lo, jhyst_up, jhyst_undefined,                                      &
       mhyst_method, mhyst_force_lo, mhyst_force_up, &
       mhyst_uporlo_jhyst, mhyst_uporlo_waterhyst, &
       mmovesect_flag1, mw_aer_mac, &
       naer, naercomp, nbin_a, nbin_a_max, &
       all_solid, all_liquid, mixed, no_aerosol
    use module_data_mosaic_asecthp, only: &
       ai_phase, &
       dens_aer, dens_water_aer, &
       hygro_aer, hyswptr_aer, &
       isize_of_ibin, itype_of_ibin, &
       massptr_aer, ncomp_aer 


    ! subr parameters
    integer, intent(in) :: imap
    integer, intent(in),    dimension(nbin_a_max) :: jaerosolstate
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(in)    :: cair_mol_cc, cair_mol_m3
    real(r8), intent(inout) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr
    real(r8), intent(inout) :: cnn(ntot_max)
    real(r8), intent(inout) :: rbox(ntot_used)
    real(r8), intent(inout), dimension(nbin_a_max) :: dp_dry_a, sigmag_a


    ! local variables
    integer :: ibin, iphase, isize, itype, jhyst_tmp
    integer :: l, ll, lunaa, mtmp, noffset
    real(r8) :: conv_diam, conv_drym, conv_gas, conv_numb, conv_watr
    real(r8) :: tmpa, tmph, tmpj, tmpr


    if ((imap < 0) .or. (imap > 1)) then
       write(*,*) &
       '*** aerchem_map_mosaic_cnn_rbox fatal error - bad imap =', imap
       stop
    end if


    ! set rbox units factors
    !
    ! note that the movesect and coag routines expect the aerosol mass units
    !    for rbox to be mass mixing ratios or mass concentrations
    ! they do not apply any molecular weight factors, so they cannot
    !    work with molar mixing ratios or molar concentrations
    !
    ! for trace gases,   rbox = [umol/mol-air = ppmv],     rbox*fact_gasmr = [molecules/cm^3-air]
    fact_gasmr = 1.0e-6_r8
    ! for aerosol mass,  rbox = [ug-AP/kg-air], rbox*fact_apmassmr = [g-AP/g-air]
    fact_apmassmr = 1.0e-9_r8
    ! for aerosol numb,  rbox = [#/kg-air],     rbox*fact_apnumbmr = [#/g-air]
    fact_apnumbmr = 1.0e-3_r8
    ! aerosol densities in mosaic aerchem and sectional routines = [g/cm^3]
    !     (mosaic densities)*fact_apdens = [g/cm^3]
    fact_apdens = 1.0_r8
    ! aerosol diameters in mosaic aerchem and sectional routines = [cm]
    !     (mosaic diameters)*fact_apdiam = [cm]
    fact_apdiam = 1.0_r8


    ! factors for converting cnn values to rbox (and dp_dry_a_ values
    tmpa = cair_mol_m3*mw_air   ! air_density in ]g/m^3]
    ! gases:             
    ! cnn is [molecules/cm^3-air],  rbox = cnn*conv_gas is [umol/mol-air]
    conv_gas  = 1.0_r8/(avogad*cair_mol_cc*fact_gasmr)
    ! aerosol number:    
    ! cnn is [#/cm^3-air],          rbox = cnn*conv_numb is [#/kg-air]  when iunits_flagaa=4
    conv_numb = 1.0e6_r8 /(fact_apnumbmr*tmpa)
    ! aerosol dry mass:  
    ! cnn is [umol/m^3-air],        rbox = cnn*conv_drym*mw_aer is [ug/kg-air]  when iunits_flagaa=4
    conv_drym = 1.0e-6_r8/(fact_apmassmr*tmpa)
    ! aerosol water:     
    ! cnn is [kg/m^3-air],          rbox = cnn*conv_watr is [ug/kg-air]  when iunits_flagaa=4
    conv_watr = 1.0e3_r8 /(fact_apmassmr*tmpa)
    ! aerosol diameter:  
    ! cnn is [um],                  dp_dry_a = cnn*conv_apdiam is [cm]
    conv_diam = 1.0e-4_r8/fact_apdiam


    if (imap == 1) goto 20000


    !
    ! imap = 0 --> map from cnn to rbox
    !
    rbox(1:ntot_used) = 0.0_r8

    do l = 1, ngas_max
       rbox(l) = cnn(l)*conv_gas
    end do

    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)
       noffset = ngas_max + (ibin-1)*naer_tot

       ! rbox may or may not be able to hold dp_dry_a and sigma_g
       ! (yes in box model, no in wrfchem)
       ! so do not used rbox for holding dp_dry_a and sigma_g
       dp_dry_a(ibin)          = cnn(kdpdry_a +noffset)*conv_diam
       sigmag_a(ibin)          = cnn(ksigmag_a+noffset)

       rbox(knum_a   +noffset) = cnn(knum_a   +noffset)*conv_numb
       rbox(kwater_a +noffset) = cnn(kwater_a +noffset)*conv_watr
       do ll = 1, naer
          l = noffset + (naer_tot-naer) + ll
          rbox(l) = cnn(l)*(conv_drym*mw_aer_mac(ll))
       end do

       if (mhyst_method == mhyst_uporlo_waterhyst) then
          ! cnn and rbox hold water_a_hyst (in water conc units)
          ! SO  convert cnn value to rbox value
          ! AND set jhyst_leg to undefined
          rbox(kjhyst_a +noffset) = cnn(kjhyst_a+noffset)*conv_watr
          jhyst_leg(ibin) = jhyst_undefined
       else
          ! cnn holds jhyst_leg, which can have three possible values
          ! rbox will eventually hold a surrogate for water_a_hyst
          ! SO  set jhyst_leg directly from cnn
          ! AND set rbox value to zero
          rbox(kjhyst_a +noffset) = 0.0_r8
          tmpa = cnn(noffset + kjhyst_a)
          if (tmpa < -0.5_r8) then
             jhyst_leg(ibin) = jhyst_undefined
          else if (tmpa < 0.5_r8) then
             jhyst_leg(ibin) = jhyst_lo
          else
             jhyst_leg(ibin) = jhyst_up
          end if
       end if

    end do

    return


20000 continue
    !
    ! imap = 1 --> map from rbox to cnn
    !
    do l = 1, ngas_max
       cnn(l) = rbox(l)/conv_gas
    end do

    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)
       noffset = ngas_max + (ibin-1)*naer_tot

       cnn(kdpdry_a +noffset) = dp_dry_a(ibin)/conv_diam
       cnn(ksigmag_a+noffset) = sigmag_a(ibin)

       cnn(knum_a   +noffset) = rbox(knum_a   +noffset)/conv_numb
       cnn(kwater_a +noffset) = rbox(kwater_a +noffset)/conv_watr
       do ll = 1, naer
          l = noffset + (naer_tot-naer) + ll
          cnn(l) = rbox(l)/(conv_drym*mw_aer_mac(ll))
       end do

       if (mhyst_method == mhyst_uporlo_waterhyst) then
          ! cnn and rbox hold water_a_hyst (in water conc units)
          ! SO  convert rbox value to cnn value
          ! AND leave jhyst_leg to unchanged
          cnn(kjhyst_a+noffset) = rbox(kjhyst_a+noffset)/conv_watr
       else
          ! cnn holds jhyst_leg, which can have three possible values
          ! rbox will eventually hold a surrogate for water_a_hyst
          ! SO  set cnn directly from jhyst_leg
          ! AND leave rbox value unchanged
          cnn(kjhyst_a+noffset) = jhyst_leg(ibin)*1.0_r8
       end if

    end do

    return
  end subroutine aerodynam_map_mosaic_species


  end module module_mosaic_aerdynam_intr
