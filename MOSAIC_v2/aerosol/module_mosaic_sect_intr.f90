module module_mosaic_sect_intr

!
! *** NOTE ***
!
! this routine currently will not work with openmp parallelization
!
! some of the variables taken from module_data_mosaic_main (cnn, pr_atm, rh, te)
! needs to be passed as subroutine arguments
!

  use module_data_mosaic_kind, only:  r8


  implicit none


  integer, parameter :: iunits_flagaa = 4
  integer, parameter :: iunits_flagbb = 1
  integer, parameter :: iunits_diagout_flagaa = 1

  !BSINGH - 05/28/2013(RCE updates)
  real(r8), save, allocatable ::   &
	  dens_aer_tmp(:,:),  &
	  volumcut_sect_tmp(:,:),  &
	  volumcen_sect_tmp(:,:),  &
	  volumlo_sect_tmp(:,:),  &
	  volumhi_sect_tmp(:,:),  &
	  dcut_sect_tmp(:,:),  &
	  dcen_sect_tmp(:,:),  &
	  dlo_sect_tmp(:,:),  &
	  dhi_sect_tmp(:,:)
  !BSINGH - 05/28/2013(RCE updates ENDS)

  contains


  !BSINGH - 05/28/2013(RCE updates)
!-----------------------------------------------------------------------
	subroutine sect_intr_allocate_memory

	use module_data_mosaic_asecthp, only:  maxd_acomp, maxd_asize, maxd_atype

	allocate( dens_aer_tmp( maxd_acomp, maxd_atype ) )
	allocate( volumcut_sect_tmp( 0:maxd_asize, maxd_atype ) )
	allocate( volumcen_sect_tmp(   maxd_asize, maxd_atype ) )
	allocate( volumlo_sect_tmp(    maxd_asize, maxd_atype ) )
	allocate( volumhi_sect_tmp(    maxd_asize, maxd_atype ) )
	allocate( dcut_sect_tmp( 0:maxd_asize, maxd_atype ) )
	allocate( dcen_sect_tmp(   maxd_asize, maxd_atype ) )
	allocate( dlo_sect_tmp(    maxd_asize, maxd_atype ) )
	allocate( dhi_sect_tmp(    maxd_asize, maxd_atype ) )

	return
	end subroutine sect_intr_allocate_memory
        !BSINGH - 05/28/2013(RCE updates ENDS)

  !-----------------------------------------------------------------------
  subroutine sectional_interface_1( dtchem, &
       rbox_sv1, rbox, &
       it_mosaic, jaerosolstate_bgn, jaerosolstate, &
       jhyst_leg, jhyst_leg_sv1, &
       dp_dry_a, dp_dry_a_sv1, &
       sigmag_a, sigmag_a_sv1, &
       mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
       fact_apmassmr, fact_apnumbmr, fact_apdens_in, fact_apdiam_in, fact_gasmr, &
       pr_atm, rh, te, cair_mol_m3, cair_mol_cc )

    use module_data_mosaic_boxmod, only: idiag_sect_coag, idiag_sect_movesect, idiag_sect_newnuc
    use module_data_mosaic_main, only: mw_air, ntot_max, ntot_used
    use module_data_mosaic_aero, only: &
       inh4_a, iso4_a, jh2o, &
       dens_aer_mac, mw_aer_mac, mw_comp_a, &
       mcoag_flag1, mmovesect_flag1, mnewnuc_flag1, &
       naer, naercomp, nbin_a, nbin_a_max
    use module_data_mosaic_asecthp, only:  maxd_asize, maxd_atype, &
       nsize_aer, ntype_aer
    use module_mosaic_movesect1d, only:  move_sections_x3
    use module_mosaic_movesect3d, only:  move_sect_3d_x1
    use module_mosaic_coag1d,     only:  mosaic_coag_1box
    use module_mosaic_coag3d,     only:  mosaic_coag_3d_1box
    use module_mosaic_newnucb,    only:  mosaic_newnuc_1box

    implicit none

    !   subr arguments
    integer, intent(in) :: it_mosaic
    integer, intent(in),    dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
    integer, intent(in),    dimension(nbin_a_max) :: jhyst_leg_sv1
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(in)    :: dtchem  ! time step [s]
    real(r8), intent(in)    :: pr_atm  ! pressure [atmospheres]
    real(r8), intent(in)    :: rh      ! relative humidity [percent]
    real(r8), intent(in)    :: te      ! temperature [K]
    real(r8), intent(in)    :: cair_mol_m3, cair_mol_cc  ! air molar densities [mol-air/m3-air] and [mol-air/cm3-air]
    real(r8), intent(in)    :: fact_apmassmr, fact_apnumbmr, fact_apdens_in, fact_apdiam_in, fact_gasmr
    real(r8), intent(in)    :: rbox_sv1(ntot_used)
    real(r8), intent(in),    dimension(nbin_a_max) :: dp_dry_a_sv1, sigmag_a_sv1
    real(r8), intent(inout) :: rbox(ntot_used)  ! gases = [ppmv],  aero mass/water = [ug-aero/kg-air],  aero numb = [#/kg-air]
    real(r8), intent(inout), dimension(nbin_a_max) :: dens_dry_a_bgn, dens_dry_a  ! [g-aero/cm3-aero]
    real(r8), intent(inout), dimension(nbin_a_max) :: mass_dry_a_bgn, mass_dry_a  ! [g-aero/cm3-air]
    real(r8), intent(inout), dimension(nbin_a_max) :: dp_dry_a, sigmag_a          ! [cm] and [--]

    !   local variables
    integer  :: idiagbb_coag, istat_coag
    integer  :: idiagbb_newnuc, istat_newnuc, itype_newnuc
    integer  :: idiagcc
    integer  :: jhyst_leg_sv2(nbin_a_max)
    integer  :: jsv
    integer  :: method_coag, method_movesect

    real(r8) :: dens_nh4so4a_newnuc
    real(r8) :: fact_apdiam, fact_apdens
    real(r8) :: rh_box, rhoair_g_cc
    real(r8) :: tmpa

    real(r8) :: cnn(ntot_used), cnn_sv1(ntot_used), cnn_sv2(ntot_used)
    real(r8) :: rbox_sv2(ntot_used)
    real(r8) :: rbox_svx(ntot_used,4)

    real(r8) :: dp_dry_a_sv2(nbin_a_max)
    real(r8) :: sigmag_a_sv2(nbin_a_max)

    real(r8) :: drydens_pregrow(maxd_asize,maxd_atype)
    real(r8) :: drydens_aftgrow(maxd_asize,maxd_atype)
    real(r8) :: drymass_pregrow(maxd_asize,maxd_atype)
    real(r8) :: drymass_aftgrow(maxd_asize,maxd_atype)
    ! drydens_xxxgrow = dry density (g/cm3) of bin
    ! drymass_xxxgrow = all-component dry-mass mixing ratio (g-AP/m3-air) of bin
    !    pregrow (aftgrow) means before (after) the last call to the
    !    gas-aerosol mass-transfer solver
    ! these arrays are used by the movesect routine to calculate transfer
    !    between bins due to particle growth/shrinkage

    real(r8) :: adrydens_box(maxd_asize,maxd_atype)
    real(r8) :: adrydpav_box(maxd_asize,maxd_atype)
    real(r8) :: adryqmas_box(maxd_asize,maxd_atype)
    real(r8) :: awetdens_box(maxd_asize,maxd_atype)
    real(r8) :: awetdpav_box(maxd_asize,maxd_atype)
    ! adrydens_box = current dry density (g/cm3) of bin
    ! adrydpav_box = current mean dry diameter (cm) of bin
    ! adryqmas_box = current all-component dry-mass mixing ratio (g-AP/m3-air) of bin
    ! awetdens_box = current wet/ambient density (g/cm3) of bin
    ! awetdpav_box = current mean wet/ambient diameter (cm) of bin
    !
    real(r8) :: adrydens_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: adrydpav_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: adryqmas_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: awetdens_box_sv(maxd_asize,maxd_atype,4)
    real(r8) :: awetdpav_box_sv(maxd_asize,maxd_atype,4)


    jhyst_leg_sv2(1:nbin_a) = jhyst_leg(1:nbin_a)
    dp_dry_a_sv2( 1:nbin_a) = dp_dry_a( 1:nbin_a)
    sigmag_a_sv2( 1:nbin_a) = sigmag_a( 1:nbin_a)
    rbox_sv2(1:ntot_used) = rbox(1:ntot_used)
    fact_apdiam = fact_apdiam_in
    fact_apdens = fact_apdens_in


    ! the movesect routine puts "initial" values into the axxxxxxx_box arrays
    !    so need to have movesect turned on (mmovesect_flag1 > 0)
    !    for any of the sectional routines to function properly
    if (mmovesect_flag1 <= 0) then
       write(*,*)   &
            '*** skipping sectional_interface_1 -- mmovesect_flag1 <= 0'
       return
    end if

    ! only iunits_flagaa = 4 has been tested, so do not allow other values
    if (iunits_flagaa /= 4) then
       write(*,*)   &
            '*** sectional_interface_1 error -- iunits_flagaa /= 4'
       stop
    end if
    if (iunits_flagbb /= 1) then
       write(*,*)   &
            '*** sectional_interface_1 error -- iunits_flagbb /= 1'
       stop
    end if
    if (iunits_diagout_flagaa /= 1) then
       write(*,*)   &
            '*** sectional_interface_1 error -- iunits_diagout_flagaa /= 1'
       stop
    end if


    ! rbox units factors
    !
    ! note that the movesect and coag routines expect the aerosol mass units
    !    for rbox to be mass mixing ratios or mass concentrations
    ! they do not apply any molecular weight factors, so they cannot
    !    work with molar mixing ratios or molar concentrations
    !
    ! for trace gases,   rbox = [umol/mol-air = ppmv],     rbox*fact_gasmr = [molecules/cm^3-air]
    !    fact_gasmr = 1.0e-6_r8
    ! for aerosol mass,  rbox = [ug-AP/kg-air], rbox*fact_apmassmr = [g-AP/g-air]
    !    fact_apmassmr = 1.0e-9_r8
    ! for aerosol numb,  rbox = [#/kg-air],     rbox*fact_apnumbmr = [#/g-air]
    !    fact_apnumbmr = 1.0e-3_r8
    ! aerosol densities in mosaic aerchem and sectional routines = [g/cm^3]
    !     (mosaic densities)*fact_apdens = [g/cm^3]
    !    fact_apdens = 1.0_r8
    ! aerosol diameters in mosaic aerchem and sectional routines = [cm]
    !     (mosaic diameters)*fact_apdiam = [cm]
    !    fact_apdiam = 1.0_r8


    ! convert densities and sizes for testing purposes
    dens_nh4so4a_newnuc = dens_aer_mac(iso4_a)
    if (iunits_flagbb == 2) call sect_iface_dens_size_testbb( 1, &
         fact_apdens, fact_apdiam, dens_nh4so4a_newnuc )


    ! do some set-up tasks
    call sect_iface_bgn_end_calcs( 0, rbox, &
         dp_dry_a, sigmag_a, &
         drydens_aftgrow, drydens_pregrow,   &
         drymass_aftgrow, drymass_pregrow,   &
         adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
         adryqmas_box,   &
         fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
         it_mosaic, jaerosolstate_bgn, jaerosolstate, jhyst_leg, &
         mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
         mw_aer_mac, cair_mol_cc, cair_mol_m3)

    jsv = 1
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    ! move_sections transfers particles between sections
    ! following condensational growth / evaporative shrinking
    !
    method_movesect = mod( max(0,mmovesect_flag1), 100 )
    if (method_movesect < 50) then
       !          call move_sections_x3( iflag, iclm, jclm, k, m, rbox,   &
       !            drydens_aftgrow, drydens_pregrow,   &
       !            drymass_aftgrow, drymass_pregrow,   &
       !            adrydens_tmp, awetdens_tmp, adrydbar_tmp, awetdbar_tmp,   &
       !            adryqmas_tmp, adryqvol_tmp )
       call move_sections_x3( 1, 1, 1, 1, 1, rbox,   &
            fact_apmassmr, fact_apnumbmr,   &
            fact_apdens, fact_apdiam,   &
            drydens_aftgrow, drydens_pregrow,   &
            drymass_aftgrow, drymass_pregrow,   &
            adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
            adryqmas_box, it_mosaic, mmovesect_flag1, idiag_sect_movesect )
    else
       call move_sect_3d_x1( 1, 1, 1, 1, 1, rbox,   &
            fact_apmassmr, fact_apnumbmr,   &
            fact_apdens, fact_apdiam,   &
            drydens_aftgrow, drydens_pregrow,   &
            drymass_aftgrow, drymass_pregrow,   &
            adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
            adryqmas_box, it_mosaic, mmovesect_flag1, idiag_sect_movesect )
    end if

    jsv = 2
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    ! do new particle nucleation
    !
    if (mnewnuc_flag1 > 0) then
       idiagbb_newnuc = +200
       if (idiag_sect_newnuc <= 0) idiagbb_newnuc = 0
       itype_newnuc = 1
       rh_box = rh*0.01
       call mosaic_newnuc_1box(   &
            istat_newnuc, idiagbb_newnuc,   &
            it_mosaic, itype_newnuc, dtchem,   &
            te, pr_atm, cair_mol_cc, rh_box,   &
            dens_nh4so4a_newnuc, mw_aer_mac(iso4_a), mw_aer_mac(inh4_a),   &
            mw_comp_a(jh2o), mw_air,  &
            fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr,   &
            rbox_sv1, rbox,   &
            adrydens_box, awetdens_box,   &
            adrydpav_box, awetdpav_box, adryqmas_box )
    end if

    jsv = 3
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    ! do particle coagulation
    !
    if (mcoag_flag1 > 0) then
       idiagbb_coag = +200
       if (idiag_sect_coag <= 0) idiagbb_coag = 0
       rhoair_g_cc = cair_mol_cc*mw_air
       method_coag = mod( max(0,mcoag_flag1), 100 )
       if (method_coag < 50) then
          call mosaic_coag_1box( istat_coag,   &
               idiagbb_coag, it_mosaic,   &
               dtchem, te, pr_atm, rhoair_g_cc, rbox,   &
               fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
               adrydens_box, awetdens_box,   &
               adrydpav_box, awetdpav_box, adryqmas_box, mcoag_flag1 )
       else
          call mosaic_coag_3d_1box( istat_coag,   &
               idiagbb_coag, it_mosaic,   &
               dtchem, te, pr_atm, rhoair_g_cc, rbox,   &
               fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam,   &
               adrydens_box, awetdens_box,   &
               adrydpav_box, awetdpav_box, adryqmas_box, mcoag_flag1 )
       end if
    end if

    jsv = 4
    rbox_svx(:,jsv) = rbox(:)
    adrydens_box_sv(:,:,jsv) = adrydens_box(:,:)
    awetdens_box_sv(:,:,jsv) = awetdens_box(:,:)
    adrydpav_box_sv(:,:,jsv) = adrydpav_box(:,:)
    awetdpav_box_sv(:,:,jsv) = awetdpav_box(:,:)
    adryqmas_box_sv(:,:,jsv) = adryqmas_box(:,:)


    ! do some clean-up tasks
    call sect_iface_bgn_end_calcs( 1, rbox, &
         dp_dry_a, sigmag_a, &
         drydens_aftgrow, drydens_pregrow,   &
         drymass_aftgrow, drymass_pregrow,   &
         adrydens_box, awetdens_box, adrydpav_box, awetdpav_box,   &
         adryqmas_box,   &
         fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
         it_mosaic, jaerosolstate_bgn, jaerosolstate, jhyst_leg, &
         mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
         mw_aer_mac, cair_mol_cc, cair_mol_m3)


    idiagcc = 1
    if (nsize_aer(1)*ntype_aer > 10000) idiagcc = 0
    if (nsize_aer(1)*ntype_aer >   100) idiagcc = 0
    ! this was needed for early testing -- should not need it for large runs
    if (idiagcc > 0) &
         call mbox_sectional_diagnostics( rbox, rbox_sv1, rbox_sv2, rbox_svx, &
         dp_dry_a, dp_dry_a_sv1, dp_dry_a_sv2, &
         sigmag_a, sigmag_a_sv1, sigmag_a_sv2, &
         drydens_aftgrow, drydens_pregrow,   &
         drymass_aftgrow, drymass_pregrow,   &
         adrydens_box_sv, awetdens_box_sv, adrydpav_box_sv,   &
         awetdpav_box_sv, adryqmas_box_sv,   &
         fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
         it_mosaic, jaerosolstate_bgn, jaerosolstate, &
         jhyst_leg, jhyst_leg_sv1, jhyst_leg_sv2, &
         cair_mol_m3, cair_mol_cc, mw_aer_mac, mw_comp_a )


    ! unconvert densities and sizes
    if (iunits_flagbb == 2) call sect_iface_dens_size_testbb( 2, &
         fact_apdens, fact_apdiam, dens_nh4so4a_newnuc )


    return
  end subroutine sectional_interface_1


  !-----------------------------------------------------------------------
  subroutine sect_iface_bgn_end_calcs( ibgn_end, rbox, &
       dp_dry_a, sigmag_a, &
       drydens_aftgrow, drydens_pregrow,   &
       drymass_aftgrow, drymass_pregrow,   &
       adrydens_box, awetdens_box, adrydpav_box,   &
       awetdpav_box, adryqmas_box,   &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       it_mosaic, jaerosolstate_bgn, jaerosolstate, jhyst_leg, &
       mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a, &
       mw_aer_mac, cair_mol_cc, cair_mol_m3)
    !
    !   does some set-up and clean-up tasks
    !
    use module_data_mosaic_boxmod, only: idiag_sect_movesect
    use module_data_mosaic_main, only: &
       avogad, &
       mw_air, naer_tot, ngas_max, ntot_max, ntot_used
    use module_data_mosaic_boxmod, only: &
       kdpdry_a, kjhyst_a, knum_a, ksigmag_a, kwater_a
    use module_data_mosaic_aero, only: &
       jhyst_lo, jhyst_up, &
       mhyst_method, mhyst_force_lo, mhyst_force_up, &
       mhyst_uporlo_jhyst, mhyst_uporlo_waterhyst, &
       mmovesect_flag1, naer, naercomp, nbin_a, nbin_a_max, &
       all_solid, all_liquid, mixed, no_aerosol
    use module_data_mosaic_asecthp, only: &
       ai_phase, isize_of_ibin, itype_of_ibin, &
       maxd_acomp, maxd_asize, maxd_atype, ncomp_aer, &
       hyswptr_aer, massptr_aer, &
       dens_aer, dens_water_aer, hygro_aer 

    use module_mosaic_movesect1d, only:  test_move_sections

    implicit none

    ! subr parameters
    integer, intent(in) :: ibgn_end, it_mosaic
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
    integer, intent(inout), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(in) :: cair_mol_cc, cair_mol_m3
    real(r8), intent(inout), dimension(nbin_a_max)  :: mass_dry_a_bgn, mass_dry_a, dens_dry_a_bgn, dens_dry_a
    real(r8), intent(inout), dimension(nbin_a_max)  :: dp_dry_a, sigmag_a
    real(r8), intent(inout), dimension(naer) :: mw_aer_mac
    real(r8), intent(inout) :: rbox(ntot_used)

    real(r8), intent(inout) :: drydens_aftgrow(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: drydens_pregrow(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: drymass_aftgrow(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: drymass_pregrow(maxd_asize,maxd_atype)

    real(r8), intent(inout) :: adrydens_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: awetdens_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: adrydpav_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: awetdpav_box(maxd_asize,maxd_atype)
    real(r8), intent(inout) :: adryqmas_box(maxd_asize,maxd_atype)

    real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, &
                            fact_apdens, fact_apdiam, fact_gasmr


    ! local variables
    integer :: ibin, iphase, isize, itype, jhyst_tmp
    integer :: l, ll, mtmp, noffset
    real(r8) :: tmpa, tmph, tmpj, tmpr



    ! on first time step, call test_move_sections
    ! note that test_move_sections only executes when mmovesect_flag1 = 80nn
    if ((it_mosaic == 1) .and. (ibgn_end == 0)) then
       call test_move_sections( 1, 1, 1, 1, 1, it_mosaic, mmovesect_flag1, idiag_sect_movesect )
    end if


    if ((ibgn_end < 0) .or. (ibgn_end > 1)) return

    if (ibgn_end == 1) goto 20000


    !
    ! ibgn_end = 0
    !     initialize the a-------_box, dry----_pregrow, and dry----aftgrow 
    !     for (mhyst_method == mhyst_uporlo_jhyst),
    !         calc/load rbox(hyswptr_aer) with "hysteresis surrogate" info
    !
    adrydens_box(:,:) = 0.0
    awetdens_box(:,:) = 0.0
    adrydpav_box(:,:) = 0.0
    awetdpav_box(:,:) = 0.0
    adryqmas_box(:,:) = 0.0

    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)
       noffset = ngas_max + (ibin-1)*naer_tot

       if (jaerosolstate_bgn(ibin) /= no_aerosol) then
          drydens_pregrow(isize,itype) = dens_dry_a_bgn(ibin)
          drymass_pregrow(isize,itype) = mass_dry_a_bgn(ibin)
       else
          drydens_pregrow(isize,itype) = -1.0
          drymass_pregrow(isize,itype) =  0.0
       end if
       if (jaerosolstate(ibin) /= no_aerosol) then
          drydens_aftgrow(isize,itype) = dens_dry_a(ibin)
          drymass_aftgrow(isize,itype) = mass_dry_a(ibin)
       else
          drydens_aftgrow(isize,itype) = -1.0
          drymass_aftgrow(isize,itype) =  0.0
       end if

       ! convert drymass_xxxgrow from mosaic aerchemistry units [g-aero/cm^3-air] 
       !    to rbox mass mixing ratio units [ug-aero/kg-air]
       tmpa = 1.0e6_r8/(cair_mol_m3*mw_air*fact_apmassmr)
       drymass_pregrow(isize,itype) = drymass_pregrow(isize,itype) * tmpa
       drymass_aftgrow(isize,itype) = drymass_aftgrow(isize,itype) * tmpa

       ! convert drydens_xxxgrow from [g-aero/cm3-aero] to [???] 
       !    (currently this is not really needed because fact_apdens = 1.0)
       tmpa = 1.0_r8/fact_apdens
       drydens_pregrow(isize,itype) = drydens_pregrow(isize,itype) * tmpa
       drydens_aftgrow(isize,itype) = drydens_aftgrow(isize,itype) * tmpa

       if (mhyst_method == mhyst_uporlo_jhyst) then
          ! cnn holds jhyst_leg
          ! rbox holds jhyst_leg*(mean_hygro*dry_volume_mixing_ratio*water_density) (see below)
          !    so it can be treated reasonably in movesect and coag routines
          ! note - the tmph (below) units are essentially the same as for mass species in rbox
          if ( jaerosolstate(ibin) == all_solid  .or. &
               jaerosolstate(ibin) == all_liquid .or. &
               jaerosolstate(ibin) == mixed      ) then
             iphase = ai_phase
             tmph = 0.0
             do ll = 1, ncomp_aer(itype)
                l = massptr_aer(ll,isize,itype,iphase)
                tmph = tmph + max(0.0_r8,rbox(l))*   &
                     (hygro_aer(ll,itype)/dens_aer(ll,itype))
             end do
             tmph = max( tmph*dens_water_aer, 1.0e-35_r8 )
             rbox(kjhyst_a +noffset) = tmph * jhyst_leg(ibin)
          else
             rbox(kjhyst_a +noffset) = 0.0_r8
          end if
       end if

    end do

    return


20000 continue
    !
    ! ibgn_end = 1
    !     load dp_dry_a with current value
    !     for (mhyst_method == mhyst_uporlo_jhyst),
    !         calc/set jhyst_leg using the info in rbox(hyswptr_aer)
    !     for (mhyst_method == other), set jhyst_leg as needed
    !
    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)
       noffset = ngas_max + (ibin-1)*naer_tot

       ! need to put final value into dp_dry_a here
       dp_dry_a(ibin) = adrydpav_box(isize,itype)

       if (mhyst_method == mhyst_uporlo_waterhyst) then
          continue
       else if (mhyst_method == mhyst_force_up) then
          jhyst_leg(ibin) = jhyst_up
       else if (mhyst_method == mhyst_force_lo) then
          jhyst_leg(ibin) = jhyst_lo
       else if (mhyst_method == mhyst_uporlo_jhyst) then
          ! convert rbox(hyswptr_aer(:,:)) 
          !    from (0/1 flag) * (mean_hygro*dry_volume_mixing_ratio*water_density)
          !    back to [0/1 flag]
          ! note - the tmph (below) units are essentially the same as for mass species in rbox
          iphase = ai_phase
          tmph = 0.0
          do ll = 1, ncomp_aer(itype)
             l = massptr_aer(ll,isize,itype,iphase)
             tmph = tmph + max(0.0_r8,rbox(l))*   &
                  (hygro_aer(ll,itype)/dens_aer(ll,itype))
          end do
          tmph = max( tmph*dens_water_aer, 1.0e-35_r8 )
          l = hyswptr_aer(isize,itype)
          tmpr = max( rbox(l), 0.0_r8 )
          if (tmpr > tmph*1.0e6) then   ! here tmpr/tmph > 1.0e6
             tmpj = 1.0e6
          else
             tmpj = tmpr/tmph
          end if
          ! *** at this point, tmpj can be anything between 0 and 1 (or even > 1)
          !        because movesect and coag can "mix" upper and lower curve
          !        particles from different bins
          !     need to decide on a reasonable "threshold" that determines how a
          !        "mixed" particles is be classified
          !        (0.99 probably NOT a good threshold)
          if (tmpj > 0.5) then
             jhyst_tmp = jhyst_up
          else
             jhyst_tmp = jhyst_lo
          end if
          jhyst_leg(ibin) = jhyst_tmp
       else
          write(*,'(/a,i10/)') &
               '*** sect_iface_bgn_end_calcs - bad mhyst_method =', mhyst_method
          stop
       end if

    end do

    return
  end subroutine sect_iface_bgn_end_calcs


  !-----------------------------------------------------------------------
  subroutine sect_diags_map_to_cnn( cair_mol_cc, cair_mol_m3, &
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

    use module_mosaic_movesect1d, only:  test_move_sections

    implicit none

    ! subr parameters
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate
    integer, intent(in), dimension(nbin_a_max) :: jhyst_leg

    real(r8), intent(in)  :: cair_mol_cc, cair_mol_m3
    real(r8), intent(in)  :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr
    real(r8), intent(in)  :: rbox(ntot_used)
    real(r8), intent(in), dimension(nbin_a_max) :: dp_dry_a, sigmag_a

    real(r8), intent(out) :: cnn(ntot_max)


    ! local variables
    integer :: ibin, iphase, isize, itype, jhyst_tmp
    integer :: l, ll, lunaa, mtmp, noffset
    real(r8) :: conv_diam, conv_drym, conv_gas, conv_numb, conv_watr
    real(r8) :: tmpa, tmph, tmpj, tmpr


    ! set factors for converting rbox (and dp_dry_a) values to cnn values
    tmpa = cair_mol_m3*mw_air   ! air_density in (g/m^3)
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
  end subroutine sect_diags_map_to_cnn


  !-----------------------------------------------------------------------
  subroutine mbox_sectional_diagnostics( rbox, rbox_sv1, rbox_sv2, rbox_svx, &
       dp_dry_a, dp_dry_a_sv1, dp_dry_a_sv2, &
       sigmag_a, sigmag_a_sv1, sigmag_a_sv2, &
       drydens_aftgrow, drydens_pregrow,   &
       drymass_aftgrow, drymass_pregrow,   &
       adrydens_box_sv, awetdens_box_sv, adrydpav_box_sv,   &
       awetdpav_box_sv, adryqmas_box_sv,   &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       it_mosaic, jaerosolstate_bgn, jaerosolstate, &
       jhyst_leg, jhyst_leg_sv1, jhyst_leg_sv2, &
       cair_mol_m3, cair_mol_cc, mw_aer_mac, mw_comp_a )
    !
    !   writes test diagnostics for the combination of movesect, newnuc, coag
    !
    use module_data_mosaic_main, only: &
       mw_air, naer_tot, ngas_max, ntot_max, ntot_used
    use module_data_mosaic_boxmod, only: &
       kdpdry_a, kjhyst_a, knum_a, ksigmag_a, kwater_a, &
       dt_sec, lun_sect_172, &
       pr_atm, rh, species, te, time_sec
    use module_data_mosaic_aero, only: &
       jh2o, mhyst_method, mhyst_uporlo_waterhyst, &
       naer, naercomp, nbin_a, nbin_a_max
    use module_data_mosaic_asecthp, only: &
       ai_phase, isize_of_ibin, itype_of_ibin, &
       maxd_asize, maxd_atype

    implicit none

    ! subr parameters
    integer, intent(in) :: it_mosaic
    integer, intent(in), dimension(nbin_a_max) :: jaerosolstate, jaerosolstate_bgn
    integer, intent(in), dimension(nbin_a_max) :: jhyst_leg, jhyst_leg_sv1, jhyst_leg_sv2

    real(r8), intent(in) :: cair_mol_m3,cair_mol_cc
    real(r8), intent(in) :: rbox(ntot_used), rbox_sv1(ntot_used), rbox_sv2(ntot_used)
    real(r8), intent(in) :: rbox_svx(ntot_used,4)
    real(r8), intent(in) :: dp_dry_a(nbin_a_max), dp_dry_a_sv1(nbin_a_max), dp_dry_a_sv2(nbin_a_max)
    real(r8), intent(in) :: sigmag_a(nbin_a_max), sigmag_a_sv1(nbin_a_max), sigmag_a_sv2(nbin_a_max)
    real(r8), intent(in) :: drydens_aftgrow(maxd_asize,maxd_atype)
    real(r8), intent(in) :: drydens_pregrow(maxd_asize,maxd_atype)
    real(r8), intent(in) :: drymass_aftgrow(maxd_asize,maxd_atype)
    real(r8), intent(in) :: drymass_pregrow(maxd_asize,maxd_atype)

    real(r8), intent(in) :: adrydens_box_sv(maxd_asize,maxd_atype,4)
    real(r8), intent(in) :: awetdens_box_sv(maxd_asize,maxd_atype,4)
    real(r8), intent(in) :: adrydpav_box_sv(maxd_asize,maxd_atype,4)
    real(r8), intent(in) :: awetdpav_box_sv(maxd_asize,maxd_atype,4)
    real(r8), intent(in) :: adryqmas_box_sv(maxd_asize,maxd_atype,4)

    real(r8), intent(in) :: fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr
    real(r8), intent(in), dimension(naer) :: mw_aer_mac
    real(r8), intent(in), dimension(naercomp) :: mw_comp_a

    ! local variables
    integer  :: ibin, iphase, isize, itype
    integer  :: jsv,  l, ll, lunaa, mtmp, noffset

    real(r8) :: cnn(ntot_used), cnn_sv1(ntot_used), cnn_sv2(ntot_used)
    real(r8) :: conv_watr
    real(r8) :: tmpa, tmpb, tmph, tmpveca(4), tmpvech(4)


! regenerate cnn array
    call sect_diags_map_to_cnn( cair_mol_cc, cair_mol_m3, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       cnn_sv1, rbox_sv1, jaerosolstate_bgn, jhyst_leg_sv1, dp_dry_a_sv1, sigmag_a_sv1 )

    call sect_diags_map_to_cnn( cair_mol_cc, cair_mol_m3, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       cnn_sv2, rbox_sv2, jaerosolstate,     jhyst_leg_sv2, dp_dry_a_sv2, sigmag_a_sv2 )

    call sect_diags_map_to_cnn( cair_mol_cc, cair_mol_m3, &
       fact_apmassmr, fact_apnumbmr, fact_apdens, fact_apdiam, fact_gasmr, &
       cnn,     rbox,     jaerosolstate,     jhyst_leg,     dp_dry_a,     sigmag_a     )


! write diagnostics
    lunaa = lun_sect_172
    if (lunaa <= 0) return

    if (it_mosaic == 1) then
       write(lunaa,'(a)') 'explanation of fort.172 contents for gas and aerosol species'
       write(lunaa,'(a)') ' '
       write(lunaa,'(a)') 'column 1 = cnn  before mosaic_box_aerchemistry calcs'
       write(lunaa,'(a)') 'column 2 = cnn  after  mosaic_box_aerchemistry calcs'
       write(lunaa,'(a)') 'column 3 = cnn  after  movesect, newnuc, and coag calcs'
       write(lunaa,'(a)') 'column 4 = rbox before movesect calcs'
       write(lunaa,'(a)') 'column 5 = rbox after  movesect and before newnuc calcs'
       write(lunaa,'(a)') 'column 6 = rbox after  newnuc and before coag calcs'
       write(lunaa,'(a)') 'column 7 = rbox after  coag calcs'
       write(lunaa,'(a)') ' '
       write(lunaa,'(a)') 'species         cnn units in    rbox units       rbox units'
       write(lunaa,'(a)') 'category        code & output   in this output   in the code'
       write(lunaa,'(a)') '--------        -------------   --------------   -----------'
       write(lunaa,'(a)') 'gas             molec/cm^3      mol/mol-air      mol/mol-air'
       write(lunaa,'(a)') 'aer dpdry       um              cm               cm'
       write(lunaa,'(a)') 'aer number      #/cm^3          #/mol-air        #/g-air'         
       write(lunaa,'(a)') 'aer water       kg/m^3          mol/mol-air      ug/g-air'
       write(lunaa,'(a)') 'aer mass        umol/m^3        mol/mol-air      ug/g-air'
       write(lunaa,'(a)') '(note:  for oin, bc, and oc, mw=1.0 so umol = ug)'
    end if


    write(lunaa,'(//a,2i5)') 'mbox_sectional_diagnostics -- it', it_mosaic
    write(lunaa,'(a,1p,4e12.4)') 'time_sec, dt_sec               ', time_sec, dt_sec
    write(lunaa,'(a,1p,4e12.4)') 'te, pr_atm, rh_pct, cair_mol_cc', te, pr_atm, rh, cair_mol_cc


! gases - h2so4, hno3, hcl, nh3
    write(lunaa,'(a)')
    do l = 1, 4
       write(lunaa,'(i4,2x,a,1p,8e12.4)') &
            l, species(l)(1:12), cnn_sv1(l), cnn_sv2(l), cnn(l),   &
            rbox_svx(l,1:4)*fact_gasmr
    end do


! aerosol parameters, 1 bin at a time
    do ibin = 1, nbin_a
       isize = isize_of_ibin(ibin)
       itype = itype_of_ibin(ibin)

       iphase = ai_phase
       do jsv = 1, 4
!         if (mhyst_method == mhyst_uporlo_waterhyst) then
!   01-oct-2014 - use this conversion for all mhyst_method
             if (iunits_diagout_flagaa > 0) then
                tmpa = cair_mol_m3*mw_air   ! air_density in (g/m^3)
                conv_watr = 1.0e3_r8 /(fact_apmassmr*tmpa)
                tmpvech(jsv) = conv_watr
             else
                tmpvech(jsv) = 1.0
             end if
!         else
!            tmph = 0.0
!            do ll = 1, ncomp_aer(itype)
!               l = massptr_aer(ll,isize,itype,iphase)
!               tmph = tmph + max(0.0_r8,rbox_svx(l,jsv))*   &
!                    (hygro_aer(ll,itype)/dens_aer(ll,itype))
!            end do
!            tmph = max( tmph*dens_water_aer, 1.0e-35_r8 )
!            tmpvech(jsv) = tmph
!         end if
       end do

       write(lunaa,'(a)')
       write(lunaa,'(4x,2x,a,1p,2i12)') &
            'jaero_state ', jaerosolstate_bgn(ibin), jaerosolstate(ibin)
       noffset = ngas_max + (ibin-1)*naer_tot
       do ll = 1, naer_tot
          !             if ((ll >= kmsa_a) .and. (ll <= klim2_a)) cycle
          l = noffset + ll
          if (ll == kjhyst_a) then
             tmpveca(1:4) = rbox_svx(l,1:4)/tmpvech(1:4)
             !              write(lunaa,'(i4,2x,a,1p,8e12.4)') &
             !              -1, 'water_jhyst ', water_a_hyst(ibin)
          else if (ll == kdpdry_a) then
             tmpveca(1:4) = dp_dry_a_sv2(ibin)*fact_apdiam
          else if (ll == ksigmag_a) then
             tmpveca(1:4) = sigmag_a_sv2(ibin)*fact_apdiam
          else
             tmpveca(1:4) = rbox_svx(l,1:4)
             if (iunits_diagout_flagaa > 0) then
                if (ll == knum_a) then
                   ! convert (1)  #/m^3 OR (2) #/m^3 OR (3) #/g-air
                   ! to #/mol-air
                   tmpa = mw_air*fact_apnumbmr
                   tmpveca(1:4) = tmpveca(1:4)*tmpa
                else if (ll >= kwater_a) then
                   ! convert (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
                   ! to mol/mol-air
                   tmpb = mw_comp_a(jh2o)
                   if (ll > kwater_a) tmpb = mw_aer_mac(ll-kwater_a)
                   tmpa = (mw_air/tmpb)*fact_apmassmr
                   tmpveca(1:4) = tmpveca(1:4)*tmpa
                end if
             end if ! (iunits_diagout_flagaa > 0)
          end if
          write(lunaa,'(i4,2x,a,1p,8e12.4)') &
               l, species(l)(1:12), cnn_sv1(l), cnn_sv2(l), cnn(l), &
               tmpveca(1:4)
       end do

       tmpa = fact_apdiam
       write(lunaa,'(4x,2x,a,1p,36x,5e12.4)') &
            'adrydpav    ', &
            adrydpav_box_sv(isize,itype,1:4)*tmpa
       write(lunaa,'(4x,2x,a,1p,36x,5e12.4)') &
            'awetdpav    ', &
            awetdpav_box_sv(isize,itype,1:4)*tmpa

       tmpa = fact_apdens
       write(lunaa,'(4x,2x,a,1p,2e12.4,12x,5e12.4)') &
            'adrydens    ', &
            drydens_pregrow(isize,itype)*tmpa, drydens_aftgrow(isize,itype)*tmpa, &
            adrydens_box_sv(isize,itype,1:4)*tmpa
       write(lunaa,'(4x,2x,a,1p,36x,5e12.4)') &
            'awetdens    ', &
            awetdens_box_sv(isize,itype,1:4)*tmpa

       if (iunits_diagout_flagaa <= 0) then
          tmpa = 1.0_r8
       else
          ! convert (1) ug/m^3 OR (2) g/m^3 OR (3) g/g-air
          ! to g/mol-air
          tmpa = mw_air*fact_apmassmr
       end if
       write(lunaa,'(4x,2x,a,1p,2e12.4,12x,6e12.4)') &
            'adryqmas    ', &
            tmpa*drymass_pregrow(isize,itype), tmpa*drymass_aftgrow(isize,itype), &
            tmpa*adryqmas_box_sv(isize,itype,1:4)
    end do

    !BSINGH - 05/28/2013(RCE updates Commented the following if construct)
    ! rce - not sure what this was for, but turn it off for allocated memory version
    !if (it_mosaic == 1) then
    !   write(lunaa,'(a)')
    !   do ll = 1, 10
    !      l = ntot_used + ll
    !      write(lunaa,'(i4,2x,a,1p,2e12.4)') &
    !           l, species(l)(1:12)
    !   end do
    !end if
    !BSINGH - 05/28/2013(RCE updates ENDS)

    return
  end subroutine mbox_sectional_diagnostics


  !-----------------------------------------------------------------------
  subroutine sect_iface_dens_size_testbb( iflagaa, &
       fact_apdens, fact_apdiam, dens_nh4so4a_newnuc )
    !
    !   when iunits_flagbb = 2, densities and sizes in module_data_mosaic_asecthp
    !   are temporarily changes from (g/cm^3) to (kg/m^3), and (cm) to (m)
    !
    use module_data_mosaic_asecthp, only: &
       maxd_acomp, &
       dens_aer, dens_mastercomp_aer, dens_water_aer, &
       dcen_sect, dcut_sect, dhi_sect, dlo_sect, &
       volumcen_sect, volumcut_sect, volumhi_sect, volumlo_sect

    implicit none

    ! subr parameters
    integer,  intent(in)    :: iflagaa
    real(r8), intent(inout) :: fact_apdens, fact_apdiam
    real(r8), intent(inout) :: dens_nh4so4a_newnuc

    ! local variables
    real(r8) :: tmpa

    real(r8), save ::   &
         dens_nh4so4a_newnuc_tmp,  &
         dens_water_aer_tmp,  & !BSINGH - 05/28/2013(RCE updates -deleted a var here)         
         dens_mastercomp_aer_tmp( maxd_acomp )

    !BSINGH - 05/28/2013(RCE updates -deleted 8 vars here)       


    if (iunits_flagbb /= 2) return
    if (iflagaa == 1) goto 10000
    if (iflagaa == 2) goto 20000
    write(*,'(2a,2(1x,i10))') &
         '*** sect_iface_dens_size_testbb fatal error - ', &
         'bad iflagaa = ', iflagaa
    stop


10000 continue
    fact_apdens = 0.001_r8   ! converts (kg/m^3) to (g/cm^3)
    fact_apdiam = 100.0_r8   ! converts (m) to (cm)

    dens_mastercomp_aer_tmp(:) = dens_mastercomp_aer(:)
    dens_aer_tmp(:,:) = dens_aer(:,:)
    dens_water_aer_tmp = dens_water_aer
    dens_nh4so4a_newnuc_tmp = dens_nh4so4a_newnuc
    dlo_sect_tmp( :,:) = dlo_sect( :,:)
    dhi_sect_tmp( :,:) = dhi_sect( :,:)
    dcen_sect_tmp(:,:) = dcen_sect(:,:)
    dcut_sect_tmp(:,:) = dcut_sect(:,:)
    volumlo_sect_tmp( :,:) = volumlo_sect( :,:)
    volumhi_sect_tmp( :,:) = volumhi_sect( :,:)
    volumcen_sect_tmp(:,:) = volumcen_sect(:,:)
    volumcut_sect_tmp(:,:) = volumcut_sect(:,:)

    tmpa = 1.0_r8/fact_apdens
    dens_mastercomp_aer(:) = dens_mastercomp_aer(:)*tmpa
    dens_aer(:,:) = dens_aer(:,:)*tmpa
    dens_water_aer = dens_water_aer*tmpa
    dens_nh4so4a_newnuc = dens_nh4so4a_newnuc*tmpa

    tmpa = 1.0_r8/fact_apdiam
    dlo_sect( :,:) = dlo_sect( :,:)*tmpa
    dhi_sect( :,:) = dhi_sect( :,:)*tmpa
    dcen_sect(:,:) = dcen_sect(:,:)*tmpa
    dcut_sect(:,:) = dcut_sect(:,:)*tmpa

    tmpa = tmpa*tmpa*tmpa
    volumlo_sect( :,:) = volumlo_sect( :,:)*tmpa
    volumhi_sect( :,:) = volumhi_sect( :,:)*tmpa
    volumcen_sect(:,:) = volumcen_sect(:,:)*tmpa
    volumcut_sect(:,:) = volumcut_sect(:,:)*tmpa

    return


20000 continue
    fact_apdens = 1.0_r8
    fact_apdiam = 1.0_r8

    dens_mastercomp_aer(:) = dens_mastercomp_aer_tmp(:)
    dens_aer(:,:) = dens_aer_tmp(:,:)
    dens_water_aer = dens_water_aer_tmp
    dens_nh4so4a_newnuc = dens_nh4so4a_newnuc_tmp

    dlo_sect( :,:) = dlo_sect_tmp( :,:)
    dhi_sect( :,:) = dhi_sect_tmp( :,:)
    dcen_sect(:,:) = dcen_sect_tmp(:,:)
    dcut_sect(:,:) = dcut_sect_tmp(:,:)

    volumlo_sect( :,:) = volumlo_sect_tmp( :,:)
    volumhi_sect( :,:) = volumhi_sect_tmp( :,:)
    volumcen_sect(:,:) = volumcen_sect_tmp(:,:)
    volumcut_sect(:,:) = volumcut_sect_tmp(:,:)

    return


  end subroutine sect_iface_dens_size_testbb


  !-----------------------------------------------------------------------


  end module module_mosaic_sect_intr

