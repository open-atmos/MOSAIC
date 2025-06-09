! Aerosol-phase Species
! subroutine for printing output at iprint time steps
!
!--------------------------------------------------------------------
!
      subroutine print_aer(idum)
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud

	use module_data_mosaic_asect, only:  &
	    dcen_sect, dlo_sect, dhi_sect,  &
	    itype_of_ibin, itype_md1_of_itype, itype_md2_of_itype,  &
	    nsize_aer, ntype_aer, ntype_md1_aer, ntype_md2_aer,  &
	    xcut_atype_md1, xcut_atype_md2

      implicit none

      integer idum, lun
      integer l, ibin, lbin, je, iv, ifix
      real(r8) :: timeofday, time_dum, cppb(ntot_max), conv1, conv2, conv3, conv4
      real(r8) :: dlogDp
      real(r8) :: SOA_tot, SOA_tot1, SOA_tot2, SOA_tot3, SOA_tot4, SOA_tot5, SOA_tot6, OC_tot, &
		  avol_tot, anum_tot, Dp_cutoff
      real(r8) :: nh4so4ratio(naerbin), no3so4ratio(naerbin), waterratio(naerbin), satratio_nh4no3(naerbin)
      real(r8) :: kg_satratio(naerbin)
      real(r8) :: aerosol_water				! mosaic func



      conv1 = 1.e15/avogad	 ! converts (molec/cc) to (nmol/m^3)
      conv2 = 1./conv1		 ! converts (nmol/m^3) to (molec/cc)
      conv3 = conv2/cair_mlc*ppb ! converts (nmol/m^3) to (ppbv)
!      conv4 = 1.e-3		 ! converts (nmol/m^3) to (umol/m^3)
      conv4 = 1.0		 ! converts (nmol/m^3) to (nmol/m^3)



!      dlogDp = 1.0/64	         ! 110 SMPS bins
!      dlogDp = 0.055243448       ! 44 FIMS bins
!      dlogDp = 1.0		 ! monodisperse
	dlogDp = 0.0236218	! 128 bins

      do l=1,ngas_max
        cppb(l) = (cnn(l)/cair_mlc)*ppb	! converting (molecules/cc) to (ppb)
      enddo

      if(idum .eq. 0)then
        time_dum = 0.0
        timeofday = time_UTC_beg
      else
        time_dum = time_hrs
        timeofday = time_UTC
      endif


      SOA_tot1 = 0.0
      SOA_tot2 = 0.0
      SOA_tot3 = 0.0
      SOA_tot4 = 0.0
      SOA_tot5 = 0.0
      SOA_tot6 = 0.0
      OC_tot   = 0.0
      avol_tot = 0.0
      anum_tot = 0.0
      do 30 ibin = 1, nbin_a
	SOA_tot1 = SOA_tot1 + aer(iapi1_a,jtotal,ibin)*1.e-3*mw_aer_mac(iapi1_a) ! ug/m^3 C5H12O5
        SOA_tot2 = SOA_tot2 + aer(iapi2_a,jtotal,ibin)*1.e-3*mw_aer_mac(iapi2_a) ! ug/m^3 C5H12O6
        SOA_tot3 = SOA_tot3 + aer(iapi3_a,jtotal,ibin)*1.e-3*mw_aer_mac(iapi3_a) ! ug/m^3 SVOC1
        SOA_tot4 = SOA_tot4 + aer(iapi4_a,jtotal,ibin)*1.e-3*mw_aer_mac(iapi4_a) ! ug/m^3 Dimer1
        SOA_tot5 = SOA_tot5 + aer(ilim1_a,jtotal,ibin)*1.e-3*mw_aer_mac(ilim1_a) ! ug/m^3 SVOC2
        SOA_tot6 = SOA_tot6 + aer(ilim2_a,jtotal,ibin)*1.e-3*mw_aer_mac(ilim2_a) ! ug/m^3 Dimer2
	OC_tot   = OC_tot   + aer(ioc_a,jtotal,ibin)  *1.e-3	   ! ug/m^3
        avol_tot = avol_tot + (num_a(ibin)*0.523598776*(Dp_dry_a(ibin)*1.e4)**3)  ! um^3/cc(air)
	anum_tot = anum_tot + num_a(ibin)
30    continue

        SOA_tot = SOA_tot1 + SOA_tot2 + SOA_tot3 + SOA_tot4 + SOA_tot5 + SOA_tot6



      do 20 ibin = 1, nbin_a

	vol_dry_a(ibin) = num_a(ibin)*0.523598776*(Dp_dry_a(ibin)*1.e4)**3

        if(mc(jc_h,ibin).gt.0.0)then
          pH(ibin) = -dlog10(mc(jc_h,ibin))
        else
          pH(ibin) = -99.0
	endif


	if(vol_dry_a(ibin) .gt. 0.0)then
	  waterratio(ibin) = (water_a(ibin)*1.e12/18.0)/vol_dry_a(ibin)
	else
	  waterratio(ibin) = -99.0
	endif



! write Dp only once
        if(idum .eq. 0)then
	  write(lun_Dpbins,*)dcen_sect(ibin,1)*1.e7	! nm
	endif



! write to "soa.aer.txt" file
!        if (iwrite_aer_species > 0)   &
!        write(lun_species,101)    &
!        time_dum, 		  &  ! hours since start of simulation
!        dcen_sect(ibin,1)*1.e7,   &  ! nm
!        num_a(ibin)/dlogDp,	  &  ! #/cc(air)
!	aer(iso4_a,jtotal,ibin),  &  ! nmol/m3
!        aer(inh4_a,jtotal,ibin),  &  ! nmol/m3
!	aer(ino3_a,jtotal,ibin),  &  ! nmol/m3
!	water_a(ibin)*1.e12/18.0     ! nmol/m3

20    continue

        if (iwrite_aer_species > 0)   &
        write(lun_numsizedist,101)(num_a(ibin)/dlogDp, ibin = 1, nbin_a)  ! #/cc

        if (iwrite_aer_species > 0)   &
        write(lun_volsizedist,101)(vol_dry_a(ibin)/dlogDp, ibin = 1, nbin_a)  ! um3/cc 

        if (iwrite_aer_species > 0)   &
        write(lun_pHsizedist,101)(pH(ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
        write(lun_nh4so4sizedist,101)(nh4so4ratio(ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
        write(lun_no3so4sizedist,101)(no3so4ratio(ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
        write(lun_kelvinsizedist,101)(kelvin(ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
        write(lun_ah2osizedist,101)(aH2O_a(ibin), ibin = 1, nbin_a)

!        if (iwrite_aer_species > 0)   &
!        write(lun_kgsatsizedist,101)(kg(ihno3_g,ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
        write(lun_kgsatsizedist,101)(sigma_soln(ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
        write(lun_satratiosizedist,101)(satratio_nh4no3(ibin), ibin = 1, nbin_a)

        if (iwrite_aer_species > 0)   &
	write(lun_species,102)time_UTC,time_dum,anum_tot,avol_tot,cppb(kh2so4),cppb(khno3),cppb(knh3), &
	                  cppb(kso2),cppb(kno),cppb(kno2),cppb(ko3),cnn(koh)


101   format(128(e13.6, 2x))
102   format(f9.4,2x,f9.4,10(2x,e13.6))


      return
      end subroutine print_aer

