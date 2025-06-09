module module_print_aer
! Aerosol-phase Species
! subroutine for printing output at iprint time steps
!
!--------------------------------------------------------------------
!
  implicit none
  private

  public:: print_aer
contains
      subroutine print_aer(idum, jaerosolstate,isteps_ASTEM,  &
           iter_MESA,aer,gas,electrolyte,mc,num_a,Dp_dry_a,Dp_wet_a, &
           area_dry_a,area_wet_a,mass_wet_a,mass_dry_a,water_a)

      use module_data_mosaic_main
      use module_data_mosaic_boxmod
      use module_data_mosaic_aero

      implicit none

      !Subroutine Arguments
      integer, intent(in) :: idum, isteps_ASTEM
      integer, intent(in), dimension(nbin_a_max) :: jaerosolstate,iter_MESA

      real(r8), intent(in), dimension(nbin_a_max) :: num_a,Dp_dry_a,Dp_wet_a,area_dry_a,water_a
      real(r8), intent(in), dimension(nbin_a_max) :: area_wet_a,mass_wet_a,mass_dry_a
      real(r8), intent(in), dimension(ngas_aerchtot) :: gas
      real(r8), intent(in), dimension(Ncation,nbin_a_max) :: mc
      real(r8), intent(in), dimension(naer,3,nbin_a_max) :: aer
      real(r8), intent(in), dimension(nelectrolyte,3,nbin_a_max) :: electrolyte

      !Local Variables
      integer l, ibin, lbin, je, iv
      real(r8) :: timeofday, time_dum, cppb(ntot_max),conv1,conv2,conv3,conv4
      real(r8) :: dlogDpdry(nbin_a), dlogDpwet(nbin_a)
      real(r8), allocatable, dimension(:) :: pH


      if(.not.allocated(pH)) allocate(pH(nbin_a_max))


      conv1 = 1.e15/avogad	 ! converts (molec/cc) to (nmol/m^3)
      conv2 = 1./conv1		 ! converts (nmol/m^3) to (molec/cc)
      conv3 = conv2/cair_mlc*ppb ! converts (nmol/m^3) to (ppbv)
!      conv4 = 1.e-3		 ! converts (nmol/m^3) to (umol/m^3)
      conv4 = 1.0		 ! converts (nmol/m^3) to (nmol/m^3)


      if(idum .eq. 0)then
        time_dum = 0.0
        timeofday = time_UTC_beg
      else
        time_dum = time_hrs
        timeofday = time_UTC
      endif




! calculate dlogDpdry and dlogDpwet
      if(nbin_a .le. 2)then

        do ibin = 1, nbin_a
          dlogDpdry(ibin) = -1.0
          dlogDpwet(ibin) = -1.0
        enddo

      else

        dlogDpdry(1) = log10(Dp_dry_a(2)/Dp_dry_a(1))
        dlogDpwet(1) = log10(Dp_wet_a(2)/Dp_wet_a(1))

        dlogDpdry(nbin_a) = log10(Dp_dry_a(nbin_a)/Dp_dry_a(nbin_a-1))
        dlogDpwet(nbin_a) = log10(Dp_wet_a(nbin_a)/Dp_wet_a(nbin_a-1))

        do ibin = 2, nbin_a-1
          dlogDpdry(ibin)=log10( (Dp_dry_a(ibin+1) + Dp_dry_a(ibin))/   &
                                 (Dp_dry_a(ibin-1) + Dp_dry_a(ibin)) )

          dlogDpwet(ibin)=log10( (Dp_wet_a(ibin+1) + Dp_wet_a(ibin))/   &
                                 (Dp_wet_a(ibin-1) + Dp_wet_a(ibin)) )
        enddo

      endif




! output species distribution
      if(iwrite_aer_species > 0)then

        if(idum .eq. 0) then
          write(lun_species,100)	! header
        end if

100   format(   &
       ' UTC(hr)    t(hr)  ibin  jstate iMESA iASTEEM  pH     ',   &
       ' Dpdry(um)     Dpwet(um)     N(#/cc)       H2SO4(g)     ',   &
       ' HNO3(g)       HCl(g)        NH3(g)        MSA(g)       ',   &
       ' ARO1(g)       ARO2(g)       ALK1(g)       OLE1(g)      ',   &
       ' API1(g)       API2(g)       LIM1(g)       LIM2(g)      ',   &
       ' SO4(t)        NO3(t)        Cl(t)         MSA(t)       ',   &
       ' NH4(t)        Na(t)         Ca(t)         OIN(t)       ',   &
       ' BC(t)         OC(t)         ARO1(t)       ARO2(t)      ',   &
       ' ALK1(t)       OLE1(t)       API1(t)       API2(t)      ',   &
       ' LIM1(t)       LIM2(t)       H2O(l)        SO4=(l)      ',   &
       ' NO3-(l)       Cl-(l)        MSA-(l)       H+(l)        ',   &
       ' NH4+(l)       Na+(l)        Ca++(l)       AmSO4(s)     ',   &
       ' Lvcite(s)     NH4HSO4(s)    NH4MSA(s)     NH4NO3(s)    ',   &
       ' NH4Cl(s)      Na2SO4(s)     Na3HSO4(s)    NaHSO4(s)    ',   &
       ' NaMSA(s)      NaNO3(s)      NaCl(s)       CaNO3(s)     ',   &
       ' CaCl2(s)      CaMSA2(s)     CaSO4(s)      CaCO3(s)     ',   &
       ' Sdry          Mdry          Swet          Mwet       ',   &
       ' dN/dlogDpdry  dN/dlogDpwet  dS/dlogDpdry  dM/dlogDpdry ',   &
       ' dS/dlogDpwet  dM/dlogDpwet  dlogDpdry     dlogDpwet')

        do ibin = 1, nbin_a

          if(mc(jc_h,ibin) .gt. 0.0)then
            pH(ibin) = -log10(mc(jc_h,ibin))
          else
            pH(ibin) = 0.0
          endif

          write(lun_species,101)   &
          timeofday,				   		&  ! time of day in UTC
          time_dum, 				   		&  ! hours since start of simulation
          ibin,					   		&  ! integer
          phasestate(jaerosolstate(ibin)),	   		&  ! 1=solid, 2=liquid, 3=mixed
          iter_mesa(ibin),			   		&  ! number of mesa iterations
          isteps_ASTEM,				   		&  ! number of astem steps
          pH(ibin),                                		&  ! pH
          Dp_dry_a(ibin)*1.e4,			   		&  ! microns (same as Dp_dry_a)
          Dp_wet_a(ibin)*1.e4,			   		&  ! microns
          num_a(ibin),				   		&  ! #/cc(air)
          gas(ih2so4_g)*conv4,			   		&  ! nmol/m^3
          gas(ihno3_g)*conv4,			   		&  ! nmol/m^3
          gas(ihcl_g)*conv4, 			   		&  ! nmol/m^3
          gas(inh3_g)*conv4,			   		&  ! nmol/m^3
          gas(imsa_g)*conv4,			   		&  ! nmol/m^3
          gas(iaro1_g)*conv4,			   		&  ! nmol/m^3
          gas(iaro2_g)*conv4,			   		&  ! nmol/m^3
          gas(ialk1_g)*conv4,			   		&  ! nmol/m^3
          gas(iole1_g)*conv4,			   		&  ! nmol/m^3
          gas(iapi1_g)*conv4,			   		&  ! nmol/m^3
          gas(iapi2_g)*conv4,			   		&  ! nmol/m^3
          gas(ilim1_g)*conv4,			   		&  ! nmol/m^3
          gas(ilim2_g)*conv4,			   		&  ! nmol/m^3
          aer(iso4_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ino3_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(icl_a, jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(imsa_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(inh4_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ina_a, jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ica_a, jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ioin_a,jtotal,ibin)*conv4,		   	&  ! ng/m^3
          aer(ibc_a,jtotal,ibin)*conv4,		   		&  ! ng/m^3
          aer(ioc_a,jtotal,ibin)*conv4,		   		&  ! ng/m^3
          aer(iaro1_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(iaro2_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ialk1_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(iole1_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(iapi1_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(iapi2_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ilim1_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ilim2_a,jtotal,ibin)*conv4,		   	&  ! nmol/m^3
          water_a(ibin)*55.555*1.e9*conv4,	   		&  ! nmol/m^3(air)
          aer(iso4_a,jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ino3_a,jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          aer(icl_a, jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          aer(imsa_a,jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          mc(jc_h,ibin)*1.e9*water_a(ibin)*conv4,	   	&  ! nmol/m^3
          aer(inh4_a,jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ina_a, jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          aer(ica_a, jliquid,ibin)*conv4,		   	&  ! nmol/m^3
          (electrolyte(je,jsolid,ibin)*conv4, je=1,nsalt),	&  ! 15 nmol/m^3
          electrolyte(jcaso4,jsolid,ibin)*conv4,     		&  ! nmol/m^3
          electrolyte(jcaco3,jsolid,ibin)*conv4,     		&  ! nmol/m^3
          area_dry_a(ibin)*1.e8,			  	&  ! area_dry_a is in cm^2/cc. output in um^2/cc
          mass_dry_a(ibin)*1.e12,			   	&  ! mass_dry_a is in g/cc. output in ug/m^3
          area_wet_a(ibin)*1.e8,			   	&  ! area_wet_a is in cm^2/cc. output in um^2/cc
          mass_wet_a(ibin)*1.e12,			   	&  ! mass_wet_a is in g/cc. output in ug/m^3
          num_a(ibin)/dlogDpdry(ibin),		  		&  ! dN/dlogDpdry (1/cc)
          num_a(ibin)/dlogDpwet(ibin),		   		&  ! dN/dlogDpwet (1/cc)
          area_dry_a(ibin)*1.e8/dlogDpdry(ibin),	   	&  ! dS/dlogDpdry in um^2/cc
          mass_dry_a(ibin)*1.e12/dlogDpdry(ibin),	   	&  ! dM/dlogDpdry in ug/m^3
          area_wet_a(ibin)*1.e8/dlogDpwet(ibin),	   	&  ! dS/dlogDpwet in um^2/cc
          mass_wet_a(ibin)*1.e12/dlogDpwet(ibin),	   	&  ! dM/dlogDpwet in ug/m^3
          dlogDpdry(ibin),   &
          dlogDpwet(ibin)
        end do

        write(lun_species,*)'  '  ! blank line


101   format(f7.3,2x,f8.3,2x,i4,2x,a6,1x,i4,2x,i4,3x,f6.3,72(2x,e28.20))

      endif ! iwrite_aer_species




      return
      end subroutine print_aer

    end module module_print_aer
