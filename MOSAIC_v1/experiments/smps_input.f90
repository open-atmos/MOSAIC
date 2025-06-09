! Idealized case conditions
! 
      Program aerinput
      implicit none

      integer jhyst, ibin, nbin_a
	real Dp_lo, Dp_up, Dp, Dp_wet, Dp_dry, Dp_nm, dlogDp, dNdlogDp_tot, Dp_cutoff
        real factor_Dp, Dpmin, Dpmax, dlogDp_max
	real DpgN_uf, DpgN_ait, DpgN_acc
	real sig_uf, sig_ait, sig_acc, sigma
	real vol_wet, vol_dry, vol_AS, vol_OC, vol_tot, vol_small, vol_large, term_AS, term_OC
        real num, num_small, num_large 
	real ratio_OC_AS, ratio_OC_SC
	real mass_AS, mass_SC, mass_OC, mass_tot, mass_OC_small, mass_OC_large
	

	real T_K, P, R, pi
	real Dg, MW, kg, Fkn, Kn, GR
	real z, speed, accom, mean_molecular_speed, freepath
	real fuchs_sutugin, sum_kg, sum_kg1, sum_kg2

	real rho_AS, rho_SC, rho_OC, rho_W, mw_AS, mw_SC, mw_OC, kappa, m0_AS, aw, fw
        real water, so4, no3, cl, nh4, msa, aro1
	real aro2, aro3, aro4, alk1, ole1
	real api1, api2, api3, api4
	real lim1, lim2, co3, na, ca, oin, oc, bc

! function
	real dNdlogDp

! input files

      open(9,  file = 'initial_sizedist.txt')	! observed initial size distribution



! output files
	open(10, file = 'aer.inp')
	open(12, file = 'so4.txt')
	open(13, file = 'nh4.txt')
	open(14, file = 'oc.txt')
	open(15, file = 'num.txt')



	pi = 3.141592654
	T_K = 298.15	! temperture (K)
	P   = 1.0		! pressure (atm)
	Dg  = 0.1		! gas-phase diffusivity (cm^2/s)
	MW  = 150		! molecular weight of diffusing species (g/mol)
	accom = 0.1


! set species parameters
	rho_AS = 1.8 	! Ammonium sulfate (g/cc)
	rho_SC = 2.2 	! NaCl (g/cc)
	rho_OC = 1.4 	! Organics (g/cc)
	rho_W  = 1.0	! water (g/cc)

	mw_AS  = 132.0	! Ammonium sulfate (g/mol)
	mw_SC  = 58.5	! NaCl (g/mol)
	mw_OC  = 1.0	! Organics (g/mol) not to be used
	
	m0_AS  = 15.6867! molality [mol/kg(water)] of AS at 50% RH
!	m0_AS  = 7.06458! molality [mol/kg(water)] of AS at 75% RH

	kappa  = 0.1	! hygroscopicity of OC
	aw     = 0.5	! water activity at 50% RH
	
! IMPORTANT: Set fw
	fw     = 0.0	! set to 1.0 if water is present at 50% RH, set to 0.0 under dry conditions
! IMPORTANT


! set up bin structure parameters
        nbin_a  = 44	! Chamber  user specified total number of bins: 10, 30, 60, 120
        vol_tot = 0.0	! initialize
        mass_tot = 0.0	! initialize

	Dp_cutoff = 100.0 ! [nm] 7_28_2014 Chamber

	
	
! Begin loop over bins	
	do ibin = 1, nbin_a

          read(9,*)Dp_nm, num	! num is in [#/cc]

	  Dp_wet = Dp_nm/1000.0		! convert nm to um




	  if(Dp_wet .lt. Dp_cutoff/1000.0)then	! small seed

	    num = num*1.0

	    vol_wet = num * 0.523599 * Dp_wet**3.0	! um^3/cc(air)
	  
            num_small = num_small + num
            vol_small = vol_small + vol_wet
	    
	    ratio_OC_AS = 6.45 ! T1: 10.012, B: 6.45					  ! small mode
!	    ratio_OC_AS = rho_OC*(Dp_nm**3.0-40.0**3.0)/(rho_AS*40.0**3.0) ! small seed: 42 nm AS core + thick OC coating


! calculate masses of AS and OC

	    term_AS = 1.0/rho_AS + fw*1000.0/(rho_W * mw_AS * m0_AS)
	    term_OC = 1.0/rho_OC + fw*(rho_W*kappa/rho_OC) * aw/(1.0 - aw)
	    mass_AS = vol_wet/(term_AS + ratio_OC_AS*term_OC)	! ug/m^3(air)
	    mass_OC = mass_AS * ratio_OC_AS			! ug/m^3(air)	    

            mass_OC_small = mass_OC_small + mass_OC
	    
	  else	! large seed
	  
	    num = num*1.0

	    vol_wet = num * 0.523599 * Dp_wet**3.0	! um^3/cc(air)

            num_large = num_large + num
            vol_large = vol_large + vol_wet
	    
	    ratio_OC_AS = 6.45 ! T1: 10.012, B: 6.45					  ! large mode
!	    ratio_OC_AS = rho_OC*(Dp_nm**3.0-50.0**3.0)/(rho_AS*50.0**3.0) ! 40 nm AS core + thick OC coating
	    

! calculate masses of AS and OC

	    term_AS = 1.0/rho_AS + fw*1000.0/(rho_W * mw_AS * m0_AS)
	    term_OC = 1.0/rho_OC + fw*(rho_W*kappa/rho_OC) * aw/(1.0 - aw)
	    mass_AS = vol_wet/(term_AS + ratio_OC_AS*term_OC)	! ug/m^3(air)
	    mass_OC = mass_AS * ratio_OC_AS			! ug/m^3(air)	    

            mass_OC_large = mass_OC_large + mass_OC
	  endif


! calculate volumes of AS and OC
	  vol_AS  = mass_AS/rho_AS
	  vol_OC  = mass_OC/rho_OC

	  vol_dry = vol_AS + vol_OC			! volume of AS and OC (excluding any water)
	  Dp_dry  = (1.909859*vol_dry/num)**(1.0/3.0)	! diameter excluding water
	  
	  
	  If(fw .eq. 0.0)then
	    Dp = Dp_wet		! no need to worry about correcting Dp
	  else
	    Dp = Dp_dry		! need to correct Dp
	  endif


! total
	  vol_tot = vol_tot + vol_wet
	  mass_tot = mass_tot + mass_AS + mass_OC	! dry mass only

	  sigma = 1.0			! dummy never used
	  jhyst = 1
	  water = 0.0
	  so4  = 1.0*mass_AS/mw_AS	! umol/m^3(air)
	  no3  = 0.0
	  cl   = 1.0*mass_SC/mw_SC	! umol/m^3(air)
	  nh4  = 2.0*mass_AS/mw_AS	! umol/m^3(air)
	  msa  = 0.0
	  aro1 = 0.0
	  aro2 = 0.0
	  aro3 = 0.0
	  aro4 = 0.0
	  alk1 = 0.0
	  ole1 = 0.0

! 7_30_2014 & 8_20_2015
	  api1 = 0.00*mass_OC/152.0	! THHP	 MW = 152
	  api2 = 0.00*mass_OC/168.0	! DHDHP  MW = 168
	  api3 = 0.00*mass_OC/136.0	! SVOC1  MW = 136
	  api4 = 0.00*mass_OC/272.0	! Dimer1 MW = 272
	  lim1 = 0.00*mass_OC/136.0	! SVOC2  MW = 136
	  lim2 = 0.00*mass_OC/272.0	! Dimer2 MW = 272
	  oc   = 1.00*mass_OC		! ug/m^3(air)
	  
	  co3  = 0.0
	  na   = 1.0*mass_SC/mw_SC	! umol/m^3(air)
	  ca   = 0.0
	  oin  = 0.0
	  bc   = 0.0

          write(10,400)ibin, num, Dp, sigma, jhyst, water, so4, no3, cl, &
		nh4, msa, aro1, aro2, aro3, aro4, alk1, ole1, &
		api1, api2, api3, api4, lim1, lim2, co3, na, ca, oin, oc, bc

	  write(12,500)'        aer_so4_bkg(',ibin,') = ',so4*1000.0,'   ! nmol/m3'
	  write(13,500)'        aer_nh4_bkg(',ibin,') = ',nh4*1000.0,'   ! nmol/m3'
	  write(14,500)'        aer_oc_bkg(',ibin,') = ',oc*1000.0,'   ! ng/m3'
          write(15,501)'        num_a_bkg(',ibin,') = ', num,'   ! #/cc'

	

! calculate overall gas-side mass transfer coefficient
          speed = mean_molecular_speed(T_K,MW)		! cm/s
	  z = MW/28.84
	  freepath = 32.*Dg/(3.*pi*z*speed)		! cm
          Kn  = 2.*freepath/(Dp_wet*1.e-4)		! Knudsen Number
          Fkn = fuchs_sutugin(Kn,accom)			! correction factor (-)
          GR = 2.*Dg*Fkn/(Dp_wet*1.e-4 * rho_OC)*1.e-12	! growth rate (cm/s)
	  kg  = 2.*pi*(Dp_wet*1.e-4)*Dg*num*Fkn		! CS - MOSAIC (1/s)
	  sum_kg = sum_kg + kg

	  if(Dp_wet .lt. Dp_cutoff/1000.0)sum_kg1 = sum_kg1 + kg	! Aitken CS (1/s)
	  if(Dp_wet .gt. Dp_cutoff/1000.0)sum_kg2 = sum_kg2 + kg	! Accumu CS (1/s)

	  write(11,411)ibin, Dp_wet, num, kg, GR

	enddo ! ibin loop
	
	

        write(6,*)'Mean molecular speed (cm/s) = ', speed
        write(6,*)'Freepath (nm) = ', freepath*1.e7
	write(6,*)'total wet volume = ', vol_tot, ' um^3/cc(air'
	write(6,*)'total dry mass   = ', mass_tot, ' ug/m^3(air)'
	write(6,*)'num_small = ', num_small, 'sum_kg1 = ', sum_kg1
	write(6,*)'num_large = ', num_large, 'sum_kg2 = ', sum_kg2
	write(6,*)'vol_small = ', vol_small
	write(6,*)'vol_large = ', vol_large
	write(6,*)'SOA_mass_small = ', mass_OC_small
	write(6,*)'SOA_mass_large = ', mass_OC_large
        write(6,*)'kg_large/kg_small   = ', sum_kg2/sum_kg1
        write(6,*)'vol_large/vol_small = ', vol_large/vol_small

400   format(i5, f12.5, f12.5, f5.1, i4, f8.2, 23(2x, e12.5))
411   format(i5, f12.5, f12.5, f12.5, 2x, e12.5)
500   format(a,i2,a,e12.5,a)
501   format(a,i2,a,f12.5,a)

      stop
      end





!----------------------------------------------------------
      function fuchs_sutugin(rkn,a)
      implicit none
      real fuchs_sutugin
! subr. arguments
      real rkn, a
! local variables
      real rnum, denom


      rnum  = 0.75*a*(1. + rkn)
      denom = rkn**2 + rkn + 0.283*rkn*a + 0.75*a
      fuchs_sutugin = rnum/denom

      return
      end function fuchs_sutugin

!----------------------------------------------------------
      function mean_molecular_speed(T, MW)	! in cm/s
      implicit none
      real mean_molecular_speed
! subr. arguments
      real T, MW	! T(K)

        mean_molecular_speed = 1.455e4 * sqrt(T/MW)

      return
      end function mean_molecular_speed

!----------------------------------------------------------
	function dNdlogDp(N_total,sigma_g, Dp, DpgN)
      implicit none

	real dNdlogDp
	real N_total, sigma_g, Dp, DpgN
	real arg


	arg = -((alog10(Dp/DpgN))**2) / (2.0*(alog10(sigma_g))**2)
	dNdlogDp = (N_total/(2.506628*alog10(sigma_g)) ) *exp(arg)

	end function dNdlogDp

