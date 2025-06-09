! Idealized case conditions
! 
      Program aerinput
      implicit none

      integer jhyst, ibin, nbin_a
	real Dp_lo, Dp_up, Dp, dlogDp, dNdlogDp_tot, dVdlogDp_tot
      real factor_Dp, Dpmin, Dpmax, dlogDp_max
	real DpgN_uf, DpgN_ait, DpgN_acc
	real sig_uf, sig_ait, sig_acc, sigma
	real num_uf, num_ait, num_acc
	real num, num1, num2, vol, vol1, vol2
	real ratio_OC_AS1, ratio_OC_AS2
	real mass_SA1, mass_SA2, mass_AB1, mass_AB2, mass_AS1, mass_AS2, &
             mass_BC1, mass_BC2, mass_OC1, mass_OC2, mass_tot1, mass_tot2
	real umol_so4, umol_nh4


	real T_K, P, R, pi
	real Dg, MW, kg, Fkn, Kn
	real z, speed, accom, mean_molecular_speed, freepath
	real fuchs_sutugin, sum_kg, sum_kg1, sum_kg2

	real rho_SA, rho_AB, rho_AS, rho_BC, rho_OC, mw_SA, mw_AS, mw_AB, mw_BC, mw_OC
	real mass_SA, mass_AB, mass_AS, mass_BC, mass_OC, vol_tot, mass_tot 
      real water, so4, no3, cl, nh4, msa, aro1
	real aro2, aro3, aro4, alk1, ole1
	real api1, api2, api3, api4
	real lim1, lim2, co3, na, ca, oin, oc, bc

! function
	real dNdlogDp

	open(10, file = 'aer.inp')
	open(11, file = 'dNdlogDp.txt')


	pi = 3.141592654
	T_K = 298.15	! temperture (K)
	P   = 1.0		! pressure (atm)
	Dg  = 0.1		! gas-phase diffusivity (cm^2/s)
	MW  = 200.0		! molecular weight of diffusing species (g/mol)
	accom = 0.1

	write(11,*)'  ibin        Dp        num     dNdlogDp       kg'


! set mode parameters
	sig_uf   = 1.3
	DpgN_uf  = 0.005	! um
	num_uf   = 10000	! #/cc(air)

	sig_ait  = 1.2
	DpgN_ait = 0.02   ! um
	num_ait  = 0.0 ! #/cc(air)

	sig_acc  = 1.3
	DpgN_acc = 0.2 ! um
	num_acc  = 500.0 ! #/cc(air)

! set species parameters
	rho_SA = 1.8 ! g/cc sulfuric acid
	rho_AS = 1.8 ! g/cc ammonium sulfate
	rho_AB = 1.8 ! g/cc ammonium bisulfate
	rho_BC = 1.8 ! g/cc
	rho_OC = 1.4 ! g/cc

	mw_SA  = 98.0   ! g/mol sulfuric acid
	mw_AS  = 132.0	! g/mol ammonium sulfate
	mw_AB  = 115.0	! g/mol ammonium bisulfate
	mw_BC  = 1.0	! g/mol
	mw_OC  = 1.0	! g/mol not to be used

	ratio_OC_AS1  = 0.0 ! 1.8472	! mass_OC/mass_AS
	ratio_OC_AS2  = 0.0 ! 1.0


! set up bin structure parameters
      nbin_a  = 128 	! 110	! user specified total number of bins: 10, 30, 60, 120
      Dpmin = 0.000973	! user specified min diameter (um) = lower boundary of the smallest bin
      Dpmax = 1.0273	! user specified max diameter (um) = upper boundary of the largest bin
        
      dlogDp_max = alog10(Dpmax/Dpmin)
      dlogDp     = dlogDp_max/float(nbin_a)
      factor_Dp  = 10.**dlogDp

	write(6,*)'dlogDp = ', dlogDp

	vol_tot = 0.0
	mass_tot = 0.0
	sum_kg = 0.0
	sum_kg1 = 0.0
	sum_kg2 = 0.0

	do ibin = 1, nbin_a

	  Dp_lo = Dpmin*factor_Dp**(ibin-1)	! lower boundary diameter for bin number 'ibin'
	  Dp_up = Dpmin*factor_Dp**ibin		! upper boundary diameter for bin number 'ibin'

	  Dp = (Dp_lo + Dp_up)/2.0		! average of lower and upper boundaries

	  num1 = dlogDp * (dNdlogDp(num_uf,sig_uf, Dp, DpgN_uf) + &	! UF + Aitken modes (#/cc(air))
		             dNdlogDp(num_ait,sig_ait, Dp, DpgN_ait))

	  num2 = dlogDp * dNdlogDp(num_acc,sig_acc, Dp, DpgN_acc)		! accumulation  mode (#/cc(air))

! compute volumes
	  vol1 = num1 * 0.523599 * Dp**3	! um^3/cc(air)
	  vol2 = num2 * 0.523599 * Dp**3	! um^3/cc(air)

! set composition of UF and Aitken modes particles in ug/m3(air)
	  mass_BC1 = 0.0 ! vol1*rho_BC
	  mass_OC1 = 0.0 ! ratio_OC_AS1*mass_AS1
	  mass_SA1 = 0.0
	  mass_AB1 = 0.0
	  mass_AS1 = vol1*rho_AS ! vol1/(1.0/rho_AS + ratio_OC_AS1/rho_OC)

! set composition of accumulation mode particles in ug/m3(air)
	  mass_BC2 = 0.0
	  mass_OC2 = 0.0 ! ratio_OC_AS2*mass_AS2
	  mass_SA2 = 0.0
	  mass_AB2 = 0.0
	  mass_AS2 = vol2*rho_AS ! vol2/(1.0/rho_AS + ratio_OC_AS2/rho_OC)

! add
	  num = num1 + num2
	  vol = vol1 + vol2
	  mass_SA = mass_SA1 + mass_SA2
	  mass_AS = mass_AS1 + mass_AS2
	  mass_AB = mass_AB1 + mass_AB2
	  mass_BC = mass_BC1 + mass_BC2
	  mass_OC = mass_OC1 + mass_OC2
	  umol_so4= mass_SA/mw_SA + mass_AB/mw_AB + mass_AS/mw_AS
	  umol_nh4= mass_AB/mw_AB + 2.*mass_AS/mw_AS
	  dNdlogDp_tot = num/dlogDp
	  dVdlogDp_tot = vol/dlogDp

! total
	  vol_tot = vol_tot + vol
	  mass_tot = mass_tot + mass_BC + mass_OC + mass_AS + mass_AB + mass_SA
	  mass_tot1 = mass_tot1 + mass_BC1 + mass_OC1 + mass_AS1 + mass_AB1 + mass_SA1
	  mass_tot2 = mass_tot2 + mass_BC2 + mass_OC2 + mass_AS2 + mass_AB2 + mass_SA2


	  sigma = 1.0	! dummy never used
	  jhyst = 1
	  water = 0.0
	  so4  = umol_so4	! umol/m^3(air)
	  no3  = 0.0
	  cl   = 0.0
	  nh4  = umol_nh4	! umol/m^3(air)
	  msa  = 0.0
	  aro1 = 0.0
	  aro2 = 0.0
	  aro3 = 0.0
	  aro4 = 0.0
	  alk1 = 0.0
	  ole1 = 0.0
	  api1 = 0.0
	  api2 = 0.0
	  api3 = 0.0
	  api4 = 0.0
	  lim1 = 0.0
	  lim2 = 0.0
	  co3  = 0.0
	  na   = 0.0
	  ca   = 0.0
	  oin  = 0.0
	  oc   = mass_OC		! ug/m^3(air)
	  bc   = mass_BC		! ug/m^3(air)

        write(10,400)ibin, num, Dp*1.0, sigma, jhyst, water, so4, no3, cl, &
		nh4, msa, aro1, aro2, aro3, aro4, alk1, ole1, &
		api1, api2, api3, api4, lim1, lim2, co3, na, ca, oin, oc, bc
	  

! calculate overall gas-side mass transfer coefficient
	  Dg = 0.05	! cm^2/s
        speed = mean_molecular_speed(T_K,MW)	! cm/s
	  z = MW/28.84
	  freepath = 32.*Dg/(3.*pi*z*speed)		! cm
        Kn  = 2.*freepath/Dp				! Knudsen Number
        Fkn = fuchs_sutugin(Kn,accom)		! correction factor (-)
	  kg  = 2.*pi*Dp*Dg*num*Fkn			! CS - MOSAIC (1/s)
	  sum_kg = sum_kg + kg

	  if(Dp .lt. 0.05)sum_kg1 = sum_kg1 + kg	! Aitken CS (1/s)
	  if(Dp .gt. 0.05)sum_kg2 = sum_kg2 + kg	! Accumu CS (1/s)

	  write(11,411)ibin, Dp, num, dNdlogDp_tot, kg


	enddo

	write(6,*)'total volume = ', vol_tot, ' um^3/cc(air'
	write(6,*)'total mass   = ', mass_tot, ' ug/m^3(air)'
	write(6,*)'mass_tot1 = ', mass_tot1
	write(6,*)'mass_tot2 = ', mass_tot2
	write(6,*)'sum_kg1 = ', sum_kg1
	write(6,*)'sum_kg2 = ', sum_kg2

400   format(i5, e12.5, f12.5, f5.1, i4, f8.2, 23(2x, e12.5))
411	format(i5, f12.5, 3(2x, e12.5))

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

