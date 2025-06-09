subroutine IntegrateChemistry( it,    & !intent-ins
     jaerosolstate, dp_wet_a,         & !intent-inouts
     aH2O_a, gam_ratio, iter_mesa     ) !intent-outs

  use module_data_mosaic_kind,  only: r8
  use module_data_mosaic_aero,  only: nbin_a_max, nsalt
  use module_data_mosaic_main,  only: mgas, maer, mcld
  use module_data_mosaic_boxmod,  only: cnn, told_sec, tcur_sec, &
                                      pr_atm, rh, te, cair_molm3

  use module_mosaic_aerdynam_intr,  only: aerosoldynamics

  implicit none

  !Subroutine arguments
  integer,  intent(in)    :: it
  integer,  intent(inout), dimension(nbin_a_max) :: jaerosolstate
  integer,  intent(out),   dimension(nbin_a_max) :: iter_MESA

  real(r8), intent(inout), dimension(nbin_a_max) :: dp_wet_a
  real(r8), intent(out),   dimension(nbin_a_max) :: aH2O_a, gam_ratio

  !Local Variables
  real(r8) :: cair_mol_m3
  real(r8) :: t_in, t_out

  t_in = told_sec
  t_out= tcur_sec
  cair_mol_m3 = cair_molm3

  if(mgas.eq.1)then
     call GasChemistry( t_in, t_out )
  endif

  if(maer.eq.1)then
     call aerosoldynamics( it, t_out, t_in,   & !intent-ins
          pr_atm, rh, te, cair_mol_m3,        &
          cnn, jaerosolstate, dp_wet_a,       & !intent-inouts
          aH2O_a, gam_ratio, iter_mesa        ) !intent-outs
     
  endif

  if(mcld.eq.1)then
     call CldChemistry( t_in, t_out )
  endif


  return
end subroutine IntegrateChemistry
