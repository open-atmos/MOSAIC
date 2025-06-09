!**********************************************************************************
! This computer software was prepared by Battelle Memorial Institute, hereinafter
! the Contractor, under Contract No. DE-AC05-76RL0 1830 with the Department of
! Energy (DOE). NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! Copyright (c) 2015 Battelle Memorial Institute
! Lead Author: Rahul A. Zaveri
!
!********************************************************************************************
! Model         : MOSAIC (Model for Simulating Aerosol Interactions & Chemistry)
!
! Last Update   : September 2015
!
! Purpose       : Simulates atmospheric aerosol formation, chemistry, dynamics, thermodynamics,
!		  and its climate-relevant properties
!
! Lead Author   : Rahul A. Zaveri, PhD
!                 Senior Research Scientist
!                 Pacific Northwest National Laboratory
!                 Atmospheric Sciences & Global Change Division
!                 1100 Dexter Ave N, Suite 400
!		  Seattle, WA 98109 USA
!                 Phone: (206) 528-3215
!                 Email: Rahul.Zaveri@pnnl.gov
!
! Terms of Use  : (1) MOSAIC and its submodules may not be included in any commercial 
!                     package, or used for any commercial applications without prior 
!                     authorization from the author.
!                 (2) The MOSAIC code may be used for educational or non-profit purposes
!                     only. Any other usage must be first approved by the author.
!                 (3) The MOSAIC code cannot be modified in any way or form or distributed
!                     without the author's prior consent.
!                 (4) No portion of the MOSAIC source code can be used in other codes
!                     without the author's prior consent.
!                 (5) The MOSAIC code is provided on an as-is basis, and the author
!                     bears no liability from its usage.
!                 (6) Publications resulting from the usage of MOSAIC must cite
!                     the references below for proper acknowledgment.
!
! Reference     : Zaveri R.A., R.C. Easter, J.D. Fast, and L.K. Peters, Model
!                 for simulating aerosol interactions and chemistry (MOSAIC),
!                 J. Geophys. Res., vol. 113, D13204, doi:10.1029/2007JD008782, 2008.
!
! Support       : Funding for the development and evaluation of MOSAIC and
!                 its sub-modules was provided by:
!                 (a) The U.S. Department of Energy (DOE) under the auspices of the
!                     Atmospheric System Research (ASR) program of the Office of Biological and
!                     Environmental Research
!                 (b) The U.S. Department of Energy (DOE) under the auspices of the
!                     Atmospheric Science Program (ASP) of the Office of Biological and
!                     Environmental Research
!                 (c) the NASA Aerosol Program and NASA Earth Science Enterprise
!                 (d) the U.S. Environmental Protection Agency (EPA) Aerosol Program
!                 (e) PNNL Laboratory Directed Research and Development (LDRD) Program
!
!
! Bugs/Problems : Please report any bugs or problems to Rahul.Zaveri@pnnl.gov
!------------------------------------------------------------------------

program main
  use module_data_mosaic_kind,  only: r8
  use module_data_mosaic_boxmod,  only: mmode
  use module_data_mosaic_aero,  only: nbin_a_max

  use module_mosaic_init_aerosol,  only: init_aerosol

  implicit none

  integer, allocatable, dimension(:) :: iter_MESA

  real(r8), allocatable,dimension(:) ::  aH2O_a,gam_ratio

  write(6,*)'   '
  write(6,*)'*****************************************************'
  write(6,*)'                     MOSAIC'
  write(6,*)'Model for Simulating Aerosol Interactions & Chemistry'
  write(6,*)'    Copyright (c) 2015 Battelle Memorial Institute'
  write(6,*)'  '
  write(6,*)'   Contact: Rahul A. Zaveri (rahul.zaveri@pnnl.gov)'
  write(6,*)'       Pacific Northwest National Laboratory'
  write(6,*)'*****************************************************'
  write(6,*)'   '
  write(6,*)'   '
  write(6,*)'simulation begins...'



  call init_data_modules   ! initializes various indices
  call SetIOfiles          ! reads inputfile

  call SetRunParameters
  call SetAirComposition
  call init_aerosol

  call LoadPeroxyParameters                     ! Aperox and Bperox

  !!      call DoMassBalance                    ! initial elemental mass balance


  if(.not.allocated(iter_MESA)) allocate(iter_MESA( nbin_a_max ))
  if(.not.allocated(aH2O_a))    allocate(aH2O_a(    nbin_a_max ))
  if(.not.allocated(gam_ratio)) allocate(gam_ratio( nbin_a_max ))

  if(mmode .eq. 1)then
     call time_integration_mode(gam_ratio, iter_mesa, aH2O_a)
  elseif(mmode .eq. 2)then
     call parametric_analysis_mode(gam_ratio, iter_mesa, aH2O_a)
  endif


  write(6,*)'   '
  write(6,*)'simulation complete.'
  write(6,*)'   '

end program main


!---------------------------------------------------------------------


