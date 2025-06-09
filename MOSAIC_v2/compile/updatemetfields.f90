      subroutine UpdateMetFields( it )
! update temperature, pressure, cloud flag, etc. here

      use module_data_mosaic_main
      use module_data_mosaic_boxmod

      implicit none

      integer, intent(in) :: it

      rh_old = rh

!      if(it .gt. 80 .and. it .le. 120)then
!
!        RH = RH + 1.0                ! case 11
!        RH = RH - 1.0                ! case 12
!
!      endif

      if (it > 1 .and. rh_updn_delta /= 0.0_r8) then
         if (rh <= rh_updn_loval) then
            rh_updn_delta = abs( rh_updn_delta )
         else if (rh >= rh_updn_hival) then
            rh_updn_delta = -abs( rh_updn_delta )
         end if
         rh = rh_old + rh_updn_delta
         rh = max( rh, 0.0_r8 )
         rh = min( rh, 100.0_r8 )
      end if

      return
      end subroutine UpdateMetFields
