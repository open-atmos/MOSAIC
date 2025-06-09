      subroutine SetIOfiles
      use module_data_mosaic_main
      use module_data_mosaic_gas
      use module_data_mosaic_aero
      use module_data_mosaic_cloud
      use module_data_mosaic_pmcmos, only: &
         lun_pmcmos_in1, lun_pmcmos_in2
      use module_data_mosaic_asect, only: lunout

      implicit none

      integer ibin, lun
      integer nchar, nbllen, nbb
      character*7 bb

! INPUT FILES
      write(6,*)'Enter gas input filename. Example: case.inp'
!      read(5,*)inputfile

      inputfile = 'case1.inp'

      lun_pmcmos_in1 = 16 ; lun_pmcmos_in2 = 17

      lun_sect_170 = 170
      lun_sect_171 = 171
      lun_sect_172 = 172
      lun_sect_180 = 180
      lun_sect_183 = 183
      lun_sect_184 = 184
      lun_sect_185 = 185
      lun_sect_186 = 186
      lun_sect_188 = 188
      lun_sect_190 = 190

      lunout = lun_sect_170

      lun_inp = 10
      open(lun_inp, file = inputfile)
      call ReadInputFile
      close(lun_inp)


!------------------------------------------------------------------
! OUTPUT FILES
!
        nchar = nbllen(inputfile) - 4

! GAS output file
      if(mgas .eq. mYES)then

        lun_gas = 19
        gas_output = 'output/'//inputfile(1:nchar)//'.gas.txt'
        open(lun_gas, file = gas_output)

      endif




! AEROSOL output files
      if(maer .eq. mYES)then	! UNCOMMENT THIS LINE


! number size distribution output file
        lun_Dpbins = 20
        species_output= 'output/'//inputfile(1:nchar)//'.Dpbins.txt'
        open(lun_Dpbins, file = species_output)

! number size distribution output file
        lun_numsizedist = 21
        species_output= 'output/'//inputfile(1:nchar)//'.numsizedist.txt'
        open(lun_numsizedist, file = species_output)

! volume size distribution output file
        lun_volsizedist = 22
        species_output= 'output/'//inputfile(1:nchar)//'.volsizedist.txt'
        open(lun_volsizedist, file = species_output)

! pH distribution output file
        lun_pHsizedist = 23
        species_output= 'output/'//inputfile(1:nchar)//'.pHsizedist.txt'
        open(lun_pHsizedist, file = species_output)

! nh4/so4 distribution output file
        lun_nh4so4sizedist = 24
        species_output= 'output/'//inputfile(1:nchar)//'.nh4so4sizedist.txt'
        open(lun_nh4so4sizedist, file = species_output)

! no3/so4 distribution output file
        lun_no3so4sizedist = 25
        species_output= 'output/'//inputfile(1:nchar)//'.no3so4sizedist.txt'
        open(lun_no3so4sizedist, file = species_output)

! kelvin distribution output file
        lun_kelvinsizedist = 26
        species_output= 'output/'//inputfile(1:nchar)//'.kelvinsizedist.txt'
        open(lun_kelvinsizedist, file = species_output)

! ah2o distribution output file
        lun_ah2osizedist = 27
        species_output= 'output/'//inputfile(1:nchar)//'.ah2osizedist.txt'
        open(lun_ah2osizedist, file = species_output)

! kg*satratio distribution output file
        lun_kgsatsizedist = 28
        species_output= 'output/'//inputfile(1:nchar)//'.kgsatsizedist.txt'
        open(lun_kgsatsizedist, file = species_output)

! Sat_ratio size distribution output file
        lun_satratiosizedist = 29
        species_output= 'output/'//inputfile(1:nchar)//'.satratiosizedist.txt'
        open(lun_satratiosizedist, file = species_output)

! species distribution output file
        lun_species = 30
        species_output= 'output/'//inputfile(1:nchar)//'.species.txt'
        open(lun_species, file = species_output)


! aerosol optical info output file
        lun_aeroptic = 31
        aeroptic_output=   &
                   'output/'//inputfile(1:nchar)//'.aeroptic.txt'
        if (maeroptic > 0) then
           open(lun_aeroptic, file = aeroptic_output)
        end if

      endif

      return
      end subroutine SetIOfiles





