module enzo_io
#ifdef ENZO

  use grafic_types
  use hdf5
  use hdf_io

#endif
contains
#ifdef ENZO

  subroutine enzo_filestat(filename, newname, filepart, n_filepart)
    character(128), intent(in)  :: filename
    character(*), intent(out) :: newname
    integer, intent(out) :: filepart, n_filepart
    character(128) :: tempstr
    character :: dimchar

    dimchar = filename(8:8)
    if (filename == "ic_deltab") then
       newname = "GridDensity"
       n_filepart = 1
       filepart = 0
    else if (filename(:7) == "ic_velb") then
#ifdef ONE_DIMENSION
       newname = "GridVelocities" // "_" // dimchar
       n_filepart = 1
       filepart = 0
#else
       newname = "GridVelocities"
       n_filepart = 3
       filepart = ichar(dimchar)-120
#endif
    else if (filename(:7) == "ic_velc") then
#ifdef ONE_DIMENSION
       newname = "ParticleVelocities" // "_" // dimchar
       n_filepart = 1
       filepart = 0
#else
       newname = "ParticleVelocities"
       n_filepart = 3
       filepart = ichar(dimchar)-120
#endif
    else ! white noise
       newname = trim(filename)
       n_filepart = 1
       filepart = 0
    endif
    newname = trim(newname)
  end subroutine enzo_filestat

#endif
end module enzo_io
