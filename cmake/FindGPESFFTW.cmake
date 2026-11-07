if(TARGET FFTW3::fftw3)
  set(GPESFFTW_FOUND TRUE)
  return()
endif()

# Some cluster FFTW config packages are incomplete; discover the library directly.
find_path(FFTW3_INCLUDE_DIR NAMES fftw3.h
  HINTS "${FFTW3_ROOT}" "$ENV{FFTW3_ROOT}" "$ENV{FFTW_ROOT}" "$ENV{EBROOTFFTW}"
  PATH_SUFFIXES include
  PATHS /opt/homebrew/include /usr/local/include /usr/include)
find_library(FFTW3_LIBRARY NAMES fftw3
  HINTS "${FFTW3_ROOT}" "$ENV{FFTW3_ROOT}" "$ENV{FFTW_ROOT}" "$ENV{EBROOTFFTW}"
  PATH_SUFFIXES lib lib64
  PATHS /opt/homebrew/lib /usr/local/lib /usr/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPESFFTW
  REQUIRED_VARS FFTW3_INCLUDE_DIR FFTW3_LIBRARY
  REASON_FAILURE_MESSAGE "Load FFTW or set FFTW3_ROOT, FFTW3_INCLUDE_DIR, and FFTW3_LIBRARY.")
if(GPESFFTW_FOUND)
  add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
  set_target_properties(FFTW3::fftw3 PROPERTIES
    IMPORTED_LOCATION "${FFTW3_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIR}")
endif()
mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)
