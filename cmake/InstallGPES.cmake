set(GPES_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/gpes")

install(TARGETS gpes EXPORT gpesTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
if(GPES_BUNDLED_SPDLOG)
  install(TARGETS gpes_logging EXPORT gpesTargets)
  install(DIRECTORY
    "${PROJECT_SOURCE_DIR}/include/gpes/third_party/spdlog"
    "${PROJECT_SOURCE_DIR}/include/gpes/third_party/fmt"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gpes/third_party")
endif()

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/gpes/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gpes"
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h" PATTERN "*.inl"
  PATTERN "source" EXCLUDE PATTERN "third_party" EXCLUDE)
# Install the Spectra headers actually used to build GPES, including any override.
install(DIRECTORY "${GPES_SPECTRA_INCLUDE_DIR}/Spectra"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gpes/third_party")
install(FILES "${PROJECT_SOURCE_DIR}/include/gpes/third_party/matplotlibcpp.h"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gpes/third_party")
install(EXPORT gpesTargets FILE gpesTargets.cmake NAMESPACE gpes::
  DESTINATION "${GPES_INSTALL_CMAKEDIR}")

if(GPES_ENABLE_PLOTTING)
  install(TARGETS gpes_plot EXPORT gpesPlotTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}")
  install(EXPORT gpesPlotTargets FILE gpesPlotTargets.cmake NAMESPACE gpes::
    DESTINATION "${GPES_INSTALL_CMAKEDIR}")
endif()

configure_package_config_file(
  "${CMAKE_CURRENT_LIST_DIR}/gpesConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/gpesConfig.cmake"
  INSTALL_DESTINATION "${GPES_INSTALL_CMAKEDIR}")
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/gpesConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}" COMPATIBILITY SameMajorVersion)
install(FILES
  "${PROJECT_BINARY_DIR}/gpesConfig.cmake"
  "${PROJECT_BINARY_DIR}/gpesConfigVersion.cmake"
  "${CMAKE_CURRENT_LIST_DIR}/FindGPESEigen3.cmake"
  "${CMAKE_CURRENT_LIST_DIR}/FindGPESFFTW.cmake"
  DESTINATION "${GPES_INSTALL_CMAKEDIR}")
