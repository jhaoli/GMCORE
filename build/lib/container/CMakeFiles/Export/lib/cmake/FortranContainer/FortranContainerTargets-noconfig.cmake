#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fortran_container" for configuration ""
set_property(TARGET fortran_container APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(fortran_container PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libfortran_container.0.0.1.dylib"
  IMPORTED_SONAME_NOCONFIG "@rpath/libfortran_container.1.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS fortran_container )
list(APPEND _IMPORT_CHECK_FILES_FOR_fortran_container "${_IMPORT_PREFIX}/lib/libfortran_container.0.0.1.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
