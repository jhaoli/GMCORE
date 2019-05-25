#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "fortran_datetime" for configuration ""
set_property(TARGET fortran_datetime APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(fortran_datetime PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libfortran_datetime.0.0.2.dylib"
  IMPORTED_SONAME_NOCONFIG "@rpath/libfortran_datetime.1.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS fortran_datetime )
list(APPEND _IMPORT_CHECK_FILES_FOR_fortran_datetime "${_IMPORT_PREFIX}/lib/libfortran_datetime.0.0.2.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
