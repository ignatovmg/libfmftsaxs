include_directories(..)
include_directories(${CHECK_INCLUDE_DIRS})
link_directories(${CHECK_LIBRARY_DIRS})
include_directories(${MOL2_INCLUDE_DIRS})
#include_directories(${GRID2_INCLUDE_DIRS})

set(tests
  saxs_test
)

# Configure data files
file(GLOB test_files "${CMAKE_CURRENT_SOURCE_DIR}/data/*")
foreach(filepath ${test_files})
  get_filename_component(filename ${filepath} NAME)
  configure_file(${filepath} ${CMAKE_CURRENT_BINARY_DIR}/${filename} COPYONLY)
endforeach()

# Add and register tests
foreach(test ${tests})
  add_executable(${test} ${test}.c)
  target_link_libraries(${test}
    fmftsaxs
    ${MOL2_LIBRARIES}
    ${CHECK_LIBRARIES}
    fftw3f
    rt
    m)
  add_test(${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
