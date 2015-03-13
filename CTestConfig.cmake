## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
##
## # The following are required to submit to the CDash dashboard:
##   ENABLE_TESTING()
##   INCLUDE(CTest)

execute_process(COMMAND "git" "rev-parse" "--abbrev-ref" "HEAD" ERROR_VARIABLE git_branch_error OUTPUT_VARIABLE git_branch OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
execute_process(COMMAND "git" "diff" OUTPUT_VARIABLE git_diff OUTPUT_STRIP_TRAILING_WHITESPACE WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
file(WRITE ${PROJECT_BINARY_DIR}/${git_branch}_diff.diff ${git_diff})
set(CTEST_EXTRA_SUBMIT_FILES ${PROJECT_BINARY_DIR}/${git_branch}_diff.diff)

set(MYNAME ${git_branch})
if(NOT git_diff STREQUAL "")
    set(MYNAME "${MYNAME}(unsynced)")
endif()
set(MYNAME "${MYNAME}-sequential-${CMAKE_BUILD_TYPE}")

set(BUILDNAME "${MYNAME}" CACHE STRING "build name variable for CDash")

set(CTEST_PROJECT_NAME "xfem")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "matrix.asc.tuwien.ac.at")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=xfem")
set(CTEST_DROP_SITE_CDASH TRUE)
