enable_language(Fortran)
# add any files that you would like to use to build a custom fluka
# e.g. add source.f and so on
set(FLUKA_FILES bfield.cpp magfld.f plasma_source.cpp source.f)
set(SRC_FILES mainFluDAG.cpp ${FLUKA_FILES})

set(CMAKE_Fortran_FLAGS "-g -cpp -O3 -fd-lines-as-comments -Wall -Waggregate-return -Wcast-align -Wline-truncation -Wno-conversion -Wno-integer-division -Wno-tabs -Wno-unused-dummy-argument -Wno-unused-function -Wno-unused-parameter -Wno-unusedvariable -Wsystem-headers -Wuninitialized -Wunused-label -mtune=generic -msse2 -mfpmath=sse -fPIC -fexpensive-optimizations -funroll-loops -fstrength-reduce -fno-automatic -finit-local-zero -ffixed-line-length-132 -fbackslash -funderscoring -frecord-marker=4 -falign-commons -fbacktrace -frange-check -fbounds-check -fdump-core -ftrapping-math -ffpe-trap=invalid,zero,overflow")

set(LINK_LIBS fludag)
set(LINK_LIBS_EXTERN_NAMES)

include_directories(${CMAKE_SOURCE_DIR}/src/fluka)
include_directories(${FLUKA_DIR}/include)

dagmc_install_exe(mainfludag)
