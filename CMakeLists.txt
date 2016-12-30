link_directories(${SUPELEC_AUDIO_LIBRARIES} ${FFTW3_LIBRARIES})
include_directories(${SUPELEC_AUDIO_INCLUDE_DIRS} ${FFTW3_INCLUDE_DIRS})

set(SOURCE_FILES melfilter.cpp mfcc.cpp dtw.cpp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -funroll-loops -Wall -O3 ${Boost_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_BUILD_TYPE Release)

add_executable(testing test.cpp ${SOURCE_FILES})
target_link_libraries(testing ${FFTW3_LIBRARIES})
target_compile_options(testing PUBLIC ${FFTW3_CFLAGS_OTHER})

add_executable(parole main.cpp ${SOURCE_FILES})
target_link_libraries(parole ${catkin_LIBRARIES} ${SUPELEC_AUDIO_LIBRARIES} ${SNDFILE_LIBRARIES} ${FFTW3_LIBRARIES})
target_compile_options(parole PUBLIC ${SUPELEC_AUDIO_CFLAGS_OTHER} ${FFTW3_CFLAGS_OTHER} ${SNDFILE_CFLAGS_OTHER} "-D_SUPELEC_AUDIO_SHORT_NAMES_")
