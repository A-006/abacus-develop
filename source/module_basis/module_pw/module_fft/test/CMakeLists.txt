remove_definitions(-D__CUDA)
remove_definitions(-D__RCOM)

AddTest(
  TARGET fft_bundle_without_fftwf_test
  LIBS parameter ${math_libs} base device FFTW3::FFTW3
  SOURCES fft_bundle_without_fftwf_test.cpp
            ../fft_bundle.cpp ../fft_cpu.cpp 
)

AddTest(
  TARGET fft_test_cpu
  LIBS parameter ${math_libs} base device FFTW3::FFTW3
  SOURCES fft_test_cpu.cpp
            ../fft_bundle.cpp ../fft_cpu.cpp 
)