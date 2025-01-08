remove_definitions(-D__CUDA)
remove_definitions(-D__RCOM)

add_compile_definitions(__ENABLE_FLOAT_FFTW)
find_package(FFTW3 REQUIRED)
AddTest(
  TARGET fft_cpu_test
  LIBS parameter ${math_libs} base device FFTW3::FFTW3 FFTW3::FFTW3_FLOAT
  SOURCES fft_test_cpu.cpp
            ../fft_bundle.cpp ../fft_cpu.cpp ../fft_cpu_float.cpp
)
AddTest(
  TARGET fft_bundle_test
  LIBS parameter ${math_libs} base device FFTW3::FFTW3 FFTW3::FFTW3_FLOAT
  SOURCES fft_bundle_test.cpp
            ../fft_bundle.cpp ../fft_cpu.cpp ../fft_cpu_float.cpp
)