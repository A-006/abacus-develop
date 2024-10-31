list(APPEND fft_files
    fft.cpp
    fft_base.cpp
    fft_cpu.cpp
    # fft_cpu_float.cpp
)

# if (USE_CUDA)
#     list(APPEND SOURCES
#             fft_gpu.cpp
#         )
# endif()

add_library(
    fft
    OBJECT
    ${fft_files}
)
target_link_libraries(fft PRIVATE -lfftw3 -lfftw3f -lfftw3l)
