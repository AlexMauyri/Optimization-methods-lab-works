add_library(numerics INTERFACE)

target_include_directories(numerics INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(numerics INTERFACE cxx_std_17)