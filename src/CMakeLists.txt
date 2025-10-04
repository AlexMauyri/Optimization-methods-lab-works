add_executable(${PROJECT_NAME} main.cpp one_dim.cpp numerics.cpp multi_dim.cpp search_result_nd.cpp search_result.cpp)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)

if(DEBUG)
    target_compile_definitions(MO_Labs PRIVATE __DEBUG__)
    message(STATUS "Debug output enabled")
endif()