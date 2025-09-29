add_executable(${PROJECT_NAME} main.cpp one_dim.cpp search_result.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE numerics)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)

if(DEBUG)
    target_compile_definitions(MO_Labs PRIVATE __DEBUG__)
    message(STATUS "Debug output enabled")
endif()