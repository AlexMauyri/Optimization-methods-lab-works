add_executable(${PROJECT_NAME} main.cpp one_dim.cpp search_result.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE numerics)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)