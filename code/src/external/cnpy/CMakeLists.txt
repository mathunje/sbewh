
find_package(ZLIB REQUIRED)
target_sources(sbewh PRIVATE
    cnpy.cpp
)
target_link_libraries(sbewh PRIVATE
    ${ZLIB_LIBRARIES}
)
# target_include_directories()
