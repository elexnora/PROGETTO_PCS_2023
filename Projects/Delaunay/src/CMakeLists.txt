list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/empty_class.hpp)

list(APPEND delaunay_sources ${CMAKE_CURRENT_SOURCE_DIR}/empty_class.cpp)

list(APPEND delaunay_includes ${CMAKE_CURRENT_SOURCE_DIR})

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/sorting.hpp)



list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/convexhull.cpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/convexhull.hpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/ipotesi.cpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/ipotesi.hpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/triangoli.cpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/funzioni.hpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/funzioni.cpp)

list(APPEND delaunay_headers ${CMAKE_CURRENT_SOURCE_DIR}/triangoli.hpp)









set(delaunay_sources ${delaunay_sources} PARENT_SCOPE)
set(delaunay_headers ${delaunay_headers} PARENT_SCOPE)
set(delaunay_includes ${delaunay_includes} PARENT_SCOPE)
