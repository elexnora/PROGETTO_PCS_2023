# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/Main_Source/googletest"
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/Main_Build/googletest"
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix"
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix/tmp"
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix/src/googletest-stamp"
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix/src"
  "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix/src/googletest-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix/src/googletest-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "C:/Users/matteo/OneDrive/Desktop/PROGETTO/externals/googletest-prefix/src/googletest-stamp${cfgdir}") # cfgdir has leading slash
endif()
