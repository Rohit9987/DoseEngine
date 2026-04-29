# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-src"
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-build"
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix"
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix/tmp"
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix/src/head_model-populate-stamp"
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix/src"
  "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix/src/head_model-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix/src/head_model-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/rohit/programming/head_model/dose_engine/build/_deps/head_model-subbuild/head_model-populate-prefix/src/head_model-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
