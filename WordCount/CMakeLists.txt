cmake_minimum_required(VERSION 3.15)
project(WordCount)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
set(CMAKE_CXX_FLAGS -lstdc++fs)
find_package(Threads REQUIRED)

add_executable(WordCount main.cpp)
add_executable(Secuencial secuencial.cpp)
target_link_libraries(WordCount Threads::Threads)
target_link_libraries(WordCount stdc++fs)
