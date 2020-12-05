cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayExample )

add_executable( sandia_decay_example sandia_decay_example.cpp )
target_link_libraries( sandia_decay_example SandiaDecay )

add_executable( batch_decay batch_decay.cpp )
target_link_libraries( batch_decay SandiaDecay )


