add_executable("parabolic" parabolic.cc bcextension.hh  bctype.hh driver.hh space_operator.hh  time_operator.hh)
target_link_dune_default_libraries("parabolic")

