# Append sources from this directory to the cam_sources list.
list(APPEND cam_sources coords_1d.F90 linear_1d_operators.F90)

sourcelist_to_parent(cam_sources)
