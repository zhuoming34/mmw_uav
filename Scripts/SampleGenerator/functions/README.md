# Functions used in Sample Generator

## remove_occlution_radAcam.m
- Modified from remove_occlution.m, added similar functions for camera.
- Selects points for camera reflectors
- Selects visible points from radar perspective

## model_point_reflecor.m
- Further generates radar reflectors from visible points

## simulate_radar_signal.m
- Simulates radar signal based on radar configuration in variable_library_radar.m

## pc2i.m
  Generates depth image from physical point cloud
