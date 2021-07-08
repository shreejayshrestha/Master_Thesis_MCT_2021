# Master_Thesis_MCT_2021

This repository contains some digital files associated with my master thesis in *Music Communication and Technology*, titled, **Study of Low Frequency Impact Sound Transmission For Auralization**. It includes the audio files of original and auralized sound pressure in the receiving room (one listening position at a corner and 5 source position on the ceiling) for the three sound sources used in my master thesis. Similarly, it contains my original scripts to read and process raw data. The scripts calls many functions owned by Peter Svensson, which are not included here, if interested, please contact me or <a href="https://www.ntnu.edu/employees/peter.svensson">Peter Svennson</a>.

Overview of the audio files:
1. *2021_06_30_mb_orig.wav* Original sound pressure for excitation with medisin ball
2. *2021_06_30_mb_aur.wav* Auralized sound pressure for excitation with medisin ball
3. *2021_06_30_with_shoes_orig.wav* Original sound pressure for excitation with footstep with shoes
4. *2021_06_30_with_shoes_aur.wav* Auralized sound pressure for excitation with footstep with shoes
5. *2021_06_30_without_shoes_orig.wav* Original sound pressure for excitation with footstep without shoes
6. *2021_06_30_without_shoes_aur.wav* Auralized sound pressure for excitation with footstep without shoes

Overview of the main scripts inside *Matlab_scripts.zip* :
1.  *acc_mb_refine_rawdata_nofilt_increased_rangem.m* This script reads and processes the raw acceleration data for medisin ball
2.  *analyze_forcesig.m* This script reads force signal of medisin ball and creates equilization filter of a modal of a force signal. This script is created by Peter Svennson.
3.  *alc_floor_naturalfreq.m*  This script calculates the natural frequency of the main floor and floating floor.
4.  *mic_footstep_without_shoes_refine_rawdata.m* This script reads and processes raw data of sound pressure recordings in the receiving room at a corner listining position and five excitation position on the floor in the sending room for the excitations with footstep without shoes on the main floor and floating floor.
5.  *mic_footstep_withshoes_refine_rawdata.m* This script reads and processes raw data of sound pressure recordings in the receiving room at a corner listining position and five excitation position on the floor in the sending room  for the excitations with footstep with shoes on the main floor and floating floor.
6.  *mic_mb_refine_rawdata_nofilt.m* This script reads and processes raw data of sound pressure recordings in the receiving room at a corner listining position and five excitation position on the floor in the sending room  for the excitations with medisin ball on the main floor and floating floor.
7.  *odsum_room_modal_analysis_updated.m* This script calculates the analytical sound pressure, room modes and modal amplitude in the receiving room for the selected listening and source positions using the modal sum theory. 
8.  *script_ir.m* This script calculates the observed room modes in the receiving room based on the measurements.
9.  *script_auralization_v7.m* Finally, this script creates the auralization of the measured sound pressure in the receiving room.
