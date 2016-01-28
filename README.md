# MotorTiming
Associated MATLAB code for Motor Timing paper

There are several branches of code for analyzing data associated with our Motor Timing paper.
VIs

"Stim_control_v17.vi" was used to collect pressure data while simultaneously stimulating the muscle via an external stimulator. Other sub-VIs in folder are called upon.

MATLAB files

"bandpass_filtfilt.m" is used in a lot of functions to filter data.

"fix_stim_artifacts_v2.m" is run on the output files from the above VI because the trigger channel often bleeds into the pressure channel.

"dummy_notmat_pressure_v2.m" segments the pressure file into each respiratory period

"compare_pressure_stims_different_trains_phase_aligned.m" extracts all of the stimulation effects from each of the patterns used in the VI, while subtracting off the mean unstimulated pressure waveform from its nearest neighbors, with alignment performed based on respiratory phase, as opposed to time relative to respiratory period onset.

"dprime_comparison_for_airsac_stim.m" takes the output of "compare_pressure_stims_different_trains_phase_aligned.m" and compares selected stimulation patterns using the d-prime metric.

"compare_targeted_spects_SE_PS_same_syl_contpitch.m" was used to quantify effects of muscle stimulation during song.

"compare_different_singing_stims_v2.m" uses the results of "compare_targeted_spects_SE_PS_same_syl_contpitch.m" across multiple (3) different stimulation patterns. This file must be edited to reflect any specifications of your specific syllable

"dprime_comparison_for_pitch_stim_v2" is used on the results of "compare_different_singing_stims_v2.m" to compare 2 different stimulation patterns using the d-prime metric.

