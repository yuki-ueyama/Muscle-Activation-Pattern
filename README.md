# Optimal-Control-of-Arm-Reaching-movement

We provide the MATLAB (MathWorks, Natick, MA, USA) codes and data sets used in this study. Additionally, it includes video files replaying the simulations. . 

# Codes
The code ‘main.m’ computes the optimal control using ILQG, and executes the simulation of the reaching movements. Note that it requires a long computation times. The code ‘ilqg_det.m’ is the function to compute ILQG and called by ‘main.m’, and originally written by Todorov, E. and Li, W 2007. The code ‘plot_result.m’ loads the data sets in the directory ‘Data’, which are precomputed results, and plots the data as replications of Figs. 3, 4, and 5. The other codes are in the directory ‘Model’. They are functions used to calculate the kinematics, Jacobian matrix, muscle dynamics, arm dynamics, and cost. 

# Data set
There are four data sets in the directory ‘Data’. The files are the precomputed results of different cost functions. 
data_p32.mat       : Case 1;
data_pv32.mat      : Case 2;
data_pf32.mat      : Case 3;
data_pvf32.mat     : Case 4;
data_stab_p32.mat  : Case 5;
data_stab_pv32.mat : Case 6.

# Videos
There are four video files in the directory ‘Video’. The files replay the simulations. 
video_p.mp4       : Case 1;
video_pv.mp4      : Case 2;
video_pf.mp4      : Case 3;
video_pvf.mp4     : Case 4;
video_stab_p.mp4  : Case 5;
video_stab_pv.mp4 : Case 6.
