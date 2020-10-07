# Optimal-Control-of-Arm-Reaching-movement

We provide the MATLAB (MathWorks, Natick, MA, USA) codes and data sets used in this study. Additionally, it includes video files replaying the simulations. . 

# Codes
The code ‘main.m’ computes the optimal control using ILQG, and executes the simulation of the reaching movements. Note that it requires a long computation times. The code ‘ilqg_det.m’ is the function to compute ILQG and called by ‘main.m’, and originally written by Todorov, E. and Li, W 49. The code ‘plot_result.m’ loads the data sets in the directory ‘Data’, which are precomputed results, and plots the data as replications of Figs. 3, 4, and 5. The other codes are in the directory ‘Model’. They are functions used to calculate the kinematics, Jacobian matrix, muscle dynamics, arm dynamics, and cost. 

# Data set
There are four data sets in the directory ‘Data’. The files ‘data1.mat’ to ‘data6.mat’ are the precomputed results in Cases 1–4, and the stabilizing controls of (11) and (12), respectively. 

# Videos
There are four video files in the directory ‘Video’. The files ‘video1.mp4’ to ‘video6.mp4’ replay the simulation in Cases 1–4, and the stabilizing controls of (11) and (12), respectively. 
