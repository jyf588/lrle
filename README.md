# lrle
Code for "Synthesis of Biologically Realistic Human Motion Using Joint Torque Actuation", SIGGRAPH 2019

https://arxiv.org/abs/1904.13041

## Environments:
To run data generation, you will need a Matlab environment with OpenSim and IPOPT packages installed:

1. Install the Opensim library and set up its Matlab scripting environment: (The current version of code uses OpenSim 3.3, which is only supported on Windows. Consider installing a Windows 7 virtual machine. Give the VM multiple cores of CPU.)
    1.1 Download and install OpenSim 3.3. Installation Guide: https://simtk-confluence.stanford.edu:8443/display/OpenSim33/Installing+OpenSim
    1.2 Configure Matlab Scripting Environment: https://simtk-confluence.stanford.edu:8443/display/OpenSim33/Scripting+with+Matlab 
2. Install IPOPT package in Matlab: Download this file: www.coin-or.org/download/binary/Ipopt/Ipopt-3.11.3-win32win64-matlabmexfiles.zip , unzip it and place the folder to some good place,  open Matlab, run addpath('D:/where-you-put-ipopt-folder') and then savepath.

To run Neural Net training, you will need a Python environment with keras (tensorflow), numpy, and scipy installed.

## Execution:
Train NN to approximate R:

1. Run MAIN_2D_R or MAIN_3D_R (MATLAB) to generate training data
2. Run train_leg_2D_R or train_leg_3D_R (Python) to train NN using the generated data
3. If to be used in Trajectory Optimization in Matlab, run keras_output_weight.py to translate trained model in Step 2 in Python to a Matlab readable NN model.

Train NN to approximate E:

1. Run rej_sample_R4E_2D or rej_sample_R4E_3D (Python) to sample valid torque vectors. (Note: this is different from the sampling method in the Siggraph paper. This new method is much faster and has better sample distribution. The effective number of samples remain similar.)
2.  Run MAIN_2D_E or MAIN_3D_E (MATLAB) to calculate the ground truth E values for the sampled vectors. They are the training data (where E are the labels, i.e. y_train).
3. Run train_leg_2D_E or train_leg_3D_E (Python) to train NN using the generated data
4. If to be used in Trajectory Optimization in Matlab, run keras_output_weight.py to translate trained model in Step 3 in Python to a Matlab readable NN model.

## Examples:
Some MATLAB examples of how to use learned R and E in trajectory optimization problems are included:

1. Jumping: run Jump_box_and_lr/JumperLeftRevoluteFoot_DC.m with addNN = false (BOX baseline) or addNN = true (Learned R); run Jump_MTU/JumperLeftRevoluteFootStatic_DC.m for the MTU baseline.
2. Max Swing: run Swing_box_and_lr/Swing_DC.m with addNN = false (BOX baseline) or addNN = true (Learned R); run Swing_MTU/SwingStaticMuscle_DC.m for the MTU baseline.
3. Fix-distance Swing: run Swing_fix_dist_lronly_and_lrle/Swing_DC.m with useMeta = false (Learned R only) or useMeta = true (Learned R and E).

The generated trajectories (joint angles in .sto file) can be visualized using OpenSim desktop application.

