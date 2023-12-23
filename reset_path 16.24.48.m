function reset_path()
% This script adds Code, Code/Common and its subfolders, and Data/Examples
% to the MATLAB path. Once the model has been selected in main.m, the
% relevant files for the chosen model and method are added.

% Set the current folder to 4YP
cd('/Users/dainjung/Year 4/4YP/BatEst');
codepath = genpath('./Code/');
addpath(codepath);
rmpath(codepath);

addpath('./Code/');
addpath(genpath('./Code/Common/'));
addpath('./Data/Examples/');

<<<<<<< HEAD
addpath('./Code/Subfunctions/');

addpath('./Classes/');
addpath(genpath('./Classes/'));

=======
>>>>>>> parent of b2b681e (refactor main_multi, add backup file for select_data_subset)
cd('..');
projectpath = genpath('./Project');
end
