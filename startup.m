% startup.m — run this once at the start of every MATLAB session.
% Adds project subfolders to the path so classes and functions are visible.
projectRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(projectRoot, 'src'));
addpath(fullfile(projectRoot, 'lib'));
addpath(fullfile(projectRoot, 'tests'));
disp('STE_project paths added.');