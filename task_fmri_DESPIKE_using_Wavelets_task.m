function task_fmri_DESPIKE_using_Wavelets_task(sub_id,task)
% DESPIKE USING WAVELETS
% 
% To use this function:
% 
% 1. cd to the BrainWavelet folder and run setup.m from the command window
%    (MEX compiler required)
% 2. run WaveletDespike as
%    'WaveletDespike('\path\to\4D_file.nii(.gz)','prefix_for_output_file');
%    which then creates
%    - prefix_wds.nii.gz    - the wavelet despiked time series
%    - prefix_noise.nii.gz  - the noise time seies removed during despike
%    - prefix_SP.txt        - the spike percentage time series

% Set path to functional images
clear all

root = '/DATA/238/yyang/workspace/973_task/preprocessing_ncoreg';
pathIMG  = fullfile(root,'FunImg_Num');
    

    
% GET LIST AND # OF SUBJECTS
% Directories of normal controls are assumed to start with 'n' (patients
% could be 'p')
% -------------------------------------------------------------------------

% subjects
num.chars    = 2;                               % # of characters to consider
subject      = struct2cell(dir(pathIMG))';  % list folder content
subject      = char(subject(:,1));              % convert to string
subject(subject(:,1)=='.',:) = [];              % find hidden folders/files (starting with '.') and delete
num.subjects = size(subject,1);                 % # of subjects
subject      = cellstr(subject);                % make cell array (for convenience)
%
% 'LimitRAM' can be omitted and WaveletDespike should be able to figure out
% the amount of RAM available. If not, set the limit manually.

for i=1:num.subjects
    cd(fullfile(pathIMG,subject{sub_id}));
	
    WaveletDespike(strcat('swafun_Num_',subject{sub_id},'.nii'),'dswa_4D'); 
end
