function task_fmri_preprocesing(sub_id,task)
% fMRI PREPROCESSING
%
% The following steps are performed (in this order):
%
% 1. Slice timing correction
% 2. Realignment        - write only the mean image. The parameters are
%                         written to file header and applied when doing the
%                         normalization
% 3. Normalization      - apply realignment and segmentation deformation to
%                         normalize to standard space)
% 4. Smoothing


%% Set paths and common options, and load files

% DIRECTORIES
% Folders containing functional and anatomical data are assumed to be at
% the same level, namely 'path.root', with subfolders for each subject, e.g.
% 'n01', 'n02' etc.
% -------------------------------------------------------------------------
path.root   = '/DATA/238/yyang/workspace/973_task/preprocessing_ncoreg';
switch task
    case 'Num'
        path.FunImg = fullfile(path.root,'FunImg_Num'); % folder containing fMRI times series
		path.T1     = fullfile(path.root,'FunImg_Num');      % folder containing T1 images
    case 'SPA'
        path.FunImg = fullfile(path.root,'FunImg_SPA'); % folder containing fMRI times series
		path.T1     = fullfile(path.root,'FunImg_Num');      % folder containing T1 images
    otherwise
        error('No such task!!!');
end

%path.T1     = fullfile(path.root,'T1Img');      % folder containing T1 images

% Get tissue probability maps from SPM. Default location (on Windows) will
% be something like 'C:\Users\username\Documents\MATLAB\spm12\tpm
path.tpm  = '/DATA/238/yyang/MatlabToolbox/spm12/tpm';


% OPTIONS
% See the full batch for detailed control of parameters.
% -------------------------------------------------------------------------
% slice timing
opt.ns  = 32;                     % no. of slices
opt.TR  = 2;                      % repetition time
opt.TA  = opt.TR-(opt.TR/opt.ns); % acquisition time, usually TR-(TR/nslices)
opt.so  = [1:2:31,2:2:32];        % slice order, e.g., [1:1:32]
opt.ref = 2;                      % reference slice

% segmentation
opt.reg = 'eastern'; % affine regularization during segmentation:
%  - 'mni'     = european brains
%  - 'eastern' = east asian brains
% normalization
opt.bbox = [-78 -112 -70;  % bounding box, [-78 -112 -70; 78 76 85] is SPM
    78   76  85]; % standard
opt.vsize_fmri = [3 3 3]; % voxel size for normalized functional images
opt.vsize_t1   = [2 2 2]; % voxel size for normalized anatomical image
opt.tpm = cellstr(fullfile(path.tpm,'TPM.nii'));
% smoothing
opt.FWHM = [6 6 6]; % smoothing kernel specified as [Sx Sy Sz]

% Subjects
num.chars    = 2;                               % # of characters to consider
subject      = struct2cell(dir(path.FunImg))';  % list folder content
subject      = char(subject(:,1));              % convert to string
subject(subject(:,1)=='.',:) = [];              % find hidden folders/files (starting with '.') and delete
num.subjects = size(subject,1);                 % # of subjects
subject      = cellstr(subject);                % make cell array (for convenience)

% Initialize SPM

spm('Defaults','fMRI');
spm_jobman('initcfg');

% PREPROCESSING
clear matlabbatch

% Load functional and structural images
FunImg = spm_select('ExtFPList',fullfile(path.FunImg,subject{sub_id}),'^fun.*.nii',7:228);
T1Img  = spm_select('FPList',fullfile(path.T1,subject{sub_id}),'^co.*.nii');
matlabbatch{1}.spm.temporal.st.scans = {cellstr(FunImg)}';
% slice timing
matlabbatch{1}.spm.temporal.st.nslices = opt.ns;
matlabbatch{1}.spm.temporal.st.tr = opt.TR;
matlabbatch{1}.spm.temporal.st.ta = opt.TA;
matlabbatch{1}.spm.temporal.st.so = opt.so;
matlabbatch{1}.spm.temporal.st.refslice = opt.ref;
matlabbatch{1}.spm.temporal.st.prefix = 'a';
% Realignment
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [0 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
% normalize
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.vol(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{3}.spm.spatial.normalise.estwrite.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.tpm = opt.tpm;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{3}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.bb = opt.bbox;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.vox = opt.vsize_fmri;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{3}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
% smoothing
matlabbatch{4}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{4}.spm.spatial.smooth.fwhm = opt.FWHM;
matlabbatch{4}.spm.spatial.smooth.dtype = 0;
matlabbatch{4}.spm.spatial.smooth.im = 0;
matlabbatch{4}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);
