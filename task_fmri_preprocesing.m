function task_fmri_preprocesing(sub_id,task)
% fMRI PREPROCESSING
%
% The following steps are performed (in this order):
%
% 1. Slice timing correction
% 2. Realignment        - write only the mean image. The parameters are
%                         written to file header and applied when doing the
%                         normalization
% 3. Co-registration    - of functional and anatomical images
% 4. Segmentation       - of anatomical image into tissue types
% 5. Normalization      - apply realignment and segmentation deformation to
%                         normalize to standard space)
% 6. Smoothing
% 7. Normalization of anatomical image

%% Set paths and common options, and load files

% DIRECTORIES
% Folders containing functional and anatomical data are assumed to be at
% the same level, namely 'path.root', with subfolders for each subject, e.g.
% 'n01', 'n02' etc.
% -------------------------------------------------------------------------
path.root   = '/DATA/238/yyang/workspace/973_task/preprocessing';
switch task
    case 'Num'
        path.FunImg = fullfile(path.root,'FunImg_Num'); % folder containing fMRI times series
        path.T1 = fullfile(path.root,'FunImg_Num'); % folder containing fMRI times series
    case 'SPA'
        path.FunImg = fullfile(path.root,'FunImg_SPA'); % folder containing fMRI times series
        path.T1 = fullfile(path.root,'FunImg_SPA'); % folder containing fMRI times series
    otherwise
        error('No such task!!!');
end


% Get tissue probability maps from SPM. Default location (on Windows) will
% be something like 'C:\Users\username\Documents\MATLAB\spm12\tpm
path.tpm  = '/DATA/238/yyang/MatlabToolbox/spm12/tpm';
tpm = cell(1,6);
for i=1:6
    tpm(i) = cellstr(fullfile(path.tpm,sprintf('TPM.nii,%d',i)));
end

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

% slice timing correction
matlabbatch{1}.spm.temporal.st.scans    = {cellstr(FunImg)};
matlabbatch{1}.spm.temporal.st.nslices  = opt.ns;
matlabbatch{1}.spm.temporal.st.tr       = opt.TR;
matlabbatch{1}.spm.temporal.st.ta       = opt.TA;
matlabbatch{1}.spm.temporal.st.so       = opt.so;
matlabbatch{1}.spm.temporal.st.refslice = opt.ref;
matlabbatch{1}.spm.temporal.st.prefix   = 'a';
% realign: estimate and write (mean img only)
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1)       = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep     = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm    = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm     = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp  = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap    = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight  = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which   = [0 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp  = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap    = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask    = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix  = 'r';
% co-registration
matlabbatch{3}.spm.spatial.coreg.estimate.ref(1)            = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{3}.spm.spatial.coreg.estimate.source            = cellstr(T1Img);
matlabbatch{3}.spm.spatial.coreg.estimate.other             = {''};
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol      = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
% segmentation
matlabbatch{4}.spm.spatial.preproc.channel.vols(1)  = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{4}.spm.spatial.preproc.channel.biasreg  = 0.001;
matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{4}.spm.spatial.preproc.channel.write    = [0 1]; % save bias corrected
matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm    = tpm(1);
matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus  = 2; % two for grey matter
matlabbatch{4}.spm.spatial.preproc.tissue(1).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(1).warped = [1 1]; % save unmodulated probability map
matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm    = tpm(2);
matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus  = 2; % two for white matter
matlabbatch{4}.spm.spatial.preproc.tissue(2).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(2).warped = [1 1];
matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm    = tpm(3);
matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus  = 2; % two for CSF
matlabbatch{4}.spm.spatial.preproc.tissue(3).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(3).warped = [1 1];
matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm    = tpm(4);
matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus  = 3; % three for bone
matlabbatch{4}.spm.spatial.preproc.tissue(4).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(4).warped = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm    = tpm(5);
matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus  = 4; % four for other tissues
matlabbatch{4}.spm.spatial.preproc.tissue(5).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(5).warped = [1 0];
matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm    = tpm(6);
matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus  = 2; % two for air(background)
matlabbatch{4}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{4}.spm.spatial.preproc.tissue(6).warped = [1 0];
matlabbatch{4}.spm.spatial.preproc.warp.mrf         = 1;
matlabbatch{4}.spm.spatial.preproc.warp.cleanup     = 1;
matlabbatch{4}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
matlabbatch{4}.spm.spatial.preproc.warp.affreg      = opt.reg;
matlabbatch{4}.spm.spatial.preproc.warp.fwhm        = 0;
matlabbatch{4}.spm.spatial.preproc.warp.samp        = 3;
matlabbatch{4}.spm.spatial.preproc.warp.write       = [0 1];
% normalize write (functional images)
matlabbatch{5}.spm.spatial.normalise.write.subj.def(1)      = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{5}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
matlabbatch{5}.spm.spatial.normalise.write.woptions.bb      = opt.bbox;
matlabbatch{5}.spm.spatial.normalise.write.woptions.vox     = opt.vsize_fmri;
matlabbatch{5}.spm.spatial.normalise.write.woptions.interp  = 4;
% smooth (functional images)
matlabbatch{6}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{6}.spm.spatial.smooth.fwhm    = opt.FWHM;
matlabbatch{6}.spm.spatial.smooth.dtype   = 0;
matlabbatch{6}.spm.spatial.smooth.im      = 0;
matlabbatch{6}.spm.spatial.smooth.prefix  = 's';
% normalize write (anatomical image)
matlabbatch{7}.spm.spatial.normalise.write.subj.def(1)      = cfg_dep('Segment: Forward Deformations', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','fordef', '()',{':'}));
matlabbatch{7}.spm.spatial.normalise.write.subj.resample(1) = cfg_dep('Coregister: Estimate: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{7}.spm.spatial.normalise.write.woptions.bb      = opt.bbox;
matlabbatch{7}.spm.spatial.normalise.write.woptions.vox     = opt.vsize_t1;
matlabbatch{7}.spm.spatial.normalise.write.woptions.interp  = 4;
spm_jobman('run',matlabbatch);
