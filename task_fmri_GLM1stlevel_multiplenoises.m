function task_fmri_GLM1stlevel_multiplenoises(sub_id,task)
% Initial 1st level analysis using only (expanded) realignment parameters
% as nuisance variable regressors. Used to generate task regressors and
% obtain estimates of task activation (i.e. statistical parametric maps)
% which is used in the PCA extraction step.

% SETUP
% =========================================================================


% Paths
path.root    = '/DATA/238/yyang/workspace/973_task/preprocessing_ncoreg';
switch task
    case 'Num'
        path.FunImg  = fullfile(path.root,'FunImg_Num'); % path to functional images
        path.spm     = fullfile(path.root,'Analysis_Num/1stlevel_swa'); % output folder
    case 'SPA'
        path.FunImg  = fullfile(path.root,'FunImg_SPA'); % path to functional images
        path.spm     = fullfile(path.root,'Analysis_SPA/1stlevel_swa'); % output folder
end
path.nvrmask = fullfile(path.root,'spatial_masks'); % path to prespecified brain mask (if used)
nvrprefix = 'nvr'; % same prefix for all NVRs!
fname = '^swafun.*.nii';
% Options
opt.PreBrainMask = 1;           % whether or not to use a previously
% generated brain mask (e.g., from
% segmentation)
opt.unit    = 'scans';          % timing units. 'scans' or 'secs'
opt.TR      = 2;                % repetition time (in seconditions)
opt.fmri_t  = 16;               % microtime resolution (default = 16). Timebins per scans when building regressors
opt.fmri_t0 = 8;                % microtime onset (default = 8). Reference timebin

opt.HighPassFilter = 128;       % highpass filter. 'Inf' to disable
opt.autocor        = 'AR(1)';   % serial correlations in time series: 'AR(1)' or 'none'

% Conditions
condition(1).name  = '0-Back';
condition(1).onset = [6; 42; 60; 132; 150; 204];
condition(1).dur   = 12;
condition(2).name  = '2-Back';
condition(2).onset = [24; 78; 96; 114; 168; 186];
condition(2).dur   = 12;
condition(3).name  = 'Cue';
condition(3).onset = [5; 23; 41; 59; 77; 95; 113; 131; 149; 167; 185; 203];
condition(3).dur   = 1;

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

% 1st level analysis

% Get functional images
FunImg = spm_select('ExtFPList',fullfile(path.FunImg,subject{sub_id}),fname);
% If chosen, get brain mask
if opt.PreBrainMask==1;         % if true, load brain mask
    path.BrainMask = fullfile(path.nvrmask,subject{sub_id},'rWholeBrain.nii');
end
% Get nuisance variable regressors
nvr_size1 = 0;
nvr_files1 = spm_select('FPList',fullfile(path.FunImg,subject{sub_id}),sprintf('^%s_rp6.*.mat',nvrprefix));
for h1=1:size(nvr_files1,1)
    t1 = struct2cell(load(nvr_files1(h1,:)));
    nvr1(:,nvr_size1+1:(nvr_size1+size(t1{:},2))) = t1{:}; % one nvr matrix
    nvr_size1 = size(nvr1,2);
end
nvr_size2 = 0;
nvr_files2 = spm_select('FPList',fullfile(path.FunImg,subject{sub_id}),sprintf('^%s_AnaPhysioPCA_5.*.mat',nvrprefix));
for h2=1:size(nvr_files2,1)
    t2 = struct2cell(load(nvr_files2(h2,:)));
    nvr2(:,nvr_size2+1:(nvr_size2+size(t2{:},2))) = t2{:}; % one nvr matrix
    nvr_size2 = size(nvr2,2);
end
nvr=cat(2,nvr1,nvr2);
nvr_size = size(nvr,2);
clear matlabbatch

%       % Get nuisance variable regressors
%     clear nvr
%     nvr_size = 0;
%     nvr_files = spm_select('FPList',fullfile(path.FunImg,subject{sub_id}),sprintf('^%s_.*.mat',nvrprefix));
%     for h=1:size(nvr_files,1)
%         t = struct2cell(load(nvr_files(h,:)));
%         nvr(:,nvr_size+1:(nvr_size+size(t{:},2))) = t{:}; % one nvr matrix
%         nvr_size = size(nvr,2);
%     end
%     clear matlabbatch

% Specification
% Experimental parameters and design
matlabbatch{1}.spm.stats.fmri_spec.dir            = cellstr(fullfile(path.spm,subject{sub_id})); % output directory
matlabbatch{1}.spm.stats.fmri_spec.timing.units   = opt.unit;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = opt.TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = opt.fmri_t;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = opt.fmri_t0;
matlabbatch{1}.spm.stats.fmri_spec.sess.scans     = cellstr(FunImg); % input scans
% Conditions
for j=1:length(condition)
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(j).name     = condition(j).name;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(j).onset    = condition(j).onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(j).duration = condition(j).dur;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(j).tmod     = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(j).pmod     = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(j).orth     = 1;
end
% Nuisance variable regressors
for k=1:nvr_size
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(k).name = sprintf('R%d',k);                    % name/prefix 'R'
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(k).val  = nvr(:,k);
end
% Highpass filtering, AR(1) modeling etc.
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = opt.HighPassFilter;                     % highpass filter
matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
if opt.PreBrainMask==1
    matlabbatch{1}.spm.stats.fmri_spec.mthresh      = -Inf; % includes all voxels (also those with negative values)
    matlabbatch{1}.spm.stats.fmri_spec.mask         = cellstr(path.BrainMask);
else
    matlabbatch{1}.spm.stats.fmri_spec.mthresh      = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask         = {''};
end
matlabbatch{1}.spm.stats.fmri_spec.cvi              = opt.autocor;      % AR(1) model
spm_jobman('run',matlabbatch);

% NVRs from SPM.xX.iC (covariates) to SPM.xX.iG (nuisance variables)
SPM_file = spm_select('FPList',fullfile(path.spm,subject{sub_id}),'SPM.mat'); % load design matrix
load(SPM_file);
SPM = nuis2iG(SPM);   % move nuisance regressors from iC to iG
save(SPM_file,'SPM'); % save SPM.mat with changes

% Estimation
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat           = cellstr(SPM_file);
matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;

% t contrasts
matlabbatch{2}.spm.stats.con.spmmat(1)               = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.con.consess{1}.tcon.name    = sprintf('Positive Effect of %s',condition(1).name);
matlabbatch{2}.spm.stats.con.consess{1}.tcon.weights = [1 0 0];
matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.name    = sprintf('Positive Effect of %s',condition(2).name);
matlabbatch{2}.spm.stats.con.consess{2}.tcon.weights = [0 1 0];
matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{3}.tcon.name    = sprintf('Positive Effect of %s',condition(3).name);
matlabbatch{2}.spm.stats.con.consess{3}.tcon.weights = [0 0 1];
matlabbatch{2}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{4}.tcon.name    = sprintf('%s > %s',condition(1).name,condition(2).name);
matlabbatch{2}.spm.stats.con.consess{4}.tcon.weights = [1 -1 0];
matlabbatch{2}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{5}.tcon.name    = sprintf('%s > %s',condition(2).name,condition(1).name);
matlabbatch{2}.spm.stats.con.consess{5}.tcon.weights = [-1 1 0];
matlabbatch{2}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
% F contrasts
matlabbatch{2}.spm.stats.con.consess{6}.fcon.name    = 'Main Effect of Task';
matlabbatch{2}.spm.stats.con.consess{6}.fcon.weights = eye(2);
matlabbatch{2}.spm.stats.con.consess{6}.fcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{7}.fcon.name    = 'Effects of Interest';
matlabbatch{2}.spm.stats.con.consess{7}.fcon.weights = eye(length(condition));
matlabbatch{2}.spm.stats.con.consess{7}.fcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.delete                  = 0;
% RUN
spm_jobman('run',matlabbatch);

