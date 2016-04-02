function task_fmri_PCA_Extraction(sub_id,task,mask)
% EXTRACT PRINCIPAL COMPONENTS FROM NOISE/WHOLE BRAIN MASK
%
% The following options are available:
%
% num.PCA       - # of principal components to be retained and saved
% mask          - which mask from which to extract PCA components. Options
%                 include
%                 NoiseROI   - mask of non-neuronal tissue (e.g., CSF and
%                              white matter)
%                 WholeBrain - whole brain mask
% orth          - whether or not to orthogonalize voxel time series wrt.
%                 task effects (i.e. remove task related variance).
% orth_method   - if orth==true, specifies how to orthogonalize wrt. task
%                 effects. Options include
%                 SPMReg  - generate task time course(s) based on SPMs
%                 RegReg  - use task regressors as task time courses
%                 VoxExcl - exclude voxels which correlate with either of
%                           the task regressors
%                 Additionally, for 'RegReg' and 'VoxExcl' it is necessary
%                 to specify which regressors to use for orthogonalization.
%                 This is done implicitly for 'SPMReg' by the choice of
%                 contrasts used.
% task          - task = Num/SPA
%
% Extracting PCAs from a whole brain mask requires some kind of design
% orthogonalization.

% SETUP
% =========================================================================


% Options
%task = 'Num';
num.PCA     = 5;
%mask        = 'NoiseROI'; % NoiseROI/WholeBrain
orth        = false;       % true/false
orth_method = 'SPMReg';   % SPMReg/RegReg/VoxExcl
if strcmp(orth_method,'RegReg')==1 || strcmp(orth_method,'VoxExcl')==1
    regs    = 1:3;        % which regressors to use for orthogonalization
end
if strcmp(orth_method,'VoxExcl')==1
    prefix  = 'reduced';  % prefix for reduced anatomical noise ROI
end
HPF.HParam  = 128;        % highpass filter cutoff in s
HPF.RT      = 2;          % TR in seconds

% Paths
path.root   = '/DATA/238/yyang/workspace/973_task/preprocessing_ncoreg';
path.mask   = fullfile(path.root,'spatial_masks'); % path to masks

switch task
    case 'Num'
        path.FunImg = fullfile(path.root,'FunImg_Num');   % path to functional images
    case 'SPA'
        path.FunImg = fullfile(path.root,'FunImg_SPA');   % path to functional images
end

path.spm    = fullfile(path.root,'ANALYSIS_Num\1stlevel_initial_dswa'); % path to 1st level analysis results


% Subjects
subject      = struct2cell(dir(path.FunImg))'; % folder content to cell
subject      = char(subject(:,1));             % convert to string
subject(subject(:,1)=='.',:) = [];             % find hidden folders/files (starting with '.') and delete
num.subjects = size(subject,1);                % get # of subjects
subject      = cellstr(subject);               % make cell array

% PRINCIPAL COMPONENET EXTRACTION
% =========================================================================
fprintf('\nPRINCIPAL COMPONENT EXTRACTION\n')
fprintf('=========================================================================\n')

fprintf('Processing subject %d of %d... ',sub_id,num.subjects);
% Get functional images
fname       = strcat('swafun_Num_',subject{sub_id},'.nii');
img.P  = spm_select('ExtFPList',fullfile(path.FunImg,subject{sub_id}),'^swafun.*.nii',7:228);

num.scans = size(img.P,1);
HPF.row   = (1:num.scans);

% Get mask and extract time series
switch mask
    case 'NoiseROI'
        aROI.P = spm_select('FPList',path.mask,'^rAnatomicalROI.nii');
    case 'WhiteMatter'
        aROI.P = spm_select('FPList',path.mask,'^rWhiteMatter.nii');
    case 'GrayMatter'
        aROI.P = spm_select('FPList',path.mask,'^rGrayMatter.nii');
    case 'CSF'
        aROI.P = spm_select('FPList',path.mask,'^rCSF.nii');
    case 'WholeBrain'
        aROI.P = spm_select('FPList',path.mask,'^rWholeBrain.nii');
end
aROI.V                 = spm_vol(aROI.P);
aROI.Y                 = spm_read_vols(aROI.V);
aROI.idx               = find(aROI.Y~=0);
[aROI.x,aROI.y,aROI.z] = ind2sub(size(aROI.Y),aROI.idx);
aROI.xyz               = [aROI.x aROI.y aROI.z]';
data                   = spm_get_data(img.P,aROI.xyz); % voxel time series
data = zscore(data);
data = spm_filter(HPF,data);
if orth==true % Task design orthogonalization
    switch orth_method
        case 'SPMReg'
            % get 1st level task SPMs
            spms.P{1} = spm_select('FPList',fullfile(path.spm,subject{sub_id}),'spmT_0001.nii');
            spms.P{2} = spm_select('FPList',fullfile(path.spm,subject{sub_id}),'spmT_0002.nii');
            % load and reshape functional data (time by voxels)
            img.V  = spm_vol(img.P);
            img.Y  = spm_read_vols(img.V);
            img.pY = permute(img.Y,[4 1 2 3]);
            img.rY = reshape(img.pY,num.scans,[]);
            if i==1
                num.spms   = length(spms.P);
                TaskSignal = zeros(num.scans,num.spms);
            end
            for j=1:num.spms
                spms.V  = spm_vol(spms.P(j));
                spms.Y  = spm_read_vols(spms.V{1,1});
                spms.rY = reshape(spms.Y,[],1);    % task (contrast) time course(s)
                TaskSignal(:,j) = spm_filter(HPF,img.rY*spms.rY); % generate task signal
            end
            TaskSignal = zscore(TaskSignal); % demean and variance normalize
            
            X = [ones(num.scans,1) TaskSignal];
            B = inv(X'*X)*X'*data;
            data_est = X*B;
            data_res = data - data_est;   % data with task variance removed
        case 'RegReg'
            if i==1
                load(spm_select('FPList',fullfile(path.spm,subject{sub_id}),'SPM.mat'));
                TaskSignal = spm_filter(HPF,SPM.xX.xKXs.X(:,regs));
            end
            TaskSignal = zscore(TaskSignal); % demean and variance normalize
            
            X = [ones(num.scans,1) TaskSignal];
            B = inv(X'*X)*X'*data;
            data_est = X*B;
            data_res = data - data_est;   % data with task variance removed
        case 'VoxExcl'
            if i==1
                load(spm_select('FPList',fullfile(path.spm,subject{sub_id}),'SPM.mat'));
                TaskSignal = spm_filter(HPF,SPM.xX.xKXs.X(:,regs));
            end
            [~,p]  = corr(TaskSignal,data); % correlate regressor and voxel time series
            [~,px] = find(p<0.2);           % find voxels correlating with design regressors
            px     = unique(px);
            data(:,px) = [];                    % remove time series of problematic voxels
            data_res   = data;
            % write (new) mask of included voxels
            aROI.Y(aROI.idx(px)) = 0;                        % Set value of problematic voxels to zero
            aROI.fname   = regexp(aROI.V.fname,'\','split'); % get original filename (with extension) from first frame
            aROI.V.fname = fullfile(path.mask,subject{sub_id},...
                sprintf('%s_%s',prefix,aROI.fname{end})); % add prefix to filename
            spm_write_vol(aROI.V,aROI.Y);
    end
else % No task orthogonalization
    data_res = data;
end

% Extract principal components
[~,PCA.comp_score,PCA.comp_latent] = princomp(data_res);
PCA.comp_score = bsxfun(@rdivide,PCA.comp_score,std(PCA.comp_score)); % variance normalization
PCA.comps      = PCA.comp_score(:,1:num.PCA); % components to keep
PCA.latent     = PCA.comp_latent(1:num.PCA);
PCA.latent_sum = sum(PCA.comp_latent);

% Save components as .mat file
switch mask
    case 'NoiseROI'
        SaveName = fullfile(path.FunImg,subject{sub_id},sprintf('nvr_AnaPhysioPCA_%d',num.PCA));
    case 'WhiteMatter'
        SaveName = fullfile(path.FunImg,subject{sub_id},sprintf('nvr_WMPhysioPCA_%d',num.PCA));
    case 'GrayMatter'
        SaveName = fullfile(path.FunImg,subject{sub_id},sprintf('nvr_GMPhysioPCA_%d',num.PCA));
    case 'CSF'
        SaveName = fullfile(path.FunImg,subject{sub_id},sprintf('nvr_CSFPhysioPCA_%d',num.PCA));
end
tmp      = PCA.comps;
save(SaveName,'tmp');
fprintf('Done!\n')
