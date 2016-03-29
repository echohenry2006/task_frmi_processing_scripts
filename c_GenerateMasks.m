% GENERATE THE FOLLOWING MASKS
% 
%   1. whole brain
%   2. (eroded) white matter
%   3. CSF
%   4. anatomical noise ROI (combined white matter and CSF)

% SETUP
% =========================================================================
clear all

% Paths
path.root = '/DATA/238/yyang/workspace/973_task/preprocessing';
    path.mask   = fullfile(path.root,'FunImg_Num');   % path to save masks
    path.T1     = fullfile(path.root,'FunImg_Num');      % path to tissue probability maps
    path.FunImg = fullfile(path.root,'FunImg_Num'); % path to functional images
    
% Subjects
subject      = struct2cell(dir(path.FunImg))'; % folder content to cell
subject      = char(subject(:,1));             % convert to string
subject(subject(:,1)=='.',:) = [];             % find hidden folders/files (starting with '.') and delete     
num.subjects = size(subject,1);                % get # of subjects
subject      = cellstr(subject);               % make cell array

% Initialise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% GENERATE INITIAL MASKS BY THRESHOLDING TISSUE PROBABILITY MAPS
% =========================================================================
clear matlabbatch
for i=1:num.subjects
    mkdir(path.mask,subject{i}); 
    matlabbatch{1}.spm.util.imcalc.input = {
        spm_select('FPList',fullfile(path.T1,subject{i}),'^wc1.*.nii') % gray matter
        spm_select('FPList',fullfile(path.T1,subject{i}),'^wc2.*.nii') % white matter
        spm_select('FPList',fullfile(path.T1,subject{i}),'^wc3.*.nii') % CSF
        };
    matlabbatch{1}.spm.util.imcalc.output         = 'WholeBrain';
    matlabbatch{1}.spm.util.imcalc.outdir         = cellstr(fullfile(path.mask,subject{i}));
    matlabbatch{1}.spm.util.imcalc.expression     = '(i1+i2+i3)>0.2';  % cutoff
    matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
    matlabbatch{2}.spm.util.imcalc.input = cellstr(spm_select('FPList',fullfile(path.T1,subject{i}),'^wc2.*.nii'));
    matlabbatch{2}.spm.util.imcalc.output         = 'WhiteMatter';
    matlabbatch{2}.spm.util.imcalc.outdir         = cellstr(fullfile(path.mask,subject{i}));
    matlabbatch{2}.spm.util.imcalc.expression     = 'i1>0.99';         % cutoff
    matlabbatch{2}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{2}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{2}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{2}.spm.util.imcalc.options.interp = 1;
    matlabbatch{2}.spm.util.imcalc.options.dtype  = 4;
    matlabbatch{3}.spm.util.imcalc.input = cellstr(spm_select('FPList',fullfile(path.T1,subject{i}),'^wc3.*.nii'));
    matlabbatch{3}.spm.util.imcalc.output         = 'CSF';
    matlabbatch{3}.spm.util.imcalc.outdir         = cellstr(fullfile(path.mask,subject{i}));
    matlabbatch{3}.spm.util.imcalc.expression     = 'i1>0.99';         % cutoff
    matlabbatch{3}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{3}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{3}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{3}.spm.util.imcalc.options.interp = 1;
    matlabbatch{3}.spm.util.imcalc.options.dtype  = 4;    
    spm_jobman('run',matlabbatch);
end

% RESLICE MASKS TO MATCH FUNCTIONAL IMAGES
% =========================================================================
clear matlabbatch
for i=1:num.subjects
    % get images
   img = spm_select('ExtFPList',fullfile(path.FunImg,subject{i}),'^swa.nii',1);
     %img = [fullfile(path.FunImg,subject{i}),'\','dswa_4D_wds.nii'];
   
    wb  = spm_select('FPList',fullfile(path.mask,subject{i}),'^WholeBrain.nii');
    wm  = spm_select('FPList',fullfile(path.mask,subject{i}),'^WhiteMatter.nii');
    csf = spm_select('FPList',fullfile(path.mask,subject{i}),'^CSF.nii');
    % reslice
    matlabbatch{1}.spm.spatial.coreg.write.ref             = cellstr(img);
    matlabbatch{1}.spm.spatial.coreg.write.source          = cellstr(wb);
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask   = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    matlabbatch{2}.spm.spatial.coreg.write.ref             = cellstr(img);
    matlabbatch{2}.spm.spatial.coreg.write.source          = cellstr(wm);
    matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
    matlabbatch{2}.spm.spatial.coreg.write.roptions.mask   = 0;
    matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';
    matlabbatch{3}.spm.spatial.coreg.write.ref             = cellstr(img);
    matlabbatch{3}.spm.spatial.coreg.write.source          = cellstr(csf);
    matlabbatch{3}.spm.spatial.coreg.write.roptions.interp = 4;
    matlabbatch{3}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
    matlabbatch{3}.spm.spatial.coreg.write.roptions.mask   = 0;
    matlabbatch{3}.spm.spatial.coreg.write.roptions.prefix = 'r';
    spm_jobman('run',matlabbatch);
end

% ERODE WHITE MATTER MASK
% =========================================================================
for i=1:num.subjects 
    rwm.P  = spm_select('FPList',fullfile(path.mask,subject{i}),'^rWhiteMatter.nii');
    rwm.V  = spm_vol(rwm.P);
    rwm.Y  = spm_read_vols(rwm.V);
    rwm.eY = spm_erode(spm_erode(rwm.Y));    % erode twice 
    rwm.pathsplit = regexp(rwm.V.fname,'\','split');
    rwm.ofname    = rwm.pathsplit{end};
    rwm.V.fname   = fullfile(path.mask,subject{i},sprintf('e%s',rwm.ofname));
    spm_write_vol(rwm.V,rwm.eY);
end

% COMBINE CSF AND WHITE MATTER MASKS TO MAKE AN 'ANATOMICAL NOISE ROI'
% =========================================================================
clear matlabbatch
for i=1:num.subjects
    matlabbatch{1}.spm.util.imcalc.input = {
        spm_select('FPList',fullfile(path.mask,subject{i}),'^rCSF.nii')
        spm_select('FPList',fullfile(path.mask,subject{i}),'^erWhiteMatter.nii')
        };
    matlabbatch{1}.spm.util.imcalc.output         = 'AnatomicalROI';
    matlabbatch{1}.spm.util.imcalc.outdir         = cellstr(fullfile(path.mask,subject{i}));
    matlabbatch{1}.spm.util.imcalc.expression     = 'i1+i2';
    matlabbatch{1}.spm.util.imcalc.var            = struct('name', {}, 'value', {});
    matlabbatch{1}.spm.util.imcalc.options.dmtx   = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask   = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype  = 4;
    spm_jobman('run',matlabbatch); 
end

% CHECK FOR MASKS WITH NO VOXELS
% =========================================================================
fprintf('\nChecking for empty masks...')
c=0;
for i=1:num.subjects % ANATOMICAL ROI
    aROI = spm_select('FPList',fullfile(path.mask,subject{i}),'^AnatomicalROI.nii');
    Y    = spm_read_vols(spm_vol(aROI));
    if max(max(max(Y)))==0
        c=c+1;
        prob{c} = aROI;
    end
end
for i=1:num.subjects % WHOLE BRAIN
    wb = spm_select('FPList',fullfile(path.mask,subject{i}),'^rWholeBrain.nii');
    Y  = spm_read_vols(spm_vol(wb));
    if max(max(max(Y)))==0
        c=c+1;
        prob{c} = wb;
    end
end
for i=1:num.subjects % WHITE MATTER
    wm = spm_select('FPList',fullfile(path.mask,subject{i}),'^erWhiteMatter.nii');
    Y  = spm_read_vols(spm_vol(wm));
    if max(max(max(Y)))==0
        c=c+1;
        prob{c} = wm;
    end
end
for i=1:num.subjects % CSF
    csf = spm_select('FPList',fullfile(path.mask,subject{i}),'^rCSF.nii');
    Y   = spm_read_vols(spm_vol(csf));
    if max(max(max(Y)))==0
        c=c+1;
        prob{c} = csf;
    end
end
fprintf(' Done!\n\n')
if exist('prob','var')
    fprintf('There were empty masks:\n')
    for i=1:length(prob)
        fprintf('%s\n',prob{i})
    end
else
    fprintf('No empty masks found!\n')
end