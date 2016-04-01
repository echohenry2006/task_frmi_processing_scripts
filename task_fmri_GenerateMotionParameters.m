function task_fmri_GenerateMotionParameters(task,method)
% PROCESS REALIGNMENT PARAMETERS
% 
% Choose method. The following options are available
% 
% original  - the original 6 realignment parameters (in effect, just
%             converting from .txt to .mat file)
% friston24 - 24 polynomial and autoregressive expansion of the 6 
%             realignment parameters
% derivative12 - 12 parameters,Use the current time point and the previous
%              time point of rigid-body 6 realign parameters. 
%              e.g., Txi, Tyi, Tzi,..., Txi-1, Tyi-1, Tzi-1...



%method = 'derivative12'; % original/friston24/derivative12

% Paths
path.root   = '/DATA/238/yyang/workspace/973_task/preprocessing_ncoreg';
switch task
case 'Num'
path.FunImg = fullfile(path.root,'FunImg_Num');
case 'SPA'
path.FunImg = fullfile(path.root,'FunImg_SPA');
end

disp(task)
disp(method)

% Subjects
subject      = struct2cell(dir(path.FunImg))'; % folder content to cell
subject      = char(subject(:,1));             % convert to string
subject(subject(:,1)=='.',:) = [];             % find hidden folders/files (starting with '.') and delete     
num.subjects = size(subject,1);                % get # of subjects
subject      = cellstr(subject);               % make cell array

switch method
    case 'original'
        for i=1:num.subjects
            rp6 = load(spm_select('FPList',fullfile(path.FunImg,subject{i}),'^rp_.*.txt')); % load realignment parameters                         % 6 parameters
%             rp6 = bsxfun(@minus,rp6,mean(rp6));                    % demean
%             rp6 = bsxfun(@rdivide,rp6,std(rp6));                   % variance normalization   
            savename = fullfile(path.FunImg,subject{i},'nvr_rp6'); % path and name to save
            save(savename,'rp6');                                  % save the resulting parameters to file
        end
    case 'friston24'
        for i=1:num.subjects
            m       = load(spm_select('FPList',fullfile(path.FunImg,subject{i}),'^rp_.*.txt')); % load realignment parameters 
            mlag    = [zeros(1,6); m(1:end-1,:)];           % calculate 1 lag
            rp24exp = [m m.^2 mlag mlag.^2];                 % compute the 24 parameters
            rp24exp = bsxfun(@minus,rp24exp,mean(rp24exp));  % demean
            rp24exp = bsxfun(@rdivide,rp24exp,std(rp24exp)); % variance normalization

            savename = fullfile(path.FunImg,subject{i},'nvr_rp24exp'); % path and name to save
            save(savename,'rp24exp');                                  % save the resulting parameters to file
        end
    case 'derivative12'
        for i=1:num.subjects
            m       = load(spm_select('FPList',fullfile(path.FunImg,subject{i}),'^rp_.*.txt')); % load realignment parameters 
            mlag    = [zeros(1,6); m(1:end-1,:)];           % calculate 1 lag
            rp12exp = [m mlag];                 % compute the 12 parameters
%             rp12exp = bsxfun(@minus,rp12exp,mean(rp12exp));  % demean
%             rp12exp = bsxfun(@rdivide,rp12exp,std(rp12exp)); % variance normalization

            savename = fullfile(path.FunImg,subject{i},'nvr_rp12exp'); % path and name to save
            save(savename,'rp12exp');                                  % save the resulting parameters to file
        end
end