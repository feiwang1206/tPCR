function recon = mreg_recon_tool_recon(how,varargin)

% function recon = mreg_recon_tool_recon(how,varargin)
%
% This function is called by the mreg_recon_tool or indirectly by
% jobs of the gridengine. It reconstructs data specified in varargin.
%
% Thimo Hugger
% 21.09.2011

tic;

if strcmp(how,'local')
    verbose_flag = 1;
    
    data = varargin{1};
    recon_details = varargin{2};
    clear varargin;
    mkdir([recon_details.pname]);
    if ~exist([recon_details.pname '/data.mat'],'file')
        save([recon_details.pname '/data.mat'],'data');
    end
%     if ~exist([recon_details.pname '/recon_details.mat'],'file')
        save([recon_details.pname '/recon_details.mat'],'recon_details');
%     end
    
    tframes = recon_details.timeframes;
    
    if ~strcmp(recon_details.recon_output_format, 'not')
        if ~exist(fullfile(recon_details.pname),'dir')
            mkdir(fullfile(recon_details.pname));
        end
        save([recon_details.pname, '/recon_details.mat'], 'recon_details','-v7.3');
    end
    

elseif strcmp(how,'sge')
    verbose_flag = 2;
    
    frame_number = varargin{1};
    wdir         = varargin{2};
    clear varargin;
    
    recon_details = mat2variable([wdir,'/recon_details.mat']);
    data = mat2variable([recon_details.pname,'/data.mat']);
    job = mat2variable([recon_details.pname,'/job_data_' num2str(frame_number),'.mat']);
    tframes = job.timeframes;

    if exist('spm_write_vol.m')~=2
        fprintf('SPM not found. No NIFTI export possible. Switching to mat-file export.\n');
        recon_details.recon_output_format = 'mat';
    end
end

Ni = recon_details.nInterleaves;
dim = recon_details.recon_resolution;
dwelltime = recon_details.dwelltime;
traj_scale = recon_details.trajectory_scaling;
Nt = recon_details.Nt;
rawdata_filename = recon_details.rawdata_filename;
traj = data.trajectory.trajectory;
traj_idx_org = data.trajectory.idx;

rawdata_header = readHeader(rawdata_filename);
if isfield(rawdata_header,'inPlaneRot') && abs(rawdata_header.inPlaneRot-pi/2)<1e-5
    for ii=1:length(traj)
        traj{ii} = traj{ii}(:,[2 1 3]);
        traj{ii}(:,2) = -traj{ii}(:,2);
    end
end


% assemble arguments for regularizedReconstruction
lengthP = 0;
for k=1:length(recon_details.penalty)
    lengthP = lengthP + length(recon_details.penalty(k).operator) + 2;
end

P = cell(1,lengthP);
counter = 1;
for k=1:length(recon_details.penalty)
    P{counter} = recon_details.penalty(k).norm;
    counter = counter + 1;
    P{counter} = recon_details.penalty(k).lambda;
    counter = counter + 1;
    for n=1:length(recon_details.penalty(k).operator)
        P{counter} = recon_details.penalty(k).operator(n).handle(recon_details.penalty(k).operator(n).args{:});
        counter = counter + 1;
    end
end
recon=[];
if strcmp(recon_details.recon_output_format, 'not')
    recon = single(zeros([dim length(tframes)]));
end
center_interleaf_idx = 1;
traj_idx = traj_idx_org;

% PC decomposition
if strcmp(recon_details.recon_type,'tPCR_step1')
    
    for n=tframes(1):1000:tframes(end)
        tframe=n:min(n+999,tframes(end));
        subfolder=[recon_details.pname '/tpcr/' num2str(n) '/'];
        if recon_details.EnableSlurm
            slurm_submit_mreg([n tframes(end)],1,subfolder,'job_slurm',2,3,'ram');
        else
            tpcr_decomposition(tframe,subfolder);
        end
    end

elseif strcmp(recon_details.recon_type,'tPCR_step2') 
    n_gap=tframes(1):1000:tframes(end);
    for n=n_gap
        tframe=n:min(n+999,tframes(end));
        subfolder=[recon_details.pname '/tpcr/' num2str(n) '/'];
        while ~exist([subfolder 'kbase/rawlevel_' num2str(floor(length(tframe)*recon_details.tPCR_acc_rate)) '.mat'],'file')
            printf('wait...')
            pause(60);
        end
        if recon_details.EnableSlurm
            slurm_submit_mreg(1:max(1,floor(length(tframe)*recon_details.tPCR_acc_rate)),1,subfolder,'job_slurm2',recon_details.frames_per_job,2,'ram');
        else
            tpcr_recon(1:max(1,floor(length(tframe)*recon_details.tPCR_acc_rate)),subfolder);
        end
    end
elseif strcmp(recon_details.recon_type,'tPCR_step3')
    for n=tframes(1):1000:tframes(end)
        tframe=n:min(n+999,tframes(end));
        subfolder=[recon_details.pname '/tpcr/' num2str(n) '/'];
        tpcr_combination(tframe,subfolder,recon_details.tPCR_acc_rate,'tpcr')
    end

    fname = fullfile([recon_details.pname '/tpcr/4D_' num2str(recon_details.tPCR_acc_rate*100) '%.nii']);
    if ~exist(fname,'file')
        for n=tframes(1):1000:tframes(end)
            subfolder=[recon_details.pname '/tpcr/' num2str(n) '/nifti/'];
            if n==tframes(1)
                temp=load_nii([subfolder '4D_' num2str(recon_details.tPCR_acc_rate*100) '%.nii']);
                delete([subfolder '4D_' num2str(recon_details.tPCR_acc_rate*100) '%.nii']);
                image4d=temp.img;
            else
                temp=load_nii([subfolder '4D_' num2str(recon_details.tPCR_acc_rate*100) '%.nii']);
                delete([subfolder '4D_' num2str(recon_details.tPCR_acc_rate*100) '%.nii']);
                image4d(:,:,:,(end+1):(end+size(temp.img,4)))=temp.img;
            end
        end
        trad=make_nii(single(image4d),[3 3 3]);
        save_nii(trad,fname);
    end
elseif strcmp(recon_details.recon_type,'standard')
    
    if recon_details.EnableSlurm
        slurm_submit_mreg(tframes,1,recon_details.pname,'job_slurm',recon_details.frames_per_job,1,'');
    else
        sr_recon(tframes,recon_details.pname)
    end

elseif strcmp(recon_details.recon_type,'standard_combine')
    tpcr_combination(tframes,recon_details.pname,[],'sr')
end

[~, tstr] = seconds2humanreadable(toc);
c = clock;
[~, host] = unix('hostname');
fprintf(['Reconstructiontime = ', tstr, '. Calculated on ', host, 'Finished on ', num2str(c(1), '%4.4i'), '/', num2str(c(2), '%2.2i'), '/', num2str(c(3), '%2.2i'), ' at ', num2str(c(4), '%2.2i'), ':', num2str(c(5), '%2.2i'), '\n']);
    
% create flag file when finished
if strcmp(how,'sge')
    unix(['touch ', recon_details.pname, '/recon_', num2str(frame_number), '.flag']);
end




end