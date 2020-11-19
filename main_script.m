%% set work time frames and path
% example setting
tframes=101:105;%the time-frame range you want to reconstruct, e.g. 101:1:2700
rawdata_path=pwd;
mreg_filename='meas_MID00242_FID18308_MREG.dat';%mreg rawdata file name, e.g. meas_MID00133_FID01338_MREG.dat
ref_filename='meas_MID00232_FID18298_gre_field_mapping.dat';%reference rawdata file name, e.g. meas_MID00100_FID01305_gre_field_mapping.dat
recon_path=pwd;%
% load your trajectory and correct gradient delay
mreg_path=which('mreg_recon_tool');
traj=loadTrajectory([mreg_path(1:(end-33)) 'grad_files/SoS_3mm_FOV_192_3_6_2_5.grad'],[],[-3.7 -3 -3]);

%% Set reconstruction parameters
[~,header]=loadData([rawdata_path '/' mreg_filename],1);
recon_details.recon_resolution=[64 64 50];%desired reconstrution spatial resolution
recon_details.voxel_size=[0.003 0.003 0.003];%mm
recon_details.DeltaT = header.te(1) + traj.idx{1}(1)*5e-6;
recon_details.rawdata_filename=header.rawdata_filename;
[~,idx_k0]=min(makesos(traj.trajectory{1}(traj.idx{1},:),2));
DORK_k0 = traj.idx{1}(idx_k0);
recon_details.DORK_k0=DORK_k0;
[recon_details.DORK_frequency,recon_details.DORK_phi_offset] = DORK_frequency...
    (recon_details.rawdata_filename,recon_details.DORK_k0,length(traj.trajectory));

recon_details.penalty.norm_string='L2-norm';
recon_details.penalty.norm=@L2Norm;
recon_details.penalty.lambda=0.2;
% recon_details.penalty.norm_string='L1-norm';
% recon_details.penalty.norm=@L1Norm;
% recon_details.penalty.lambda=5e-6;

recon_details.penalty.operator_string='identity';
recon_details.penalty.operator.handle = @identityOperator;
recon_details.penalty.operator.args = {};
recon_details.pname=[recon_path '/' recon_details.penalty.norm_string '_' ...
    num2str(recon_details.penalty.lambda)];

% recon_details.penalty.operator_string='TV';
% recon_details.penalty.operator(1).handle = @finiteDifferenceOperator;
% recon_details.penalty.operator(1).args = {1};
% recon_details.penalty.operator(2).handle = @finiteDifferenceOperator;
% recon_details.penalty.operator(2).args = {2};
% recon_details.penalty.operator(3).handle = @finiteDifferenceOperator;
% recon_details.penalty.operator(3).args = {3};
% recon_details.pname=[pathname recon_details.penalty.norm_string '_' ...
%     num2str(recon_details.penalty.lambda)];

recon_details.offresonance_correction_flag=1;%off-resonance correction
recon_details.max_iterations=20;
recon_details.recon_output_format='nifti';
recon_details.global_frequency_shift=0;
recon_details.tolerance=1e-5;
recon_details.timeframes=tframes;
recon_details.coil_compression=0;%coil array compression
recon_details.num_coil=32;
recon_details.acceleration_rate=0.1;%tPCR acceleration rate, 0.1 means 10 times acceleration,..
%0.01 means 100 times acceleration

mkdir(recon_details.pname);
save([recon_details.pname '/recon_details.mat'],'recon_details');

% calculate the sensitivity map and field map
if ~exist([recon_details.pname '/data.mat'],'file')
    data=sensitivity_field_map([rawdata_path '/' mreg_filename],[rawdata_path '/' ref_filename],recon_details);
    save([recon_details.pname '/data'],'data');
end


%% Standard Reconstruction
%%
sr_recon(recon_details.timeframes,recon_details.pname);

%when you want to use slurm to accelerate
% slurm_submit_mreg(recon_details.timeframes,1,recon_details.pname,'job_slurm',50,1,'');

%% combine 3d images to 4d
tpcr_combination(recon_details.timeframes,recon_details.pname,1,'sr');



%% tPCR

%% step 1: decomposition
for n=recon_details.timeframes(1):1000:recon_details.timeframes(end)
    tframe=n:1:min(n+999,recon_details.timeframes(end));
    subfolder1=[recon_details.pname '/tpcr/' num2str(n) '/'];
    if ~exist([subfolder1 'SV/Vt.mat'],'file')
        tpcr_decomposition(tframe,subfolder1);
      %when you want to use slurm to accelerate
%         slurm_submit_mreg([n recon_details.timeframes(end)],1,subfolder1,'job_slurm',2,3,'ram');
    end
end

%% step 2: reconstruct principal components - segmented
nn=recon_details.timeframes(1):1000:recon_details.timeframes(end);
for n=nn
    tframe=n:1:min(n+999,recon_details.timeframes(end));
    subfolder1=[recon_details.pname '/tpcr/' num2str(n) '/'];
    while ~exist([subfolder1 'kbase/rawlevel_' num2str(floor(length(tframe)*recon_details.acceleration_rate)+1) '.mat'],'file')
        printf('wait...')
        pause(180);
    end
    tpcr_recon(1:max([1,length(tframe)*recon_details.acceleration_rate]),subfolder1);
  %when you want to use slurm to accelerate
%     slurm_submit_mreg(1:max([1,length(tframe)*accelerate]),1,subfolder1,'job_slurm2',2,2,'ram');
end

%% step 3: recombination
 
fname = fullfile([recon_details.pname '/tpcr/4D_' num2str(recon_details.acceleration_rate*100) '%.nii']);
if ~exist(fname,'file')
    for n=recon_details.timeframes(1):1000:recon_details.timeframes(end)
        tframe=n:min(n+999,recon_details.timeframes(end));
        subfolder1=[recon_details.pname '/tpcr/' num2str(n) '/'];
        tpcr_combination(tframe,subfolder1,recon_details.acceleration_rate,'tpcr')
    end

    for n=recon_details.timeframes(1):1000:recon_details.timeframes(end)
        subfolder1=[recon_details.pname '/tpcr/' num2str(n) '/nifti/'];
        if n==recon_details.timeframes(1)
            temp=load_nii([subfolder1 '4D_' num2str(recon_details.acceleration_rate*100) '%.nii']);
            delete([subfolder1 '4D_' num2str(recon_details.acceleration_rate*100) '%.nii']);
            image4d=temp.img;
        else
            temp=load_nii([subfolder1 '4D_' num2str(recon_details.acceleration_rate*100) '%.nii']);
            delete([subfolder1 '4D_' num2str(recon_details.acceleration_rate*100) '%.nii']);
            image4d(:,:,:,(end+1):(end+size(temp.img,4)))=temp.img;
        end
    end
    trad=make_nii(single(image4d),[3 3 3]);
    save_nii(trad,fname);
end

