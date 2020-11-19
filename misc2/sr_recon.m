function recon=sr_recon(tframes,pathname)
% reconstruction with SR
%     tframes: the time frames to reconstruct
%     pathname: the folder where the rawdatas store and the reconstructed images save
% 19.11.2020
% Fei Wang

    
    
    % load data which includes the sensitivity map, field map and
    % trajectory
    data=load([pathname '/data.mat']);
    data=data.data;
    recon_details=load([pathname '/recon_details.mat']);
    recon_details=recon_details.recon_details;
    traj = data.trajectory.trajectory;
    traj_idx = data.trajectory.idx;

    if recon_details.coil_compression==1
        %% coil compression
        smaps_mean=mean(data.smaps,4);
        data.smaps=data.smaps-repmat(smaps_mean,[1 1 1 size(data.smaps,4)]);
        smaps2d=reshape(data.smaps,[64*64*50 size(data.smaps,4)]);
        [Vt,~]=eig(smaps2d'*smaps2d);
        Vt=Vt(:,end:-1:(end-recon_details.num_coil+1));
        Ut=smaps2d*Vt;    
        data.smaps=reshape(Ut,[64 64 50 recon_details.num_coil])+repmat(smaps_mean,[1 1 1 recon_details.num_coil]);
    end
    % make folders for storation
    filename_save=[recon_details.pname '/sr'];
    if strcmp(recon_details.recon_output_format,'mat')
        mkdir([filename_save '/mat']);
    end
    if strcmp(recon_details.recon_output_format,'nifti')
        mkdir([filename_save '/nifti']);
    end
    
    for ii=1:length(traj)
        K{ii}=traj{ii}(traj_idx{ii},:);
    end
    
    % make forward operator
    if recon_details.offresonance_correction_flag==1
        Fg=orc_segm_nuFTOperator_multi(K,size(data.wmap),data.smaps,data.wmap,5e-6,10,recon_details.DeltaT,'minmax:kb');
    else
        Fg=nuFTOperator_multi(K,size(data.wmap),data.smaps,'','','minmax:kb');
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
    

    % reconstruction  
    for k=tframes
        if mod(k-1,length(traj))==0
            for ii=1:length(traj)
                interleaf_idx=ii;
                kk=k+ii-1;
                [rawdata, header] = loadData([recon_details.rawdata_filename],kk);
                shift_raw_data
                %% coil compression
                if recon_details.coil_compression==1
                    rawdata_mean=mean(rawdata,2);
                    rawdata=(rawdata-repmat(rawdata_mean,[1 size(rawdata,2)]))*Vt;
                    rawdata=rawdata+repmat(rawdata_mean,[1 recon_details.num_coil]);
                end
                if ii==1
                    rawdata2d=rawdata;
                else
                    rawdata2d((end+1):(end+length(rawdata)),:)=rawdata;
                end
            end

            if strcmp(recon_details.recon_output_format,'mat')
                fname0=[filename_save '/mat/' num2str(k) '.mat'];
                if ~exist(fname0,'file')
                    k
                    [recon] = regularizedReconstruction(Fg,rawdata2d(:),P{:}, ...
                        'tol',recon_details.tolerance, ...
                        'maxit',recon_details.max_iterations, ...
                        'verbose_flag', 1);

                    fname = fullfile([filename_save '/mat/' num2str(k) '.mat']);
                    save(fname,'recon');
                end
            else
                fname0=[filename_save '/nifti/real/' num2str(k) '.nii'];
                if ~exist(fname0,'file')
                    k
                    [recon] = regularizedReconstruction(Fg,rawdata2d(:),P{:}, ...
                        'tol',recon_details.tolerance, ...
                        'maxit',recon_details.max_iterations, ...
                        'verbose_flag', 1);

                    if ~exist([filename_save '/nifti/real'],'dir')
                        mkdir([filename_save '/nifti/real']);
                    end
                    if ~exist([filename_save '/nifti/imag'],'dir')
                        mkdir([filename_save '/nifti/imag']);
                    end
                    fname = fullfile([filename_save '/nifti/real/' num2str(k) '.nii']);
                    trad=make_nii(single(real(recon)),[3 3 3]);
                    save_nii(trad,fname);
                    fname = fullfile([filename_save '/nifti/imag/' num2str(k) '.nii']);
                    trad=make_nii(single(imag(recon)),[3 3 3]);
                    save_nii(trad,fname);
                end
            end
        end
    end

    function shift_raw_data

    %% DORK and off-resonance correction
    t = header.te(1)*ones(size(rawdata,1),1) + (0:size(rawdata,1)-1)'*header.dwelltime;
    % which is caused by heating of the system if without frequency adjustment,
    %here 60 is estimated by eyes
    if isempty(recon_details.DORK_frequency)
        freq_offset = recon_details.global_frequency_shift;
        phi_offset = 0;
    else
        freq_offset = recon_details.DORK_frequency(kk) + recon_details.global_frequency_shift;
        phi_offset = recon_details.DORK_phi_offset(kk);
    end
    rawdata = rawdata .* repmat(exp(-1i*(phi_offset+freq_offset.*t)), [1 size(rawdata, 2)]);

    %% Adding additional Phase to the data for data shifting
    rawdata = rawdata(traj_idx{interleaf_idx},:) .* repmat(exp(1i*traj{interleaf_idx}(traj_idx{interleaf_idx},[2 1 3])...
        *(data.shift+[0.5 -0.5 0.5])'), [1 size(rawdata, 2)]);%from 0 to 0.5
    end
end
