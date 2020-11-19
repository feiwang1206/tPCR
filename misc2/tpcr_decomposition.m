function tpcr_decomposition(tframes,subfolder)

    
    % make folders for storation
    if ~exist([subfolder '/kbase'],'dir')
        mkdir([subfolder '/kbase']);
        mkdir([subfolder '/SV']);
        mkdir([subfolder '/rawdata_mean']);
    end
    shift=0;
    
    % load data which includes the sensitivity map, field map and
    % trajectory
    num=0;
    for i=length(subfolder):-1:1
        cha=subfolder(i);
        if strcmp(cha,'/')
            pathname=subfolder(1:i);
            num=num+1;
            if num==3
                break;
            end
        end
    end

    data=load([pathname 'data.mat']);
    data=data.data;
    recon_details=load([pathname 'recon_details.mat']);
    recon_details=recon_details.recon_details;
    traj = data.trajectory.trajectory;
    traj_idx = data.trajectory.idx;
    if recon_details.coil_compression==1
        %% coil compression with demean
        smaps_mean=mean(data.smaps,4);
        data.smaps=data.smaps-repmat(smaps_mean,[1 1 1 size(data.smaps,4)]);
        smaps2d=reshape(data.smaps,[64*64*50 size(data.smaps,4)]);
        [Vt,~]=eig(smaps2d'*smaps2d);
        Vt=Vt(:,end:-1:1);
        num_coil=recon_details.num_coil;
        Vt=Vt(:,1:num_coil);
        Ut=smaps2d*Vt;    
        data.smaps=reshape(Ut,[64 64 50 num_coil])+repmat(smaps_mean,[1 1 1 num_coil]);
    end
        
    % load the kspace rawdata
    tic;
%     [rawdata, header] = loadData([recon_details.rawdata_filename],tframes);
%     shift_raw_data1
%     if length(traj)==1
%         interleaf_idx=1;
%         shift_raw_data2
%         rawdata2d1=reshape(rawdata,[size(rawdata,1)*size(rawdata,2) size(rawdata,3)]);
%         rawdata_mean=mean(rawdata2d1,2);
%     else
%         rawdata3d0=rawdata;
%         tframes_idx=mod(tframes-1,length(traj))+1;
%         for ii=1:length(traj)
%             tframes_ii=tframes(tframes_idx==ii);
%             rawdata=rawdata3d0(:,:,tframes_ii-tframes(1)+1);
%             interleaf_idx=ii;
%             shift_raw_data2
%             if ii==1
%                 rawdata3d(:,:,1:length(tframes_ii))=rawdata;
%             else
%                 rawdata3d(:,:,(end+1):(end+length(tframes_ii)))=rawdata;
%             end
%         end
%         rawdata2d1=reshape(rawdata3d,[size(rawdata3d,1)*size(rawdata3d,2) size(rawdata3d,3)]);
%         rawdata_mean=mean(rawdata2d1,2);
%     end
%     save([subfolder '/rawdata_mean/rawdata_mean.mat'],'rawdata_mean');

    for n=1:length(tframes)       
        for ii=1:length(traj)
            interleaf_idx=ii;
            kk=tframes(n)+ii-1;
            [rawdata, header] = loadData([recon_details.rawdata_filename],kk);
            shift_raw_data
            %% coil compression
            if recon_details.coil_compression==1
                rawdata_mean_coil=mean(rawdata,2);
                rawdata=(rawdata-repmat(rawdata_mean_coil,[1 size(rawdata,2)]))*Vt;
                rawdata=rawdata+repmat(rawdata_mean_coil,[1 num_coil]);
            end
            if ii==1
                rawdata2d=rawdata;
            else
                rawdata2d((end+1):(end+length(rawdata)),:)=rawdata;
            end
        end
        
        rawdata2d1(:,n)=rawdata2d(:);
    end
    rawdata_mean=mean(rawdata2d1,2);
    save([subfolder '/rawdata_mean/rawdata_mean.mat'],'rawdata_mean');
toc
    % singluar value decomposition
%     [Ut1,St1,Vt1]=svd(rawdata2d1,0);
    [Vt,St]=eig(rawdata2d1'*rawdata2d1);
    St=sqrt(St(end:-1:1,end:-1:1));
    Ut=rawdata2d1*Vt;
    Ut=Ut(:,end:-1:1);Vt=Vt(:,end:-1:1);
    fname=[subfolder '/SV/St.mat'];
    save(fname,'St');
    fname=[subfolder '/SV/Vt.mat'];
    save(fname,'Vt');
toc
    % save the scaled spatial part
    for n=1:length(tframes)
        fname=[subfolder '/kbase/rawlevel_',num2str(n),'.mat'];
        rawlevel=Ut(:,n);
        rawlevel=reshape(rawlevel,[],size(data.smaps,4));
        save(fname,'rawlevel');
    end
%     Ut=reshape(Ut,size(rawdata));
%     trad=make_nii(single(real(Ut)),[3 3 3]/3);
%     save_nii(trad,[subfolder '/SV/Ut_r.nii']);
%     trad=make_nii(single(imag(Ut)),[3 3 3]/3);
%     save_nii(trad,[subfolder '/SV/Ut_i.nii']);
    toc;

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

    function shift_raw_data1
        %% DORK and off-resonance correction
        t = header.te(1)*ones(size(rawdata,1),1) + (0:size(rawdata,1)-1)'*header.dwelltime;
        recon_details.global_frequency_shift=shift;
        % which is caused by heating of the system if without frequency adjustment,
        %here 60 is estimated by eyes
        if isempty(recon_details.DORK_frequency)
            freq_offset = recon_details.global_frequency_shift;
            phi_offset = 0;
        else
            freq_offset = recon_details.DORK_frequency(tframes) + recon_details.global_frequency_shift;
            phi_offset = recon_details.DORK_phi_offset(tframes);
        end
        rawdata = rawdata .* exp(-1i*(repmat(permute(phi_offset,[2 3 1]),[size(rawdata, 1) size(rawdata, 2) 1])+...
            repmat(permute(freq_offset,[2 3 1]),[size(rawdata, 1) size(rawdata, 2) 1]).*repmat(t,[1 size(rawdata, 2) size(rawdata, 3)])));
    end

    function shift_raw_data2
        %% Adding additional Phase to the data for data shifting
        rawdata = rawdata(traj_idx{interleaf_idx},:,:) .* repmat(exp(1i*traj{interleaf_idx}...
            (traj_idx{interleaf_idx},[2 1 3])*(data.shift+[0.5 -0.5 0.5])'), ...
            [1 size(rawdata, 2) size(rawdata, 3)]);%from 0 to 0.5
    end
end
