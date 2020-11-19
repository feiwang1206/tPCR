function recon=tpcr_recon(tframes,pathname)


    % make folders for storation
    filename_save=[pathname];
    if ~exist([filename_save '/ibase'],'dir')
        mkdir([filename_save '/ibase']);
    end
    
    % load data which includes the sensitivity map, field map and
    % trajectory
    num=0;
    for i=length(pathname):-1:1
        cha=pathname(i);
        if strcmp(cha,'/')
            upfolder=pathname(1:i);
            num=num+1;
            if num==3
                break;
            end
        end
    end
    data=load([upfolder 'data.mat']);
    data=data.data;
    recon_details=load([upfolder 'recon_details.mat']);
    recon_details=recon_details.recon_details;
    traj = data.trajectory.trajectory;
    traj_idx = data.trajectory.idx;
    %% coil compression
    if recon_details.coil_compression==1
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
    
%     % reconstruct the mean image
%     fname0=[filename_save '/mean.mat'];
%     if exist(fname0,'file')
%         load(fname0);
%     else
%         fname1=[pathname 'rawdata_mean/rawdata_mean.mat'];
%         load(fname1);
%         [recon,recon_iter] = regularizedReconstruction(Fg,rawdata_mean(:),P{:}, ...
%             'tol',recon_details.tolerance, ...
%             'maxit',20, ...
%             'verbose_flag', 1);
%         save(fname0,'recon_iter');
%     end
%     
%     
%     %determine mean number of iterations and the mean convergence speed
%     iter0=length(recon_iter{3});
%     y=log10(recon_iter{3}(1:iter0));
%     ll=length(y);
%     clear A;
%     A(:,1)=ones(ll,1);
%     A(:,2)=1:ll;
%     cor=A\y';
%     iter_step=-1/cor(2);
    iter_step=67;
    
    %calculate the distribution of number of iterations of each components
    load([pathname 'SV/St.mat']);
    st=diag(St);
    load([pathname 'rawdata_mean/rawdata_mean']);
    norm_mean=l2norm(rawdata_mean);
%     rawlevel_r=load_nii([pathname 'SV/Ut_r.nii']);
%     rawlevel_i=load_nii([pathname 'SV/Ut_i.nii']);
%     rawlevel_multi=rawlevel_r.img(:,:,tframes(n))+1i*rawlevel_i.img(:,:,tframes(n));
    

%     fname0=[pathname 'kbase/rawlevel_1.mat'];
%     load(fname0);
    iter=zeros(length(st),1);
    iter1=1;
    while(mean(iter(:,1))<recon_details.max_iterations)
        iter1=iter1+1;
        for n=1:length(st)
            iter(n,1)=max(0,min(1000,iter1+floor(iter_step*(log10(st(n)/st(1)))))); 
        end
    end
    % reconstruct each principal component
    
    for n=1:length(tframes)
        fname=[pathname 'kbase/rawlevel_',num2str(tframes(n)),'.mat'];tframes(n)
        load(fname);
%         rawlevel=rawlevel_multi(:,:,n);
        norm1=l2norm(double(rawlevel(:)));

        %scaling the regularization parameter according to the
        %magnitude if using l1-norm regularization
        if strcmp(recon_details.penalty.norm_string,'L1-norm')
            P{2} = double(lambda*norm1/norm_mean);
        end
            
        fname0=[filename_save 'ibase/',num2str(tframes(n)),'.mat'];
        if ~exist(fname0)
            if iter(tframes(n),1)>0
                [recon,recon_iter] = regularizedReconstruction_redis(Fg,rawlevel(:),P{:}, ...
                    'maxit',iter(tframes(n),:), ...
                    'verbose_flag', 1);
                recon_iter{1}{1}=[];
                save(fname0,'recon_iter');
            
                niifolder=[filename_save '/nifti/real'];
                if ~exist(niifolder,'dir')
                    mkdir(niifolder);
                end
                fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
                trad=make_nii(single(real(recon)),[3 3 3]);
                save_nii(trad,fname);

                niifolder=[filename_save '/nifti/imag'];
                if ~exist(niifolder,'dir')
                    mkdir(niifolder);
                end
                fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
                trad=make_nii(single(imag(recon)),[3 3 3]);
                save_nii(trad,fname);
            end
%         else
%             if iter(tframes(n),1)>0
%                 recon0_r=load_nii([filename_save '/nifti/real/' num2str(tframes(n)) '.nii']);
%                 recon0_i=load_nii([filename_save '/nifti/imag/' num2str(tframes(n)) '.nii']);
%                 recon0=recon0_r.img+1i*recon0_i.img;
%                 [recon,recon_iter] = regularizedReconstruction_redis(Fg,rawlevel(:),P{:}, ...
%                     'maxit',10, ...
%                     'verbose_flag', 1,...
%                     'z0',double(recon0));
%                 recon_iter{1}{1}=[];
%                 recon_iter0=load([filename_save '/ibase/' num2str(tframes(n)) '.mat']);
%                 recon_iter{1}{2}(1)=recon_iter0.recon_iter{1}{2}(1)+recon_iter{1}{2}(1);
%                 recon_iter{3}=[recon_iter0.recon_iter{3} recon_iter{3}];
%                 recon_iter{5}=[recon_iter0.recon_iter{5} recon_iter{5}];
%                 save(fname0,'recon_iter');
%             
%                 niifolder=[filename_save '/nifti20/real'];
%                 if ~exist(niifolder,'dir')
%                     mkdir(niifolder);
%                 end
%                 fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
%                 trad=make_nii(single(real(recon)),[3 3 3]);
%                 save_nii(trad,fname);
% 
%                 niifolder=[filename_save '/nifti20/imag'];
%                 if ~exist(niifolder,'dir')
%                     mkdir(niifolder);
%                 end
%                 fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
%                 trad=make_nii(single(imag(recon)),[3 3 3]);
%                 save_nii(trad,fname);
%             end
        end
    end
end
