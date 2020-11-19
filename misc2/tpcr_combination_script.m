function tpcr_combination_script(start,send,step,pathname)
rate=[0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1];
load([pathname '/recon_details.mat']);
for n=recon_details.tframes(1):1000:recon_details.tframes(end)
    tframe=n:min(n+999,recon_details.tframes(end));
    subfolder=[recon_details.pname '/tpcr/' num2str(n) '/'];
    tpcr_combination(tframe,subfolder,rate(start:step:send),'tpcr')
end
for iirate=start:step:send
    fname = fullfile([recon_details.pname '/tpcr/4D_' num2str(rate(iirate)*100) '%.nii']);
    if ~exist(fname,'file')

        for n=recon_details.tframes(1):1000:recon_details.tframes(end)
            subfolder=[recon_details.pname '/tpcr/' num2str(n) '/nifti/'];
            if n==recon_details.tframes(1)
                temp=load_nii([subfolder '4D_' num2str(rate(iirate)*100) '%.nii']);
                delete([subfolder '4D_' num2str(rate(iirate)*100) '%.nii']);
                image4d=temp.img;
            else
                temp=load_nii([subfolder '4D_' num2str(rate(iirate)*100) '%.nii']);
                delete([subfolder '4D_' num2str(rate(iirate)*100) '%.nii']);
                image4d(:,:,:,(end+1):(end+size(temp.img,4)))=temp.img;
            end
        end
        trad=make_nii(single(image4d),[3 3 3]);
        save_nii(trad,fname);
    end
end
end