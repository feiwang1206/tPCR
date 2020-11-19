function tpcr_combination(tframes,pathname,rate,method)

matwrite=0;
nii4dwrite=1;

if strcmp(method,'sr')
    for i=1:length(tframes)
        image_r=load_nii([pathname '/sr/nifti/real/' num2str(tframes(i)) '.nii']);
        image_i=load_nii([pathname '/sr/nifti/imag/' num2str(tframes(i)) '.nii']);
        image=abs(image_r.img+image_i.img*1i);
        image4d(:,:,:,i)=image;
    end
    trad=make_nii(abs(single(image4d)),[3 3 3]);
    save_nii(trad,[pathname '/sr/4D.nii']);
elseif strcmp(method,'tpcr')
    clear recon3d_tpcr recon reconbase
    pname=[pathname '/nifti'];
    fname=[pathname '/SV/Vt.mat'];load(fname);
    baseorder=[1:floor(length(tframes))]';
    iter2=max(1,floor(length(baseorder)*rate(end))+1);
    for m=baseorder'
        fname1=[pname '/real/' num2str(m) '.nii'];
        fname2=[pname '/imag/' num2str(m) '.nii'];
        if m<=iter2
            if exist(fname1,'file') && exist(fname2,'file')
                ibase_r=load_nii(fname1);
                ibase_i=load_nii(fname2);
                reconbase(:,:,:,m)=ibase_r.img+1i*ibase_i.img;
            else
                reconbase(:,:,:,m)=0;
            end
        end
    end
    reconbase2d=reshape(reconbase,[size(reconbase,1)*size(reconbase,2)*size(reconbase,3) size(reconbase,4)]);
    for iirate=1:length(rate)
        fname4D=[pname '/4D_' num2str(rate(iirate)*100) '%.nii'];
%         if ~exist(fname4D,'file')
            iter2=max(1,floor(length(baseorder)*rate(iirate))+1);
            reconbase2d_sub=reconbase2d(:,1:iter2);
            recon4d=single(reshape(reconbase2d_sub*(Vt(:,baseorder(1:iter2)))',...
                [size(reconbase,1) size(reconbase,2) size(reconbase,3) length(tframes)]));
            clear reconbase2d_sub
            if matwrite==1
                save([pname '/recon4d.mat'],'recon4d');
            end
            if nii4dwrite==1
                trad=make_nii(abs(single(recon4d)),[3 3 3]);
                save_nii(trad,fname4D);
            end
%         end
    end
end