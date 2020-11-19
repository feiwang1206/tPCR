function data_estimation(mreg_fname,reference_name)

    % load rawdata
    [rawdata, header] = loadData(mreg_fname,1);
    save('header','header');

    % estimate the sensitivity map and field map
    reference = loadReference(reference_name);
    reference = resizeReference(reference, [0.003 0.003 0.003].*[64 64 50], [64 64 50],header);
    data.smaps      = reference.smaps;
    data.wmap       = reference.wmap;
    data.anatomical = reference.anatomical;
    data.shift      = reference.shift;
    data.Cor_angle  = reference.Cor_angle;
    data.ref_header = reference.header;
    data.sensmode   = reference.sensmode;

    % load and save trajectory
    tra=loadTrajectory(['/raid/home/extern/wangfe/software/github/mreg_recon_tool/grad_files/' ...
        header.trajectory],[],[-3.7 -3 -3]);%rotation correction
    data.trajectory=tra;
    save('data','data');

    % make and save the reconstruction details
    recon_details.recon_resolution=size(data.wmap);
    recon_details.DeltaT = header.te(1) + data.trajectory.idx{1}(1)*5e-6;
    recon_details.rawdata_filename=mreg_fname;
    [~,idx_k0]=min(makesos(data.trajectory.trajectory{1}(data.trajectory.idx{1},:),2));
    DORK_k0 = data.trajectory.idx{1}(idx_k0);
    recon_details.DORK_k0=DORK_k0;
    [recon_details.DORK_frequency,recon_details.DORK_phi_offset] = DORK_frequency(recon_details.rawdata_filename,recon_details.DORK_k0,length(data.trajectory.trajectory));
    save('recon_details','recon_details'); 
end