function data=sensitivity_field_map(mreg_fname,ref_fname,recon_details)
    
    [~, header] = loadData(mreg_fname,1);

        reference = loadReference(ref_fname);
        reference = resizeReference(reference, recon_details.voxel_size.*recon_details.recon_resolution, recon_details.recon_resolution,header);
        data.smaps      = reference.smaps;
        data.wmap       = reference.wmap;
        data.anatomical = reference.anatomical;
        data.shift      = reference.shift;
        data.Cor_angle  = reference.Cor_angle;
        data.ref_header = reference.header;
        data.sensmode   = reference.sensmode;

        tra=loadTrajectory(['/raid/groupspace/ovoc/wang/old/gradient/' header.trajectory],[],[-3.7 -3 -3]);
        data.trajectory=tra;
   
end