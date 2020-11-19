function slurm_submit_mreg(timeframes,step,pathname,pname,frames_per_job,mark,nodelist)

% function rpath = mreg_recon_tool_slurm_init(data,recon_details,frames_per_job)
% 
% Submits jobs to the Slurm.
% 
% arguments:
% data = struct with necessary data
% recon_details = struct recon parameters
% frames_per_job = number of timeframes for each job
% 
% Jakob Asslaender
% jakob.asslaender@uniklinik-freiburg.de
% 07.03.2013
% Slurm edit: Bruno Riemenschneider,Fei Wang
% 05.10.2017

% This sets whether the grid engine uses the compiled version. If this
% doesn't work for some reason: set it to zero and you  are back to the old
% version. If you change it: Change it in mreg_recon_tool_resubmit as well!

cdir = pwd;

all_tpts = length(timeframes);
job_start_number = 0;

mcmd = 'matlab -nodisplay -singleCompThread -r';
if mark==1
    methodname=[pathname '/sr/' pname];
    if ~exist(methodname,'dir')
        mkdir(methodname);
    end % don't overwrite any jobs that are already in this folder
    pathname2=['''' pathname ''''];
    ustr = [mcmd, ' "sr_recon_script(%i,%i,%i,',pathname2, ');" '];
elseif mark==2
    methodname=[pathname '/tpcr_' pname];
    if ~exist(methodname,'dir')
        mkdir(methodname);
    end % don't overwrite any jobs that are already in this folder
    pathname2=['''' pathname ''''];
    ustr = [mcmd, ' "tpcr_recon_script(%i,%i,%i,',pathname2, ');" '];
elseif mark==3
    methodname=[pathname '/tpcr_' pname];
    if ~exist(methodname,'dir')
        mkdir(methodname);
    end % don't overwrite any jobs that are already in this folder
    pathname2=['''' pathname ''''];
    ustr = [mcmd, ' "tpcr_decomposition_script(%i,%i,%i,',pathname2, ');" '];
end

job.frames_per_job = frames_per_job;

for k=0:(floor(all_tpts/frames_per_job))
    
    if k == floor(all_tpts/frames_per_job)
        pts = k*frames_per_job+1:length(timeframes);
    else
        pts = k*frames_per_job+1:(k+1)*frames_per_job;
    end

    if ~isempty(pts)
    job.args = {k+job_start_number};
    job.timeframes = timeframes(pts);
    save(fullfile(methodname, ['/job_data_', num2str(k+job_start_number, '%i')]),'job');

    outfname = fullfile(methodname, ['/job_' num2str(k+job_start_number, '%i') '.out']);
    errfname = fullfile(methodname, ['/job_' num2str(k+job_start_number, '%i') '.err']);
    fname = fullfile(methodname, ['/job_' num2str(k+job_start_number, '%i') '.sbatch']);
    fid = fopen(fname,'w');
    %...do the slurmy stuff
    fprintf(fid,'%s\n\n','#!/bin/bash');
    fprintf(fid,'%s\n','#SBATCH --job-name=MREG_reco');
    fprintf(fid,'%s\n',['#SBATCH --output=' outfname]);
    fprintf(fid,'%s\n',['#SBATCH --error=' errfname]);
%    fprintf(fid,'%s\n','#SBATCH --nodes=1');
    fprintf(fid,'%s\n','#SBATCH --partition=engine');
    fprintf(fid,'%s\n',['#SBATCH --nodelist=' nodelist ' # Options: [yeti,nessie,ram]']);
    fprintf(fid,'%s\n','#SBATCH --ntasks=1');
    fprintf(fid,'%s\n','#SBATCH --cpus-per-task=1');
    fprintf(fid,'%s\n','#SBATCH --time=01-00:00:00');
%    fprintf(fid,'%s\n','#SBATCH --mem=100');
%    fprintf(fid,'%s\n','#SBATCH --mail-type=ALL');
%    fprintf(fid,'%s\n','#SBATCH --mail-user=');

    fprintf(fid,'\n%s',sprintf(ustr,job.timeframes(1),job.timeframes(end),step));
    fclose(fid);

%    unix(['chmod 777 ' fname]);
    unix(['sbatch ' fname]);
    end
end

cd(cdir);
