spm('defaults', 'FMRI');
spm_jobman('initcfg');

IDCell=load('subjectIDs.mat');
IDs=IDCell.IDs;

groupLevelPath='D:\GraduationStudySPM\PPIAnalysis';

% Con_0001.nii => Faces vs. Baseline
contrasts=cell(3,100);
for idx=1:100
    contrasts{1,idx}=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',num2str(IDs(1,idx)),'\con_0001.nii,1');
    contrasts{2,idx}=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',num2str(IDs(1,idx)),'\con_0002.nii,1');
    contrasts{3,idx}=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',num2str(IDs(1,idx)),'\con_0003.nii,1');
end


mkdir("PPIAnalysis\Faces2nd");
matlabbatch = [];
matlabbatch{1}.spm.stats.factorial_design.dir = {strcat(groupLevelPath, '\Faces2nd')};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(contrasts(1,:))';
matlabbatch{1}.spm.stats.factorial_design.cov = struct([]); % No covariates for 1-sample t-test
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Run the design specification
spm_jobman('run', matlabbatch);

% Load the SPM.mat file from the second-level analysis directory
spmMatPath = fullfile(strcat(groupLevelPath, '\Faces2nd'), 'SPM.mat');

% Set up the batch job for model estimation
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmMatPath};

% Run the model estimation
spm_jobman('run', matlabbatch);

% Set up the batch job for contrast manager
matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {spmMatPath};


% Define the T-contrast
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Faces vs. Baseline';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1; % Already generated a contrast in 1st-level
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% Run the contrast manager
spm_jobman('run', matlabbatch);
