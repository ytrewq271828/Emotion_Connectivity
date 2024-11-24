spm('defaults', 'FMRI');
spm_jobman('initcfg');

IDCell=load('subjectIDs.mat');
IDs=IDCell.IDs;

groupLevelPath='D:\GraduationStudySPM\Emotion_Task_2nd';
if(~exist('csvData.mat','file'))
    csvPath='D:\GraduationStudySPM\HCP_HYA_Demographics.csv';
    csvData=readcell(csvPath);
    save csvData.mat csvData
else
    csvData=load('csvData.mat');
    csvData=csvData.csvData;
end

if(~exist('contrasts.mat','file'))
    contrasts=cell(3,100);
    for idx=1:100
        contrasts{1,idx}=strcat('D:\GraduationStudySPM\Emotion_Task_1st\',num2str(IDs(1,idx)),'\con_0001.nii,1');
        %contrasts{2,idx}=strcat('D:\GraduationStudySPM\Emotion_Task_1st\',num2str(IDs(1,idx)),'\con_0002.nii,1');
        %contrasts{3,idx}=strcat('D:\GraduationStudySPM\Emotion_Task_1st\',num2str(IDs(1,idx)),'\con_0003.nii,1');
    end
    save contrasts.mat contrasts
else
    contrasts=load('contrasts.mat');
    contrasts=contrasts.contrasts;
end

marsbar('on');

% ROI that showed significantly different activation in faces vs nonfaees in Hariri 2002
% Talairach [16, -4, -16] => MNI [17, -1, -23]
center = [17, -1, -23];  
radius = 7;              % Example radius in mm

%roi = maroi_sphere(struct('centre', center, 'radius', radius));


%roiFile ='D:\GraduationStudySPM\sphere_7-17_-1_-23_roi.mat';
%roiNii='D:\GraduationStudySPM\sphere_7-17_-1_-23_roi.nii';
%Roi = maroi('load', roiFile);
 
%saveroi(Roi, roiNii);


roiGUI='D:\GraduationStudySPM\sphere_10-17_-1_-23.nii';

matlabbatch = [];
matlabbatch{1}.spm.stats.factorial_design.dir = {strcat(groupLevelPath,'\roianalysis')};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(contrasts(1,:))';
matlabbatch{1}.spm.stats.factorial_design.cov = struct([]); % No covariates for 1-sample t-test
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {roiGUI};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Run the design specification
spm_jobman('run', matlabbatch);

% Load the SPM.mat file from the second-level analysis directory
spmMatPath = fullfile(strcat(groupLevelPath,'\roianalysis'), 'SPM.mat');

% Set up the batch job for model estimation
matlabbatch = [];
matlabbatch{1}.spm.stats.fmri_est.spmmat = {spmMatPath};

% Run the model estimation
spm_jobman('run', matlabbatch);

% Set up the batch job for contrast manager
matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {spmMatPath};


% Define the T-contrast
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Faces vs. Shapes';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = 1; % Already generated a contrast in 1st-level
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% Run the contrast manager
spm_jobman('run', matlabbatch);




% spm_path = 'E:\GraduationStudySPM\Emotion_Task_2nd\roianalysis\SPM.mat';
% D = mardo(spm_path);
% 
% Y = get_marsy(ROI, D, 'mean');
% 
% activation_values = summary_data(Y);
% 
% disp(activation_values);