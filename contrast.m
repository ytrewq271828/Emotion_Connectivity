% Initialize SPM
spm('defaults', 'FMRI');
spm_jobman('initcfg');

IDCell=load('subjectIDs.mat');
IDs=IDCell.IDs;

for idx=1:100
    subjectID=num2str(IDs(1,idx));
    if exist(strcat('E:\GraduationStudySPM\Emotion_Task_1st\',subjectID,'\marsS.mat'),'file')
        continue
    end
    % Define the batch job for contrast manager
    matlabbatch = [];
    

    % Specify the SPM.mat file path from the first-level analysis
    spmMatPath = strcat('E:\GraduationStudySPM\Emotion_Task_1st\',subjectID,'\SPM.mat');
    load(spmMatPath);
    %% Atlas load
    
%     if ~exist('D:\GraduationStudySPM\Neuromorphometrics_mask.nii','file')
%         xA=spm_atlas('load');
%         S=spm_atlas('select',xA,'Right Amygdala');
%         mask_file=spm_atlas('mask',xA,S);
%         mask_file = spm_write_vol(mask_file,mask_file.dat);
%         mask_file = rmfield(mask_file,'dat');
%     %     atlas = spm_atlas('load', 'Neuromorphometrics');  
%     %     roi_label = {'Right Amygdala'};  
%     %     roi_indices = spm_atlas('select', atlas, roi_label);
%     %     mask_file = spm_atlas('mask', atlas, roi_indices);
%     %     mask_file = spm_write_vol(mask_file,mask_file.dat);
%     %     mask_file = rmfield(mask_file,'dat');
%     else
%         mask_file=spm_atlas('load',"Neuromorphometrics_mask.nii");
%     end
    %xA=spm_atlas('load','Right Amygdala');
    %S=spm_atlas('select',xA);
     
    % for i = 1:size(S,2)
    %     fname=strcat(S{i},'.nii');
    %     VM=spm_atlas('mask',xA,S{i});
    %     VM.fname=fname;
    %     spm_write_vol(VM,spm_read_vols(VM));
    % end
     
    
    
    %% Contrast Manager
    
    contrastZeros=zeros(1,12);
    contrast1.name='Faces vs. Shapes';
    contrast1.weights=[0.5 -0.5 contrastZeros 0.5 -0.5 contrastZeros 0 0];
    contrast1.sessrep='none';
    
    contrast2.name='Faces vs. Baseline';
    contrast2.weights=[0.5 0 contrastZeros 0.5 0 contrastZeros 0 0];
    contrast2.sessrep='none';
    
    contrast3.name='Shapes vs. Baseline';
    contrast3.weights=[0 0.5 contrastZeros 0 0.5 contrastZeros 0 0];
    contrast3.sessrep='none';
    
    matlabbatch{1}.spm.stats.con.spmmat = {spmMatPath};
    matlabbatch{1}.spm.stats.con.consess{1}.tcon=contrast1;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon=contrast2;
    matlabbatch{1}.spm.stats.con.consess{3}.tcon=contrast3;
    matlabbatch{1}.spm.stats.con.delete = 1;
    %matlabbatch{1}.spm.stats.con.mask = {mask_file};
    spm_jobman('run', matlabbatch);
    
    
    roiFile = 'E:\GraduationStudySPM\marsbar-aal-rois-0.3\marsbar-aal-rois-0.3\MNI_Amygdala_R_roi.mat';
    extantPath='E:\GraduationStudySPM\';
    % Make marsbar design object
    D  = mardo(spmMatPath);
    D = cd_images(D, extantPath);
    save_spm(D);
    % Make marsbar ROI object
    R  = maroi(roiFile);
    % Fetch data into marsbar data object
    Y  = get_marsy(R, D, 'mean');
    % Get contrasts from original design
    xCon = get_contrasts(D);
    % Estimate design on ROI data
    E = estimate(D, Y);
    % Put contrasts from original design back into design object
    E = set_contrasts(E, xCon);
    % get design betas
    b = betas(E);
    % get stats and stuff for all contrasts into statistics structure
    marsS = compute_contrasts(E, 1:length(xCon));
    
    save(strcat('E:\GraduationStudySPM\Emotion_Task_1st\',subjectID,'\marsS.mat'),'marsS');
end
% % Define the contrast number
% contrast_number = 1;  % Replace with your contrast number
% 
% % Load the mask
% %mask = spm_read_vols(spm_vol('D:\GraduationStudySPM\Right_Amygdala_Mask.nii'));
% mask=spm_read_vols(spm_vol('Neuromorphometrics_mask.nii'));
% % Apply the mask to the SPM results
% 
% conPath=strcat('D:\GraduationStudySPM\Emotion_Task_1st\',num2str(subjectID));
% cd(conPath)
% 
% SPM.xCon(contrast_number).Vcon = spm_vol(SPM.xCon(contrast_number).Vcon.fname);
% SPM.xCon(contrast_number).Vspm = spm_vol(SPM.xCon(contrast_number).Vspm.fname);
% con_data = spm_read_vols(SPM.xCon(contrast_number).Vcon);
% spm_data = spm_read_vols(SPM.xCon(contrast_number).Vspm);
% 
% % Apply the mask to the data
% masked_con_data = con_data .* mask;
% masked_spm_data = spm_data .* mask;
% 
% % Save the masked data
% V = SPM.xCon(contrast_number).Vspm;
% V.fname = strcat('D:\GraduationStudySPM\Emotion_Task_1st\',subjectID,'masked_contrast.nii');
% spm_write_vol(V, masked_spm_data);
% 
% % Visualize the masked results
% spm_image('Display', V.fname);
