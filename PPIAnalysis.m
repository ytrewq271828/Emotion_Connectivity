% Step 1: Initialize SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

IDCell=load('subjectIDs.mat');
IDs=IDCell.IDs;
numScans=176;
duration=2.*ones(1,36);

for idx=1:100
    ID=num2str(IDs(1,idx));
    if exist(strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\SPM.mat'),'file')
        continue
    end

    disp(strcat("Subject #",ID));
    mkdir(strcat("PPIAnalysis\combined\",ID));
    matlabbatch=cell(1);

    niiPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_LR_Files');
    niiPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_RL_Files');
    filesLR=spm_select('ExtFPList', niiPathLR, ['^', ID], 1:numScans);
    filesRL=spm_select('ExtFPList', niiPathRL, ['^', ID], 1:numScans);
    filesPPI=cellstr([filesLR;filesRL]);

    regPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_Original_Files\',ID,'_3T_tfMRI_Emotion_preproc\',ID,'\MNINonLinear\Results\tfMRI_EMOTION_LR');
    regPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_Original_Files\',ID,'_3T_tfMRI_Emotion_preproc\',ID,'\MNINonLinear\Results\tfMRI_EMOTION_RL');
    regLR = strcat(regPathLR,'\Movement_Regressors.txt');
    regRL = strcat(regPathRL,'\Movement_Regressors.txt');
    regLRMat=readmatrix(regLR);
    regRLMat=readmatrix(regRL);
    writematrix(cat(1,regLRMat,regRLMat),'All_Regressors.txt', 'Delimiter', '\t');

    runDuration=0.72*176;
    onsetPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_Original_Files\',ID,'_3T_tfMRI_EMOTION_preproc\',ID,'\MNINonLinear\Results\tfMRI_EMOTION_RL\EVs');
    onsetPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_Original_Files\',ID,'_3T_tfMRI_EMOTION_preproc\',ID,'\MNINonLinear\Results\tfMRI_EMOTION_RL\EVs');
    fearIDLR=fopen(strcat(onsetPathLR,'\fear.txt'));
    fearCellLR=textscan(fearIDLR,'%f');
    fearMatLR=fearCellLR{1,1}; 
    fearOnsetLR=zeros(1,18);
    for i=0:2
        veryOnset=fearMatLR(3*i+1);
        for j=1:6
            fearOnsetLR(1,6*i+j)=veryOnset+3*(j-1)+3;
        end
    end    
    fearIDRL=fopen(strcat(onsetPathRL,'\fear.txt'));
    fearCellRL=textscan(fearIDRL,'%f');
    fearMatRL=fearCellRL{1,1};
    fearOnsetRL=zeros(1,18);
    for i=0:2
        veryOnset=fearMatRL(3*i+1);
        for j=1:6
            fearOnsetRL(1,6*i+j)=veryOnset+3*(j-1)+3;
        end
    end
    fearOnsetAll=[fearOnsetLR fearOnsetRL];

    neutIDLR=fopen(strcat(onsetPathLR,'\neut.txt'));
    neutCellLR=textscan(neutIDLR,'%f');
    neutMatLR=neutCellLR{1,1};
    neutOnsetLR=zeros(1,18);
    for i=0:2
        veryOnset=neutMatLR(3*i+1);
        for j=1:6
            neutOnsetLR(1,6*i+j)=veryOnset+3*(j-1)+3;
        end
    end
    neutIDRL=fopen(strcat(onsetPathRL,'\neut.txt'));
    neutCellRL=textscan(neutIDRL,'%f');
    neutMatRL=neutCellRL{1,1};
    neutOnsetRL=zeros(1,18);
    for i=0:2
        veryOnset=neutMatRL(3*i+1);
        for j=1:6
            neutOnsetRL(1,6*i+j)=veryOnset+3*(j-1)+3;
        end
    end
    neutOnsetAll=[neutOnsetLR neutOnsetRL];


    matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = filesPPI;

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = 'Faces';
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = fearOnsetAll; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = duration; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).name = 'Shapes';
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).onset = neutOnsetAll;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).duration = duration; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;

    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {'D:\GraduationStudySPM\All_Regressors.txt'};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;

    zipMaskPath=strcat(regPathRL, '\brainmask_fs.2.nii.gz');
    gunzip(zipMaskPath);
    matlabbatch{1}.spm.stats.fmri_spec.mask={strcat(regPathRL,'\brainmask_fs.2.nii')};


    matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


    contrast1.name='Faces';
    contrast1.weights=[1 0];
    contrast1.sessrep='none';

    contrast2.name='Shapes';
    contrast2.weights=[0 1];
    contrast2.sessrep='none';

    contrast3.name='Faces - Shapes';
    contrast3.weights=[1 -1];
    contrast3.sessrep='none';

    contrast4.name='Effect of Interest';
    contrast4.weights=[1 0; 0 1];

    matlabbatch{3}.spm.stats.con.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon=contrast1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon=contrast2;
    matlabbatch{3}.spm.stats.con.consess{3}.tcon=contrast3;
    matlabbatch{3}.spm.stats.con.consess{4}.fcon=contrast4;
    matlabbatch{3}.spm.stats.con.delete = 1;


    spm_jobman('run', matlabbatch);
end




for idx=1:100

    ID=num2str(IDs(1,idx));
    %mkdir(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID))
    spm_mat_file=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\SPM.mat');
    voi_file=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\VOI_Amygdala_L_1.mat');
    if exist(voi_file,'file')
        continue
    end

    disp(strcat("Subject #",ID));
    matlabbatch=cell(1);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr(spm_mat_file);
    matlabbatch{1}.spm.util.voi.adjust = 4; % Effects of interest contrast number 
    matlabbatch{1}.spm.util.voi.session = 1; % Session index 
    matlabbatch{1}.spm.util.voi.name = 'Amygdala_L'; % VOI name

    matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''}; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1; % Index of contrast for choosing voxels 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none'; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.01; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});

    % Define large fixed outer sphere 
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [-26 -6 -16]; 
    % Set coordinates here 
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10; % Radius (mm) 
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;

    % Define smaller inner sphere which jumps to the peak of the outer sphere 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = [0 0 0]; % Leave this at zero 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 6; % Set radius here (mm) 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1; % Index of SPM within the batch 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2'; % Index of the outer sphere within the batch

    % Include voxels in the thresholded SPM (i1) and the mobile inner sphere (i3) 
    matlabbatch{1}.spm.util.voi.expression = 'i1 & i3';

    % Run the batch 
    spm_jobman('run',matlabbatch); 

end

%% Amygdala PPI File Generation
cd D:\GraduationStudySPM
mkdir(strcat('PPIAnalysis\','PPI_Faces'))
marsbar('on')
for i=1:100
    ID=num2str(IDs(1,i));
    if exist(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Amygdala_Faces.mat'),'file')
        continue
    end

    mkdir(strcat('PPIAnalysis\PPI_Faces\',ID));
    saveDir=strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\', ID);
    analysis_dir = strcat('D:\GraduationStudySPM\PPIAnalysis\combined\', ID);
    if ~exist(strcat(analysis_dir,'\VOI_Amygdala_L_1.mat'),"file")
        continue
    end
    disp(ID)


    weights=[1 1 1; 2 1 0]; % Referring to SPM.Sess.U - Faces vs. Baseline

    SPM_VOI=load(fullfile(analysis_dir, 'SPM.mat'));
    VOIFile=load(fullfile(analysis_dir, 'VOI_Amygdala_L_1.mat'));

    PPI = spm_peb_ppi(SPM_VOI.SPM, 'ppi', VOIFile.xY, weights, 'PPI_Amygdala_Faces', 1);
    save(fullfile(saveDir, 'PPI_Amygdala_Faces.mat'), 'PPI')

    niiPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_LR_Files');
    niiPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_RL_Files');
    filesLR=spm_select('ExtFPList', niiPathLR, ['^', ID], 1:numScans);
    filesRL=spm_select('ExtFPList', niiPathRL, ['^', ID], 1:numScans);
    filesPPI=cellstr([filesLR;filesRL]);

    matlabbatch = cell(1);
    matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; 
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = filesPPI;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name=PPI.name;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val=PPI.ppi;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name='Amygdala_L_BOLD';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val=PPI.Y;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name='Faces only';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val=PPI.P;

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {'D:\GraduationStudySPM\All_Regressors.txt'};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

    matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    contrast.name='PPI_Amygdala_Faces';
    contrast.weights=[1 zeros(1,14)];
    contrast.sessrep='none';

    matlabbatch{3}.spm.stats.con.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon=contrast;
    matlabbatch{3}.spm.stats.con.delete = 1;

    spm_jobman('run', matlabbatch);

    movefile(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat'), strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Amygdala_Faces.mat'));
end



%% 2nd-level analysis with amygdala VOI PPI results
mkdir PPIAnalysis\PPI_Faces\2nd
cd PPIAnalysis\PPI_Faces\2nd
if(~exist('contrasts.mat','file'))
    contrasts=cell(1,100);
    idx=1;
    cIdx=0;
    while idx<101
        if ~exist(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',num2str(IDs(1,idx)),'\SPM.mat'),'file')
            continue
        else
            cIdx=cIdx+1;
            contrasts{1,cIdx}=strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',num2str(IDs(1,idx)),'\con_0001.nii,1');
        end
        idx=idx+1;
    end
    save contrasts.mat contrasts
else
    contrasts=load('contrasts.mat');
    contrasts=contrasts.contrasts;
end


matlabbatch = [];
matlabbatch{1}.spm.stats.factorial_design.dir = {'D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\2nd\'};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = cellstr(contrasts(1,:))';
matlabbatch{1}.spm.stats.factorial_design.cov = struct([]); % No covariates for 1-sample t-test
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;


% Load the SPM.mat file from the second-level analysis directory
spmMatPath = fullfile('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\2nd\', 'SPM.mat');

matlabbatch{2}.spm.stats.fmri_est.spmmat = {spmMatPath};


matlabbatch{3}.spm.stats.con.spmmat = {spmMatPath};


% Define the T-contrast
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Faces only';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1; % Already generated a contrast in 1st-level
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% Run the contrast manager
spm_jobman('run', matlabbatch);


%% Generating a VOI that survives the 2nd-level PPI Analysis
for idx=1:100

    ID=num2str(IDs(1,idx));
    %mkdir(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID))
    spm_mat_file=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\SPM.mat');
    voi_file_INS=strcat('D:\GraduationStudySPM\PPIAnalysis\combined\',ID,'\VOI_Insula_L_1.mat');
    if exist(voi_file_INS,'file')
        continue
    end

    disp(strcat("Subject #",ID));
    matlabbatch=cell(1);
    matlabbatch{1}.spm.util.voi.spmmat = cellstr(spm_mat_file);
    matlabbatch{1}.spm.util.voi.adjust = 4; % Effects of interest contrast number 
    matlabbatch{1}.spm.util.voi.session = 1; % Session index 
    matlabbatch{1}.spm.util.voi.name = 'Insula_L'; % VOI name

    matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat = {''}; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast = 1; % Index of contrast for choosing voxels 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc = 'none'; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh = 0.01; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.extent = 0; 
    matlabbatch{1}.spm.util.voi.roi{1}.spm.mask = struct('contrast', {}, 'thresh', {}, 'mtype', {});

    % Define large fixed outer sphere 
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre = [-34 18 -2]; 
    % Set coordinates here 
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius = 10; % Radius (mm) 
    matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;

    % Define smaller inner sphere which jumps to the peak of the outer sphere 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre = [0 0 0]; % Leave this at zero 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius = 6; % Set radius here (mm) 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm = 1; % Index of SPM within the batch 
    matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2'; % Index of the outer sphere within the batch

    % Include voxels in the thresholded SPM (i1) and the mobile inner sphere (i3) 
    matlabbatch{1}.spm.util.voi.expression = 'i1 & i3';

    % Run the batch 
    spm_jobman('run',matlabbatch); 

end

%% Creating 4 PPIs for Plotting - V2/V5 & Faces/Baseline

% Amygdala - Baseline
cd D:\GraduationStudySPM
%mkdir(strcat('PPIAnalysis\','PPI_Faces'))
marsbar('on')
for i=1:100
    ID=num2str(IDs(1,i));
    if exist(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Amygdala_Baseline.mat'),'file')
        continue
    end

    %mkdir(strcat('PPIAnalysis\PPI_Faces\',ID));
    saveDir=strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\', ID);
    analysis_dir = strcat('D:\GraduationStudySPM\PPIAnalysis\combined\', ID);
    if ~exist(strcat(analysis_dir,'\VOI_Amygdala_L_1.mat'),"file")
        continue
    end
    disp(ID)


    weights=[1 1 0; 2 1 0]; % Baseline effect - 3rd column is all zero

    SPM_VOI=load(fullfile(analysis_dir, 'SPM.mat'));
    VOIFile=load(fullfile(analysis_dir, 'VOI_Amygdala_L_1.mat'));

    PPI = spm_peb_ppi(SPM_VOI.SPM, 'ppi', VOIFile.xY, weights, 'PPI_Amygdala_Baseline', 1);
    save(fullfile(saveDir, 'PPI_Amygdala_Baseline.mat'), 'PPI')

    niiPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_LR_Files');
    niiPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_RL_Files');
    filesLR=spm_select('ExtFPList', niiPathLR, ['^', ID], 1:numScans);
    filesRL=spm_select('ExtFPList', niiPathRL, ['^', ID], 1:numScans);
    filesPPI=cellstr([filesLR;filesRL]);

    matlabbatch = cell(1);
    matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; 
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = filesPPI;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name=PPI.name;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val=PPI.ppi;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name='Amygdala_L_BOLD';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val=PPI.Y;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name='Baseline';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val=PPI.P;

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {'D:\GraduationStudySPM\All_Regressors.txt'};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

    matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    contrast.name='PPI_Amygdala_Baseline';
    contrast.weights=[0 zeros(1,14)];
    contrast.sessrep='none';

    matlabbatch{3}.spm.stats.con.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon=contrast;
    matlabbatch{3}.spm.stats.con.delete = 1;

    spm_jobman('run', matlabbatch);

    movefile(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat'), strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Amygdala_Baseline.mat'));
end

% Insula - Faces
for i=1:100
    ID=num2str(IDs(1,i));
    if exist(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Insula_Faces.mat'),'file')
        continue
    end

    %mkdir(strcat('PPIAnalysis\PPI_Faces\',ID));
    saveDir=strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\', ID);
    analysis_dir = strcat('D:\GraduationStudySPM\PPIAnalysis\combined\', ID);
    if ~exist(strcat(analysis_dir,'\VOI_Insula_L_1.mat'),"file")
        continue
    end
    disp(ID)


    weights=[1 1 1; 2 1 0]; % Faces Effect

    SPM_VOI=load(fullfile(analysis_dir, 'SPM.mat'));
    VOIFile=load(fullfile(analysis_dir, 'VOI_Insula_L_1.mat'));

    PPI = spm_peb_ppi(SPM_VOI.SPM, 'ppi', VOIFile.xY, weights, 'PPI_Insula_Faces', 1);
    save(fullfile(saveDir, 'PPI_Insula_Faces.mat'), 'PPI')

    niiPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_LR_Files');
    niiPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_RL_Files');
    filesLR=spm_select('ExtFPList', niiPathLR, ['^', ID], 1:numScans);
    filesRL=spm_select('ExtFPList', niiPathRL, ['^', ID], 1:numScans);
    filesPPI=cellstr([filesLR;filesRL]);

    matlabbatch = cell(1);
    matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; 
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = filesPPI;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name=PPI.name;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val=PPI.ppi;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name='Insula_L_BOLD';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val=PPI.Y;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name='Faces only';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val=PPI.P;

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {'D:\GraduationStudySPM\All_Regressors.txt'};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

    matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    contrast.name='PPI_Insula_Faces';
    contrast.weights=[1 zeros(1,14)];
    contrast.sessrep='none';

    matlabbatch{3}.spm.stats.con.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon=contrast;
    matlabbatch{3}.spm.stats.con.delete = 1;

    spm_jobman('run', matlabbatch);

    movefile(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat'), strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Insula_Faces.mat'));

end

% Insula - Baseline
for i=1:100
    ID=num2str(IDs(1,i));
    if exist(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Insula_Baseline.mat'),'file')
        continue
    end

    %mkdir(strcat('PPIAnalysis\PPI_Faces\',ID));
    saveDir=strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\', ID);
    analysis_dir = strcat('D:\GraduationStudySPM\PPIAnalysis\combined\', ID);
    if ~exist(strcat(analysis_dir,'\VOI_Insula_L_1.mat'),"file")
        continue
    end
    disp(ID)


    weights=[1 1 0; 2 1 0]; % Baseline effect - 3rd column is all zero

    SPM_VOI=load(fullfile(analysis_dir, 'SPM.mat'));
    VOIFile=load(fullfile(analysis_dir, 'VOI_Insula_L_1.mat'));

    PPI = spm_peb_ppi(SPM_VOI.SPM, 'ppi', VOIFile.xY, weights, 'PPI_Insula_Baseline', 1);
    save(fullfile(saveDir, 'PPI_Insula_Baseline.mat'), 'PPI')

    niiPathLR=strcat('D:\GraduationStudySPM\Emotion_Task_LR_Files');
    niiPathRL=strcat('D:\GraduationStudySPM\Emotion_Task_RL_Files');
    filesLR=spm_select('ExtFPList', niiPathLR, ['^', ID], 1:numScans);
    filesRL=spm_select('ExtFPList', niiPathRL, ['^', ID], 1:numScans);
    filesPPI=cellstr([filesLR;filesRL]);

    matlabbatch = cell(1);
    matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72; 
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = filesPPI;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name=PPI.name;
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val=PPI.ppi;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name='Insula_L_BOLD';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val=PPI.Y;

    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name='Faces only';
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val=PPI.P;

    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {'D:\GraduationStudySPM\All_Regressors.txt'};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

    matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

    contrast.name='PPI_Insula_Baseline';
    contrast.weights=[0 zeros(1,14)];
    contrast.sessrep='none';

    matlabbatch{3}.spm.stats.con.spmmat = {strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat')};
    matlabbatch{3}.spm.stats.con.consess{1}.tcon=contrast;
    matlabbatch{3}.spm.stats.con.delete = 1;

    spm_jobman('run', matlabbatch);

    movefile(strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM.mat'), strcat('D:\GraduationStudySPM\PPIAnalysis\PPI_Faces\',ID,'\SPM_Insula_Baseline.mat'));

end