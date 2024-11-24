 path='E:\GraduationStudySPM\Emotion_Task_Original_Files';
 fileNames=dir(path);
 fileCell=struct2cell(fileNames);
 names=fileCell(1,3:end);
 IDs=zeros(size(names));
 for idx=1:length(names)
     name=names(1,idx);
     splits=split(name, '_');
     IDs(1,idx)=str2double(splits{1,1});
 end
 save subjectIDs.mat IDs
 
 for idx=1:length(names)
      ID=num2str(IDs(1,idx));
      subDir='E:\GraduationStudySPM\Emotion_Task_LR_Files\';
      if exist(strcat(subDir,ID,'LR','.nii'),'file')==2
          %disp("ddd");
          continue
      else
          pathName=strcat('Emotion_Task_Original_Files','\',ID,'_3T_tfMRI_EMOTION_preproc','\',ID,'\','MNINonLinear\Results\tfMRI_EMOTION_LR');
          oldPath=cd(pathName);
          gunzip('tfMRI_EMOTION_LR.nii.gz');
          movefile('tfMRI_EMOTION_LR.nii', strcat(ID,'LR','.nii'))
          movefile(strcat(ID,'LR','.nii'), 'E:\GraduationStudySPM\Emotion_Task_LR_Files')
          cd(oldPath);
      end
  end
 
  for idx=1:length(names)
      ID=num2str(IDs(1,idx));
      subDir='E:\GraduationStudySPM\Emotion_Task_RL_Files\';
      if exist(strcat(subDir,ID,'RL','.nii'),'file')==2
          %disp("ddd");
          continue
      else
          pathName=strcat('Emotion_Task_Original_Files','\',ID,'_3T_tfMRI_EMOTION_preproc','\',ID,'\','MNINonLinear\Results\tfMRI_EMOTION_RL');
          oldPath=cd(pathName);
          gunzip('tfMRI_EMOTION_RL.nii.gz');
          movefile('tfMRI_EMOTION_RL.nii', strcat(ID,'RL','.nii'))
          movefile(strcat(ID,'RL','.nii'), 'E:\GraduationStudySPM\Emotion_Task_RL_Files')
          cd(oldPath);
      end
  end
 
 spm('defaults', 'FMRI');
 spm_jobman('initcfg');
 
for idx=1:100
    subjectID=num2str(IDs(1,idx));
    if exist(strcat('E:\GraduationStudySPM\Emotion_Task_1st\',subjectID,'\SPM.mat'),'file')
        continue
    end

    disp(strcat("Subject #",subjectID));
    mkdir(strcat("Emotion_Task_1st\",subjectID));
    matlabbatch=cell(1);
    
    numScans=176;
    onsetIdx=1:3:7;
    %duration=2.*ones(1,18);
    duration=zeros(1,18);
    %% SPM 1st level analysis - Model Specification

    matlabbatch{1}.spm.stats.fmri_spec.dir = {strcat('E:\GraduationStudySPM\Emotion_Task_1st\',subjectID)};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 0.72;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    
    
    %% Model Specification - session (scan) #1 - LR
    niiPathLR=strcat('E:\GraduationStudySPM\Emotion_Task_LR_Files\');
    onsetPathLR=strcat('E:\GraduationStudySPM\Emotion_Task_Original_Files\',subjectID,'_3T_tfMRI_EMOTION_preproc\',subjectID,'\MNINonLinear\Results\tfMRI_EMOTION_LR\EVs');
    regPathLR=strcat('E:\GraduationStudySPM\Emotion_Task_Original_Files\',subjectID,'_3T_tfMRI_Emotion_preproc\',subjectID,'\MNINonLinear\Results\tfMRI_EMOTION_LR');
    
    %% Session #1 - retrieving scan files
    filesLR=spm_select('ExtFPList', niiPathLR, ['^', subjectID],1:numScans);
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = cellstr(filesLR);
    
    
    %% Session #1 - Condition #1 - Faces
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
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).name = 'Faces';
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).onset = fearOnsetLR; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).duration = duration; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(1).orth = 1;
    
    %% Session #1 - Condition #2 - Shapes
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
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).name = 'Shapes';
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).onset = neutOnsetLR;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).duration = duration; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond(2).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
    
    %% Session #1 - multiple regressors
    
    %regFileLR=fopen(strcat(regPathLR,'\Movement_Regressors.txt'));
    %regFileIDLR=textscan(regFileLR, '%f');
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {strcat(regPathLR,'\Movement_Regressors.txt')};
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
    
    
    %% Model Specification - session (scan) #2 - RL
    niiPathRL=strcat('E:\GraduationStudySPM\Emotion_Task_RL_Files\');
    onsetPathRL=strcat('E:\GraduationStudySPM\Emotion_Task_Original_Files\',subjectID,'_3T_tfMRI_EMOTION_preproc\',subjectID,'\MNINonLinear\Results\tfMRI_EMOTION_RL\EVs');
    regPathRL=strcat('E:\GraduationStudySPM\Emotion_Task_Original_Files\',subjectID,'_3T_tfMRI_Emotion_preproc\',subjectID,'\MNINonLinear\Results\tfMRI_EMOTION_RL');
    
    %% Session #2 - retrieving scan files
    filesRL=spm_select('ExtFPList', niiPathRL, ['^', subjectID],1:numScans);
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans = cellstr(filesRL);
    
    %% Session #2 - Condition #1 - Faces
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
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).name = 'Faces';
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).onset = fearOnsetRL; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).duration = duration; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(1).orth = 1;
    
    
    %% Session #2 - Condition #2 - Shapes
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
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).name = 'Shapes';
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).onset = neutOnsetRL;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).duration = duration; 
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).tmod = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond(2).orth = 1;
    
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
    
    %% Session #1 - multiple regressors
    %regFileIDRL=textscan(regFileRL, '%f');
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {strcat(regPathRL,'\Movement_Regressors.txt')};
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;
    
    zipMaskPath=strcat(regPathRL,'\brainmask_fs.2.nii.gz');
    gunzip(zipMaskPath);
    matlabbatch{1}.spm.stats.fmri_spec.mask={strcat(regPathRL,'\brainmask_fs.2.nii')};
    
    %% Model Estimation
    matlabbatch{2}.spm.stats.fmri_est.spmmat = {strcat('E:\GraduationStudySPM\Emotion_Task_1st\',subjectID,'\SPM.mat')};
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
    
    spm_jobman('run', matlabbatch);
end

