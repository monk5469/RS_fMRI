%==========================================================================
% Script for fMRI analysis using SPM12
%
% Based on:
% Lee, S., et al. (2022). Effects of face repetition on ventral visual stream connectivity 
% using dynamic causal modelling of fMRI data. *NeuroImage*.
%
% This version includes the implementation of an F-contrast across all 
% experimental conditions to enhance the identification of task-related 
% activation. The aim is to facilitate more robust ROI definition for 
% subsequent effective connectivity analysis.
%==========================================================================

clear;

%% -------- Load workspace --------
data_dir = uigetdir(pwd, 'Select project directory');
load(fullfile(fileparts(data_dir), 'workspace_data.mat'));

conds = {'F_Init','F_Im','F_L','U_Init','U_Im','U_L','S_Init','S_Im','S_L'};
conditions = length(conds);
%% -------- 1st Level Analysis --------
for s = 1:length(objects_dir)

    ses_dir = sessions_dir{s};
    subj_target_dir = target_dir{s};
    output_dir = fullfile(subj_target_dir, 'bold', 'Imm_Del');

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    load(fullfile(subj_target_dir, 'bold', 'trial_onsets'));

    % Initialize batch
    matlabbatch = [];

    % Design specification
    matlabbatch{1}.spm.stats.fmri_spec.dir = {output_dir};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    onsets = cell(1, length(conds));
    R = []; scans = {}; time_so_far = 0; nscan = zeros(1, length(ses_dir));

    for i = 1:length(ses_dir)

        ses_path = fullfile(subj_target_dir, ses_dir(i).name);
        scans{i} = spm_select('FPList', ses_path, '^swar.*\.nii');
        rp_file = spm_select('FPList', ses_path, '^rp.*\.txt');

        R = [R; load(rp_file)];

        for c = 1:conditions
            onsets{c} = [onsets{c}; trial_onsets{1,i}{1,c} + time_so_far];
        end

        nscan(i) = length(scans{i});
        time_so_far = time_so_far + nscan(i) * 2;
    end

    % Save realignment regressors
    move_file = fullfile(output_dir, 'rp_all_sess.mat');
    save(move_file, 'R');

    % Set scans
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(strvcat(scans{:}));

    % Define conditions
    for c = 1:conditions
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name = conds{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset = onsets{c};
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod = struct('name', {}, 'param', {}, 'poly', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth = 1;
    end

    % Regressors and other params
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {move_file};
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

    % Run design
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);

    % Concatenate sessions
    spm_fmri_concatenate(fullfile(output_dir, 'SPM.mat'), nscan);
end

%% -------- Estimate (1st Level) --------
est_dir = 'Imm_Del'; sessrep = 'none';

for s = 1:length(objects_dir)
    cd(fullfile(target_dir{s}, 'bold', est_dir));
    load SPM;
    spm_spm(SPM);
end

%% -------- Define Contrasts (1st Level) --------
contrastsF.names = {'effects vs baseline', 'Imm_Del - repetition'};
contrastsF.weights = {
    eye(9)', ...
    [0 1 -1 0 0 0 0 0 0;
     0 0 0 0 1 -1 0 0 0;
     0 0 0 0 0 0 0 1 -1]
};

contrastsT.names = {'F_Init','F_Im','F_L','U_Init','U_Im','U_L','S_Init','S_Im','S_L', ...
                   'faces (F+U) > scrambled (S)', 'main effects of RS of faces'};
contrastsT.weights = {
    [1 0 0 0 0 0 0 0 0], [0 1 0 0 0 0 0 0 0], [0 0 1 0 0 0 0 0 0],
    [0 0 0 1 0 0 0 0 0], [0 0 0 0 1 0 0 0 0], [0 0 0 0 0 1 0 0 0],
    [0 0 0 0 0 0 1 0 0], [0 0 0 0 0 0 0 1 0], [0 0 0 0 0 0 0 0 1],
    [0.25 0.125 0.125 0.25 0.125 0.125 -0.5 -0.25 -0.25],
    [0.5 -0.25 -0.25 0.5 -0.25 -0.25 0 0 0]
};

for s = 1:length(objects_dir)
    spm_mat = spm_select('FPList', fullfile(target_dir{s}, 'bold', est_dir), '^SPM\.mat$');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr(spm_mat);

    for i = 1:length(contrastsF.names)
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.name = contrastsF.names{i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.convec = contrastsF.weights{i};
        matlabbatch{1}.spm.stats.con.consess{i}.fcon.sessrep = sessrep;
    end

    for i = 1:length(contrastsT.names)
        idx = i + length(contrastsF.names);
        matlabbatch{1}.spm.stats.con.consess{idx}.tcon.name = contrastsT.names{i};
        matlabbatch{1}.spm.stats.con.consess{idx}.tcon.convec = contrastsT.weights{i};
        matlabbatch{1}.spm.stats.con.consess{idx}.tcon.sessrep = sessrep;
    end

    matlabbatch{1}.spm.stats.con.delete = 1;
    spm_jobman('run', matlabbatch);
end

clear matlabbatch;
