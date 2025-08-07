%==========================================================================
% Script for fMRI DCM analysis using SPM12
%
% Based on:
% Lee, S., et al. (2022). Effects of face repetition on ventral visual 
% stream connectivity using dynamic causal modelling of fMRI data. *NeuroImage*.
%
% This version includes a key modification to the stimulus input matrix (B-matrix)
% to examine modulatory effects across multiple repetition and familiarity conditions.
%==========================================================================

% addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest; 
% addpath /imaging/henson/users/rh01/SPM/Rik_Updates; % for updated version of spm_regions that handles NaN voxels, and for version of spm_dcm_peb_fit that allows parfor
% addpath /imaging/henson/users/sl01/FacesData/Scripts;
clear
data_dir = uigetdir(pwd, 'Select project directory');
load([fileparts(data_dir) filesep 'workspace_data.mat']);
% rwd     = '/imaging/henson/users/sl01/FacesData';
% owd     = '/imaging/henson/users/sl01/FacesData';
% dosubs  = [1 2 3 5 6 8 9 10 12 14 15 16 17 18 19 23 24 25]; % all except s11 for debrief
objects_dir([4:6,11,12]) = [];
target_dir([4:6,11,12]) = [];
% stat_dir = 'stat_concat_deb';
% DCM_dir  = 'DCM';
% ROI_dir  = 'ROI_time_series';

stat_dir = 'Imm_Del';
DCM_dir  = 'Imm_Del_DCM';
ROI_dir  = 'Imm_Del_ROI_adjust_eye(9)_subs_11';
Nses = 1;
sed_dir = sessions_dir{1};
%% DCM
ROInames = {'rEVC';'rOFA';'rFFA';'lEVC';'lOFA';'lFFA'};
VOInames = {};

for v = 1:length(ROInames)
    for ses = 1:Nses
        VOInames{v,ses} = ['VOI_' ROInames{v} sprintf('_%d', ses)];
    end
end

scaleU  = 0; % Whether to scale (Z-score) B inputs - suggestion of Peter
% Only contrasts that show sig effect (2-tailed) below
% contrasts = [1 1 1  1 1 1  1 1 1;
%              1 1 1  1 1 1  0 0 0; % faces
%              0 1 0  0 1 0  0 1 0; % immediate repetition
%              0 0 1  0 0 1  0 0 1; % delayed repetition            
%              1 1 1  0 0 0  0 0 0]; % fame

contrasts = [1 1 1  1 1 1  1 1 1;
             0 1 0  0 0 0  0 0 0;  % familiar-immediate
             0 0 1  0 0 0  0 0 0;  % familiar-delayed
             0 0 0  0 1 0  0 0 0;  % unfamiliar-immediate
             0 0 0  0 0 1  0 0 0;  % unfamiliar-delayed
             0 0 0  0 0 0  0 1 0;  % scrambled-immediate
             0 0 0  0 0 0  0 0 1]; % scrambled-delayed
Ncon = size(contrasts,1);

%% 2 ROI - right OFA & FFA
VOIs = [2:3], Nvoi = length(VOIs), dirname = sprintf('Right%d_adjust_eye(9)_subs_11_group6_12R',Nvoi)
% VOIs = [5:6], Nvoi = length(VOIs), dirname = sprintf('Left%d',Nvoi)
Model(1).a = ones(Nvoi); % fully connected
Model(1).b(:,:,1) = zeros(Nvoi); % no modulation for the first input
for n = 2:Ncon
    Model(1).b(:,:,n) = ones(Nvoi); % modulations for the four remaining inputs/contrasts of interest
end
Model(1).c = zeros(Nvoi,Ncon);
Model(1).c(:,1) = 1; % inputs to OFA, FFA
%% 3 ROI - right EVC, OFA & FFA
% VOIs = [1:3], Nvoi = length(VOIs), dirname = sprintf('Right%d',Nvoi)
% % VOIs = [4:6], Nvoi = length(VOIs), dirname = sprintf('Left%d',Nvoi)
% Model(1).a = ones(Nvoi); % fully connected
% Model(1).b(:,:,1) = zeros(Nvoi); % no modulation for the first input
% for n = 2:Ncon
%     Model(1).b(:,:,n) = ones(Nvoi); % modulations for the four remaining inputs/contrasts of interest
% end
% Model(1).c = zeros(Nvoi,Ncon);
% Model(1).c(1,1) = 1; % input to EVC only

%% 6 ROI - bilateral EVC, OFA & FFA
% VOIs = [1:6], Nvoi = length(VOIs), dirname = sprintf('Bilateral%d',Nvoi)
% Model(1).a = [1 1 1 0 0 0; % rEVC, rOFA, rFFA, lEVC, lOFA, lFFA
%               1 1 1 0 1 0;
%               1 1 1 0 0 1;
%               0 0 0 1 1 1;
%               0 1 0 1 1 1;
%               0 0 1 1 1 1];
% Model(1).b = zeros(Nvoi); % no modulation for the first input
% for b = 2:Ncon
%     Model(1).b(:,:,b) = Model(1).a; % modulations for the four remaining inputs/contrasts of interest
% end
% Model(1).c = zeros(Nvoi,Ncon);
% Model(1).c(1,1) = 1; Model(1).c(4,1) = 1; % inputs to rEVC & lEVC

%% Collate DCMs into GCM
GCM = {}; nsubses = 0;
for s=1:length(objects_dir)

    % sub = dosubs(s);
    % swd = sprintf('subject_%02d',sub);
    data_path = fullfile(target_dir{s},'bold',stat_dir);
    out_path  = fullfile(target_dir{s},'bold',DCM_dir);
    if ~exist(out_path) 
        mkdir(out_path);
    end
    out_path  = fullfile(out_path,dirname);
    if ~exist(out_path) 
        mkdir(out_path);
    end

    load(fullfile(data_path,'SPM.mat'));

    for ses=1:Nses

        for v=1:Nvoi
            check_path = fullfile(target_dir{s},'bold',ROI_dir);
            file_check = fullfile(check_path,[VOInames{VOIs(v),ses} '.mat']);
            if ~exist(file_check,'file')
                continue;
            end
            load(fullfile(target_dir{s},'bold',ROI_dir,VOInames{VOIs(v),ses}),'xY');
            DCM_allM.xY(ses,v) = xY;
        end

        nsubses = nsubses+1;

        DCM = [];

        for v=1:Nvoi
            DCM.xY(v) = DCM_allM.xY(ses,v);
        end

        DCM.n = length(DCM.xY);      % number of regions
        DCM.v = length(DCM.xY(1).u); % number of time points

        DCM.Y.dt  = SPM.xY.RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        max_len = max(arrayfun(@(x) length(x.u), DCM.xY));  % 找出所有 u 向量的最大长度
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end

        DCM.Y.Q    = spm_Ce(ones(1,DCM.n)*DCM.v); % Models autocorrelation

        DCM.U.dt   = SPM.Sess(ses).U(1).dt;

        for ui = 1:size(contrasts,1)
            ind = find(contrasts(ui,:)~=0);

            tmp = sparse(length(SPM.Sess(ses).U(1).u)-32,1);
            for c = 1:length(ind)
                tmp = tmp + contrasts(ui,ind(c)) * SPM.Sess(ses).U(ind(c)).u(33:end);
            end

            if scaleU & ui>1 % Don't scale first modulation (ui=1) because in C? (Peter)
                tmp = tmp/std(tmp);
            end

            DCM.U.u(:,ui) = tmp;

            tmp = SPM.Sess(ses).U(ind(1)).name{1};
            for c = 2:length(ind)
                 if contrasts(ui,ind(c)) == 1, sym = '+'; else sym = '-'; end
                 tmp = [tmp sym SPM.Sess(ses).U(ind(c)).name{1}];
            end
            DCM.U.name{ui} = tmp;  

            %DCM.U.idx(ui,:) = [ind zeros(Ncon-length(ind))]; % not sure this necessary
        end

        DCM.delays = repmat(SPM.xY.RT,Nvoi,1)/2;
        DCM.TE     = 0.03;  % 30ms on CBU

        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic = 0;
        DCM.options.centre     = 1; % centre regressors (default = 1)
        DCM.options.nograph    = 1;
        DCM.options.maxnodes   = 8;
        DCM.options.maxit      = 128; %32;
        DCM.options.hidden     = [];
        DCM.options.induced    = 0;

        DCM.M.options = struct(); % needed else crashes in estimation

        for m = 1:length(Model)        
            DCM.a        = Model(m).a;
            DCM.b        = Model(m).b;
            DCM.c        = Model(m).c;
            DCM.d        = zeros(DCM.n,DCM.n,0); % needed else crashes in estimation

            filename = fullfile(out_path,sprintf('DCM_mod%d_ses%d.mat',m,ses));
             save(filename,'DCM');

            GCM{nsubses,m} = filename;            
        end
    end
end

% cd(owd)
PEBwd = fullfile(data_dir,'Imm_Del_PEB',dirname)
try mkdir(PEBwd); end
cd(PEBwd)
save('GCM_defined','GCM','-v7.3');

%% Estimate each subject's model and group PEB model
M=[]; M.X = ones(size(GCM,1),1);
M.Q = 'single';
use_parfor = true; 
% [GCM,PEB] = region_spm_dcm_peb_fit(GCM(:,1),M,{'A','B','C'},use_parfor); % need to have parfor-edited version of this function on path...
% save('GCM_peb_fitABC','GCM','PEB')
[GCM,PEB] = region_spm_dcm_peb_fit(GCM(:,1),M,{'A','B'},use_parfor); % need to have parfor-edited version of this function on path...
save('GCM_peb_fitAB','GCM','PEB')

%% BMR
BMA = spm_dcm_peb_bmc(PEB);
save('BMA_AB','BMA');

%% BMC
vals = [0 0; 0 1; 1 0; 1 1];
Nval = size(vals,1);

bGCM = rmfield(GCM{1,1},'M');

%% 2 ROIs
% 2 family comparisons: self, OFA<->FFA
Nfam = 6;
FamP = nan(Nfam,Ncon-1);
final_PF = [];
BMA_all = cell(1, Ncon-1); 
for con = 2:Ncon
    GCMs = cell(1,Nval^2);
    m = 1; % First model is full model
    GCMs{1,m} = bGCM;
    
    family = nan(Nfam,Nval^2);
    family(:,1) = ones(Nfam,1)*2;

    for self = 1:Nval
        for OFAFFA = 1:Nval
            if m < Nval^2 % Ensuring last model is not repeat of first (full) model
                m=m+1;
                GCMs{1,m} = bGCM;
                GCMs{1,m}.b(1,1,con) = vals(self,1); GCMs{1,m}.b(2,2,con) = vals(self,2);
                GCMs{1,m}.b(2,1,con) = vals(OFAFFA,1); GCMs{1,m}.b(1,2,con) = vals(OFAFFA,2);
                
                family(1,m) = max(vals(self,:))+1; % Any self
                family(2,m) = max(vals(OFAFFA,:))+1; % Any OFA<->FFA
                family(3,m) = vals(self,1)+1; % OFA self
                family(4,m) = vals(self,2)+1; % FFA self
                family(5,m) = vals(OFAFFA,1)+1; % OFA->FFA
                family(6,m) = vals(OFAFFA,2)+1; % FFA->OFA
            end
        end
    end

    [BMA,BMR] = spm_dcm_peb_bmc(PEB, GCMs);
    BMA_all{con-1}   = BMA;
    for f = 1:Nfam
        [BMA_fam_task,fam_task,final_PF] = spm_dcm_peb_bmc_fam(BMA, BMR, family(f,:), 'ALL',final_PF);
        fprintf('\ncontrast %d, family %d\n',con-1,f)
        fam_task.family.post
        FamP(f,con-1) = fam_task.family.post(2);
        %pause
    end
end
round(FamP,2)

% % 2 family comparisons: input OFA, input FFA
% vals = [1 1; 1 0; 0 1; 0 0];
% FamP = nan(1,4);
% GCMs = {};
% GCMs{1,1} = bGCM; % First model is full model
% 
% family = nan(2,4);
% family(:,1) = [1 1]*2;
% for m = 2:Nval
%     GCMs{1,m} = bGCM;
%     GCMs{1,m}.c(1,1) = vals(m,1); GCMs{1,m}.c(2,1) = vals(m,2);
%     %GCMs{1,m}.b(1,2,2) = vals(m,1); GCMs{1,m}.b(2,1,2) = vals(m,2);
% 
%     family(1,m) = vals(m,1)+1; % Input OFA
%     family(2,m) = vals(m,2)+1; % Input FFA
% end
% 
% [BMA,BMR] = spm_dcm_peb_bmc(PEB, GCMs);
% 
% FamP = nan(1,2);
% for f = 1:size(family,1)
%     [BMA_fam_task,fam_task] = spm_dcm_peb_bmc_fam(BMA, BMR, family(f,:), 'ALL');
%     fprintf('family %d\n',f)
%     FamP(f) = fam_task.family.post(2);
% end
% round(FamP,2)

%% 3 ROIs
% 4 family comparisons: self, EVC->OFA/FFA, OFA<->FFA, EVC<-OFA/FFA
% Nfam = 12;
% FamP = nan(Nfam,Ncon-1);
% for con = 2:Ncon
%     GCMs = cell(1,Nval^4);
%     m = 1; % First model is full model
%     GCMs{1,m} = bGCM;
% 
%     family = nan(Nfam,Nval^4);
%     family(:,1) = ones(Nfam,1)*2;
% 
%     for self = 1:Nval
%         for forw = 1:Nval
%             for OFAFFA = 1:Nval
%                 for backw = 1:Nval
%                         if m < Nval^4 % Ensuring last model is not repeat of first (full) model
%                             m=m+1;
%                             GCMs{1,m} = bGCM;
%                             GCMs{1,m}.b(2,2,con) = vals(self,1); GCMs{1,m}.b(3,3,con) = vals(self,2);
%                             GCMs{1,m}.b(2,1,con) = vals(forw,1); GCMs{1,m}.b(3,1,con) = vals(forw,2);
%                             GCMs{1,m}.b(1,2,con) = vals(backw,1); GCMs{1,m}.b(1,3,con) = vals(backw,2);
%                             GCMs{1,m}.b(3,2,con) = vals(OFAFFA,1); GCMs{1,m}.b(2,3,con) = vals(OFAFFA,2);
% 
%                             family(1,m) = max(vals(self,:))+1; % Any self
%                             family(2,m) = max(vals(forw,:))+1; % Any forward from EVC
%                             family(3,m) = max(vals(OFAFFA,:))+1; % Any OFA<->FFA
%                             family(4,m) = max(vals(backw,:))+1; % Any backward to EVC
%                             family(5,m) = vals(self,1)+1; % OFA self
%                             family(6,m) = vals(self,2)+1; % FFA self
%                             family(7,m) = vals(OFAFFA,1)+1; % OFA->FFA
%                             family(8,m) = vals(OFAFFA,2)+1; % OFA<-FFA
%                             family(9,m) = vals(forw,1)+1; % EVC->OFA
%                             family(10,m) = vals(forw,2)+1; % EVC->FFA
%                             family(11,m) = vals(backw,1)+1; % EVC<-OFA
%                             family(12,m) = vals(backw,2)+1; % EVC<-FFA
%                         end
%                 end
%             end
%         end
%     end
% 
%     [BMA,BMR] = spm_dcm_peb_bmc(PEB, GCMs);
% 
%     for f = 1:Nfam
%         [BMA_fam_task,fam_task] = spm_dcm_peb_bmc_fam(BMA, BMR, family(f,:), 'ALL');
%         fprintf('\ncontrast %d, family %d\n',con-1,f)
%         fam_task.family.post
%         FamP(f,con-1) = fam_task.family.post(2);
%     end
% end
% round(FamP,2)
% 
% %% 6 ROIs
% % 5 family comparisons: self, EVC->OFA/FFA, OFA<->FFA, EVC<-OFA/FFA, bilat
% Nfam = 15;
% FamP = nan(Nfam,Ncon-1);
% for con = 2:Ncon
%     GCMs = cell(1,Nval^5);
%     m = 1; % First model is full model
%     GCMs{1,m} = bGCM;
% 
%     family = nan(Nfam,Nval^5);
%     family(:,1) = ones(Nfam,1)*2;
% 
%     for self = 1:Nval
%         for forw = 1:Nval
%             for OFAFFA = 1:Nval
%                 for backw = 1:Nval
%                     for bila = 1:Nval
%                         if m < Nval^5 % Ensuring last model is not repeat of first (full) model
%                             m=m+1;
%                             GCMs{1,m} = bGCM;
%                             GCMs{1,m}.b(2,2,con) = vals(self,1); GCMs{1,m}.b(3,3,con) = vals(self,2);
%                             GCMs{1,m}.b(2,1,con) = vals(forw,1); GCMs{1,m}.b(3,1,con) = vals(forw,2);
%                             GCMs{1,m}.b(3,2,con) = vals(OFAFFA,1); GCMs{1,m}.b(2,3,con) = vals(OFAFFA,2);
%                             GCMs{1,m}.b(1,2,con) = vals(backw,1); GCMs{1,m}.b(1,3,con) = vals(backw,2);
%                             GCMs{1,m}.b(5,2,con) = vals(bila,1); GCMs{1,m}.b(6,3,con) = vals(bila,1); 
% 
%                             GCMs{1,m}.b(5,5,con) = vals(self,1); GCMs{1,m}.b(6,6,con) = vals(self,2);
%                             GCMs{1,m}.b(5,4,con) = vals(forw,1); GCMs{1,m}.b(6,4,con) = vals(forw,2);
%                             GCMs{1,m}.b(6,5,con) = vals(OFAFFA,1); GCMs{1,m}.b(5,6,con) = vals(OFAFFA,2);
%                             GCMs{1,m}.b(4,5,con) = vals(backw,1); GCMs{1,m}.b(4,6,con) = vals(backw,2);
%                             GCMs{1,m}.b(2,5,con) = vals(bila,2); GCMs{1,m}.b(3,6,con) = vals(bila,2);
% 
%                             family(1,m) = max(vals(self,:))+1; % Any self
%                             family(2,m) = max(vals(OFAFFA,:))+1; % Any OFA<->FFA
%                             family(3,m) = max(vals(forw,:))+1; % Any forward from EVC
%                             family(4,m) = max(vals(backw,:))+1; % Any backward to EVC
%                             family(5,m) = max(vals(bila,:))+1; % Any OFA<->OFA FFA<->FFA
% 
%                             family(6,m) = vals(self,1)+1; % L+R OFA self
%                             family(7,m) = vals(self,2)+1; % L+R FFA self
%                             family(8,m) = vals(forw,1)+1; % L+R forw OFA
%                             family(9,m) = vals(forw,2)+1; % L+R forw FFA
%                             family(10,m) = vals(backw,1)+1; % L+R back OFA
%                             family(11,m) = vals(backw,2)+1; % L+R back FFA
%                             family(12,m) = vals(OFAFFA,1)+1; % L+R OFA->FFA
%                             family(13,m) = vals(OFAFFA,2)+1; % L+R FFA->OFA
%                             family(14,m) = vals(bila,1)+1; % R->L OFA+FFA
%                             family(15,m) = vals(bila,2)+1; % L->R OFA+FFA
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
%     [BMA,BMR] = spm_dcm_peb_bmc(PEB, GCMs);
% 
%     for f = 1:Nfam
%         [BMA_fam_task,fam_task] = spm_dcm_peb_bmc_fam(BMA, BMR, family(f,:), 'ALL');
%         fprintf('\ncontrast %d, family %d\n',con-1,f)
%         fam_task.family.post
%         FamP(f,con-1) = fam_task.family.post(2);
%     end
% end
% round(FamP,2)


%% Examine predicted univariate results
load('GCM_peb_fitAB');

Nvoi = 2; H = 1; %1 for right; 2 for left
cols = [1:9]; %9con
betas = []; dcm_betas = []; pvar = [];

for s = 1:length(objects_dir)
    for v = 1:Nvoi
        % sub = dosubs(s);
        % swd = sprintf('subject_%02d',sub);
        
        %% Get GLM betas (Might already have betas from above!)
        load(fullfile(target_dir{s},'bold',stat_dir,'SPM.mat'));
        if Nvoi == 2 && H == 1
            VOI = load(fullfile(target_dir{s},'bold',ROI_dir,VOInames{v+1}));
        elseif Nvoi == 3 && H == 1
            VOI = load(fullfile(target_dir{s},'bold',ROI_dir,VOInames{v}));
        elseif Nvoi == 2 && H == 2
            VOI = load(fullfile(target_dir{s},'bold',ROI_dir,VOInames{v+4}));
        elseif Nvoi == 3 && H == 2
            VOI = load(fullfile(target_dir{s},'bold',ROI_dir,VOInames{v+3}));
        elseif Nvoi == 6
            VOI = load(fullfile(target_dir{s},'bold',ROI_dir,VOInames{v}));
        end
        tmp = pinv(SPM.xX.xKXs.X)*VOI.Y;
        betas(s,v,:) = tmp(cols);
        
%        cb(r,s,:) = pinv([detrend(contrasts(2:end,:)',0) ones(9,1)])*bv;
        
        %% Get DCM betas
        Y2 = GCM{s,1}.y(:,v); % fitted data
        Y2 = Y2/GCM{s,1}.Y.scale;
%%        Y = Y2 + GCM{s,1}.R(:,v); % still not same as VOI.Y above?
%         X0 = GCM{s,1}.Y.X0;
%         Y2 = Y2 - X0*spm_inv(X0'*X0)*(X0'*Y2); % remove constant and drifts from DCM fit (why in there anyway?)
%         tmp = pinv(SPM.xX.X)*Y2;
        tmp = pinv(SPM.xX.xKXs.X)*Y2;
        dcm_betas(s,v,:) = tmp(cols);
        
        PSS   = sum(GCM{s,1}.y(:,v).^2);
        RSS   = sum(GCM{s,1}.R(:,v).^2);
        pvar(s,v) = PSS/(PSS + RSS);        
    end
    allF(s) = GCM{s,1}.F;
end
figure,boxplot(pvar)
figure,boxplot(allF)
mean(pvar)
mean(allF)


cw = [
     1/6 1/6 1/6         1/6 1/6 1/6         -1/3 -1/3 -1/3; % Perception - faces > scrambled
     1 -1 0              1 -1 0               1 -1 0; % Imm Rep
     1 0 -1              1 0 -1               1 0 -1; % Del Rep
     1/3 1/3 1/3         -1/3 -1/3 -1/3       0 0 0; % Recognition - famous > unfamous     
     1 -1/2 -1/2         1 -1/2 -1/2          0 0 0; % Imm vs Del Rep
     1/2 -1/2 0          1/2 -1/2 0           -1 1 0; % Perception X Imm Rep
     1/2 0 -1/2          1/2 0 -1/2           -1 0 1; % Perception X Del Rep
     1 -1 0              -1 1 0                0 0 0; % Recognition X Imm Rep
     1 0 -1              -1 0 1                0 0 0; % Recognition X Del Rep
    ];

T=[]; p=[]; figure; clf
sp = [2 4 6 1 3 5];
for v=1:Nvoi
    for c=1:size(cw,1)
         dc(:,c) = squeeze(betas(:,v,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    if Nvoi == 2 && H == 1
        VOItitle = VOInames{v+1};
    elseif Nvoi == 3 && H == 1
        VOItitle = VOInames{v};
    elseif Nvoi == 2 && H == 2
        VOItitle = VOInames{v+4};
    elseif Nvoi == 3 && H == 2
        VOItitle = VOInames{v+3};
    elseif Nvoi == 6
        VOItitle = VOInames{v};
    end
    fprintf('%s\n',VOItitle)
    disp([cw T' p'])
    
    for c=1:size(cw,1)
         dc(:,c) = squeeze(dcm_betas(:,v,:))*cw(c,:)';
         T(c) = mean(dc(:,c))/(std(dc(:,c))/sqrt(size(dc,1)));
    end  
    p = t2p(T,size(dc,1)-1,2);
    
    fprintf('DCM: %s\n',VOItitle)
    disp([cw T' p'])
    
    subplot(3,2,sp(v)),
    mbeta = mean(squeeze(betas(:,v,:)))';
    mbeta = mbeta/mean(mbeta)
    mdcm_beta = mean(squeeze(dcm_betas(:,v,:)))';  
    mdcm_beta = mdcm_beta/mean(mdcm_beta)
    bar([mbeta mdcm_beta])
    legend('real data','DCM fitted');
    Title = {'rEVC','rOFA','rFFA','lEVC','lOFA','lFFA'};
    if Nvoi == 2 && H == 1
        title(Title{v+1},'Interpreter','none')
    elseif Nvoi == 3 && H == 1
        title(Title{v},'Interpreter','none')
    elseif Nvoi == 2 && H == 2
        title(Title{v+4},'Interpreter','none')
    elseif Nvoi == 3 && H == 2
        title(Title{v+3},'Interpreter','none')
    elseif Nvoi == 6
        title(Title{v},'Interpreter','none')
    end
    xticklabels({'IniFF','ImmFF','DelFF','IniNF','ImmNF','DelNF','IniSF','ImmSF','DelSF'});
end

b=9;
for v = 1:Nvoi
    corr(betas(:,v,b),dcm_betas(:,v,b))
end


