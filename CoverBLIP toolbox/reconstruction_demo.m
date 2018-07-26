% MRF reconstruction Demo
%
% features:
% Brainweb digital phantom used for simulations
% IR-BSSFP fingerprint dictionary (with a set of pseudo-random excitations) 
% Multi-shot Echo Planar Imaging (EPI) Cartesian k-space acquisition
%
% Reconstruction algorithms:
% 1-CoverBLIP [Golbabaee, Chen, Wiaux, Davies; arxiv'18]
% accelerated iterative matched-filtering using cover tree's searches
%
% 2- TM [Dan Ma et al.'12]
% non-iterative template matching using brute-force searches  
%
% 3- BLIP [Davies et al.'13]
% itrative brute-force matched filtering
%
%
% If you find this package useful, please kindly refer to our article:
% ** M Golbabaee, Z Chen, Y Wiaux, M Davies, "Cover tree compressed sensing for fast MR fingerprint recovery", 
%    IEEE 27th International Workshop on Machine Learning for Signal, 2017 
%
% (c) Mohammad Golbabaee, Zhouye Chen 2017
%%

clc;
clear;
close all;
clear path

%% Set up initial parameters
disp('initialization...');
addpath './covertree_matlab_wrapper';
addpath('./tools');
addpath('./data');
addpath('./data/Bloch');


%% set hyperparameters
input.ALGO = 'CoverBLIP';
input.SNR = 50;
input.CS = [1/16]; % subsampling ratio
%---(Low-rank SVD) subspace reconstruction --%
input.svd_k = [200]; %number of svd components

%% load dictionary and data
if exist('./data/Bloch/brain_phantom_bSSFP_10tr.mat','file')==2
    load ./data/Bloch/brain_phantom_bSSFP_10tr.mat
    load ./data/Bloch/brain_dict_bSSFP_10tr.mat
else
    disp('Dictionary/Phantom data construction: if you run for the first time this step might take some time...');
    [dict,dict_norm,lut,phantom,density,T1_phantom,T2_phantom,df_phantom]=...
        generate_dict_phantom();
end

[N,M,L]=size(phantom);
% original dict---------------------
dict= dict.*( ones(size(dict,1),1)*dict_norm );
dict_norm = sqrt(sum(abs(dict).^2,1));
dict = dict./( ones(size(dict,1),1)*dict_norm );

%% k-space measurements
CS = input.CS;
% Create k-space undersampling (Multi-shot EPI)
H = sampling_mask_FASTsingle('deterministic comb',N,M,L, CS);
disp('initiate acquisition...');
% make measruements
data.y = H.forward(phantom);
SNR = input.SNR;
% Add measurement noise
data.y = data.y + randn(size(data.y))*std(data.y(:))*10^(-SNR/20)*1/sqrt(numel(data.y));

%% Pre-processing
data.Dnorm = dict_norm;
data.N = N;
data.M = M;
data.L = L;
disp('pre-processing step...');
% calculate SVD subspace bases
SVDk = input.svd_k;
param.k = SVDk;
if(param.k<data.L)
    [Vkt,~,~] = svd(dict*dict');
    Vkt = Vkt(:,1:param.k);
    param.svd = 1;
else
    Vkt = 1;
    param.svd = 0;
end
data.Vkt = @(z) Vkt*z;
data.Vk = @(z) Vkt'*z;
Dk = data.Vk(dict);
tmp = sqrt(sum(abs(Dk).^2,1 ));
Dk=Dk./(ones(param.k,1)*tmp);
tmp=tmp.*dict_norm;
D.times = @(X) ctranspose(Dk)*X; %complex conjugate transpose for matched filter
D.atoms = @(ind) Dk(:,ind);
D.scale = @(ind) tmp(ind);

% pre-process covertree structure on (normalized) dictionary
if(any(strcmp(input.ALGO,'CoverBLIP')))
    CTs = mexCovertree(transpose(Dk));
end

%% Reconstruction
disp('commence reconstruction...');

method = input.ALGO;
data.Vk = @(z) Vkt'*z;
data.Vkt = @(z) Vkt*z;
data.F = sampling_mask_FASTsingle('deterministic comb',N,M,L, CS,data.Vk,data.Vkt);
data.D = D;
param.step = 1/CS;
param.tol = 1e-4;

switch method
    case 'TM' %single step template matching [Ma et al. Nature'13]
        disp('TM...');
        param.algo = 'BLIP';
        param.iter = 1;
        data.D.match = @(z, ~, ~) find_nearest_MRF_mat(D, z, data.N, data.M);
        
        res= algo_IPA(data, param);       
        res.param = param;
        
    case 'BLIP' % iterative exact matched-filtering BLIP [Davies et al. SIAM'13]
        disp('BLIP...');
        param.algo = 'BLIP';
        param.iter = 50;
        data.D.match = @(z, ~, ~) find_nearest_MRF_mat(D, z, data.N, data.M);
        
        res = algo_IPA(data, param);
        res.param = param;
        
    case 'CoverBLIP'% iterative approximate matched-filtering CoverBLIP [Golbabaee et a. arxiv'18]
        data.D.match = @(z, epsilon, curr_ind) find_nearest_MRF_prox_ct(CTs, z, epsilon, D, curr_ind);
        disp('CoverBLIP...');
        param.algo = 'BLIP';
        param.iter = 50;
        param.current_est = 0;
        
        param.epsilon = 0.4; % covertree search inaccuracy level
        res= algo_IPA(data, param);
        res.param = param;
end


par.recon.mask_thresh = 61; % threshold on PD map to generate mask
par.ind.mtx_reco = N;
par.ind.Nslices = 1;
sig = get_mask(res.pd,par);

T1_IPA = lut(res.dm(:),1); T1_IPA= reshape(T1_IPA, N, M);
T2_IPA = lut(res.dm(:),2); T2_IPA= reshape(T2_IPA, N, M);
df_IPA = lut(res.dm(:),3); df_IPA= reshape(df_IPA, N,M);
res.qmap(:,:,1)=T1_IPA.*sig;
res.qmap(:,:,2)=T2_IPA.*sig;
res.qmap(:,:,3)=df_IPA.*sig;

%% visualize reconstructed maps
figure(1); 
subplot(141);imagesc(T1_phantom);colormap jet; set(gca,'xtick',[],'ytick',[]);caxis([0 2000]); axis image; title('T1 map'); colorbar('southoutside')
ylabel('ground truth'); 
subplot(142);imagesc(T2_phantom);colormap jet; set(gca,'xtick',[],'ytick',[]);caxis([0 300]); axis image; title('T2 map');colorbar('southoutside')
subplot(143);imagesc(df_phantom);colormap jet;set(gca,'xtick',[],'ytick',[]);caxis([-.050 .050]); axis image;title('B0 map');colorbar('southoutside')
subplot(144);imagesc(density/max(density(:)));colormap jet;set(gca,'xtick',[],'ytick',[]);title('PD map');axis image; colorbar('southoutside')



figure(2); title('ground truth maps')
subplot(141);imagesc(res.qmap(:,:,1));colormap jet; set(gca,'xtick',[],'ytick',[]);caxis([0 2000]); axis image; title('T1 map'); colorbar('southoutside')
ylabel(method)
subplot(142);imagesc(res.qmap(:,:,2));colormap jet; set(gca,'xtick',[],'ytick',[]);caxis([0 300]); axis image; title('T2 map'); colorbar('southoutside')
subplot(143);imagesc(res.qmap(:,:,3));colormap jet;set(gca,'xtick',[],'ytick',[]);caxis([-.050 .050]); axis image;title('B0 map'); colorbar('southoutside')
subplot(144);imagesc(abs(res.pd)/max(abs(res.pd(:))));colormap jet;set(gca,'xtick',[],'ytick',[]);title('PD map');axis image;colorbar('southoutside')


