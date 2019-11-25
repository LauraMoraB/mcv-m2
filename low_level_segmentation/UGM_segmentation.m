clear all;
close all;
clc;

rng(42);
%
%im_name='3_12_s.bmp';
%im_name='2_1_s.bmp';
im_name='7_9_s.bmp';

% TODO: Update library path
% Add  library paths
basedir='UGM/';
addpath(genpath(basedir));
output = 'results/7_9_s_2_labels50_icm';


%Set model parameters
%cluster color
K=2; % Number of color clusters (=number of states of hidden variables)

%Pair-wise parameters
smooth_term=[0.0 10]; % Potts Model
gamma=50; % Contrast-sensitive Model

%Load images
im_full = double(imread(im_name));
im = imresize(im_full, 0.3);

nRows = size(im,1);
nCols = size(im,2);
nNodes = nRows*nCols;

%Convert to LAB colors space
% TODO: Uncomment if you want to work in the LAB space
%
im = RGB2Lab(im);



%Preparing data for GMM fiting
%
% TODO: define the unary energy term: data_term
% nodePot = P( color at pixel 'x' | Cluster color 'c' )  


im = double(im);
x=reshape(im, [size(im,1)*size(im,2) size(im,3)]);
gmm_color = fitgmdist(x, K);
data_term = gmm_color.posterior(x);
mu_color = gmm_color.mu;

nodePot = data_term;



%Building 4-grid
%Build UGM Model for 4-connected segmentation
disp('create UGM model');

% Create UGM data
%[edgePot,edgeStruct] = CreateGridUGMModel(nRows, nCols, K, smooth_term);
[edgePot,edgeStruct] = CreateGridUGMModelContrast(im, K, gamma, false);


if ~isempty(edgePot)

    % color clustering
    [~,c] = max(reshape(data_term,[nRows*nCols K]),[],2);
    im_c= reshape(mu_color(c,:),size(im));
    
    % Call different UGM inference algorithms
    display('Loopy Belief Propagation'); tic;
    [nodeBelLBP,edgeBelLBP,logZLBP] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);toc;
    [~,c] = max(nodeBelLBP,[],2);
    im_lbp = reshape(mu_color(c,:),size(im));
    
    % Max-sum
    display('Max-sum'); tic;
    decodeLBP = UGM_Decode_LBP(nodePot,edgePot,edgeStruct);
    im_bp= reshape(mu_color(decodeLBP,:),size(im));
    toc;
    
    % TODO: apply other inference algorithms and compare their performance
    fprintf('Running ICM decoding - ');tic;
    ICMDecoding = UGM_Decode_ICM(nodePot,edgePot,edgeStruct);
    im_icm= reshape(mu_color(ICMDecoding,:),size(im));
    toc;
    
    if K<=2
        % Linear Programming
        fprintf('Lin Prog - '); tic;
        linProgDecoding = UGM_Decode_LinProg(nodePot,edgePot,edgeStruct);
        im_lp= reshape(mu_color(linProgDecoding,:),size(im));
        toc;

        fprintf('Running Graph Cut decoding - ');tic;
        optimal_GC = UGM_Decode_GraphCut(nodePot,edgePot,edgeStruct);
        im_graphCut= reshape(mu_color(optimal_GC,:),size(im));
        toc;
    end
    
    
    fig=figure;
    if K<=2
        subplot(3,3,1),imshow(Lab2RGB(im));xlabel('Original');
        subplot(3,3,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
        subplot(3,3,3),imshow(Lab2RGB(im_bp),[]);xlabel('Max-Sum');
        subplot(3,3,4),imshow(Lab2RGB(im_lbp),[]);xlabel('Loopy Belief Propagation');
        subplot(3,3,5),imshow(Lab2RGB(im_graphCut),[]);xlabel('Graph Cut');
        subplot(3,3,6),imshow(Lab2RGB(im_lp),[]);xlabel('Linear Programming');
        subplot(3,3,7),imshow(Lab2RGB(im_icm),[]);xlabel('ICM');

        
    else
        subplot(2,3,1),imshow(Lab2RGB(im));xlabel('Original');
        subplot(2,3,2),imshow(Lab2RGB(im_c),[]);xlabel('Clustering without GM');
        subplot(2,3,3),imshow(Lab2RGB(im_bp),[]);xlabel('Max-Sum');
        subplot(2,3,4),imshow(Lab2RGB(im_lbp),[]);xlabel('Loopy Belief Propagation');
        subplot(2,3,5),imshow(Lab2RGB(im_icm),[]);xlabel('ICM');
       
    end
    
else
   
    error('You have to implement the CreateGridUGMModel.m function');

end

print(fig,output,'-dpng')
