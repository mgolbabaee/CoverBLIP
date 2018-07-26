function H = sampling_mask_FASTsingle(sampling, N,M,L, percentage, Vk,Vkt)
% This function generates forward and adjoint sub-sampled FFT operators
% corresponding to the Multi-shot Echo Planar Imaging (EPI) MRF aqcuisition. 
% Inputs: 
%                 sampling: 'deterministic comb' i.e. the Multi-shot EPI
%                 acquisition
%                 N: number of rows
%                 M: number of columns
%                 percentage: sub-sampling rate
%                 vk,vkt: compression/decompression bases (when SVD subspace option is on)

% Outputs: Forward model
%                 H.forward
%                 H.adjoint
%
% (c) Mohammad Golbabaee, 2017
%%
if nargin<6
    Vk=@(x) x; Vkt=Vk;
end
switch sampling %'random comb'; %'spiral'; % 'deterministic'; %  'random line'; %    
    case 'deterministic comb' % generates uniform comb undersampling that cycles 1 point in k-space each time
        step = round(1/percentage);
        no_of_steps = floor(N/step);
        nb_meas = no_of_steps*M;
        mask = zeros(nb_meas*L,1);
        comb = zeros(N,1); comb(1:step:step*nb_meas/M) = 1;
        for ind = 1:L
            comb = comb([N,1:N-1]);  %cyclic shift
            template = comb*ones(1,M);
            mask(nb_meas*(ind-1)+1:nb_meas*(ind)) = N*M*(ind-1)+find(template(:)==1);            
        end   
end % end switch

H.forward = @(x) single(FW_helper(x, mask,N,M,L,Vkt));
H.adjoint = @(y) single(Adj_helper(y, mask,N,M,L,Vk));
end

function x = FW_helper(x, mask,N,M,L,Vkt)
x=reshape(fft2(x),N*M,[]);
x=transpose( Vkt(transpose(x)) );
x=x(mask)/sqrt(N*M);
end

function tmp = Adj_helper(y, mask,N,M,L,Vk)
tmp = zeros(N*M*L,1);
tmp(mask)=y;
tmp = transpose(Vk( transpose( reshape(tmp, [], L) )));
%y = full(sparse(mask,1,y,N*M*L,1));
tmp = ifft2(reshape(tmp,N, M,[]))*sqrt(N*M);
end