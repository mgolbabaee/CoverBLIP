function mask = get_mask(data,par)
% % Get mask 
% =========================================================================
% mask = get_mask(data,par)
%   Input
%       data         [Nx,Ny,(Nz),...]  data            
%
%       par.        Anything that somehow controls the recon pipeline goes into this struct
%
% %
%   Output
%       mask          [Nx,Ny,(Nz)]     mask  
%
% %   
%   Function
% =========================================================================
% (c) P. Gomez TUM.de
% =========================================================================
if par.ind.Nslices == 1
    data = data(:,:,:); % concatenate all extra dims
    mask = squeeze(sum(abs(data),3)).^2; %sum over extra dims
    mask = mask/max(mask(:));
    [x,y] = gradient(mask);
    gradnorm = sqrt(x.^2 + y.^2);
    mask = mask + gradnorm;
else
    data = data(:,:,:,:);
    mask = squeeze(sum(abs(data),4)).^2; %sum over extra dims
       mask = mask/max(mask(:));
    [x,y,z] = gradient(mask);
    gradnorm = sqrt(x.^2 + y.^2 + z.^2);
    mask = mask + gradnorm;
end
mask = mask./max(mask(:));
h = fspecial('gaussian',par.ind.mtx_reco,par.ind.mtx_reco/4); %add gaussian prior to the mask: data is probably in the center of the img
h = h/max(h(:));
mask = mask.*h;
mask = mask>prctile(mask(:),par.recon.mask_thresh); %mask all data bellow thresh of X percentile
mask = imfill(mask,'holes');

end

