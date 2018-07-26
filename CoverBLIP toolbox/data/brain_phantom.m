function [phantom, density,T1_phantom,T2_phantom,df_phantom] = brain_phantom(RFpulses, TR)
% This function builds an MRF phantom using the Brainweb anatomical phantom
% along with selected T1, T2 and pd values based on existing tables 
% The anatomical model was described in:
% Design and Construction of a Realistic Digital Brain Phantom,
% DL Collins, AP Zijdenbos, V Kollokian, JG Sled, NJ Kabani, CJ Holmes and
% AC Evans. 
% IEEE Transactions on Medical Imaging, vol.17, No.3, p.463--468, June
% 1998.  
% 
% The function takes in RFpulses and TR sequences and returns the phantom
% and density image. 
% 
% Mohammad Golbabaee, 2017

% load in MNI segmented brain phantom
% MNI segmented brain phantom is segmented into:
%(0=Background, 1=CSF, 2=Grey Matter, 3=White Matter, 4=Fat, 5=Muscle/Skin, 6=Skin, 7=Skull, 8=Glial Matter, 9=Connective)

 [imaVol,scaninfo] = loadminc('phantom_1.0mm_normal_crisp.mnc'); % volume is  217   181   181

 %... and select a slice
 
 Slice = zeros(256,256);
 Slice(21:237,37:217) = imaVol(:,:,80);
 Slice = Slice.*(Slice<7);  % remove skull Glial matter and connective
 
% Tissue table chosen to match chemicals from Basics of MRI webbook:
% NOTE the pd are chosen to be particularly challenging.

% 8 is empty!
% From Basics of MRI (Web book)
% Tissue        T1 (s)      T2 (ms)     r*            phantom index
% Background    -           -           0                   0
% CSF           0.8 - 20    110 - 2000  70-230              1
% Gray          1.09 - 2.15 61 - 109    85 - 125            2
% White         0.76 - 1.08 61-100      70-90               3    
% Adipose       0.2 - 0.75  53 - 94     50 - 100            4  
% Muscle        0.95 - 1.82 20 - 67     45 - 90             5/6
 
 T1 = [0,5012,1545,811,530,1425,1425];
 T2 = [0, 512,83,77,77,41,41];
 pd = [0, 100,100,80,80,80,80];
 df = [0, -20, -40, -30, 50, 250, 250]/1000;
 
 L = length(TR);
 phantom = zeros(256,256,L);
 density = zeros(256,256);
 T1_phantom = zeros(256,256);
 T2_phantom = zeros(256,256);
 df_phantom = zeros(256,256);
 for i = 0:6
     
     density = density+(Slice==i).*pd(i+1);  
     T1_phantom = T1_phantom+(Slice==i).*T1(i+1);          
     T2_phantom = T2_phantom+(Slice==i).*T2(i+1);
     df_phantom = df_phantom+(Slice==i).*df(i+1);
     
     
     [Ix,Jx] = find(Slice==i);
     response = fastMRFdictionary_Grisword(RFpulses, TR ,T1(i+1), T2(i+1), df(i+1));
     for j = 1:length(Ix)
         phantom(Ix(j),Jx(j),:) = density(Ix(j),Jx(j)).*reshape(response,1,1,L);
     end
     
     

 end

end