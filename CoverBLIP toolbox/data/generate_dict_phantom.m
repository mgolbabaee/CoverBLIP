function [dict,dict_norm,lut,phantom,density,T1_phantom,T2_phantom,df_phantom]= generate_dict_phantom()
% this function generates IR-BSSFP dictionary and synthetic MRF data from
% Brainweb phantom
% Outputs: 
%                 dict: IR-BSSFP dictionaty (noormalized)
%                 dict_norm: fingerprints norm
%                 lut: look-up table for all atoms
%
%                 phantom: time-series of MRF images
%                 density: ground-truth proton density map
%                 T1_phantom: ground-truth T1 density map (msec)
%                 T2_phantom: ground-truth T2 density map (msec)
%                 df_phantom: ground-truth off-resonance density map (kHz)
%          
%
% (c) Mohammad Golbabaee, 2017
%
load FAnTR_bSSFP.mat
RFpulses = RFpulses(1:1000);
TR = TR(1:1000)*1; % Arnold's 15msec
%TR = TR(1:1000); % 10msec

[phantom, density, T1_phantom, T2_phantom, df_phantom] = brain_phantom(RFpulses, TR);
save ./data/Bloch/brain_phantom_bSSFP_10tr.mat phantom T1_phantom T2_phantom density df_phantom;

T1 = [100:40:2000, 2200:200:6000];
T2 = [20:2:100, 110:4:200, 220:20:600];

off = [-250:40:-190, -50:2:50, 190:40:250]/1000;

[dict, dict_norm, lookup_table] = brain_dict_true(RFpulses, TR, T1, T2, off);
%%
lut=zeros(size(lookup_table,1),3);
for i = 1:size(lookup_table,1)
    lut(i,3) = lookup_table{i}(1);
    lut(i,1) = lookup_table{i}(2);
    lut(i,2) = lookup_table{i}(3);
end

save -v7.3 ./data/Bloch/brain_dict_bSSFP_10tr.mat RFpulses T1 T2 TR dict dict_norm lut

