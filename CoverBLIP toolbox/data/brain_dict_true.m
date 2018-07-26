function [dict, dict_norm, lookup_table] = brain_dict_true(RFpulses, TR, T1, T2, off)
% this function generates IR-BSSFP dictionary 
% Inputs: 
%                 RFpulses: sequence of L Flip Angles (RAD)
%                 TR: sequqnce of L repetition times (msec)
%                 T1: set of T1 values (msec)
%                 T2: set of T2 values (msec)
%                 off: set of off-resonance frequencies (kHz)

% Outputs: 
%                 dict: IR-BSSFP dictionaty (noormalized)
%                 dict_norm: fingerprints norm
%                 lut: look-up table for all atoms

% (c) Mohammad Golbabaee, 2017

L = size(RFpulses,1);

d = numel(T1)*numel(T2)*numel(off);
T1_set = zeros(numel(T1)*numel(T2), 1);
T2_set = zeros(numel(T1)*numel(T2), 1);
counter = 1;
for i = T1
    for j = T2
        T1_set(counter) = i;
        T2_set(counter) = j;
        counter = counter + 1;
    end
end
% Build lookup table
lookup_table = cell(d, 1);
counter = 1;
for k = off
    for i = T1
        for j = T2
            lookup_table{counter} = [k i j];
            counter = counter + 1;
        end
    end
end

dict = zeros(d, L);
Dprime = numel(T1)*numel(T2);
for k = 1:numel(off)
    dict((k-1)*Dprime+1:k*Dprime, :) = ...
        fastMRFdictionary_Grisword(RFpulses, TR, T1_set, T2_set, off(k));
end

dict = single(dict);
dict = transpose(dict);
dict_norm = single(sqrt(sum(abs(dict).^2, 1)));
dict = dict./(ones(L,1)*dict_norm); %NOTE dict is now normalized