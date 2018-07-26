function [dict,dictz] = fastMRFdictionary_Grisword(RFpulses, TR ,T1, T2, df,initial_state)
%
% This function calculates the MRF signal evolution for an IR-bSSFP-based
% readout for a single material (T1, T2, df) for the input RF pulse train
% and TR.
%
% This code is based on the single isochromat method popularized by Klaus
% Scheffler. This works for bSSFP as long as one can assume that both the
% flip angle and off-resonance frequency are essentially homogeneous inside
% a voxel.
%
%       INPUTS:
%               RFpulses        Complex RF pulse train
%                                   (1 x #timepoints)
%
%               TR              TR  (1 x #timepoints)
%
%               T1, T2, df      T1, T2, and off-resonance frequency
%
%                               Make sure these are all in the same
%                               time units!
%
%               initial_state   Defines whether we use inversion recovery:
%                              initial_state = -1 (default) or not:
%                              initial_state = 1;
%
%       OUTPUTS:
%
%               dict            Calculated dictionary
%
%
% Based on a version by "Nov. 15, 2012     Dan Ma, Mark Griswold"
% This is a fast vectorized version written by MED 22-04-13
% 
%CORRECTION DONE FOR TE=TR/2 sampling by Mo Golobabaee 01.04.2016

if nargin <6
    initial_state = -1;
elseif abs(initial_state)~=1
    error('fastMRFdictionary', 'initial state can only be +/- 1')
end

if nargin<5
    df = 0;
end

if nargin<4
    T2 = 50;
end

if nargin<3
    T1 = 1000;
end

if nargin<2
    TR = rand(1,1000)*4+10;
end

if nargin<1
    RFpulses = 10*pi./180*randn(1,1000);
end

N=length(TR);
D = length(T1);  % We assume that T2 and df are the same size

rf=abs(RFpulses);
rph=angle(RFpulses);

dict=zeros(D,N);
dictz = dict;

% Create and initialize the states matrix based on whether we invert or not

M1 = repmat([0 0 initial_state].',1,D);


%m1=[0 0 -1].';  % This assumes a perfect inversion pulse with no delay

for i=1:N
  
    % Define the linear dynamic operators for pulse i
    
    rx =    [   1.0 0.0 0.0;        % rotation matrix for pulse
                0.0 cos(rf(i)) sin(rf(i));
                0.0 -sin(rf(i)) cos(rf(i))];
            
    rdzp =  [   cos(rph(i)) sin(rph(i)) 0.0; % RF phase
                -sin(rph(i)) cos(rph(i)) 0.0;
                0.0 0.0 1.0];
    
    rdzm =  [   cos(-rph(i)) sin(-rph(i)) 0.0;
                -sin(-rph(i)) cos(-rph(i)) 0.0;
                0.0 0.0 1.0];
                 
    M1=rdzp*rx*rdzm*M1;    % do RF pulse for all states
    
    % relaxation terms
    E12 = [exp(-TR(i)./(2*T2)');exp(-TR(i)./(2*T2)');exp(-TR(i)./(2*T1)')]; %write the relaxation terms as a long matrix 
    
    % relax 1st 1/2 of TR
    M1 = E12.*M1;
    M1(3,:) = M1(3,:)+(1-exp(-TR(i)./(2*T1)'));
    
    beta=df.'.*TR(i)*2*pi; % beta is now a Dx1 vector associated with df

    % Off-resonance for first 1/2 of TR
    M1(1:2,:) = M1(1:2,:).*repmat([cos(beta./2);cos(beta/2)], 1, size(M1, 2)) + ...
        M1(2:-1:1,:).*repmat([sin(beta./2);-sin(beta./2)], 1, size(M1, 2));
    
    % update dictionary
    %Sample assuming TE=TR/2     
    dict(:,i) = M1(1,:)+1i.*M1(2,:);           
    dictz(:,i)= M1(3,:);
    
    % 2nd Half relax times
    M1 = E12.*M1;
    M1(3,:) = M1(3,:)+(1-exp(-TR(i)/2./T1'));   
    %dictz(:,i)= M1(3,:);
   
    % Off-resonance for 2nd 1/2 of TR
    M1(1:2,:) = M1(1:2,:).*repmat([cos(beta./2);cos(beta/2)], 1, size(M1, 2)) + ...
       M1(2:-1:1,:).*repmat([sin(beta./2);-sin(beta./2)], 1, size(M1, 2));
    
end