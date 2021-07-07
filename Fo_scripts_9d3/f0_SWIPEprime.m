function [p, s, t]=f0_SWIPEprime(x,param)
% Calculates pitch and pitch strength using the SWIPE' algorithm.
%
% Outputs:
% p - Pitch value vector (Hz)
% s - Pitch strength vector
% t - Times at which pitch is estimated (s)
%
% Inputs:
% x - Speech signal vector
% param - Algorithm parameters struct. with following fields. Any number of
%         fields can be specified. Undefined fields will be filled with
%         default values. If param is not supplied then all the parameters
%         are set to default values.
%
%         --------------------- Parameter fields --------------------------
%         param.fs - Sampling frequency of audio data in Hz.
%                   (Default: 44.1e3 Hz)
%         param.plim - Pitch search range (Hz). (Default: [30 2000])
%         param.dt - Estimates pitch every dt seconds. (Default: 0.001)
%         param.dlog2p - Pitch candidates are distributed every dlog2p
%                        units on a base-2 logarithmic scale of Hz.
%                        (Default: 1/48)
%         param.dERBs - Spacing in ERBs. (Default: 0.1)
%         param.overlap - Overlap ratio for windows. (Default: 0.5)
%         param.sTHR - Pitch estimates with a strength lower than sTHR are
%                      treated as undefined. (Default: 0.2)
%         param.medFiltOrder - Order of median filter for filtering the
%                              pitch contour. (Default: 1 no filtering)
%
% References: Camacho, A., & Harris, J. G. (2008). "A sawtooth waveform
%             inspired pitch estimator for speech and music" The Journal of
%             the Acoustical Society of America, 124(3), 1638–52.
%             doi:10.1121/1.2951592
%
% Savyasachi Singh (University of Florida)          (Updated) March 8, 2013

error(nargchk(1, 2, nargin)) % Sanity Check
% Default Parameters
if nargin<2, param=struct; end % Create empty structure

% Process param struct
param=process_param(param,{'fs'; 'plim'; 'dt'; 'dlog2p'; 'dERBs'; 'overlap'; 'sTHR'; 'medFiltOrder'}...
    ,{44.1e3; [30 2e3]; 0.001; 1/48; 0.1; 0.5; 0.2; 1});

if strcmp(computer('arch'), 'win64')
    [p, s, t]=SWIPEPmex(x, param.fs, param);
elseif strcmp(computer('arch'), 'win32')
    [p, s, t]=SWIPEPrimeMex(x, param);
else
    error('Unsupported machine.')
end


%--------------------------------------------------------------------------
function out=process_param(in,fields,defaults)
% This function processes input parameter struct so that undefined fields
% are filled with default values.
%
% Outputs:
% out - Output parameter struct
%
% Inputs:
% in - Input parameter struct with missing fields or no fields at all.
% fields - Cell array of field name strings of size # fields -by- 1.
% defaults - Cell array of default values of fields. These values could be
%            of any valid MATLAB data type including cells and structs.
%            Should be of same size as fields.
%
% Savyasachi Singh (University of Florida)

if size(fields,2)~=1 || size(defaults,2)~=1
    error('fields and defaults should be cell arrays of size # fields -by- 1.');
end
if length(fields)~=length(defaults)
    error('fields and defaults should be cell arrays of same size.');
end

% Find already existing fields
curr_fields=fieldnames(in); % Already defined fields
curr_values=struct2cell(in); % Already defined field values

temp_fields=[curr_fields; fields]; % Concatenate fields to curr_fields in the end
temp_values=[curr_values; defaults]; % Concatenate defaults to curr_values in the end

[temp, indx]=unique(temp_fields,'first'); % Get names and indices of already defined fields and those which are not defined
out=cell2struct(temp_values(indx),temp,1); % Finally create the output struct

% EOF

