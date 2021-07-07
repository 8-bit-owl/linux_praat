function [dataout] = MakeSameLength(datain,datatm,datatm0)
%MAKESAMELENGTH pulls in two vectors. then it interpolates the one to match
%the spacing of the other.
%
%   INPUT:
%   datain - data array that corresponds to the datatm
%   datatm - is the time array that goes with datain
%   datatm0 - is a time array to fit the datain and datatm to. 
%
%   OUTPUT:
%   dataout - it an interpolated datain array that now is fit to the time
%   array presented in data0
%
%   ejh - 20190418

tmp=find(datatm0<datatm(1)); 
if isempty(tmp)
   tmp2=datatm;
else
   tmp2=datatm-(datatm(1)-datatm0(tmp(end)));
end
[Lia,~] = ismember(round(100*((datatm0))),round(100*tmp2));
dataout=NaN(length(datatm0),1);
dataout(Lia)=datain;
dataout=single(dataout);


end

