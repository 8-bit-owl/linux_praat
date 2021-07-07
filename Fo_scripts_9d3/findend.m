function [lngth_]=findend(in_file_)
% FINDEND scrolls through a text file to get the length of it
% 
%  Eric Hunter NCVS


%clear
%in_file_='0S1_V1_sr125_ns13_ln3_fb1_a0.out';
numl_=2;

cnvg_=-1;
disp('finding length...')
leave_=0;     

cntr_=0;
fid=fopen(in_file_);
while 1
cntr_=cntr_+1;
tline = fgetl(fid);
if ~ischar(tline), break, end
%disp(tline)
end
fclose(fid);
lngth_=cntr_;
