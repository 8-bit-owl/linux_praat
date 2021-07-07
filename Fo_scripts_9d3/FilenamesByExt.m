function filename=FilenamesByExt(loadstr,filedirectory)
%FilenamesByExt
%   returns a list of filenames with a particular extension.  this function
%   was originally based on the following...  (eh, 20041108)
%Kelvin50GetFilename
%	
%	send a character string to be loaded (i.e. '.txt')
% 
% 	returns a filename chosen from directory
%
%		Eric Hunter 10/16/2000
%

% loadstr='.bmp'

tmp2=length(loadstr);
if nargin > 1
d = dir(filedirectory);
else
d = dir;
end
str = {d.name};
cntr=0;
for n=1:length(str)
  tmp=str{n};
  if length(tmp)>tmp2,
    if tmp(end-tmp2+1:end)==loadstr
      cntr=cntr+1;
    end
  end
end
str2=cell(cntr,1);
cntr=0;
for n=1:length(str)
  tmp=str{n};
  if length(tmp)>tmp2,
    if tmp(end-tmp2+1:end)==loadstr
      cntr=cntr+1;
      str2{cntr}=tmp;    
    end
  end
end
str=str2;
clear str2 tmp cntr d
filename=str;