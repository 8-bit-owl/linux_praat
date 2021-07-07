function r = fftE_eh(x,frame_sz,frame_step,Fs)
% computes the spectral engergy 
%
% r = rms(x,[frame_sz]) computes the root mean square of
% signal x. if frame_sz is included, the rms value of
% each frame is computed depending on frame_sz (frame size).
%
%   Eric Hunter, NCVS@UTAH, 20120411
%

if nargin ==1
   r = sqrt(mean(x.^2));
end

if nargin >2
   num_frames=floor((length(x)/frame_step))+1;
   r = zeros(1,num_frames);
   wind=1:frame_sz;
   
   fu = min([8000 floor(Fs/2)-1]);
   fl = 100;
   c = fir1(250,[fl fu]/(Fs/2),'high');  % filter
   x = filtfilt(c,1,x);
   
   for n = 1:num_frames,
      index = (n-1)*floor(frame_step)+1:(n-1)*floor(frame_step)+frame_sz;
      if(index(end)>length(x))
          index = index(1):length(x);
      end
      if isempty(index),
          r(n)=0;
      else
          tmp=abs(fft(x(index)));
          tmp=tmp(1:floor(length(index)/2));
          r(n) = sum(tmp);
      end
   end
elseif nargin==2
   num_frames=floor(length(x)/frame_sz);
   r = zeros(1,num_frames);
   for n = 1:num_frames,
      index = round((n-1)*frame_sz+1:(n)*frame_sz);
      tmp=abs(fft(x(index)));
      tmp=tmp(1:floor(length(index)/2));
      r(n) = sum(tmp);
%       r(n) = sqrt(mean(x(index).^2));
   end
end

r(r <= 0) = .000000000000000001;


return
