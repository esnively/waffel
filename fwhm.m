function width = fwhm(x,y)

% function width = fwhm(x,y)
%
% Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% and its polarity.
% The FWHM result in 'width' will be in units of 'x'
%
%
% Rev 1.2, April 2006 (Patrick Egan)

y = smooth( y / max(y));
N = length(y);

if (max(y) < 2*mean(y))
    width = (x(end)-x(1))/2;
    return
end


lev50 = 0.5;
if y(1) < lev50                  % find index of center (max or min) of pulse
    [garbage,centerindex]=max(y);
    Pol = +1;
%    disp('Pulse Polarity = Positive')
else
    [garbage,centerindex]=min(y);
    Pol = -1;
    disp('Pulse Polarity = Negative')
end
i = 2;
while ((i <= N-5)& sign(y(i)-lev50) == sign(y(i-1)-lev50) )
    i = i+1;
end                                   %first crossing is between v(i-1) & v(i)
interp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + interp*(x(i)-x(i-1));
i = centerindex+1;                    %start search for next crossing at center
while ( (i <= N-5) & (sign(y(i)-lev50) == sign(y(i-1)-lev50)) )
    i = i+1;
end
if i ~= N
    Ptype = 1;  
%    disp('Pulse is Impulse or Rectangular with 2 edges')
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
else
    Ptype = 2; 
    disp('Step-Like Pulse, no second edge')
    ttrail = NaN;
    width = NaN;
end


