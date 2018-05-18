function [x,y,z]=fcnhexcolor(a,b,c)
%RGBCONV   Convert hex color to or from matlab rgb vector.
%   RGB = RGBCONV(HEX)  where HEX is a string of hexadecimal numbers
%      'RRGGBB' where the parts RR, GG, BB are two digit hexadecimal
%      numberes corresponding to the strengths of the colors red,
%      green and blue. RGB is the matlab rgb-vector.
%   RGB = RGBCONV(HEX)  where HEX is a string matrix of the form
%      ['RR';'GG';'BB'], or a cell string {'RR','GG','BB'}, or
%   RGB = RGBCONV('RR','GG','BB')
%   [R,G,B] = RGBCONV(...)  gives the corresponding matlab rgb-values
%   as scalar values in R, G and B correspondingly.
%
%   HEX = RGBCONV(RGB)  where RGB is a matlab rgb-vector with
%      three elements each corresponding to strengs of
%      red, green and blue and ranging between 0 and 1.
%   HEX = RGBCONV(R,G,B)  where R, G and B are scalar values
%      ranging between 0 and 1 each corresponding to red, green
%      and blue.
%      HEX is a string containing the hex values and is of length 6.
%   ['RR','GG','BB'] = RGBCONV(...)  gives the corresponding hex
%      values as two character strings 'RR', 'GG' and 'BB'
%
%   Examples:
%      rgbconv([.1 .2 .3])
%      rgbconv(.55,.6,.12)
%      [r,g,b]=rgbconv([.1 .2 .3])
%      rgbconv('1a334d')
%      rgbconv('1a','33','4d')
%      [r,g,b]=rgbconv({'F1','2A','55'})
%
%   See also COLORMAP.

% Copyright (c) 2001-10-09, B. Rasmus Anthin.
% Revisited 2002-03-06.
% Renewed 2003-12-16.

flag=0;           %meaning hex->rgb
switch(nargin)
case 1
   if iscell(a)
      arg=[a{1} a{2} a{3}];
   elseif ischar(a) & all(size(a)==[3 2])
      arg=[a(1,:) a(2,:) a(3,:)];
   elseif ischar(a) & all(size(a)==[6 1])
      arg=a';
   elseif ischar(a) & all(size(a)==[1 6])
      arg=a;
   elseif isnumeric(a) & all(size(a)==[1 3]) | all(size(a)==[3 1])
      arg=a;
      flag=1;
   else
      error('Syntax error.')
   end
case 3
   if ischar(a) & ischar(b) & ischar(c) & ...
         length(a)==2 & length(b)==2 & length(c)==2
      arg=[a(:)' b(:)' c(:)'];
   elseif isnumeric(a) & isnumeric(b) & isnumeric(c) & ...
         length(a)==1 & length(b)==1 & length(c)==1
      arg=[a b c];
      flag=1;
   end
otherwise
   error('Must be 1 or 3 input arguments.')
end
if flag        % rgb->hex
   out=dec2hex(round(255*arg));
   out=reshape(out',1,6);
   switch nargout
   case {0,1}
      x=out;
   case 3
      x=out(1:2);
      y=out(3:4);
      z=out(5:6);
   otherwise
      error('Must be 1 or 3 output arguments.')
   end
else           % hex->rgb
   out=[hex2dec(arg(1:2)) hex2dec(arg(3:4)) hex2dec(arg(5:6))]/255;
   switch nargout
   case {0,1}
      x=out;
   case 3
      x=out(1);
      y=out(2);
      z=out(3);
   otherwise
      error('Must be 1 or 3 output arguments.')
   end
end