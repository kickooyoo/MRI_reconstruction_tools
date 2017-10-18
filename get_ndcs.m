function n = get_ndcs(N, varargin)
%function n = get_ndcs(N, varargin)
% 
% returns array of integer indices of length N
% zero is centered for odd N, at ndx N/2+1 for even N (FFT convention)
%
% varargin: index of zero
if nargin == 1
        zndx = ceil((N + 1)/2);
elseif nargin == 2 && mod(varargin{1}, 1) == 0
        zndx = varargin{1};
else
        display('invalid varargin for get_ndcs, centered at zero');
        zndx = ceil((N + 1)/2);
end

if mod(N, 2) == 0
        n = [- N/2 : N/2 - 1].';
        shift = (N/2 + 1) - zndx; 
else
        n = [- (N - 1)/2 : (N - 1)/2].'
        (N - 1)/2
        shift = (N - 1)/2 + 1 - zndx;
end
n = n + shift;

end

