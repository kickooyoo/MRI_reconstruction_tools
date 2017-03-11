function out = interp1nd(in, xq, interp_dim, varargin)
%function out = interp1nd(in, xq, interp_dim, varargin)
% for when you want to interpolated only along one dim of an nd volume
% unlike interpn, don't need to hog memory with nd grid points

arg.pf = 0; % to do: parallelize?
arg.extrap = '';
arg.method = 'linear'; % for interpn
arg = vararg_pair(arg, varargin);

if ndims(in) > 5
	display('hard coded for up to 5 dims only');
	keyboard
end

maxN = ndims(in);
[N(1), N(2), N(3), N(4), N(5)] = size(in);
Nid = N(interp_dim);
Nother = prod(N)/Nid;

permute_order = [setdiff(1:maxN, interp_dim) interp_dim];

% reshape so in [prod(Nother) Ninterp_dim]
a = reshape(permute(in, permute_order), Nother, Nid); 

x = 1:Nid;
% want to vectorize to make it fast but not sure how to handle edge conditions...
for ii = 1:Nother
	if ~isempty(arg.extrap)
		b(ii,:) = interp1(x, a(ii,:), xq, arg.method, arg.extrap);
	else
		b(ii,:) = interp1(x, a(ii,:), xq, arg.method);
	end
end

M = N;
M(interp_dim) = length(xq);

out = ipermute(reshape(b, permute(M, permute_order)), permute_order);

