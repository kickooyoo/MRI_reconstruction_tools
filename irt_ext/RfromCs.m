  function R = RfromCs(varargin)
%|function R = RfromCs(varargin)
%|
%| Create regularization object from "analysis" matrix (or fatrix2) C1.
%|
%| R(x) = sum_k wt_k pot([C1*x]_k)
%|
%| option
%|	'C1'		from Godwt1 or such.  default: 1
%|			must be capable of providing abs(C1)
%|                      now capable of handling multiple {C1}
%|	'pot'		from potential_fun. default: quadratic
%|                      now capable of handling multiple {pot}
%|	'wt'		size = sum_i size(C1_i*x), including beta. default: 1
%|      'mask'          logical support of x, default: true(size(x))
%|
%| out
%|	R		strum with methods:
%|			R.cgrad(R,x) column gradient
%|			R.penal(R,x) penalty value
%|			R.denom(R,x) denominator for SQS
%|			R.C(R,x) - of limited use.  C = diag(wt) * C1
%|
%| Copyright 2011-1-21, Jeff Fessler, University of Michigan
%| generalized to multiple C and pot, Mai Le, 2016-12-16

if nargin < 1, help(mfilename), error(' '), end
if nargin == 1 && streq(varargin{1}, 'test'), RfromCs_test, return, end

arg.C1 = 1;
arg.Cs = []; % use Cs for Cdiffs
arg.pot = {potential_fun('quad')};
arg.wt = {1};
arg.mask = []; % can give explicitly or read from Cdiffs
arg.chat = 1;
arg = vararg_pair(arg, varargin);

% convert Cdiffs to {C1}
if ~isempty(arg.Cs) && (arg.C1 == 1)
        if length(arg.Cs) == 1
        arg.C1 = arg.Cs.arg.Cc;
        elseif iscell(arg.Cs)
                arg.C1 = arg.Cs;
        else
                display('unk')
                keyboard
        end
                
end
% provided mask overrides Cdiffs mask
if ~isempty(arg.Cs) && isempty(arg.mask)
	if iscell(arg.Cs)
		Cs_mask = arg.Cs{1}.arg.mask;
		if ndims(Cs_mask) == 5
			arg.mask = Cs_mask(:,:,1,1,1);
		elseif ndims(Cs_mask) == 2
			arg.mask = Cs_mask;
		else
			display('weird Cs mask dims')
			keyboard
		end
	else
		arg.mask = arg.Cs.arg.mask;
	end
end
% for backwards compatibility with RfromC (single C1)
if ~iscell(arg.C1)
        arg.C1 = {arg.C1};
end
arg.NC1 = length(arg.C1);
if ~iscell(arg.pot) 
        arg.pot = {arg.pot};
end
% TODO: error checking on {wt} size
% if one pot provided for multiple C, apply for all
if (length(arg.pot) == 1) && (arg.NC1 > 1)
     arg.pot = repmat(arg.pot, arg.NC1, 1);
end

meth = {
	'cgrad', @RfromCs_cgrad, '(R,x)';
	'penal', @RfromCs_penal, '(R,x)';
	'denom', @RfromCs_denom, '(R,x)';
	'C', @RfromCs_C, '(R)';
};

R = strum(arg, meth);

% TODO different C and pot for each dim R from Cs
% quadratic for time, respiratory dimensions, leave huber for spatial
% dimensions appealing to have faster quadratic penalty on resp for MMI

% RfromCs_penal()
function penal = RfromCs_penal(dummy, arg, x)
penal = [];
if any(col(~arg.mask))
	x = embed(x, arg.mask);
end
for ii = 1:arg.NC1
        C1 = arg.C1{ii};
        tmp = C1 * x;
        potfun = arg.pot{ii};
        pot = potfun.potk(tmp);
        pot_mask = repmat(arg.mask, [ones(1, ndims(arg.mask)) size(pot, ndims(pot))]);
        pot = pot .* pot_mask;
        curr_penal = sum(arg.wt{ii}(:) .* pot(:));
        penal = [penal; curr_penal];
end

% RfromCs_cgrad()
function cgrad = RfromCs_cgrad(dummy, arg, x)
cgrad = 0;
if any(col(~arg.mask))
	x = embed(x, arg.mask);
end
for ii = 1:arg.NC1
	if arg.chat, display(sprintf('line 113 %d/%d in R cgrad at %s', ii, arg.NC1, datestr(now))), end
        tmp = reshape(arg.C1{ii} * x, arg.C1{ii}.odim);
        potfun = arg.pot{ii};
        wpot = potfun.wpot(tmp);
        if numel(arg.wt{ii}) == 1
                wt = wpot * arg.wt{ii}; % scalar wt
	elseif numel(arg.wt{ii}) == numel(wpot)/numel(x) % vector wt
		wt = wpot .* reshape(kron(arg.wt{ii}, ones(1, numel(x))), size(wpot));
	elseif numel(arg.wt{ii}) == numel(wpot) % for each elem
                wt = wpot .* reshape(arg.wt{ii}, size(wpot));
	else	
		display('bad choice of wt dims')
		keyboard
        end
	if arg.chat, display(sprintf('line 127 %d/%d in R cgrad at %s', ii, arg.NC1, datestr(now))), end
        curr_cgrad = arg.C1{ii}' * (wt .* tmp);
	if arg.chat, display(sprintf('line 129 %d/%d in R cgrad at %s', ii, arg.NC1, datestr(now))), end
	if any(col(~arg.mask))
		curr_cgrad = masker(curr_cgrad, arg.mask);
	end
        cgrad = cgrad + curr_cgrad;
end

% RfromCs_denom()
function denom = RfromCs_denom(dummy, arg, x)
denom = 0;
if any(col(~arg.mask))
	x = embed(x, arg.mask);
end
for ii = 1:arg.NC1
        Ca = abs(arg.C1{ii});
	if arg.chat, display(sprintf('line 140 %d/%d in R denom at %s', ii, arg.NC1, datestr(now))), end
        tmp = Ca * x;
	if arg.chat, display(sprintf('line 142 %d/%d in R denom at %s', ii, arg.NC1, datestr(now))), end
        potfun = arg.pot{ii};
        wpot = potfun.wpot(tmp);
        if numel(arg.wt{ii}) == 1
                wt = wpot * arg.wt{ii}; % scalar wt
	elseif numel(arg.wt{ii}) == numel(wpot)/numel(x) % vector wt
		wt = wpot .* reshape(kron(arg.wt{ii}, ones(1, numel(x))), size(wpot));
	elseif numel(arg.wt{ii}) == numel(wpot) % for each elem
                wt = wpot .* reshape(arg.wt{ii}, size(wpot));
	else	
		display('bad choice of wt dims')
		keyboard
        end
	%if arg.chat, display(sprintf('line 139 %d/%d in R denom at %s', ii, arg.NC1, datestr(now))), end
        c1 = Ca * ones(Ca.idim);
	if arg.chat, display(sprintf('line 157 %d/%d in R denom at %s', ii, arg.NC1, datestr(now))), end
        curr_denom = Ca' * (wt .* c1);
	if arg.chat, display(sprintf('line 159 %d/%d in R denom at %s', ii, arg.NC1, datestr(now))), end
	if any(col(~arg.mask))
        	curr_denom = masker(curr_denom, arg.mask);
	end
        denom = denom + curr_denom;
end

% RfromCs_C() %% buggy
function C = RfromCs_C(dummy, arg)
% cannot vertcat fatrix2 with empty matrix
% cannot change imask of C1 by multiplication with D
curr_wt = arg.wt{1} .* arg.mask; 
D = Gdiag(sqrt(curr_wt(:)), 'mask', true(size(curr_wt)));
C = D * arg.C1{1};
for ii = 2:arg.NC1
        curr_wt = arg.wt{ii} .* arg.mask; 
        D = Gdiag(sqrt(curr_wt(:)), 'mask', true(size(curr_wt)));
        curr_C = D * arg.C1{ii};
        C = [C; curr_C];
end

function RfromCs_test
close all
Nx = 2^5; 
Ny = 2^4;
Nz = 2^6;
Nt = 1;
Nresp = 2;
mask = true(Nx, Ny, Nz, Nt, Nresp); 
mask(2,:,:,:,:) = false; % stress
offsets = [1 Nx Nx*Ny];
x = zeros(size(mask)); 
x(end/2+1,end/2+1,end/2+1,:,1) = 1;
x(end/2+1,end/2+1,end/2-1,:,1) = 1;
xm = masker(x, mask); 
Cs = Cdiffs([Nx Ny Nz Nt Nresp], 'offsets', offsets, 'type_diff', 'circshift', 'mask', mask);
Cr = Crespdiff(squeeze(x), 'z_search', 3, 'mask_2D', mask(:,:,1,1,1), 'pf', 0, 'doweight', 0);
Cs = {Cs; Cr};
pot = {potential_fun('lange1', 0.01)};%
%         potential_fun('huber', 0.01)};
%         potential_fun('quad')};
wt = {3; 1};
R = RfromCs('pot', pot, 'Cs', Cs, 'wt', wt, 'mask', mask);

% cx = [R.C1{1} * x];% R.C1{2} * x];
% CfromR = R.C(R);
% cx = CfromR * x;
% cx = Cs * x;
im plc 2 2
cgrad = R.cgrad(R, xm);
denom = R.denom(R, xm);
penal = R.penal(R, xm);
im(1, x), cbar
% im(2, cx), cbar
im(3, cgrad), cbar
im(4, denom), cbar
