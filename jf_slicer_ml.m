 function jf_slicer_ml(data, varargin)
%function jf_slicer_ml(data, [options])
%|
%| slice 3d data interactively (along 3rd dimension).
%| Uses scroll wheel to sweep through all slices.
%|
%| in
%|	data	[nx ny nz]
%|
%| options
%|	clim	[1 2]		clim arg to im()
%|	iz	[1]		initial slice (default: nz/2+1)
%|
%| Jeff Fessler, University of Michigan
%| edited by Mai Le to traverse 3rd dim of 4D objects
%| if 'mid3', traverse 4th dim

if ~nargin, help(mfilename), error(mfilename), end
if streq(data, 'test'), jf_slicer_test, return, end

arg.clim = [];
arg.iz = [];
arg.mid3 = [];
arg.colormap = 'gray';
arg.auto = 0;
arg.aspect = []; % how???
arg = vararg_pair(arg, varargin);
if ~isreal(data)
	printm 'warning: taking abs of complex data'
	data = abs(data);
end

if isempty(arg.clim)
	arg.clim = double(minmax(data)');
end

data = squeeze(data); % ML
if isempty(arg.mid3)
	nz = size(data, 3);
else
	nz = size(data,4);
end

if isempty(arg.iz)
	iz = ceil((nz+1)/2);
else
	iz = arg.iz;
end

clamp = @(iz) max(min(iz, nz), 1);

stuff.data = data;
stuff.arg = arg;

%jf_slicer_call(iz, stuff)
%drawnow
jf_slicer_show

if arg.auto
	scroll = true;
	set(gcf, 'KeyPressFcn', 'scroll = false;')
	while scroll
		iz = iz + 1;
		iz = mod(iz - 1, nz) + 1;
		pause(0.2)
		jf_slicer_show
	end
else
	set(gcf, 'WindowScrollWheelFcn', @jf_slicer_scroll)
end

%h = jf_add_slider('callback', @jf_slicer_call, 'data', {stuff}, ...
%	'min', 1, 'max', nz, 'sliderstep', [1 1]/(nz-1), 'value', iz);


function jf_slicer_scroll(src, event)
	iz = iz + event.VerticalScrollCount;
	iz = clamp(iz);
	jf_slicer_show
end % jf_slicer_scroll


function jf_slicer_show
	if ~isempty(arg.mid3)
		if isscalar(arg.mid3) && (arg.mid3 == 1)
			im('mid3', squeeze(data(:,:,:,iz)), arg.clim), cbar %ML
		else
			mid3 = [squeeze(data(:,:,arg.mid3(3),iz)) squeeze(data(:,arg.mid3(2),:,iz)); squeeze(data(arg.mid3(1),:,:,iz)).' zeros(size(data,3))];
			im(mid3, arg.clim);
		end
		colormap(gca, arg.colormap);
	else
		im(squeeze(data(:,:,iz,:)), arg.clim), cbar %ML
	end
	xlabelf('%d / %d', iz, nz)
end % jf_slicer_show

end % jf_slicer


function jf_slicer_call(iz, stuff)
persistent iz_save
if isempty(iz_save)
	iz_save = -1;
end
iz = round(iz);
if iz ~= iz_save
	iz_save = iz;
	arg = stuff.arg;
	im(stuff.data(:,:,iz), arg.clim), cbar
end

end % jf_slicer_call


function jf_slicer_test
data = reshape(1:7*8*9, 7, 8, 9);
im plc 1 2
im subplot 2
if im
	jf_slicer_ml(data, 'clim', [0 500])
end
end % jf_slicer_test
