 function ml_slicer(data, varargin)
%function ml_slicer(data, [options])
%|
%| slice 3d data interactively (along 3rd dimension).
%| Uses scroll wheel to sweep through all slices.
%|
%| in
%|	data	[nx ny nz nc]
%|
%| options
%|	clim	[1 2]		clim arg to im()
%|	iz	[1]		initial slice (default: nz/2+1)
%|
%| Jeff Fessler, University of Michigan
%| edited by Mai Le to traverse 3rd dim of 4D objects
%| if 'mid3', traverse 4th dim

if ~nargin, help(mfilename), error(mfilename), end
if streq(data, 'test'), ml_slicer_test, return, end

arg.clim = [];
arg.iz = [];
arg.mid3 = [];
arg.colormap = 'gray';
arg.auto = 0;
arg.aspect = [1 1 1];
arg.line = []; % Ndim x 2, Ndim = 2 or 3
arg.sliders_on = 0;
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

if ~isempty(arg.mid3)
% mid3 plotting constants
border = 30; % pixels?
x0 = 10;
y0 = 10;
[Nx, Ny, Nz, N4, N5] = size(data);
Ax = Nx * arg.aspect(1);
Ay = Ny * arg.aspect(2);
Az = Nz * arg.aspect(3);
% arg.line f= diag(arg.aspect) * arg.line; % (3x3) x (3x2)
if isscalar(arg.mid3) && (arg.mid3 == 1)
        arg.mid3 = round([Nx Ny Nz]/2);
end

sxy = subplot(221); 
sxz = subplot(223); 
syz = subplot(222); 
full_width = (x0 + Ay + 2*border + Az + 10 + 100); % 10 for colobar
full_height = (y0 + Az + 2*border + Az + 100);
pos = get(gcf, 'Position');
P = max(full_width/pos(3), full_height/pos(4));
% left bottom width height
set(sxy, 'Units', 'pixels', 'Position', [x0/P, (y0 + Az + border)/P, Ay/P, Ax/P]);
set(sxz, 'Units', 'pixels', 'Position', [x0/P, y0/P, Ay/P, Az/P]);
set(syz, 'Units', 'pixels', 'Position', [(x0 + Ay + border)/P, (y0 + Az + border)/P, Az/P, Ax/P]);
end

clamp = @(iz) max(min(iz, nz), 1);

stuff.data = data;
stuff.arg = arg;

%ml_slicer_call(iz, stuff)
%drawnow
ml_slicer_show

if arg.auto
	scroll = true;
	set(gcf, 'KeyPressFcn', 'scroll = false;')
	while scroll
		iz = iz + 1;
		iz = mod(iz - 1, nz) + 1;
		pause(0.2)
		ml_slicer_show
	end
else
	set(gcf, 'WindowScrollWheelFcn', @ml_slicer_scroll)
end
set(gcf,'ResizeFcn',@Resize_clbk);

if arg.sliders_on && ~isempty(arg.mid3)
        x_slider = uicontrol('Style', 'slider','Min', 1, 'Max', Nx, ...
                'Position', [(x0 + Ay + border)/P (y0 + 2*border)/P Nx 20], ...
                'Value', arg.mid3(1), 'Callback', @refresh_mid3); 
        y_slider = uicontrol('Style', 'slider','Min', 1, 'Max', Ny, ...
                'Position', [(x0 + Ay + border)/P (y0 + border)/P Ny 20], ...
                'Value', arg.mid3(2), 'Callback', @refresh_mid3); 
        z_slider = uicontrol('Style', 'slider','Min', 1, 'Max', Nz, ...
                'Position', [(x0 + Ay + border)/P y0/P Nz 20], ...
                'Value', arg.mid3(3), 'Callback', @refresh_mid3); 
end

%h = jf_add_slider('callback', @ml_slicer_call, 'data', {stuff}, ...
%	'min', 1, 'max', nz, 'sliderstep', [1 1]/(nz-1), 'value', iz);


function ml_slicer_scroll(src, event)
	iz = iz + event.VerticalScrollCount;
	iz = clamp(iz);
	ml_slicer_show
end % ml_slicer_scroll

function ml_slicer_show
        if ~isempty(arg.mid3)
                xy = squeeze(data(:,:,arg.mid3(3),iz));
                xz = squeeze(data(:,arg.mid3(2),:,iz));
                yz = squeeze(data(arg.mid3(1),:,:,iz)).';
                % imagesc to fill subplot space, .' to match im orientation
                subplot(sxz), imagesc(xz.', arg.clim); 
                colormap(gca, arg.colormap); 
                title(sprintf('%d/%d x-z plane', arg.mid3(2), Ny))
                subplot(sxy), imagesc(xy.', arg.clim); 
                colormap(gca, arg.colormap); 
                title(sprintf('%d/%d x-y plane', arg.mid3(3), Nz))
                subplot(syz), imagesc(yz.', arg.clim); 
                colormap(gca, arg.colormap); cbar
                title(sprintf('%d/%d y-z plane', arg.mid3(1), Nx))
                % restore position after colorbar
                set(syz, 'Position', [(x0 + Ay + border)/P, (y0 + Az + border)/P, Az/P, Ax/P]);
                % cbar?
                if ~isempty(arg.line)
                        in_x_view = arg.mid3(1) <= max(arg.line(1,:)) && arg.mid3(1) >= min(arg.line(1,:));
                        in_y_view = arg.mid3(2) <= max(arg.line(2,:)) && arg.mid3(2) >= min(arg.line(2,:));
                        in_z_view = arg.mid3(3) <= max(arg.line(3,:)) && arg.mid3(3) >= min(arg.line(3,:));
                        if in_y_view
                                subplot(sxz), hold on, plot(arg.line(1,:), arg.line(3,:), 'r');
                        end
                        if in_z_view
                                subplot(sxy), hold on, plot(arg.line(1,:), arg.line(2,:), 'r');
                        end
                        if in_x_view
                                % yz plot axes flipped
                                subplot(syz), hold on, plot(arg.line(3,:), arg.line(2,:), 'r');
                        end
                end
        else
                im(squeeze(data(:,:,iz,:)), arg.clim), cbar %ML
                if ~isempty(arg.line)
                end
        end
        xlabelf('%d / %d', iz, nz)
end % ml_slicer_show


function Resize_clbk(src, event)
        
pos = get(gcf, 'Position');
P = max(full_width/pos(3), full_height/pos(4));
% left bottom width height
set(sxy, 'Units', 'pixels', 'Position', [x0/P, (y0 + Az + border)/P, Ay/P, Ax/P]);
set(sxz, 'Units', 'pixels', 'Position', [x0/P, y0/P, Ay/P, Az/P]);
set(syz, 'Units', 'pixels', 'Position', [(x0 + Ay + border)/P, (y0 + Az + border)/P, Az/P, Ax/P]);
ml_slicer_show
end

function refresh_mid3(src, event)
        arg.mid3(1) = round(get(x_slider,'Value')); 
        arg.mid3(2) = round(get(y_slider,'Value')); 
        arg.mid3(3) = round(get(z_slider,'Value')); 
        ml_slicer_show;
end

end % ml_slicer



function ml_slicer_call(iz, stuff)              
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

end % ml_slicer_call


function ml_slicer_test
data = reshape(1:7*8*9, 7, 8, 9);
im plc 1 2
im subplot 2
if im
	ml_slicer(data, 'clim', [0 500])
end
end % ml_slicer_test
