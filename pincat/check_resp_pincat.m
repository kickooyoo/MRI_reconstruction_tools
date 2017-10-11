% loads binary pincat outputs into a matrix
% run in directory with binary files
% wrapper for load_pincat
% 
% Mai Le 05/24/16
% cleaned up 05/18/17

%run in /y/mtle/resp_pincat
cd '/y/mtle/resp_pincat/Nt20/'
pincat_dir = './';
Nframes = 20; ;
resp_vals = 0:0.1:0.9;
for resp_ndx = 1:length(resp_vals)
	for frame = 0:Nframes - 1
		p(:,:,:,frame + 1,resp_ndx) = load_pincat(pincat_dir, 'prefix', sprintf('pincat_resp%.1f_frame%d_act', resp_vals(resp_ndx), frame));
	end
end
