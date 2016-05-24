pincat_dir = '.';

resp_vals = [0.1 0.3 0.5 0.7];
for resp_ndx = 1:length(resp_vals)
	for frame = 1:Nframes
		p(:,:,:,:,resp_ndx) = load_pincat(pincat_dir, 'prefix', sprintf('pincat_resp%.1f_out_frame%d', resp_vals(resp_ndx), frame));
	end
end
return
for ii = 1:10
	p1_shift = p(:,:,:,1 + ii:end, 1);
	p3_crop = p(:,:,:,1:end - ii, 2);
	err(ii) = calc_NRMSE_over_mask(p1_shift, p3_crop, true(size(p1_shift)));
end

% write_pincat_time_curves((1:50)*1e3, 'time_curves/liver_curve.txt')
