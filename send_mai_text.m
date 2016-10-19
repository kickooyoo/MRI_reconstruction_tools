function send_mai_text(msg);
%function send_mai_text(msg);
% if not working may need to turn on access for less secure apps in gmail account

[~, machine] = system('hostname');
if ~isempty(findstr(machine, 'eecs')) || ~isempty(findstr(machine, 'iv1'))
	send_text_message('2149126246','tmobile','Matlab',msg);
end
