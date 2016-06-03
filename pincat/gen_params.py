#!/usr/bin/env python
import fileinput
import subprocess
import array
import sys
import os
import os.path

organ_line = 'liver_activity = '
resp_line = 'resp_start_ph_index = '

liver_vals = range(0, 150, 5)

# resp_start_ph_index 
resp_vals = [0.1, 0.3, 0.5, 0.7]


base_file = '/home/mtle/Documents/xcat-prog/general.samp.mri.par'

for resp_val in resp_vals :
	for frame in range(0, len(liver_vals)) :
		write_file = 'par/pincat_resp' + str(resp_val) + '_frame' + str(frame) + '.par'

		filedata = None
		f = open(base_file, 'r') 
		filedata = f.read()
		f.close()

		newdata = filedata.replace(organ_line + str(120), organ_line + str(liver_vals[frame]))

		newdata = newdata.replace(resp_line + str(0.7), resp_line + str(resp_val))

		f = open(write_file, 'w') 
		f.write(newdata)
		f.close()

subprocess.call(['chmod', '-R', '+w', './par'])
