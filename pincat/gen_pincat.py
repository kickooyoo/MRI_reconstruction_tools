#!/usr/bin/env python
import subprocess
import array
import smtplib
import time
import os.path

smtpObj = smtplib.SMTP_SSL('smtp.gmail.com', 465)
smtpObj.login('051364@gmail.com', 'M3a1t4L1a5b8')
#smtpObj.sendmail('051364@gmail.com', '2149126246@tmomail.net', 'test text from Python')


Nframes = 20
resp_vals = [x / 10.0 for x in range(0, 11, 1)]
#[0.1, 0.3, 0.5, 0.7]
for resp_val in resp_vals :
	for frame in range(0, Nframes) :
                start = time.time()
        	resp_file = '/y/mtle/resp_pincat/pincat_resp' + str(resp_val) + '_frame' + str(frame)
		par_file = '/y/mtle/resp_pincat/par/pincat_resp' + str(resp_val) + '_frame' + str(frame) + '.par'
                if not os.path.isfile(par_file) :
                        print(par_file, 'not found')
                call = './dxcat2', par_file, resp_file
                print call
	        subprocess.call(call) 
                print time.time() - start, 'seconds to generate one resp phase'

smtpObj.sendmail('051364@gmail.com', '2149126246@tmomail.net', 'done with PINCAT')
smtpObj.quit()
