-------------------------------------------
Files Necessary to run PINCAT
-------------------------------------------
- executable 
	(uses old XCAT with license bypass by MM)
	iv1/home/mtle/Documents/xcat-prog/dxcat2

- parameter file
	default CT contrast:
	
	iv1/home/mtle/Documents/xcat-prog/general.samp.par 
	
	quasi MR contrast 
	(based on gen_cine_256x100z50f.par from Behzad Sharif's PINCAT package):
	
	iv1/home/mtle/Documents/xcat-prog/general.samp.mri.par 

-------------------------------------------
Generating XCAT/PINCAT images
-------------------------------------------

Examples:

1) Single volume:
./dxcat <params_file> <resp_file_prefix>

Creates .bin files with image values and log file

2) Generate volumes over a range of respiratory/contrast vals:
./gen_pincat.py

note: must run in dxcat2 's original dir, not following soft link

3) Generate volumes over custom respiratory/contrast vals:
./gen_params.py
chmod 777 *.par % so dxcat can read new par files
./gen_pincat.py


-------------------------------------------
Reading XCAT/PINCAT images
-------------------------------------------
 
1) Read in range of respiratory/contrast volumes:
check_resp_pincat

2) Read in a single volume:
img = load_pincat(<dir_name>);


-------------------------------------------
Notes
-------------------------------------------

Behzad Sharif's PINCAT package, shared over email 04/18/16:
iv1/home/mtle/Documents/PINCAT
