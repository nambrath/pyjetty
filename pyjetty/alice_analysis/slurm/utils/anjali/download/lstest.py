import os

train_output_dir = "/alice/data/2016/LHC16q/000265525/pass2_CENT_woSDD/AOD244/PWGHF/HF_TreeCreator/776_20230504-0006_child_1"
temp_filelist_name = "subdirs_temp_000265525.txt"
cmd = 'alien_ls {} > {}'.format(train_output_dir, temp_filelist_name)
os.system(cmd)
with open(temp_filelist_name) as f:
	subdirs_all = f.read().splitlines()
	for x in subdirs_all:
		print(x, x.strip(" /").isdigit())
	subdirs = [ x for x in subdirs_all if x.isdigit() ]
	print(temp_filelist_name, subdirs)
