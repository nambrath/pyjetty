import sys
import ROOT
import ctypes
from tqdm import tqdm

def debug(line, fout):
	file = ROOT.TFile(line)
	h = file.Get("h_matched_jet_ENC_RL2_bb_JetPt_R0.4_1.0")
	
	if h:
		h.GetXaxis().SetRange(61,80)
		if h.GetMaximum() > 0.1:
			maxbin = h.GetMaximumBin()
			xx, yy, zz = ctypes.c_int(0), ctypes.c_int(0), ctypes.c_int(0)
			h.GetBinXYZ(maxbin, xx, yy, zz)
		
			print(line[94:], h.GetMaximum(), xx.value, yy.value);
			output = line[94:] + " " + f'{h.GetMaximum():.5f} ' + int(xx.value) + " " + int(yy.value)
			fout.write(output)

# /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/anjali/25665164/slurm-25665164_191.out

if __name__ == "__main__":
	
	fout = open("out_25665164.txt","a")
	filelist = 'files_25665164.txt' # 'files_23791406.txt'
	
	with open(filelist) as f:
		for line in f:
			line = line.strip()
			debug(line, fout)
	
	fout.write("===============================================")
	fout.close()
			
