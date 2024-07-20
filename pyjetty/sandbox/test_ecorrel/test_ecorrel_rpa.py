#!/usr/bin/env python3

import array
import numpy as np
import matplotlib.pyplot as plt
from pyjetty.mputils import BoltzmannEvent
from pyjetty.mputils.mputils import logbins, linbins
import operator as op
import itertools as it
import sys
import os
import argparse
from tqdm import tqdm
from heppy.pythiautils import configuration as pyconf
import pythiaext
import pythiafjext
import pythia8
import fjtools
import ecorrel
import fjcontrib
import fjext
import fastjet as fj
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
# ROOT.gSystem.AddDynamicPath('$HEPPY_DIR/external/roounfold/roounfold-current/lib')
# ROOT.gSystem.Load('libRooUnfold.dylib')
# ROOT.gSystem.AddDynamicPath('$PYJETTY_DIR/cpptools/lib')
# ROOT.gSystem.Load('libpyjetty_rutilext')
# _test = ROOT.RUtilExt.Test()
from tqdm import tqdm
import argparse
import os


def get_args_from_settings(ssettings):
	sys.argv = sys.argv + ssettings.split()
	parser = argparse.ArgumentParser(description='pythia8 fastjet on the fly')
	pyconf.add_standard_pythia_args(parser)
	parser.add_argument('--ncorrel', help='max n correlator', type=int, default=2)
	parser.add_argument('--output', default="test_ecorrel_randtrk.root", type=str)
	parser.add_argument('--user-seed', help='pythia seed',
											default=1111, type=int)
	args = parser.parse_args()
	print(args)
	return args

def hadron_efficiency(cms_rpa, hpt):
	b = cms_rpa.FindBin(hpt)
	prob =  cms_rpa.GetBinContent(b)
	if prob < 1:
		keep = np.random.random() <= prob
		if keep:
			return 1
		return 0
	else:
		add = np.random.random() <= (prob-1)
		if add:
			return 2
		return 1
	# 0: drop particle
	# 1: keep particle
	# 2: add particle

def randomize_added(parts, rpa_mask, jet=None):
	out = []
	n = len(parts)
	for i in range(n):
		r = rpa_mask[i]
		p = parts[i]
		if r:
			out.append(p)
			if r==2:
				trk = fj.PseudoJet()
				if jet:
					_r = np.random.uniform(0,1)*0.4
					_t = np.random.uniform(0,2*np.pi)
					eta = jet.rapidity() + _r*np.sin(_t)
					phi = jet.phi() + _r*np.cos(_t)
				else:
					eta = np.random.uniform(-1.3,1.3)
					phi = np.random.uniform(0,2*np.pi)
				trk.reset_PtYPhiM(p.perp(), eta, phi)
				trk.set_user_index(888)
				out.append(trk)
	return out

def main():
	mycfg = []
	ssettings = "--py-ecm 5000 --py-pthatmin 20"
	args = get_args_from_settings(ssettings)
	pythia_hard = pyconf.create_and_init_pythia_from_args(args, mycfg)
	jet_R0 = 0.4
	max_eta_jet = 0.9
	max_eta_hadron = max_eta_jet + jet_R0
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	jet_selector = fj.SelectorPtMin(20.0) & fj.SelectorAbsEtaMax(max_eta_jet) & (~fj.SelectorIsPureGhost())
	pfc_selector0 = fj.SelectorPtMin(0.)
	pfc_selector1 = fj.SelectorPtMin(1.)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)

	# jet_def_sub = fj.JetDefinition(fj.kt_algorithm, jet_R0)
	# jet_selector_sub = fj.SelectorAbsEtaMax(max_eta_jet) & (~fj.SelectorNHardest(2)) & (~fj.SelectorIsPureGhost()) # used to calculate rho
	# median_subtractor = fj.JetMedianBackgroundEstimator(jet_selector_sub, jet_def_sub, fj.AreaDefinition(fj.active_area_explicit_ghosts))
	# jet_selector_csa = fj.SelectorAbsEtaMax(max_eta_jet) & (~fj.SelectorIsPureGhost()) # used to get jets for C_area

	cms_rpa_fin = ROOT.TFile("cms_rpa.root")
	cms_rpa = cms_rpa_fin.Get("Table 16/Hist1D_y1")

	nbinsy = int(36.)
	lbinsy = logbins(1.e-3, 1., nbinsy)
	nbinsx = int(60.)
	lbinsx = linbins(20., 80., nbinsx)
	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()

	hjpt = ROOT.TH1F('hjpt', 'hjpt', 60, 20, 80)
	hjpt_rpa_br = ROOT.TH1F('hjpt_rpa_br', 'hjpt_rpa_br', 60, 20, 80)
	hjpt_rpa_bc = ROOT.TH1F('hjpt_rpa_bc', 'hjpt_rpa_bc', 60, 20, 80)
	hjpt_rpa_ar = ROOT.TH1F('hjpt_rpa_ar', 'hjpt_rpa_ar', 60, 20, 80)
	hjpt_rpa_ac = ROOT.TH1F('hjpt_rpa_ac', 'hjpt_rpa_ac', 60, 20, 80)

	yieldnbins = int(33)
	yieldbins = logbins(0.3, 200, yieldnbins)
	hyield = ROOT.TH1F('hyield', 'hyield', yieldnbins, yieldbins)
	hyield_rpa_br = ROOT.TH1F('hyield_rpa_br', 'hyield_rpa_br', yieldnbins, yieldbins)
	hyield_rpa_bc = ROOT.TH1F('hyield_rpa_bc', 'hyield_rpa_bc', yieldnbins, yieldbins)
	
	hNtracks = ROOT.TH1F('hNtracks', 'hNtracks', 50, 0, 100)
	hNtracks_rpa_br = ROOT.TH1F('hNtracks_rpa_br', 'hNtracks_rpa_br', 50, 0, 100)
	hNtracks_rpa_bc = ROOT.TH1F('hNtracks_rpa_bc', 'hNtracks_rpa_bc', 50, 0, 100)
	
	hZ = ROOT.TH1F('z_ch','z_ch', 10, 0, 1)
	hZ_rpa_br = ROOT.TH1F('z_ch_rpa_br','z_ch_rpa_br', 10, 0, 1)
	hZ_rpa_bc = ROOT.TH1F('z_ch_rpa_bc','z_ch_rpa_bc', 10, 0, 1)
	hZ_rpa_ar = ROOT.TH1F('z_ch_rpa_ar','z_ch_rpa_ar', 10, 0, 1)
	hZ_rpa_ac = ROOT.TH1F('z_ch_rpa_ac','z_ch_rpa_ac', 10, 0, 1)
		
	hec0 = []
	hec1 = []
	hec1_rpa_br = [] 
	hec1_rpa_bc = [] 
	hec1_rpa_ar = [] 
	hec1_rpa_ac = [] 
	if args.ncorrel < 2:
		args.ncorrel = 2
	if args.ncorrel > 5:
		args.ncorrel = 5
	print('[i] n correl up to', args.ncorrel)
	for i in range(args.ncorrel - 1):
		h = ROOT.TH2F('hec0_{}'.format(i+2), 'hec0_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec0.append(h)
		h = ROOT.TH2F('hec1_{}'.format(i+2), 'hec1_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1.append(h)
		h = ROOT.TH2F('hec1_rpa_br_{}'.format(i+2), 'hec1_rpa_br_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_rpa_br.append(h)
		h = ROOT.TH2F('hec1_rpa_bc_{}'.format(i+2), 'hec1_rpa_bc_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_rpa_bc.append(h)
		h = ROOT.TH2F('hec1_rpa_ar_{}'.format(i+2), 'hec1_rpa_ar_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_rpa_ar.append(h)
		h = ROOT.TH2F('hec1_rpa_ac_{}'.format(i+2), 'hec1_rpa_ac_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_rpa_ac.append(h)

	for n in tqdm(range(args.nev)):
		if not pythia_hard.next():
				continue
		parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		# parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)
		# parts_pythia_h_selected_charged = [p for p in parts_pythia_h_selected if pythiafjext.getPythia8Particle(p).isCharged()]
		mult = len(parts_pythia_h_selected)
		hNtracks.Fill(mult)
		if mult != 0:
			hpt = [h.pt() for h in parts_pythia_h_selected]
		else:
			mult = 1
			hpt = [-1]
		hyield.FillN(mult, array.array('d',hpt), array.array('d',[1]*mult))

		rpa_mask_before = [hadron_efficiency(cms_rpa, h.pt()) for h in parts_pythia_h_selected]
		parts_pythia_h_rpa_bc = np.repeat(parts_pythia_h_selected, rpa_mask_before)
		parts_pythia_h_rpa_br = randomize_added(parts_pythia_h_selected, rpa_mask_before)

		mult_rpa_bc = len(parts_pythia_h_rpa_bc)
		hNtracks_rpa_bc.Fill(mult_rpa_bc)
		if mult_rpa_bc != 0:
			hpt_rpa_bc = [h.pt() for h in parts_pythia_h_rpa_bc]
		else:
			hpt_rpa_bc = [-1]
			mult_rpa_bc = 1
		hyield_rpa_bc.FillN(mult_rpa_bc, array.array('d',hpt_rpa_bc), array.array('d',[1]*mult_rpa_bc))
		
		mult_rpa_br = len(parts_pythia_h_rpa_br)
		hNtracks_rpa_br.Fill(mult_rpa_br)
		if mult_rpa_br != 0:
			hpt_rpa_br = [h.pt() for h in parts_pythia_h_rpa_br]
		else:
			hpt_rpa_bc = [-1]
			mult_rpa_br = 1
		hyield_rpa_br.FillN(mult_rpa_br, array.array('d',hpt_rpa_br), array.array('d',[1]*mult_rpa_br))
		if mult_rpa_br != mult_rpa_bc:
			print(hpt)
			print(hpt_rpa_br)
			print(hpt_rpa_bc)

		jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected))) 
		if len(jets_h) < 1:
				continue
		
		jets_rpa_bc = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_rpa_bc))) 
		jets_rpa_br = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_rpa_br))) 
		
		for j in jets_h:
			hjpt.Fill(j.perp())

			# alternative: push constutents to a vector in python
			_vc = fj.vectorPJ()
			_ = [_vc.push_back(c) for c in j.constituents() if c.user_index() != -1]

			n = len(_vc)
			hZ.FillN(n, array.array('d', [c.perp()/j.perp() for c in _vc]), 
						array.array('d', [1]*n))
			# n-point correlator with all charged particles
			cb = ecorrel.CorrelatorBuilder(_vc, j.perp(), args.ncorrel, 1, -9999, -9999)
   
			# select only charged constituents with 1 GeV cut
			_vc1 = fj.vectorPJ()
			_ = [_vc1.push_back(c) for c in pfc_selector1(_vc)]
			# n-point correlator with charged particles pt > 1
			cb1 = ecorrel.CorrelatorBuilder(_vc1, j.perp(), args.ncorrel, 1, -9999, -9999)
   
			for i in range(args.ncorrel - 1):
				if cb.correlator(i+2).rs().size() > 0:
					n = cb.correlator(i+2).rs().size()
					hec0[i].FillN(	n, array.array('d', [j.perp()]*n),
									array.array('d', cb.correlator(i+2).rs()), 
									array.array('d', cb.correlator(i+2).weights()))
				if cb1.correlator(i+2).rs().size() > 0:
					n = cb1.correlator(i+2).rs().size()
					hec1[i].FillN(	n, array.array('d', [j.perp()]*n),
									array.array('d', cb1.correlator(i+2).rs()), 
									array.array('d', cb1.correlator(i+2).weights()))

			rpa_eec(j, _vc, hjpt_rpa_ar, hZ_rpa_ar, hec1_rpa_ar, cms_rpa, False)
			rpa_eec(j, _vc, hjpt_rpa_ac, hZ_rpa_ac, hec1_rpa_ac, cms_rpa, True)
					
		for j in jets_rpa_br:
			fill_eec(j, hjpt_rpa_br, hZ_rpa_br, hec1_rpa_br)
		for j in jets_rpa_bc:
			fill_eec(j, hjpt_rpa_bc, hZ_rpa_bc, hec1_rpa_bc)
		
	njets = hjpt.Integral()
	if njets == 0:
		njets = 1.

	fout.cd()

	pythia_hard.stat()

	fout.Write()
	fout.Close()
	print('# jets:', njets)
	print('[i] written ', fout.GetName())

def fill_eec(j, hjpt, hZ, hec1):
	pfc_selector1 = fj.SelectorPtMin(1.)
	hjpt.Fill(j.perp())

	# alternative: push constutents to a vector in python
	_vc = fj.vectorPJ()
	_ = [_vc.push_back(c) for c in j.constituents() if c.user_index() != -1]
	
	n = len(_vc)
	hZ.FillN(n, array.array('d', [c.perp()/j.perp() for c in _vc]), 
				array.array('d', [1]*n))

	# select only charged constituents with 1 GeV cut
	_vc1 = fj.vectorPJ()
	_ = [_vc1.push_back(c) for c in pfc_selector1(_vc)]
	# n-point correlator with charged particles pt > 1
	cb1 = ecorrel.CorrelatorBuilder(_vc1, j.perp(), 2, 1, -9999, -9999)

	for i in range(2 - 1):
		if cb1.correlator(i+2).rs().size() > 0:
			n = cb1.correlator(i+2).rs().size()
			hec1[i].FillN(	n, array.array('d', [j.perp()]*n),
							array.array('d', cb1.correlator(i+2).rs()), 
							array.array('d', cb1.correlator(i+2).weights()))
			
def rpa_eec(jet, vc, hjpt, hZ, hec1, cms_rpa, coll):
	_vc = vc
	pfc_selector1 = fj.SelectorPtMin(1.)

	rpa_mask = [hadron_efficiency(cms_rpa, h.pt()) for h in _vc]
	if coll:
		_vc_rpa = np.repeat(_vc, rpa_mask)
	else: 
		_vc_rpa = randomize_added(_vc, rpa_mask, jet=jet)
	jpt = np.sum([c.pt() for c in _vc_rpa])

	hjpt.Fill(jpt)

	n = len(_vc_rpa)
	hZ.FillN(n, array.array('d', [c.perp()/jpt for c in _vc_rpa]), 
				array.array('d', [1]*n))
	
	# select only charged constituents with 1 GeV cut
	_vc1_rpa = fj.vectorPJ()
	_ = [_vc1_rpa.push_back(c) for c in pfc_selector1(_vc_rpa)]
	# n-point correlator with charged particles pt > 1
	cb1_rpa = ecorrel.CorrelatorBuilder(_vc1_rpa, jpt, 2, 1, -9999, -9999)

	for i in range(2- 1):
		if cb1_rpa.correlator(i+2).rs().size() > 0:
			n = cb1_rpa.correlator(i+2).rs().size()
			hec1[i].FillN(	n, array.array('d', [jpt]*n),
							array.array('d', cb1_rpa.correlator(i+2).rs()), 
							array.array('d', cb1_rpa.correlator(i+2).weights()))

if __name__ == '__main__':
	main()
