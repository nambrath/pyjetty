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

def reweight(hz, z1, z2):
	b1 = int(z1*10)
	if b1 == 10: b1 = 9
	b2 = int(z2*10)
	if b2 == 10: b2 = 9
	zr1 =  hz[b1]
	zr2 =  hz[b2]
	return zr1*zr2

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

	hz_20 = [1.1475, 1.025, 0.9175,	0.87, 0.845, 0.7925, 0.8275, 0.725, 0.7025, 0.5325]
	hz_30 = [1.0966, 0.98495, 0.96796, 0.94126, 0.90971, 0.89272, 0.91942, 0.885437, 0.849029, 0.81262]

	nbinsy = int(36.)
	lbinsy = logbins(1.e-3, 1., nbinsy)
	nbinsx = int(60.)
	lbinsx = linbins(20., 80., nbinsx)
	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()

	hjpt = ROOT.TH1F('hjpt', 'hjpt', 60, 20, 80)

	hNtracks = ROOT.TH1F('hNtracks', 'hNtracks', 50, 0, 100)
	
	hec0 = []
	hec1 = []
	hec1_z = [] 
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
		h = ROOT.TH2F('hec1_z_{}'.format(i+2), 'hec1_z_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_z.append(h)

	for n in tqdm(range(args.nev)):
		if not pythia_hard.next():
				continue
		parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		# parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)
		# parts_pythia_h_selected_charged = [p for p in parts_pythia_h_selected if pythiafjext.getPythia8Particle(p).isCharged()]
		mult = len(parts_pythia_h_selected)
		hNtracks.Fill(mult)

		jets_h = fj.sorted_by_pt(jet_selector(jet_def(parts_pythia_h_selected))) 
		if len(jets_h) < 1:
				continue
		
		for j in jets_h:
			hjpt.Fill(j.perp())

			# alternative: push constutents to a vector in python
			_vc = fj.vectorPJ()
			_ = [_vc.push_back(c) for c in j.constituents() if c.user_index() != -1]

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
					
					if 20 < j.perp() < 40:
						for index in range(n):
							part1 = cb1.correlator(i+2).indices1()[index]
							part2 = cb1.correlator(i+2).indices2()[index]
							z1 = _vc1[part1].perp() / j.perp()
							z2 = _vc1[part2].perp() / j.perp()
							if j.perp() < 30: hz = hz_20
							if 30 <= j.perp(): hz = hz_30
							hec1_z[i].Fill(j.perp(), cb1.correlator(i+2).rs()[index], 
											cb1.correlator(2).weights()[index] * reweight(hz, z1, z2))	
		
	njets = hjpt.Integral()
	if njets == 0:
		njets = 1.

	fout.cd()

	pythia_hard.stat()

	fout.Write()
	fout.Close()
	print('# jets:', njets)
	print('[i] written ', fout.GetName())


if __name__ == '__main__':
	main()
