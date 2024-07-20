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
	parser.add_argument('--background', help='add background', action='store_true')
	parser.add_argument('--rhosub', help='do rho subtraction', action='store_true')
	parser.add_argument('--perpcone', help='do perpcone subtraction', action='store_true')
	parser.add_argument('--bgrho', help='do thermal rho subtraction', action='store_true')
	parser.add_argument('--bgperp', help='do thermal perpcone subtraction', action='store_true')
	parser.add_argument('--output', default="test_ecorrel_rw.root", type=str)
	parser.add_argument('--user-seed', help='pythia seed',
											default=1111, type=int)
	args = parser.parse_args()
	print(args)
	return args


def main():
	mycfg = []
	ssettings = "--py-ecm 5000 --py-pthatmin 60"
	args = get_args_from_settings(ssettings)
	pythia_hard = pyconf.create_and_init_pythia_from_args(args, mycfg)
	jet_R0 = 0.4
	max_eta_jet = 0.9
	max_eta_hadron = max_eta_jet + jet_R0
	parts_selector_h = fj.SelectorAbsEtaMax(max_eta_hadron)
	jet_selector = fj.SelectorPtMin(60.0) & fj.SelectorAbsEtaMax(max_eta_jet) & (~fj.SelectorIsPureGhost())
	pfc_selector0 = fj.SelectorPtMin(0.)
	pfc_selector1 = fj.SelectorPtMin(1.)

	# print the banner first
	fj.ClusterSequence.print_banner()
	print()
	# set up our jet definition and a jet selector
	jet_def = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def)
	jet_def_bg = fj.JetDefinition(fj.antikt_algorithm, jet_R0)
	print(jet_def_bg)
	jet_selector_bg = fj.SelectorAbsEtaMax(max_eta_jet)

	jet_def_sub = fj.JetDefinition(fj.kt_algorithm, jet_R0)
	jet_selector_sub = fj.SelectorAbsEtaMax(max_eta_jet) & (~fj.SelectorNHardest(2)) & (~fj.SelectorIsPureGhost()) # used to calculate rho
	median_subtractor = fj.JetMedianBackgroundEstimator(jet_selector_sub, jet_def_sub, fj.AreaDefinition(fj.active_area_explicit_ghosts))
	jet_selector_csa = fj.SelectorAbsEtaMax(max_eta_jet) & (~fj.SelectorIsPureGhost()) # used to get jets for C_area

	nbinsy = int(36.)
	lbinsy = logbins(1.e-3, 1., nbinsy)
	nbinsx = int(60.)
	lbinsx = linbins(20., 80., nbinsx)
	fout = ROOT.TFile(args.output, 'recreate')
	fout.cd()

	hjpt = ROOT.TH1F('hjpt', 'hjpt', 60, 20, 80)
	hjpt_sub = ROOT.TH1F('hjpt_sub', 'hjpt_sub', 60, 20, 80)
	hjpt_bg = ROOT.TH1F('hjpt_bg', 'hjpt_bg', 50, 0, 100)
	hjpt_m = ROOT.TH1F('hjpt_m', 'hjpt_m', 50, 0, 100)
	hjpt_bg_sub = ROOT.TH1F('hjpt_bg_sub', 'hjpt_bg_sub', 50, 0, 100)
	
	hNtracks = ROOT.TH1F('hNtracks', 'hNtracks', 50, 0, 100)
	hjpt_nperp = ROOT.TH2F('hjpt_nperp', 'hjpt_nperp', 60, 20, 80, 20, 0, 20)
	hnjet_nperp_2040 = ROOT.TH2F('hnjet_nperp_2040', 'hnjet_nperp_2040', 20, 0, 20, 15, 0, 15)
	hnjet_nperp_4060 = ROOT.TH2F('hnjet_nperp_4060', 'hnjet_nperp_4060', 20, 0, 20, 15, 0, 15)
	hnjet_nperp_6080 = ROOT.TH2F('hnjet_nperp_6080', 'hnjet_nperp_6080', 20, 0, 20, 15, 0, 15)
	
	do_rho = False
	do_perp = False
	if args.rhosub:
		do_rho = True
	if args.perpcone:
		do_perp = True
		
	hec0 = []
	hec1 = []
	hec1_sub = []
	hec1_perp = []
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
		h = ROOT.TH2F('hec1_sub_{}'.format(i+2), 'hec1_sub_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_sub.append(h)
		h = ROOT.TH2F('hec1_perp_{}'.format(i+2), 'hec1_perp_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
		hec1_perp.append(h)
 
	be = None
	if args.background:
		be = BoltzmannEvent(mean_pt=0.7, multiplicity=40 * max_eta_hadron * 2, max_eta=max_eta_hadron, max_pt=100)
		print(be)
	do_bgrho = False
	do_bgperp = False
	if args.bgrho:
		do_bgrho = True
	if args.bgperp:
		do_bgperp = True

	hbgtrks = ROOT.TH1F('hbgtrks', 'hbgtrks', 10, -4, 6)
	
	hecm_bg = []
	hec1_bg = []
	hec1_bg_sub = []
	hec1_bg_perp = []
	if be:
		for i in range(args.ncorrel - 1):
			h = ROOT.TH2F('hecm_bg_{}'.format(i+2), 'hecm_bg_{}'.format(i+2), 50, linbins(0.,100.,50), nbinsy, lbinsy)
			hecm_bg.append(h)
			h = ROOT.TH2F('hec1_bg_{}'.format(i+2), 'hec1_bg_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
			hec1_bg.append(h)
			h = ROOT.TH2F('hec1_bg_sub_{}'.format(i+2), 'hec1_bg_sub_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
			hec1_bg_sub.append(h)
			h = ROOT.TH2F('hec1_bg_perp_{}'.format(i+2), 'hec1_bg_perp_{}'.format(i+2), nbinsx, lbinsx, nbinsy, lbinsy)
			hec1_bg_perp.append(h)

	for n in tqdm(range(args.nev)):
		if not pythia_hard.next():
				continue
		# parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal, pythiafjext.kCharged], 0, True)
		parts_pythia_h = pythiafjext.vectorize_select(pythia_hard, [pythiafjext.kFinal], 0, True)
		parts_pythia_h_selected = parts_selector_h(parts_pythia_h)
		parts_pythia_h_selected_charged = [p for p in parts_pythia_h_selected if pythiafjext.getPythia8Particle(p).isCharged()]
		hNtracks.Fill(len(parts_pythia_h_selected_charged))

		csa = fj.ClusterSequenceArea(parts_pythia_h_selected, jet_def, fj.AreaDefinition(fj.active_area_explicit_ghosts))
		jets_h = fj.sorted_by_pt(jet_selector(csa.inclusive_jets())) 
		if len(jets_h) < 1:
				continue
		
		if do_rho:
			csa_kt = fj.ClusterSequenceArea(parts_pythia_h_selected, jet_def_sub, fj.AreaDefinition(fj.active_area_explicit_ghosts))
			median_subtractor.set_cluster_sequence(csa_kt)
			jets_kt = fj.sorted_by_pt(jet_selector_csa(csa_kt.inclusive_jets()))		
			rho = median_subtractor.rho() * np.sum(np.array([j.area() for j in jets_kt])) / (2*np.pi*1.8)

		# print('njets:', hjpt.Integral())
		
		for j in jets_h:
			hjpt.Fill(j.perp())
			# note: the EEC,... takes vector<PseudoJet> while PseudoJet::constituents() returns a tuple in python
			# so we use a helper function (in SWIG only basic types typles handled easily...)
			# vconstits = ecorrel.constituents_as_vector(j)
			# eecs_alt = ecorrel.EEC(vconstits, scale=j.perp())   
			# e3cs_alt = ecorrel.E3C(vconstits, scale=j.perp())
			# e4cs_alt = ecorrel.E4C(vconstits, scale=j.perp())

			# alternative: push constutents to a vector in python
			_v = fj.vectorPJ()
			_ = [_v.push_back(c) for c in j.constituents() if c.user_index() != -1]

			# select only charged constituents
			_vc = fj.vectorPJ()
			_ = [_vc.push_back(c) for c in _v
							if pythiafjext.getPythia8Particle(c).isCharged()]
			# n-point correlator with all charged particles
			cb = ecorrel.CorrelatorBuilder(_v, j.perp(), args.ncorrel, 1, -9999, -9999)
   
			# select only charged constituents with 1 GeV cut
			_vc1 = fj.vectorPJ()
			_ = [_vc1.push_back(c) for c in pfc_selector1(_v)
							if pythiafjext.getPythia8Particle(c).isCharged()]
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
			
			# median subtraction
			if do_rho:
				jpt_corr = j.perp() - rho*j.area()
				hjpt_sub.Fill(jpt_corr)
				cb1_sub = ecorrel.CorrelatorBuilder(_vc1, jpt_corr, args.ncorrel, 1, -9999, -9999)
				for i in range(args.ncorrel - 1):
					if cb1_sub.correlator(i+2).rs().size() > 0:
						n = cb1_sub.correlator(i+2).rs().size()
						hec1_sub[i].FillN(	n, array.array('d', [jpt_corr]*n),
											array.array('d', cb1_sub.correlator(i+2).rs()), 
											array.array('d', cb1_sub.correlator(i+2).weights()))
				
				if do_perp:
					h_aux = [hjpt_nperp, hnjet_nperp_2040, hnjet_nperp_4060, hnjet_nperp_6080]
					analyze_perp_cones(parts_pythia_h_selected, j, _vc1, jpt_corr, jet_R0, hec1_perp, h_aux=h_aux)
			
			if be:
				# print('add background')
				bg_parts = be.generate(offset=10000)
				full_event = fj.vectorPJ()
				_tmp = [full_event.push_back(c) for c in pfc_selector1(j.constituents())
											if pythiafjext.getPythia8Particle(c).isCharged()]
				_tmp = [full_event.push_back(psj) for psj in bg_parts]

				# print('nparts in the jet', len(j.constituents()), 'nparts in the hybrid event', full_event.size())

				csa_bg = fj.ClusterSequenceArea(full_event, jet_def_bg, fj.AreaDefinition(fj.active_area_explicit_ghosts))
				jets_hbg = fj.sorted_by_pt(jet_selector_bg(csa_bg.inclusive_jets())) 

				jbg_selected = None
				for jbg in jets_hbg:
					mpt = fjtools.matched_pt(jbg, j)
					if mpt < 0.8:
						continue
					jbg_selected = jbg
				
				if jbg_selected is None:
					# print('match eta={} phi={} - no match found'.format(j.eta(), j.phi()))
					continue
				else:
					# print('match eta={} phi={} - matched to eta={} phi={}'.format(j.eta(), j.phi(), jbg_selected.eta(), jbg_selected.phi()))
					_vc1_bg = fj.vectorPJ()
					_ = [_vc1_bg.push_back(c) for c in pfc_selector1(jbg_selected.constituents())]
					
					hbgtrks.Fill(len(_vc1_bg)-len(_vc1))
					hjpt_bg.Fill(jbg_selected.perp())
					hjpt_m.Fill(j.perp())

					# n-point correlator with charged particles pt > 1
					cb1_bg = ecorrel.CorrelatorBuilder(_vc1_bg, jbg_selected.perp(), args.ncorrel, 1, -9999, -9999)
					cb1_m = ecorrel.CorrelatorBuilder(_vc1, j.perp(), args.ncorrel, 1, -9999, -9999)

					for i in range(args.ncorrel - 1):
						if cb1_bg.correlator(i+2).rs().size() > 0:
							n = cb1_bg.correlator(i+2).rs().size()
							hec1_bg[i].FillN(	n, array.array('d', [jbg_selected.perp()]*n),
												array.array('d', cb1_bg.correlator(i+2).rs()), 
												array.array('d', cb1_bg.correlator(i+2).weights()))
						if cb1_m.correlator(i+2).rs().size() > 0:
							n = cb1_m.correlator(i+2).rs().size()
							hecm_bg[i].FillN(	n, array.array('d', [j.perp()]*n),
												array.array('d', cb1_m.correlator(i+2).rs()), 
												array.array('d', cb1_m.correlator(i+2).weights()))
					
					if do_bgrho:
						combined_event = fj.vectorPJ()
						_ = [combined_event.push_back(c) for c in parts_pythia_h_selected]
						_ = [combined_event.push_back(c) for c in bg_parts]
						# csa_bg_kt = fj.ClusterSequenceArea(combined_event, jet_def_sub, fj.AreaDefinition(fj.active_area_explicit_ghosts))
						# median_subtractor.set_cluster_sequence(csa_bg_kt)
						# jets_bg_kt = fj.sorted_by_pt(jet_selector_csa(csa_bg_kt.inclusive_jets()))		
						# rho = median_subtractor.rho() * np.sum(np.array([j.area() for j in jets_bg_kt])) / (2*np.pi*1.8)

						# jpt_bg_corr = jbg_selected.perp() - rho*jbg_selected.area()
						# hjpt_bg_sub.Fill(jpt_bg_corr)
						
						# cb1_bg_sub = ecorrel.CorrelatorBuilder(_vc1_bg, jpt_bg_corr, args.ncorrel, 1, -9999, -9999)
						# for i in range(args.ncorrel - 1):
						# 	if cb1_bg_sub.correlator(i+2).rs().size() > 0:
						# 		n = cb1_bg_sub.correlator(i+2).rs().size()
						# 		hec1_bg_sub[i].FillN(	n, array.array('d', [jpt_bg_corr]*n),
						# 							array.array('d', cb1_bg_sub.correlator(i+2).rs()), 
						# 							array.array('d', cb1_bg_sub.correlator(i+2).weights()))
						
						if do_bgperp:
							analyze_perp_cones(combined_event, j, _vc1_bg, jbg_selected.perp(), jet_R0, hec1_bg_perp)
		

	njets = hjpt.Integral()
	if njets == 0:
		njets = 1.
	njets_bg = hjpt_bg.Integral()

	fout.cd()

	for hg in [hec0, hec1, hec1_sub, hec1_perp, hec1_bg, hecm_bg, hec1_bg_sub, hec1_bg_perp]:
		for i in range(args.ncorrel - 1):
			hg[i].Sumw2()
			# intg = hg[i].Integral()
			# if intg > 0:
			# 	hg[i].Scale(1./intg)
			# if i > 0:
			# 	fout.cd()
			# 	hc = hg[i].Clone(hg[i].GetName() + '_ratio_to_EEC')
			# 	hc.Sumw2()
			# 	hc.Divide(hg[0])  
			# 	hc.Write()

	pythia_hard.stat()
	if be:
		be.histogram_pt.Write()
		be.histogram_phi.Write()
		be.histogram_eta.Write()

	fout.Write()
	fout.Close()
	print('# jets:', njets)
	print('# jets with background:', njets_bg)
	print('[i] written ', fout.GetName())

def analyze_perp_cones(parts, jet, _vc1, jpt_corr, jetR, h, h_aux=None):

	perp_jet1 = fj.PseudoJet()
	perp_jet1.reset_PtYPhiM(jet.pt(), jet.rapidity(), jet.phi() + np.pi/2, jet.m())
	perp_jet2 = fj.PseudoJet()
	perp_jet2.reset_PtYPhiM(jet.pt(), jet.rapidity(), jet.phi() - np.pi/2, jet.m())

	constituents = fj.vectorPJ()
	for c in _vc1:
		constituents.push_back(c)
	
	parts_in_perpcone1 = find_parts_around_jet(parts, perp_jet1, jetR)
	parts_in_perpcone1 = rotate_parts(parts_in_perpcone1, -np.pi/2)
	
	parts_in_perpcone2 = find_parts_around_jet(parts, perp_jet2, jetR)
	parts_in_perpcone2 = rotate_parts(parts_in_perpcone2, +np.pi/2)

	parts_in_cone1 = fj.vectorPJ()
	for part in constituents:
		part.set_user_index(1)
		parts_in_cone1.append(part)
	for part in parts_in_perpcone1:
		part.set_user_index(-1)
		parts_in_cone1.append(part)

	parts_in_cone2 = fj.vectorPJ()
	for part in constituents:
		part.set_user_index(1)
		parts_in_cone2.append(part)
	for part in parts_in_perpcone2:
		part.set_user_index(-1)
		parts_in_cone2.append(part)

	analyze_accepted_cone(parts_in_cone1, jet, jetR, jpt_corr, h, h_aux)
	analyze_accepted_cone(parts_in_cone2, jet, jetR, jpt_corr, h, h_aux)

def analyze_accepted_cone(cone_parts, jet, jetR, jpt_corr, h, h_aux=None):
	c_select = fj.vectorPJ()
	c_select_perp = fj.vectorPJ()
	cone_parts_sorted = fj.sorted_by_pt(cone_parts)
	for part in cone_parts_sorted:
		if part.pt() < 1:
			break
		c_select.append(part) # NB: use the break statement since constituents are already sorted
		if part.user_index() < 0:
			c_select_perp.append(part)

	if h_aux:
		nperp = len(c_select_perp)
		njet = len([c for c in jet.constituents() if c.perp() >= 1])
		h_aux[0].Fill(jpt_corr, nperp)
		if 20 <= jpt_corr < 40:
			h_aux[1].Fill(njet, nperp)
		elif 40 <= jpt_corr < 60:
			h_aux[2].Fill(njet, nperp)
		elif 60 <= jpt_corr < 80:
			h_aux[3].Fill(njet, nperp)

	cbp = ecorrel.CorrelatorBuilder(c_select, jpt_corr, 2, 1, -9999, -9999)		
	for i in range(2 - 1): #args.ncorrel = 2
		for index in range(cbp.correlator(2).rs().size()):
			pair_type_label = ''
			pair_type = check_pair_type(cbp, 2, c_select, index)
			pair_type = check_pair_type(cbp, 2, c_select, index)
			w = pair_type_weight(pair_type)
			h[i].Fill(jpt_corr, cbp.correlator(2).rs()[index], w * cbp.correlator(2).weights()[index])

def find_parts_around_jet(parts, jet, jetR):
	cone_parts = fj.vectorPJ()
	for part in parts:
		if jet.delta_R(part) <= jetR:
			cone_parts.push_back(part)

	cone_select = [c for c in cone_parts if c.pt()>1]
	if len(cone_select) >= 10:
		print("BIG PERP CONE", jet.pt())
		generate_event_display(parts, jet)
		for part in parts:
			if jet.delta_R(part) <= jetR:
				if part.pt() > 1 and pythiafjext.getPythia8Particle(part):
					print(part.pt(), pythiafjext.getPythia8Particle(part).id())
	return cone_parts

def rotate_parts(parts, rotate_phi):
	parts_rotated = fj.vectorPJ()
	for part in parts:
		pt_new = part.pt()
		y_new = part.rapidity()
		phi_new = part.phi() + rotate_phi
		m_new = part.m()
		part.reset_PtYPhiM(pt_new, y_new, phi_new, m_new)
		parts_rotated.push_back(part)
	return parts_rotated

def check_pair_type(corr_builder, ipoint, constituents, index):
	part1 = corr_builder.correlator(ipoint).indices1()[index]
	part2 = corr_builder.correlator(ipoint).indices2()[index]
	type1 = constituents[part1].user_index()
	type2 = constituents[part2].user_index()

	if type1 < 0 and type2 < 0:
		# print('bkg-bkg (',type1,type2,') pt1',constituents[part1].perp()
		return 0 # means bkg-bkg
	if type1 < 0 and type2 >= 0:
		# print('sig-bkg (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
		return 1 # means sig-bkg
	if type1 >= 0 and type2 < 0:
		# print('sig-bkg (',type1,type2,') pt1',constituents[part1].perp(),'pt2',constituents[part2].perp())
		return 1 # means sig-bkg
	if type1 >= 0 and type2 >= 0:
		# print('sig-sig (',type1,type2,') pt1',constituents[part1].perp()
		return 2 # means sig-sig

def pair_type_weight(pair_type):
	if pair_type % 2 == 0: # bb or ss
		return 1
	else: # sb
		return -1


def generate_event_display(event, jet):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	for particle in event:
		particle = pythiafjext.getPythia8Particle(particle)
		if particle:
			if particle.isCharged():
				color = 'b'  # Charged particles in blue
			else:
				color = 'r'  # Neutral particles in red
			ax.plot([0, particle.px()], [0, particle.py()], [0, particle.pz()], color)

	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.set_title('Event Display')

	plt.savefig('events/event_{:.2f}.png'.format(jet.pt()))

if __name__ == '__main__':
	main()
