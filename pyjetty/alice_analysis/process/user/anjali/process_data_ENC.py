#!/usr/bin/env python3

"""
	Analysis class to read a ROOT TTree of track information
	and do jet-finding, and save basic histograms.
	
	Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse
import sys

# Data analysis and plotting
import ROOT
import yaml
import numpy as np
import array 
import math

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjtools
import ecorrel

# Base class
from pyjetty.alice_analysis.process.user.substructure import process_data_base

def linbins(xmin, xmax, nbins):
	lspace = np.linspace(xmin, xmax, nbins+1)
	arr = array.array('f', lspace)
	return arr

def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr

################################################################
class ProcessData_ENC(process_data_base.ProcessDataBase):

	#---------------------------------------------------------------
	# Constructor
	#---------------------------------------------------------------
	def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
	
		# print(sys.path)
		# print(sys.modules)

		# Initialize base class
		super(ProcessData_ENC, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
		
		self.observable = self.observable_list[0]


	#---------------------------------------------------------------
	# Initialize histograms
	#---------------------------------------------------------------
	def initialize_user_output_objects(self):

		for jetR in self.jetR_list:
			for observable in self.observable_list:
				for trk_thrd in self.obs_settings[observable]:

					obs_label = self.utils.obs_label(trk_thrd, None) 
					if self.is_pp or self.is_pA:
							# Init ENC histograms
							if 'ENC' in observable:
								for ipoint in range(2, 3):
										name = 'h_{}_JetPt_R{}_{}'.format(observable + str(ipoint), jetR, trk_thrd)
										pt_bins = linbins(0,200,200)
										RL_bins = logbins(1E-4,1,50)
										h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
										h.GetXaxis().SetTitle('p_{T,ch jet}')
										h.GetYaxis().SetTitle('R_{L}')
										setattr(self, name, h)

										name = 'h_{}Pt_JetPt_R{}_{}'.format(observable + str(ipoint), jetR, trk_thrd)
										pt_bins = linbins(0,200,200)
										ptRL_bins = logbins(1E-3,1E2,60)
										h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
										h.GetXaxis().SetTitle('p_{T,ch jet}')
										h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
										setattr(self, name, h)

							name = 'hSigJetpt_Partpt_R{}_{}'.format(jetR,obs_label)
							h = ROOT.TH2F(name, name, 100, 0, 100., 200, 0, 50.)
							setattr(self, name, h)

							self.mult_labels = ['']
							if self.mult_threshold != 0:
								self.mult_labels = ['_lm','_hm']
								if type(self.mult_threshold) is list:
									self.mult_labels.append('_mm')


							if 'EEC_noweight' in observable:
								name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)

							if 'EEC_weight2' in observable: # NB: weight power = 2
								name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)

							if 'jet_pt' in observable:
								name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								h = ROOT.TH1D(name, name, 200, pt_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('Counts')
								setattr(self, name, h)

								for mult_label in self.mult_labels:
									name = 'h_{}_JetPt_R{}_{}'.format(observable + mult_label, jetR, obs_label)
									pt_bins = linbins(0,200,200)
									h = ROOT.TH1D(name, name, 200, pt_bins)
									h.GetXaxis().SetTitle('p_{T,ch jet}')
									h.GetYaxis().SetTitle('Counts')
									setattr(self, name, h)
							
							if self.do_perpendicular_cone:
								name = 'h_perpEEC_JetPt_R{}_{}'.format(jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)
								
								# background subtraction histograms
								perpcone_R_list = [jetR]
								self.pair_type_labels = ['']
								if self.do_median_subtraction:
									self.pair_type_labels = ['_bb','_sb','_ss']

								for perpcone_R in perpcone_R_list:
									if 'jet_pt' in observable:
										name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, 'pt', jetR, trk_thrd)
										pt_bins = linbins(-200,200,400)
										h = ROOT.TH1D(name, name, 200, pt_bins)
										h.GetXaxis().SetTitle('p_{T,perp cone}')
										h.GetYaxis().SetTitle('Counts')
										setattr(self, name, h)

										name = 'h_perpcone{}_Nconst_JetPt_R{}_{}'.format(perpcone_R, jetR, trk_thrd)
										pt_bins = linbins(-200,200,400)
										Nconst_bins = linbins(0,50,50)
										h = ROOT.TH2D(name, name, 200, pt_bins, 50, Nconst_bins)
										h.GetXaxis().SetTitle('p_{T,ch jet}')
										h.GetYaxis().SetTitle('N_{const}')
										setattr(self, name, h)

									for pair_type_label in self.pair_type_labels:

										for mult_label in self.mult_labels:

											if 'ENC' in observable:
												for ipoint in range(2, 3):
													name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label + mult_label, jetR, trk_thrd)
													pt_bins = linbins(0,200,200)
													RL_bins = logbins(1E-4,1,50)
													h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
													h.GetXaxis().SetTitle('p_{T,ch jet}')
													h.GetYaxis().SetTitle('R_{L}')
													setattr(self, name, h)

													name = 'h_perpcone{}_{}Pt_JetPt_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label + mult_label, jetR, trk_thrd)
													pt_bins = linbins(0,200,200)
													ptRL_bins = logbins(1E-3,1E2,60)
													h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
													h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
													setattr(self, name, h)

											if 'EEC_noweight' in observable or 'EEC_weight2' in observable:
												name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, observable + pair_type_label + mult_label, jetR, obs_label)
												pt_bins = linbins(0,200,200)
												RL_bins = logbins(1E-4,1,50)
												h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
												h.GetXaxis().SetTitle('p_{T,ch jet}')
												h.GetYaxis().SetTitle('R_{L}')
												setattr(self, name, h)

							if self.do_reshuffle:
								name = 'h_reshuffle_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)
										
	#---------------------------------------------------------------
	# Calculate pair distance of two fastjet particles
	#---------------------------------------------------------------
	def calculate_distance(self, p0, p1):   
		dphiabs = math.fabs(p0.phi() - p1.phi())
		dphi = dphiabs

		if dphiabs > math.pi:
			dphi = 2*math.pi - dphiabs

		deta = p0.eta() - p1.eta()
		return math.sqrt(deta*deta + dphi*dphi)

	#---------------------------------------------------------------
	# Rotate tracks inside jet
	#---------------------------------------------------------------
	def reshuffle_parts(self, parts, jet, jetR):
		# rotate parts in azimuthal direction
		parts_reshuffled = fj.vectorPJ()
		for part in parts:
			pt_new = part.pt()
			m_new = part.m()
			R_new = 1 # initialize to a big radius
			while R_new > jetR:
				phi_new = np.random.uniform(-jetR,+jetR)
				y_new = np.random.uniform(-jetR,+jetR)
				R_new = math.sqrt(phi_new*phi_new+y_new*y_new)
			phi_new = part.phi() + phi_new
			y_new = part.rapidity() + y_new
			
			# print('before',part.phi())
			part.reset_PtYPhiM(pt_new, y_new, phi_new, m_new)
			# print('after',part.phi())
			parts_reshuffled.push_back(part)
		
		return parts_reshuffled
	
	def check_pair_type(self, corr_builder, ipoint, constituents, index):
		part1 = corr_builder.correlator(ipoint).indices1()[index]
		part2 = corr_builder.correlator(ipoint).indices2()[index]
		type1 = constituents[part1].user_index()
		type2 = constituents[part2].user_index()

		# NB: match the strings in self.pair_type_label = ['bb','sb','ss']
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


	#---------------------------------------------------------------
	# This function is called once for each jet subconfiguration
	#---------------------------------------------------------------
	def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
													obs_label, jet_pt_subtracted, suffix):

		constituents = fj.sorted_by_pt(jet.constituents())
		c_select = fj.vectorPJ()
		trk_thrd = obs_setting

		for c in constituents:
			if c.pt() < trk_thrd:
				break
			getattr(self, 'hSigJetpt_Partpt_R{}_{}'.format(jetR,obs_label)).Fill(jet.perp(), c.pt())
			c_select.append(c) # NB: use the break statement since constituents are already sorted

		if self.ENC_pair_cut:
			dphi_cut = -9999 # means no dphi cut
			deta_cut = 0.008
		else:
			dphi_cut = -9999
			deta_cut = -9999

		jet_pt = jet_pt_subtracted

		hname = 'h_{}_JetPt_R{}_{}'
		new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
		for observable in self.observable_list:
			if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
				for ipoint in range(2, 3):
					for index in range(new_corr.correlator(ipoint).rs().size()):

						if 'ENC' in observable:
							getattr(self, hname.format(observable + str(ipoint), jetR, obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
							getattr(self, hname.format(observable + str(ipoint) + 'Pt', jetR, obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: fill pt*RL

						if ipoint==2 and 'EEC_noweight' in observable:
							getattr(self, hname.format(observable, jetR, obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

						if ipoint==2 and 'EEC_weight2' in observable:
							getattr(self, hname.format(observable, jetR, obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

			if 'jet_pt' in observable:
				getattr(self, hname.format(observable, jetR, obs_label)).Fill(jet_pt)  
				if self.mult_threshold != 0:
					if self.isHighMult: 
						mult_label = self.mult_labels[1]
					elif self.isLowMult: 
						mult_label = self.mult_labels[0]
					else:
						mult_label = self.mult_labels[2]
					getattr(self, hname.format(observable + mult_label, jetR, obs_label)).Fill(jet_pt)  
			
			if self.do_reshuffle:
				_c_reshuffle = self.reshuffle_parts(c_select, jet, jetR)
				cb_reshuffle = ecorrel.CorrelatorBuilder(_c_reshuffle, jet_pt, 2, 1, dphi_cut, deta_cut)

				for ipoint in range(2, 3):
						for index in range(cb_reshuffle.correlator(ipoint).rs().size()):
										getattr(self, 'h_reshuffle_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)).Fill(jet_pt, cb_reshuffle.correlator(ipoint).rs()[index], cb_reshuffle.correlator(ipoint).weights()[index])


	#---------------------------------------------------------------
	# This function is called once for each jet subconfiguration
	#---------------------------------------------------------------
	def fill_perpcone_histograms(self, jet, jetR, perp_cone_particles, obs_setting, obs_label, suffix, rho_bge=0):

		trk_thrd = obs_setting
		c_select = fj.vectorPJ()

		for p in fj.sorted_by_pt(perp_cone_particles):
			if p.pt() < trk_thrd:
				break
			c_select.append(p) # NB: use the break statement since constituents are already sorted

		if self.ENC_pair_cut:
			dphi_cut = -9999 # means no dphi cut
			deta_cut = 0.008
		else:
			dphi_cut = -9999
			deta_cut = -9999

		hname = 'h_perpEEC_JetPt_R{}_{}'
		new_corr = ecorrel.CorrelatorBuilder(c_select, jet.perp()-rho_bge*jet.area(), 2, 1, dphi_cut, deta_cut)
		for ipoint in range(2, 3):
			for index in range(new_corr.correlator(ipoint).rs().size()):
				getattr(self, hname.format(jetR, obs_label)).Fill(jet.perp()-rho_bge*jet.area(), new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])   

	#---------------------------------------------------------------
	# Perp cone background pair subtraction histograms
	#---------------------------------------------------------------
	def fill_perp_cone_histograms(self, cone_parts, cone_R, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed, suffix, rho_bge = 0):

		# calculate perp cone pt after subtraction. Notice that the perp cone already contain the particles from signal "jet". Signal and background can be identified using user_index()
		cone_px = 0
		cone_py = 0
		cone_npart = 0
		for part in cone_parts:
			if part.user_index() < 0:
				cone_px = cone_px + part.px()
				cone_py = cone_py + part.py()
				cone_npart = cone_npart + 1
		cone_pt = math.sqrt(cone_px*cone_px + cone_py*cone_py)
		# print('cone pt', cone_pt-rho_bge*jet.area(), '(', cone_pt, ')')
		cone_pt = cone_pt-rho_bge*jet.area() # ideally this should fluctuate around 0
		# print('jet pt', jet_pt_ungroomed, '(', jet.perp(), ')')

		# combine sig jet and perp cone with trk threshold cut
		trk_thrd = obs_setting
		c_select = fj.vectorPJ()
		c_select_perp = fj.vectorPJ()

		cone_parts_sorted = fj.sorted_by_pt(cone_parts)
		# print('perp cone nconst:',len(cone_parts_sorted))
		for part in cone_parts_sorted:
			if part.pt() < trk_thrd:
				break
			c_select.append(part) # NB: use the break statement since constituents are already sorted
			if part.user_index() < 0:
				c_select_perp.append(part)

		nconst_perp = len(c_select_perp)
		# print('cone R',cone_R)
		# print('total cone nconst (with thrd cut):',len(c_select))
		# print('perp cone nconst (with thrd cut):',nconst_perp)

		if self.ENC_pair_cut:
			dphi_cut = -9999 # means no dphi cut
			deta_cut = 0.008
		else:
			dphi_cut = -9999
			deta_cut = -9999

		hname = 'h_perpcone{}_{}_JetPt_R{}_{}{}'
		if self.do_median_subtraction:
			jet_pt = jet_pt_ungroomed # jet_pt_ungroomed stores subtracted jet pt for energy weight calculation and pt selection for there is a non-zero UE energy density
		else:
			jet_pt = jet.perp()

		new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
		for observable in self.observable_list:

			if 'jet_pt' in observable:
				getattr(self, hname.format(cone_R, 'pt', jetR, obs_label, suffix)).Fill(cone_pt)
				getattr(self, hname.format(cone_R, 'Nconst', jetR, obs_label, suffix)).Fill(jet_pt, nconst_perp)

			if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
				for ipoint in range(2, 3):
					for index in range(new_corr.correlator(ipoint).rs().size()):

						# # processing only like-sign pairs when self.ENC_pair_like is on
						# if self.ENC_pair_like and (not self.is_same_charge(new_corr, ipoint, c_select, index)):
						# 	continue

						# # processing only unlike-sign pairs when self.ENC_pair_unlike is on
						# if self.ENC_pair_unlike and self.is_same_charge(new_corr, ipoint, c_select, index):
						# 	continue

						# separate out sig-sig, sig-bkg, bkg-bkg correlations for EEC pairs
						pair_type_label = ''
						mult_label = ''
						if self.do_median_subtraction:
							pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
							pair_type_label = self.pair_type_labels[pair_type]
							if self.mult_threshold != 0:	
								if self.isHighMult: 
									mult_label = self.mult_labels[1]
								elif self.isLowMult: 
									mult_label = self.mult_labels[0]
								else:
									mult_label = self.mult_labels[2]
								

						if 'ENC' in observable:
							# print('hname is',hname.format(cone_R, observable + str(ipoint) + pair_type_label, jetR, obs_label))
							getattr(self, hname.format(cone_R, observable + str(ipoint) + pair_type_label + mult_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
							getattr(self, hname.format(cone_R, observable + str(ipoint) + pair_type_label + mult_label + 'Pt', jetR, obs_label, suffix)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: fill pt*RL

						if ipoint==2 and 'EEC_noweight' in observable:
							getattr(self, hname.format(cone_R, observable + pair_type_label + mult_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

						if ipoint==2 and 'EEC_weight2' in observable:
							getattr(self, hname.format(cone_R, observable + pair_type_label + mult_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))


##################################################################
if __name__ == '__main__':
	# Define arguments
	parser = argparse.ArgumentParser(description='Process data')
	parser.add_argument('-f', '--inputFile', action='store',
											type=str, metavar='inputFile',
											default='AnalysisResults.root',
											help='Path of ROOT file containing TTrees')
	parser.add_argument('-c', '--configFile', action='store',
											type=str, metavar='configFile',
											default='config/analysis_config.yaml',
											help="Path of config file for analysis")
	parser.add_argument('-o', '--outputDir', action='store',
											type=str, metavar='outputDir',
											default='./TestOutput',
											help='Output directory for output to be written to')
	
	# Parse the arguments
	args = parser.parse_args()
	
	print('Configuring...')
	print('inputFile: \'{0}\''.format(args.inputFile))
	print('configFile: \'{0}\''.format(args.configFile))
	print('ouputDir: \'{0}\"'.format(args.outputDir))
	print('----------------------------------------------------------------')
	
	# If invalid inputFile is given, exit
	if not os.path.exists(args.inputFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
		sys.exit(0)
	
	# If invalid configFile is given, exit
	if not os.path.exists(args.configFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
		sys.exit(0)

	analysis = ProcessData_ENC(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
	analysis.process_data()
