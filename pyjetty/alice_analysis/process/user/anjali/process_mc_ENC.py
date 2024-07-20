#!/usr/bin/env python3

"""
	Analysis class to read a ROOT TTree of MC track information
	and do jet-finding, and save response histograms.
	
	Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import os
import sys
import argparse

# Data analysis and plotting
import numpy as np
import ROOT
import yaml
import array
import math
# from array import *

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib
import fjtools
import ecorrel

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.user.substructure import process_mc_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils.csubtractor import CEventSubtractor

def linbins(xmin, xmax, nbins):
	lspace = np.linspace(xmin, xmax, nbins+1)
	arr = array.array('f', lspace)
	return arr

def logbins(xmin, xmax, nbins):
	lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
	arr = array.array('f', lspace)
	return arr

################################################################
class EEC_pair:
	def __init__(self, _index1, _index2, _weight, _r, _pt):
		self.index1 = _index1
		self.index2 = _index2
		self.weight = _weight
		self.r = _r
		self.pt = _pt

	def is_equal(self, pair2):
		return (self.index1 == pair2.index1 and self.index2 == pair2.index2) \
			or (self.index1 == pair2.index2 and self.index2 == pair2.index1)
	
	def __str__(self):
		return "EEC pair with (index1, index2, weight, RL, pt) = (" + \
			str(self.index1) + ", " + str(self.index2) + ", " + str(self.weight) + \
			", " + str(self.r) + ", " + str(self.pt) + ")"

################################################################
class ProcessMC_ENC(process_mc_base.ProcessMCBase):

	#---------------------------------------------------------------
	# Constructor
	#---------------------------------------------------------------
	def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
	
		# Initialize base class
		super(ProcessMC_ENC, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
		
		self.observable = self.observable_list[0]

		self.pair_eff_file = ROOT.TFile.Open("/global/cfs/cdirs/alice/wenqing/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/wenqing/PairEff.root","READ")
		# self.dpbin = 5
		# self.dp_lo = [0, 0.1, 0.2, 0.4, 1]
		# self.dp_hi = [0.1, 0.2, 0.4, 1, 2]
		self.dpbin = 10
		self.dp_lo = [0, 0.02, 0.06, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6, 1]
		self.dp_hi = [0.02, 0.06, 0.1, 0.14, 0.2, 0.3, 0.4, 0.6, 1, 2]
		self.h1d_eff_vs_dR_in_dq_over_p = []
		for idp in range(self.dpbin):
				hname = 'h1d_eff_vs_dR_in_dq_over_p_{}'.format(idp)
				self.h1d_eff_vs_dR_in_dq_over_p.append( ROOT.TH1D(self.pair_eff_file.Get(hname)) )

	#---------------------------------------------------------------
	# Determine pair efficiency with the pair
	# property and input histograms
	#---------------------------------------------------------------
	def get_pair_eff(self, dist, dq_over_p):
		# return pair efficiency (from 0 to 1)
		idpbin = -9999
		for idp in range(self.dpbin):
				if math.fabs(dq_over_p)>=self.dp_lo[idp] and math.fabs(dq_over_p)<self.dp_hi[idp]:
						idpbin = idp
		
		pair_eff = 1 # set pair efficeincy to 1 if dq_over_p>=2
		if idpbin>=0:
			if math.log10(dist)<0 and math.log10(dist)>-3:
					ibin = self.h1d_eff_vs_dR_in_dq_over_p[idpbin].FindBin(math.log10(dist))
					pair_eff = self.h1d_eff_vs_dR_in_dq_over_p[idpbin].GetBinContent(ibin)
			elif math.log10(dist)>=0:
					pair_eff = 1 # overflow
			else:
					pair_eff = 0 # NB: underflow set to 0 efficiency. Maybe too aggressive but should be fine since we plan to measure down to dist ~1E-2
		
		return pair_eff

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
	# Initialize histograms
	#---------------------------------------------------------------
	def initialize_user_output_objects_R(self, jetR):

		# delta-eta-delta-phi distribution for truth and matched reco particles
		name = 'h2d_matched_part_deta_pt'
		eta_bins = linbins(-0.05, 0.05, 100)
		pt_bins = linbins(0, 10, 100)
		h = ROOT.TH2D(name, name, 100, eta_bins, 100, pt_bins)
		h.GetXaxis().SetTitle('#Delta#eta')
		h.GetYaxis().SetTitle('p_{T}')
		setattr(self, name, h)

		name = 'h2d_matched_part_dphi_pt'
		phi_bins = linbins(-0.02, 0.02, 100)
		pt_bins = linbins(0, 10, 100)
		h = ROOT.TH2D(name, name, 100, phi_bins, 100, pt_bins)
		h.GetXaxis().SetTitle('#Delta#phi')
		h.GetYaxis().SetTitle('p_{T}')
		setattr(self, name, h)

		# matching efficiency of reco tracks wrt truth tracks, as a function of pt
		name = 'h1d_truth_part_pt'
		pt_bins = linbins(0,10,200)
		h = ROOT.TH1D(name, name, 200, pt_bins)
		h.GetXaxis().SetTitle('p_{T}')
		setattr(self, name, h)
		
		name = 'h1d_det_part_pt'
		pt_bins = linbins(0,10,200)
		h = ROOT.TH1D(name, name, 200, pt_bins)
		h.GetXaxis().SetTitle('p_{T}')
		setattr(self, name, h)

		name = 'h1d_truth_matched_part_pt'
		pt_bins = linbins(0,10,200)
		h = ROOT.TH1D(name, name, 200, pt_bins)
		h.GetXaxis().SetTitle('p_{T}')
		setattr(self, name, h)
		
		name = 'h1d_det_matched_part_pt'
		pt_bins = linbins(0,10,200)
		h = ROOT.TH1D(name, name, 200, pt_bins)
		h.GetXaxis().SetTitle('p_{T}')
		setattr(self, name, h)
		
		for observable in self.observable_list:

			for trk_thrd in self.obs_settings[observable]:

				obs_label = self.utils.obs_label(trk_thrd, None) 

				self.pair_type_labels = ['']
				if self.do_median_subtraction and not self.do_feedin_check:
					self.pair_type_labels = ['_bb','_sb','_ss']

				# Init ENC histograms (both det and truth level)
				for pair_type_label in self.pair_type_labels:
					if 'ENC' in observable:
						for ipoint in range(2, 3):
							name = 'h_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							name = 'h_{}{}{}Pt_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) # pt scaled histograms (currently only for unmatched jets)
							pt_bins = linbins(0,200,200)
							ptRL_bins = logbins(1E-3,1E2,60)
							h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
							setattr(self, name, h)

							# Truth histograms
							name = 'h_{}{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							name = 'h_{}{}{}Pt_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label) # pt scaled histograms (currently only for unmatched jets)
							pt_bins = linbins(0,200,200)
							ptRL_bins = logbins(1E-3,1E2,60)
							h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
							setattr(self, name, h)

							# Matched det histograms
							name = 'h_matched_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							# Matched det histograms (with matched truth jet pT filled to the other axis)
							name = 'h_matched_extra_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							name = 'h_matched_extramig_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							name = 'h_matched_extrawt_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							name = 'h_matched_extratru_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							# Matched truth histograms
							name = 'h_matched_{}{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)
							
							if self.thermal_model:
								name = 'h_{}{}{}_JetPt_R{}_{}'.format(observable, ipoint, pair_type_label, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)

					if 'EEC_noweight' in observable or 'EEC_weight2' in observable:
						name = 'h_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						# Truth histograms
						name = 'h_{}{}_JetPt_Truth_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						# Matched det histograms
						name = 'h_matched_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						# Matched det histograms (with matched truth jet pT filled to the other axis)
						name = 'h_matched_extra_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						name = 'h_matched_extramig_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						name = 'h_matched_extrawt_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						name = 'h_matched_extratru_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}^{truth}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)

						# Matched truth histograms
						name = 'h_matched_{}{}_JetPt_Truth_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
						pt_bins = linbins(0,200,200)
						RL_bins = logbins(1E-4,1,50)
						h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
						h.GetXaxis().SetTitle('p_{T,ch jet}')
						h.GetYaxis().SetTitle('R_{L}')
						setattr(self, name, h)
						
						if self.thermal_model:
							name = 'h_{}{}_JetPt_R{}_{}'.format(observable, pair_type_label, jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

				self.jetptslices = [(20,40),(40,60),(60,80)]
				self.feedintypes = ['FF','FIa','FIb','FOa','FOb']
				if self.do_feedin_check:
					if 'ENC' in observable:
						for ipoint in range(2, 3):
							for ptlo, pthi in self.jetptslices:
								for fi in self.feedintypes:
									name = 'h_matched{}_{}{}_JetPt{}{}_R{}_{}'.format(fi, observable, ipoint, ptlo, pthi, jetR, obs_label)
									pt_bins = linbins(0,200,200)
									RL_bins = logbins(1E-4,1,50)
									h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
									h.GetXaxis().SetTitle('p_{T,ch jet}')
									h.GetYaxis().SetTitle('R_{L}')
									setattr(self, name, h)

									name = 'h_matched{}_{}{}_JetPt{}{}_Truth_R{}_{}'.format(fi, observable, ipoint, ptlo, pthi, jetR, obs_label)
									h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
									h.GetXaxis().SetTitle('p_{T,ch jet}')
									h.GetYaxis().SetTitle('R_{L}')
									setattr(self, name, h)

									name = 'h_matched{}_{}_JetPt{}{}_Truth_vs_Det_R{}_{}'.format(fi, observable, ptlo, pthi, jetR, obs_label)
									pt_bins = linbins(0,200,200)
									h = ROOT.TH2D(name, name, 200, pt_bins, 200, pt_bins)
									h.GetXaxis().SetTitle('p_{T,ch jet}^{det}')
									h.GetYaxis().SetTitle('p_{T,ch jet}^{truth}')
									setattr(self, name, h)
				
				
				name = 'hSigJetpt_Partpt_R{}_{}'.format(jetR,obs_label)
				h = ROOT.TH2F(name, name, 100, 0, 100., 200, 0, 50.)
				setattr(self, name, h)

				name = 'hSigJetpt_Partpt_MC_R{}_{}'.format(jetR,obs_label)
				h = ROOT.TH2F(name, name, 100, 0, 100., 200, 0, 50.)
				setattr(self, name, h)

				name = 'hSigJetpt_Partpt_UE_R{}_{}'.format(jetR,obs_label)
				h = ROOT.TH2F(name, name, 100, 0, 100., 200, 0, 50.)
				setattr(self, name, h)

				for ptlo, pthi in [(20,40),(40,60),(60,80)]:
					name = 'hJetConstituents_{}{}_R{}_{}'.format(ptlo,pthi,jetR,obs_label)
					h = ROOT.TH1F(name, name, 50, 0, 50)
					setattr(self, name, h)
					name = 'hJetConstituents_MC_{}{}_R{}_{}'.format(ptlo,pthi,jetR,obs_label)
					h = ROOT.TH1F(name, name, 50, 0, 50)
					setattr(self, name, h)
					name = 'hJetConstituents_UE_{}{}_R{}_{}'.format(ptlo,pthi,jetR,obs_label)
					h = ROOT.TH1F(name, name, 50, 0, 50)
					setattr(self, name, h)
				
				if 'jet_pt' in observable:
					name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					h = ROOT.TH1D(name, name, 200, pt_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}')
					h.GetYaxis().SetTitle('Counts')
					setattr(self, name, h)

					name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					h = ROOT.TH1D(name, name, 200, pt_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}')
					h.GetYaxis().SetTitle('Counts')
					setattr(self, name, h)

					# Matched det histograms
					name = 'h_matched_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					h = ROOT.TH1D(name, name, 200, pt_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}')
					h.GetYaxis().SetTitle('Counts')
					setattr(self, name, h)

					# Matched truth histograms
					name = 'h_matched_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					h = ROOT.TH1D(name, name, 200, pt_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}')
					h.GetYaxis().SetTitle('Counts')
					setattr(self, name, h)

					# Correlation between matched det and truth
					name = 'h_matched_{}_JetPt_Truth_vs_Det_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					h = ROOT.TH2D(name, name, 200, pt_bins, 200, pt_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}^{det}')
					h.GetYaxis().SetTitle('p_{T,ch jet}^{truth}')
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

							name = 'h_perpcone{}_{}_JetPt_Truth_R{}_{}'.format(perpcone_R, 'pt', jetR, trk_thrd)
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

							name = 'h_perpcone{}_Nconst_JetPt_Truth_R{}_{}'.format(perpcone_R, jetR, trk_thrd)
							pt_bins = linbins(-200,200,400)
							Nconst_bins = linbins(0,50,50)
							h = ROOT.TH2D(name, name, 200, pt_bins, 50, Nconst_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('N_{const}')
							setattr(self, name, h)

						for pair_type_label in self.pair_type_labels:

							if 'ENC' in observable:
								for ipoint in range(2, 3):
									name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label, jetR, trk_thrd)
									pt_bins = linbins(0,200,200)
									RL_bins = logbins(1E-4,1,50)
									h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
									h.GetXaxis().SetTitle('p_{T,ch jet}')
									h.GetYaxis().SetTitle('R_{L}')
									setattr(self, name, h)

									name = 'h_perpcone{}_{}Pt_JetPt_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label, jetR, trk_thrd)
									pt_bins = linbins(0,200,200)
									ptRL_bins = logbins(1E-3,1E2,60)
									h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
									h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
									setattr(self, name, h)

									name = 'h_perpcone{}_{}_JetPt_Truth_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label, jetR, trk_thrd)
									pt_bins = linbins(0,200,200)
									RL_bins = logbins(1E-4,1,50)
									h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
									h.GetXaxis().SetTitle('p_{T,ch jet}')
									h.GetYaxis().SetTitle('R_{L}')
									setattr(self, name, h)

									name = 'h_perpcone{}_{}Pt_JetPt_Truth_R{}_{}'.format(perpcone_R, observable + str(ipoint) + pair_type_label, jetR, trk_thrd)
									pt_bins = linbins(0,200,200)
									ptRL_bins = logbins(1E-3,1E2,60)
									h = ROOT.TH2D(name, name, 200, pt_bins, 60, ptRL_bins)
									h.GetYaxis().SetTitle('p_{T,ch jet}R_{L}') # NB: y axis scaled by jet pt (applied jet by jet)
									setattr(self, name, h)

							if 'EEC_noweight' in observable or 'EEC_weight2' in observable:
								name = 'h_perpcone{}_{}_JetPt_R{}_{}'.format(perpcone_R, observable + pair_type_label, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)

								name = 'h_perpcone{}_{}_JetPt_Truth_R{}_{}'.format(perpcone_R, observable + pair_type_label, jetR, obs_label)
								pt_bins = linbins(0,200,200)
								RL_bins = logbins(1E-4,1,50)
								h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
								h.GetXaxis().SetTitle('p_{T,ch jet}')
								h.GetYaxis().SetTitle('R_{L}')
								setattr(self, name, h)

				# # Diagnostic
				# if 'jet_diag' in observable:
				#   name = 'h_{}_JetEta_R{}_{}'.format(observable, jetR, obs_label)
				#   pt_bins = linbins(0,200,200)
				#   eta_bins = linbins(-10,10,200)
				#   h = ROOT.TH2D(name, name, 200, pt_bins, 200, eta_bins)
				#   h.GetXaxis().SetTitle('p_{T,ch jet}')
				#   h.GetYaxis().SetTitle('#eta_{ch jet}')
				#   setattr(self, name, h)

				#   name = 'h_{}_JetEta_Truth_R{}_{}'.format(observable, jetR, obs_label)
				#   pt_bins = linbins(0,200,200)
				#   eta_bins = linbins(-10,10,200)
				#   h = ROOT.TH2D(name, name, 200, pt_bins, 200, eta_bins)
				#   h.GetXaxis().SetTitle('p_{T,ch jet}')
				#   h.GetYaxis().SetTitle('#eta_{ch jet}')
				#   setattr(self, name, h)

				# Init pair distance histograms (both det and truth level)
				# average track pt bins
				self.trk_pt_lo = [0, 1, 2, 3, 5, 7, 10]
				self.trk_pt_hi = [1, 2, 3, 5, 7, 10, 100]
				# track pt asymmetry bins: (pt_trk1-pt_trk2)/(pt_trk1+pt_trk2)
				self.trk_alpha_lo = [0, 0.2, 0.4, 0.6, 0.8]
				self.trk_alpha_hi = [0.2, 0.4, 0.6, 0.8, 1]
				if 'EEC_detail' in observable:
					# inclusive
					name = 'h_{}_JetPt_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					RL_bins = logbins(1E-4,1,50)
					h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}')
					h.GetYaxis().SetTitle('R_{L}')
					setattr(self, name, h)

					name = 'h_{}_JetPt_Truth_R{}_{}'.format(observable, jetR, obs_label)
					pt_bins = linbins(0,200,200)
					RL_bins = logbins(1E-4,1,50)
					h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
					h.GetXaxis().SetTitle('p_{T,ch jet}')
					h.GetYaxis().SetTitle('R_{L}')
					setattr(self, name, h)

					# fine bins
					for ipt in range( len(self.trk_pt_lo) ):
						for ialpha in range( len(self.trk_alpha_lo) ):
							name = 'h_{}{}{}_{:.1f}{:.1f}_JetPt_R{}_{}'.format(observable, self.trk_pt_lo[ipt], self.trk_pt_hi[ipt], self.trk_alpha_lo[ialpha], self.trk_alpha_hi[ialpha], jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

							name = 'h_{}{}{}_{:.1f}{:.1f}_JetPt_Truth_R{}_{}'.format(observable, self.trk_pt_lo[ipt], self.trk_pt_hi[ipt], self.trk_alpha_lo[ialpha], self.trk_alpha_hi[ialpha], jetR, obs_label)
							pt_bins = linbins(0,200,200)
							RL_bins = logbins(1E-4,1,50)
							h = ROOT.TH2D(name, name, 50, pt_bins, 50, RL_bins)
							h.GetXaxis().SetTitle('p_{T,ch jet}')
							h.GetYaxis().SetTitle('R_{L}')
							setattr(self, name, h)

				# Residuals and responses (currently not filled or used)
				for trk_thrd in self.obs_settings[observable]:
				
					for ipoint in range(2, 3):
						if not self.is_pp and not self.is_pA:
							for R_max in self.max_distance:
								self.create_response_histograms(observable, ipoint, jetR, trk_thrd, R_max)          
						else:
							self.create_response_histograms(observable, ipoint, jetR, trk_thrd)
					

	#---------------------------------------------------------------
	# This function is called once for each jet subconfiguration
	# Fill 2D histogram of (pt, obs)
	#---------------------------------------------------------------
	def create_response_histograms(self, observable, ipoint, jetR, trk_thrd, R_max = None):
	
		if R_max:
			suffix = '_Rmax{}'.format(R_max)
		else:
			suffix = ''

		# Create THn of response for ENC
		dim = 4;
		title = ['p_{T,det}', 'p_{T,truth}', 'R_{L,det}', 'R_{L,truth}']
		nbins = [30, 20, 100, 100]
		min = [0., 0., 0., 0.]
		max = [150., 200., 1., 1.]
		name = 'hResponse_JetPt_{}{}_R{}_{}{}'.format(observable, ipoint, jetR, trk_thrd, suffix)
		self.create_thn(name, title, dim, nbins, min, max)
		
		name = 'hResidual_JetPt_{}{}_R{}_{}{}'.format(observable, ipoint, jetR, trk_thrd, suffix)
		h = ROOT.TH3F(name, name, 20, 0, 200, 100, 0., 1., 200, -2., 2.)
		h.GetXaxis().SetTitle('p_{T,truth}')
		h.GetYaxis().SetTitle('R_{L}')
		h.GetZaxis().SetTitle('#frac{R_{L,det}-R_{L,truth}}{R_{L,truth}}')
		setattr(self, name, h)

	def get_pair_eff_weights(self, corr_builder, ipoint, constituents):
		# NB: currently applying the pair eff weight to both 2 point correlator and higher point correlators. Need to check if the same pair efficiency effect still work well for higher point correlators
		weights_pair = []
		for index in range(corr_builder.correlator(ipoint).rs().size()):
			part1 = corr_builder.correlator(ipoint).indices1()[index]
			part2 = corr_builder.correlator(ipoint).indices2()[index]
			if part1!=part2: # FIX ME: not sure, but for now only apply pair efficiency for non auto-correlations
				# Need to find the associated truth information for each pair (charge and momentum)
				part1_truth = constituents[part1].python_info().particle_truth
				part2_truth = constituents[part2].python_info().particle_truth
				q1 = constituents[part1].python_info().charge
				q2 = constituents[part2].python_info().charge
				dist = corr_builder.correlator(ipoint).rs()[index] # NB: use reconstructed distance since it's faster and should be equivalent to true distance because there is no angular smearing on the track momentum. To switch back to the true distance, use: self.calculate_distance(part1_truth, part2_truth)
				dq_over_p = q1/part1_truth.pt()-q2/part2_truth.pt()
				# calculate pair efficeincy and apply it as an additional weight
				weights_pair.append( self.get_pair_eff(dist, dq_over_p) )
			else:
				weights_pair.append( 1 )
		return weights_pair

	def is_same_charge(self, corr_builder, ipoint, constituents, index):
		part1 = corr_builder.correlator(ipoint).indices1()[index]
		part2 = corr_builder.correlator(ipoint).indices2()[index]
		q1 = constituents[part1].python_info().charge
		q2 = constituents[part2].python_info().charge

		if q1*q2 > 0:
			return True
		else:
			return False

	def check_pair_type(self, corr_builder, ipoint, constituents, index):
		# print("TESTING:",type(corr_builder.correlator(ipoint)))
		# print("TESTING:",type(corr_builder.correlator(ipoint).indices1()))
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
	# Fill 2D histogram of (pt, obs)
	#---------------------------------------------------------------
	def fill_observable_histograms(self, hname, jet, jet_groomed_lund, jetR, obs_setting,
																 grooming_setting, obs_label, jet_pt_ungroomed):
		
		constituents = fj.sorted_by_pt(jet.constituents())
		c_select = fj.vectorPJ()
		trk_thrd = obs_setting

		# NB: use jet_pt_ungroomed instead of jet.perp() for PbPb, which include the UE subtraction
		if self.do_median_subtraction:
			jet_pt = jet_pt_ungroomed
		else:
			jet_pt = jet.perp()
		
		
		for c in constituents:
			if c.pt() < trk_thrd:
				break
			c_select.append(c) # NB: use the break statement since constituents are already sorted
		
		if self.ENC_pair_cut and (not 'Truth' in hname):
			dphi_cut = -9999 # means no dphi cut
			deta_cut = 0.008
		else:
			dphi_cut = -9999
			deta_cut = -9999

		new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
		for observable in self.observable_list:
			if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
				for ipoint in range(2, 3):
					if self.ENC_fastsim and (not 'Truth' in hname): # NB: only apply pair efficiency effect for fast sim and det level distributions
						weights_pair = self.get_pair_eff_weights(new_corr, ipoint, c_select)

					for index in range(new_corr.correlator(ipoint).rs().size()):

						# processing only like-sign pairs when self.ENC_pair_like is on
						if self.ENC_pair_like and (not self.is_same_charge(new_corr, ipoint, c_select, index)):
							continue

						# processing only unlike-sign pairs when self.ENC_pair_unlike is on
						if self.ENC_pair_unlike and self.is_same_charge(new_corr, ipoint, c_select, index):
							continue
						
						# separate out sig-sig, sig-bkg, bkg-bkg correlations for EEC pairs
						pair_type_label = ''
						if self.do_median_subtraction and not self.do_feedin_check:
							pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
							pair_type_label = self.pair_type_labels[pair_type]

						if 'ENC' in observable:
							if self.ENC_fastsim and (not 'Truth' in hname):
								getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index])
								getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index]) # NB: fill pt*RL

							else:
								getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
								getattr(self, hname.format(observable + str(ipoint) + pair_type_label + 'Pt',obs_label)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])

						if ipoint==2 and 'EEC_noweight' in observable:
							if self.ENC_fastsim and (not 'Truth' in hname):
								getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], weights_pair[index])
							else:
								getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

						if ipoint==2 and 'EEC_weight2' in observable:
							if self.ENC_fastsim and (not 'Truth' in hname):
								getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index]*weights_pair[index],2))
							else:
								getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

			if 'jet_pt' in observable:
				getattr(self, hname.format(observable,obs_label)).Fill(jet_pt)
			
			# NB: for now, only perform this check on data and full sim
			if 'EEC_detail' in observable and self.ENC_fastsim==False: 
				ipoint = 2 # EEC is 2 point correlator
				for index in range(new_corr.correlator(ipoint).rs().size()):
					part1 = new_corr.correlator(ipoint).indices1()[index]
					part2 = new_corr.correlator(ipoint).indices2()[index]
					pt1 = c_select[part1].perp()
					pt2 = c_select[part2].perp()
					pt_avg = (pt1+pt2)/2
					alpha = math.fabs(pt1-pt2)
					
					for _, (pt_lo, pt_hi) in enumerate(zip(self.trk_pt_lo,self.trk_pt_hi)):
						for _, (alpha_lo, alpha_hi) in enumerate(zip(self.trk_alpha_lo,self.trk_alpha_hi)):
							if pt_avg >= pt_lo and pt_avg < pt_hi and alpha >= alpha_lo and alpha < alpha_hi:
								if 'noweight' in observable:
									getattr(self, hname.format(observable + str(pt_lo) + str(pt_hi) + '_' + '{:.1f}'.format(alpha_lo) + '{:.1f}'.format(alpha_hi),obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])
								else:
									getattr(self, hname.format(observable + str(pt_lo) + str(pt_hi) + '_' + '{:.1f}'.format(alpha_lo) + '{:.1f}'.format(alpha_hi),obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
								break

					# fill inclusively
					if 'noweight' in observable:
						getattr(self, hname.format(observable,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])
					else:
						getattr(self, hname.format(observable,obs_label)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
				
	#---------------------------------------------------------------
	# This function is called per observable per jet subconfigration 
	# used in fill_matched_jet_histograms
	# This function is created because we cannot use fill_observable_histograms 
	# directly because observable list loop inside that function
	#---------------------------------------------------------------
	def fill_matched_observable_histograms(self, hname, observable, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed, jet_pt_for_weight=-999, fillhist=False):
		
		constituents = fj.sorted_by_pt(jet.constituents())
		c_select = fj.vectorPJ()
		trk_thrd = obs_setting

		if self.do_median_subtraction:
			jet_pt = jet_pt_ungroomed
		else:
			jet_pt = jet.perp()

		if jet_pt_for_weight > 0:
			jet_pt = jet_pt_for_weight

		nconstituents = 0
		for c in constituents:
			if c.pt() < trk_thrd:
				break
			c_select.append(c) # NB: use the break statement since constituents are already sorted
			nconstituents += 1
			if fillhist:
				getattr(self, 'hSigJetpt_Partpt_R{}_{}'.format(jetR,obs_label)).Fill(jet_pt, c.pt())
		
		if fillhist:
			c_select_MC = [c for c in c_select if c.user_index()>=0]
			nconstituents_MC = len(c_select_MC)
			for c in c_select_MC:
				getattr(self, 'hSigJetpt_Partpt_MC_R{}_{}'.format(jetR,obs_label)).Fill(jet_pt, c.pt())
			
			c_select_UE = [c for c in c_select if c.user_index()<0]
			nconstituents_UE = len(c_select_UE)
			for c in c_select_UE:
				getattr(self, 'hSigJetpt_Partpt_UE_R{}_{}'.format(jetR,obs_label)).Fill(jet_pt, c.pt())
			
			for ptlo, pthi in [(20,40),(40,60),(60,80)]:
				if jet_pt < pthi and jet_pt > ptlo:
					getattr(self, 'hJetConstituents_{}{}_R{}_{}'.format(ptlo,pthi,jetR,obs_label)).Fill(nconstituents)
					getattr(self, 'hJetConstituents_MC_{}{}_R{}_{}'.format(ptlo,pthi,jetR,obs_label)).Fill(nconstituents_MC)
					getattr(self, 'hJetConstituents_UE_{}{}_R{}_{}'.format(ptlo,pthi,jetR,obs_label)).Fill(nconstituents_UE)
					break
		
			if 60 < jet_pt < 80 and nconstituents_UE>0 and trk_thrd==1 and 'ENC' in observable and 'extra' not in hname:
				print("NEW JET")
				print(jet_pt, nconstituents, nconstituents_MC, nconstituents_UE)
				print([(c.pt(),c.user_index()) for c in c_select if c.user_index()<0])
		
		if self.ENC_pair_cut and (not 'Truth' in hname):
			dphi_cut = -9999 # means no dphi cut
			deta_cut = 0.008
		else:
			dphi_cut = -9999
			deta_cut = -9999

		new_corr = ecorrel.CorrelatorBuilder(c_select, jet_pt, 2, 1, dphi_cut, deta_cut)
		if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
			for ipoint in range(2, 3):
				if self.ENC_fastsim and (not 'Truth' in hname): # NB: only apply pair efficiency effect for fast sim and det level distributions
					weights_pair = self.get_pair_eff_weights(new_corr, ipoint, c_select)

				for index in range(new_corr.correlator(ipoint).rs().size()):
					pair_type_label = ''
					if self.do_median_subtraction and not self.do_feedin_check:
						pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
						pair_type_label = self.pair_type_labels[pair_type]
					
					if 'ENC' in observable:
						if self.ENC_fastsim and (not 'Truth' in hname):
							getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt_ungroomed, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]*weights_pair[index])
						else:
							getattr(self, hname.format(observable + str(ipoint) + pair_type_label,obs_label)).Fill(jet_pt_ungroomed, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: use jet_pt_ungroomed instead of jet_pt so if jet_pt_ungroomed is different from jet_pt, it will be used. This is mainly for matched jets study

					if ipoint==2 and 'EEC_noweight' in observable:
						if self.ENC_fastsim and (not 'Truth' in hname):
							getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_ungroomed, new_corr.correlator(ipoint).rs()[index], weights_pair[index])
						else:
							getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_ungroomed, new_corr.correlator(ipoint).rs()[index])

					if ipoint==2 and 'EEC_weight2' in observable:
						if self.ENC_fastsim and (not 'Truth' in hname):
							getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_ungroomed, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index]*weights_pair[index],2))
						else:
							getattr(self, hname.format(observable + pair_type_label,obs_label)).Fill(jet_pt_ungroomed, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))

		if 'jet_pt' in observable:
			getattr(self, hname.format(observable,obs_label)).Fill(jet_pt)

	def analyze_matched_pairs(self, fj_particles_det, fj_particles_truth, jetR=0.4):
		
		#handle case with no truth particles
		if type(fj_particles_truth) is float:
				print("EVENT WITH NO TRUTH PARTICLES!!!")
				return
		
		#handle case with no det particles
		if type(fj_particles_det) is float:
				print("EVENT WITH NO DET PARTICLES!!!")
				fj_particles_det = []
				
		for det_part in fj_particles_det:
				hname = 'h1d_det_part_pt'
				getattr(self, hname).Fill(det_part.perp())

		############################# TRACK MATCHING ################################
		# set all indices to dummy index
		dummy_index = -1
		for i in range(len(fj_particles_truth)):
			fj_particles_truth[i].set_user_index(dummy_index)
		for i in range(len(fj_particles_det)):
			fj_particles_det[i].set_user_index(dummy_index)
		
		fj_particles_det_sel = [t for t in fj_particles_det if np.abs(t.eta()) < 0.9 ]
		fj_particles_truth_sel = [t for t in fj_particles_truth if np.abs(t.eta()) < 0.9 ]
		
		# perform matching, give matches the same user_index
		index = 0
		det_used = []
		# note: CANNOT loop like this: <for truth_part in fj_particles_truth:>
		for itruth in range(len(fj_particles_truth_sel)):
			truth_part = fj_particles_truth_sel[itruth]
		
			truth_part.set_user_index(index)
				
			candidates = []
			candidates_R = []
				
			for idet in range(len(fj_particles_det_sel)):
				det_part = fj_particles_det_sel[idet]
				
				delta_R = self.calculate_distance(truth_part, det_part)
				if delta_R < 0.05 and abs((det_part.perp() - truth_part.perp()) / truth_part.perp()) < 0.1 \
								and det_part not in det_used:
						candidates.append(det_part)
						candidates_R.append(delta_R)

			# if match found
			if len(candidates) > 0:
				det_match = candidates[np.argmin(candidates_R)]
				det_match.set_user_index(index)
				det_used.append(det_match)

				deta = det_match.eta() - truth_part.eta()
				dphi = det_match.phi() - truth_part.phi()

				hname = 'h2d_matched_part_deta_pt'
				getattr(self, hname).Fill(deta, truth_part.perp())

				hname = 'h2d_matched_part_dphi_pt'
				getattr(self, hname).Fill(dphi, truth_part.perp())

				hname = 'h1d_truth_matched_part_pt'
				getattr(self, hname).Fill(truth_part.perp())
				
				hname = 'h1d_det_matched_part_pt'
				getattr(self, hname).Fill(det_match.perp())


			hname = 'h1d_truth_part_pt'
			getattr(self, hname).Fill(truth_part.perp())

			index += 1

		# handle unmatched particles, give them all different user_index s
		for i in range(len(fj_particles_truth_sel)):
			part = fj_particles_truth_sel[i]
			if part.user_index() == dummy_index:
				print("DUMMY INDEX TRUTH PARTICLE")
				part.set_user_index(index)
				index += 1
		for i in range(len(fj_particles_det_sel)):
			part = fj_particles_det_sel[i]
			if part.user_index() == dummy_index:
				part.set_user_index(index)
				index += 1

		# print(len(det_used)/len(fj_particles_truth_sel))

		############################# JET RECO ################################
	#   # Set jet definition and a jet selector
	#   jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
	#   jet_selector_det = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)

	#   cs_truth = fj.ClusterSequence(fj_particles_truth, jet_def)
	#   truth_jets = fj.sorted_by_pt(jet_selector_det(cs_truth.inclusive_jets()))
		
	#   cs_det = fj.ClusterSequence(fj_particles_det, jet_def)
	#   det_jets = fj.sorted_by_pt(jet_selector_det(cs_det.inclusive_jets()))

	# # truth level EEC pairs
	#   truth_pairs = []
	#   for jet in truth_jets:
	#     truth_pairs += self.get_EEC_pairs(jet, ipoint=2)

	#   # det level EEC pairs
	#   det_pairs = []
	#   for jet in det_jets:
	#     det_pairs += self.get_EEC_pairs(jet, ipoint=2)


		########################## TTree output generation #########################
		# # composite of truth and smeared pairs, fill the TTree preprocessed
		# dummyval = -9999

		# # pair mathcing
		# for t_pair in truth_pairs:

		#     gen_energy_weight = t_pair.weight
		#     gen_R_L = t_pair.r
		#     gen_jet_pt = t_pair.pt
		#     obs_thrown = 0

		#     match_found = False
		#     for d_pair in det_pairs:
		#         if d_pair.is_equal(t_pair):
		#             obs_energy_weight = d_pair.weight
		#             obs_R_L = d_pair.r
		#             obs_jet_pt = d_pair.pt
		#             match_found = True
		#             break
		#     if not match_found:
		#         obs_energy_weight = dummyval
		#         obs_R_L = dummyval
		#         obs_jet_pt = dummyval
		#         obs_thrown = 1
				
		#     # find TTree TODO
		#     name = 'preprocessed_np_mc'
		#     getattr(self, name).append([gen_energy_weight, gen_R_L, gen_jet_pt, obs_energy_weight, obs_R_L, obs_jet_pt, self.pt_hat, self.event_number])
				
				
		"""
		line = ""
		for part in fj_particles_truth:
				line += str(part.user_index()) + " "
		print(line)
		
		line = ""
		for part in fj_particles_det:
				line += str(part.user_index()) + " "
		print(line)
		"""
				


	def get_EEC_pairs(self, jet, ipoint=2):
		pairs = []

		jet_pt = jet.perp()

		#push constutents to a vector in python
		_v = fj.vectorPJ()
		_ = [_v.push_back(c) for c in jet.constituents()]

		# n-point correlator with all charged particles
		max_npoint = 2
		weight_power = 1
		dphi_cut = -9999
		deta_cut = -9999
		cb = ecorrel.CorrelatorBuilder(_v, jet_pt, max_npoint, weight_power, dphi_cut, deta_cut)

		EEC_cb = cb.correlator(ipoint)

		EEC_weights = EEC_cb.weights() # cb.correlator(npoint).weights() constains list of weights
		EEC_rs = EEC_cb.rs() # cb.correlator(npoint).rs() contains list of RL
		EEC_indicies1 = EEC_cb.indices1() # contains list of 1st track in the pair (index should be based on the indices in c_select)
		EEC_indicies2 = EEC_cb.indices2() # contains list of 2nd track in the pair

		for i in range(len(EEC_rs)):
			event_index1 = _v[EEC_indicies1[i]].user_index()
			event_index2 = _v[EEC_indicies2[i]].user_index()
			pairs.append(EEC_pair(event_index1, event_index2, EEC_weights[i], EEC_rs[i], jet_pt))

		return pairs
	
	#---------------------------------------------------------------
	# This function is called per jet subconfigration 
	# Fill matched jet histograms
	#---------------------------------------------------------------
	def fill_matched_jet_histograms(self, jet_det, jet_det_groomed_lund, jet_truth,
																	jet_truth_groomed_lund, jet_pp_det, jetR,
																	obs_setting, grooming_setting, obs_label,
																	jet_pt_det_ungroomed, jet_pt_truth_ungroomed, R_max, suffix, **kwargs):
		# If jetscape, we will need to correct substructure observable for holes (pt is corrected in base class)
		# For ENC in PbPb, jet_pt_det_ungroomed stores the corrected jet pT
		if self.jetscape:
				holes_in_det_jet = kwargs['holes_in_det_jet']
				holes_in_truth_jet = kwargs['holes_in_truth_jet']

		cone_parts_in_det_jet = kwargs['cone_parts_in_det_jet']
		if cone_parts_in_det_jet != None:
			parts_in_perpcone1 = cone_parts_in_det_jet[0]
			parts_in_perpcone2 = cone_parts_in_det_jet[1]
			
			parts_in_perpcone1_truth = cone_parts_in_det_jet[2]
			parts_in_perpcone2_truth = cone_parts_in_det_jet[3]

		# Todo: add additonal weight for jet pT spectrum
		# if self.rewight_pt:
		#   w_pt = 1+pow(jet_truth,0.2)
		# else:
		#   w_pt = 1
		
		if self.do_median_subtraction:
			# print('evt #',self.event_number)
			jet_pt_det = jet_pt_det_ungroomed
			jet_pt_truth = jet_pt_truth_ungroomed
			rho_bg = kwargs['rho']
			if len(rho_bg) == 2:
				rho_bg_truth = rho_bg[1]
				rho_bg = rho_bg[0]

			# print('Det: pT',jet_det.perp(),'(',jet_pt_det,')','phi',jet_det.phi(),'eta',jet_det.eta())
			# print('Truth: pT',jet_truth.perp(),'phi',jet_truth.phi(),'eta',jet_truth.eta())
		else:
			jet_pt_det = jet_det.pt()
			jet_pt_truth = jet_truth.pt()


		if cone_parts_in_det_jet == None:

			for observable in self.observable_list:

				hname = 'h_matched_{{}}_JetPt_R{}_{{}}'.format(jetR)
				self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, fillhist=True)

				hname = 'h_matched_{{}}_JetPt_Truth_R{}_{{}}'.format(jetR)				
				self.fill_matched_observable_histograms(hname, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth)

				# fill RL vs matched truth jet pT for det jets (only fill these extra histograms for ENC or pair distributions)
				if 'ENC' in observable or 'EEC_noweight' in observable or 'EEC_weight2' in observable:
					hname = 'h_matched_extra_{{}}_JetPt_R{}_{{}}'.format(jetR)
					self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth) # NB: use the truth jet pt so the reco jets histograms are comparable to matched truth jets. However this also means that two identical histograms will be filled fot jet_pt observable
					# print("NEW JET")
					# print("DT:", len(jet_det.constituents()), jet_pt_truth)
					hname = 'h_matched_extramig_{{}}_JetPt_R{}_{{}}'.format(jetR)
					self.fill_matched_observable_histograms(hname, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det) # NB: use the truth jet pt so the reco jets histograms are comparable to matched truth jets. However this also means that two identical histograms will be filled fot jet_pt observable
					# print("TD:", len(jet_truth.constituents()), jet_pt_det)
					hname = 'h_matched_extrawt_{{}}_JetPt_R{}_{{}}'.format(jetR)
					self.fill_matched_observable_histograms(hname, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det, jet_pt_for_weight=jet_pt_truth) 
					# print("DD:", len(jet_det.constituents()), jet_pt_det)
					hname = 'h_matched_extratru_{{}}_JetPt_R{}_{{}}'.format(jetR)
					self.fill_matched_observable_histograms(hname, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth, jet_pt_for_weight=jet_pt_truth) 
					# print("TT:", len(jet_truth.constituents()), jet_pt_truth)
					# NB: use the truth jet pt so the reco jets histograms are comparable to matched truth jets. However this also means that two identical histograms will be filled fot jet_pt observable

					if self.do_feedin_check:
						if 'ENC' in observable:
							for ipoint in range(2, 3):
								for ptlo, pthi in self.jetptslices:
									fi = 0
									if ptlo <= jet_pt_det <= pthi and ptlo <= jet_pt_truth <= pthi:
										fi = 'FF'
									elif ptlo <= jet_pt_truth <= pthi:
										if jet_pt_det > pthi:
											fi = 'FIa'
										if jet_pt_det < ptlo:
											fi = 'FIb'
									elif ptlo <= jet_pt_det <= pthi:
										if jet_pt_truth > pthi:
											fi = 'FOa'
										if jet_pt_truth < ptlo:
											fi = 'FOb'
									if fi != 0:
										hname = 'h_matched{}_{{}}_JetPt{}{}_R{}_{{}}'.format(fi, ptlo, pthi, jetR)
										hname_truth = 'h_matched{}_{{}}_JetPt{}{}_Truth_R{}_{{}}'.format(fi, ptlo, pthi, jetR)
										ptname = 'h_matched{}_{}_JetPt{}{}_Truth_vs_Det_R{}_{}'.format(fi, observable, ptlo, pthi, jetR, obs_label)

										self.fill_matched_observable_histograms(hname, observable, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det)
										self.fill_matched_observable_histograms(hname_truth, observable, jet_truth, jet_truth_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth)
										getattr(self, ptname).Fill(jet_pt_det, jet_pt_truth)


				# Fill correlation between matched det and truth jets
				if 'jet_pt' in observable:
					hname = 'h_matched_{}_JetPt_Truth_vs_Det_R{}_{}'.format(observable, jetR, obs_label)
					getattr(self, hname).Fill(jet_pt_det, jet_pt_truth) 

		else:
			if self.do_perpendicular_cone:
				self.fill_perp_cone_histograms(parts_in_perpcone1, jetR, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det_ungroomed, suffix, rho_bg)
				self.fill_perp_cone_histograms(parts_in_perpcone2, jetR, jet_det, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_det_ungroomed, suffix, rho_bg)
				self.fill_perp_cone_histograms(parts_in_perpcone1_truth, jetR, jet_truth, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth_ungroomed, suffix, rho_bg_truth, truth=True)
				self.fill_perp_cone_histograms(parts_in_perpcone2_truth, jetR, jet_truth, jet_det_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_truth_ungroomed, suffix, rho_bg_truth, truth=True)
			 
		# # Find all subjets
		# trk_thrd = obs_setting
		# cs_subjet_det = fj.ClusterSequence(jet_det.constituents(), self.subjet_def[trk_thrd])
		# subjets_det = fj.sorted_by_pt(cs_subjet_det.inclusive_jets())

		# cs_subjet_truth = fj.ClusterSequence(jet_truth.constituents(), self.subjet_def[trk_thrd])
		# subjets_truth = fj.sorted_by_pt(cs_subjet_truth.inclusive_jets())
		
		# if not self.is_pp:
		#   cs_subjet_det_pp = fj.ClusterSequence(jet_pp_det.constituents(), self.subjet_def[trk_thrd])
		#   subjets_det_pp = fj.sorted_by_pt(cs_subjet_det_pp.inclusive_jets())

		# # Loop through subjets and set subjet matching candidates for each subjet in user_info
		# if self.is_pp:
		#     [[self.set_matching_candidates(subjet_det, subjet_truth, subjetR, 'hDeltaR_ppdet_pptrue_ENC_R{}_{}'.format(jetR, subjetR)) for subjet_truth in subjets_truth] for subjet_det in subjets_det]
		# else:
		#     # First fill the combined-to-pp matches, then the pp-to-pp matches
		#     [[self.set_matching_candidates(subjet_det_combined, subjet_det_pp, subjetR, 'hDeltaR_combined_ppdet_ENC_R{}_{}_Rmax{}'.format(jetR, subjetR, R_max), fill_jet1_matches_only=True) for subjet_det_pp in subjets_det_pp] for subjet_det_combined in subjets_det]
		#     [[self.set_matching_candidates(subjet_det_pp, subjet_truth, subjetR, 'hDeltaR_ppdet_pptrue_ENC_R{}_{}_Rmax{}'.format(jetR, subjetR, R_max)) for subjet_truth in subjets_truth] for subjet_det_pp in subjets_det_pp]
			
		# # Loop through subjets and set accepted matches
		# if self.is_pp:
		#     [self.set_matches_pp(subjet_det, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det in subjets_det]
		# else:
		#     [self.set_matches_AA(subjet_det_combined, subjetR, 'hSubjetMatchingQA_R{}_{}'.format(jetR, subjetR)) for subjet_det_combined in subjets_det]

		# # Loop through matches and fill histograms
		# for observable in self.observable_list:
		
		#   # Fill inclusive subjets
		#   if 'inclusive' in observable:

		#     for subjet_det in subjets_det:
				
		#       z_det = subjet_det.pt() / jet_det.pt()
					
		#       # If z=1, it will be default be placed in overflow bin -- prevent this
		#       if np.isclose(z_det, 1.):
		#         z_det = 0.999
					
		#       successful_match = False

		#       if subjet_det.has_user_info():
		#         subjet_truth = subjet_det.python_info().match
					
		#         if subjet_truth:
						
		#           successful_match = True

		#           # For subjet matching radius systematic, check distance between subjets
		#           if self.matching_systematic:
		#             if subjet_det.delta_R(subjet_truth) > 0.5 * self.jet_matching_distance * subjetR:
		#               continue
							
		#           z_truth = subjet_truth.pt() / jet_truth.pt()
							
		#           # If z=1, it will be default be placed in overflow bin -- prevent this
		#           if np.isclose(z_truth, 1.):
		#             z_truth = 0.999
							
		#           # In Pb-Pb case, fill matched pt fraction
		#           if not self.is_pp:
		#             self.fill_subjet_matched_pt_histograms(observable,
		#                                                    subjet_det, subjet_truth,
		#                                                    z_det, z_truth,
		#                                                    jet_truth.pt(), jetR, subjetR, R_max)
							
		#           # Fill histograms
		#           # Note that we don't fill 'matched' histograms here, since that is only
		#           # meaningful for leading subjets
		#           self.fill_response(observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
		#                              z_det, z_truth, obs_label, R_max, prong_match=False)
																 
		#       # Fill number of subjets with/without unique match, as a function of zr
		#       if self.is_pp:
		#         name = 'h_match_fraction_{}_R{}_{}'.format(observable, jetR, subjetR)
		#         getattr(self, name).Fill(jet_truth.pt(), z_det, successful_match)
															 
		#   # Get leading subjet and fill histograms
		#   if 'leading' in observable:
			
		#     leading_subjet_det = self.utils.leading_jet(subjets_det)
		#     leading_subjet_truth = self.utils.leading_jet(subjets_truth)
				
		#     # Note that we don't want to check whether they are geometrical matches
		#     # We rather want to correct the measured leading subjet to the true leading subjet
		#     if leading_subjet_det and leading_subjet_truth:
					
		#       z_leading_det = leading_subjet_det.pt() / jet_det.pt()
		#       z_leading_truth = leading_subjet_truth.pt() / jet_truth.pt()
					
		#       # If z=1, it will be default be placed in overflow bin -- prevent this
		#       if np.isclose(z_leading_det, 1.):
		#         z_leading_det = 0.999
		#       if np.isclose(z_leading_truth, 1.):
		#         z_leading_truth = 0.999
					
		#       # In Pb-Pb case, fill matched pt fraction
		#       if not self.is_pp:
		#         match = self.fill_subjet_matched_pt_histograms(observable,
		#                                                        leading_subjet_det, leading_subjet_truth,
		#                                                        z_leading_det, z_leading_truth,
		#                                                        jet_truth.pt(), jetR, subjetR, R_max)
		#       else:
		#         match = False
					
		#       # Fill histograms
		#       self.fill_response(observable, jetR, jet_pt_det_ungroomed, jet_pt_truth_ungroomed,
		#                          z_leading_det, z_leading_truth, obs_label, R_max, prong_match=match)
					
		#       # Plot deltaR distribution between the detector and truth leading subjets
		#       # (since they are not matched geometrically, the true leading may not be the measured leading
		#       deltaR = leading_subjet_det.delta_R(leading_subjet_truth)
		#       name = 'hDeltaR_det_truth_{}_R{}_{}'.format(observable, jetR, subjetR)
		#       if not self.is_pp:
		#         name += '_Rmax{}'.format(R_max)
		#       getattr(self, name).Fill(jet_truth.pt(), z_leading_truth, deltaR)
				
	#---------------------------------------------------------------
	# Do prong-matching
	#---------------------------------------------------------------
	# def fill_subjet_matched_pt_histograms(self, observable, subjet_det, subjet_truth,
	#                                       z_det, z_truth, jet_pt_truth, jetR, subjetR, R_max):
		
		# # Get pp det-level subjet
		# # Inclusive case: This is matched to the combined subjet (and its pp truth-level subjet)
		# # Leading case: This is matched only to the pp truth-level leading subjet
		# subjet_pp_det = None
		# if subjet_truth.has_user_info():
		#   subjet_pp_det = subjet_truth.python_info().match
		# if not subjet_pp_det:
		#   return
																		 
		# matched_pt = fjtools.matched_pt(subjet_det, subjet_pp_det)
		# name = 'h_{}_matched_pt_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
		# getattr(self, name).Fill(jet_pt_truth, z_det, matched_pt)
		
		# # Plot dz between det and truth subjets
		# deltaZ = z_det - z_truth
		# name = 'h_{}_matched_pt_deltaZ_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
		# getattr(self, name).Fill(jet_pt_truth, matched_pt, deltaZ)

		# # Plot dR between det and truth subjets
		# deltaR = subjet_det.delta_R(subjet_truth)
		# name = 'h_{}_matched_pt_deltaR_JetPt_R{}_{}_Rmax{}'.format(observable, jetR, subjetR, R_max)
		# getattr(self, name).Fill(jet_pt_truth, matched_pt, deltaR)

		# match = (matched_pt > 0.5)
		# return match

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
	def fill_perp_cone_histograms(self, cone_parts, cone_R, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting, obs_label, jet_pt_ungroomed, suffix, rho_bge = 0, truth=False):

		if self.is_pA:
			suffix = ''

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
		if truth:
			hname = 'h_perpcone{}_{}_JetPt_Truth_R{}_{}{}'

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
						if self.do_median_subtraction:
							pair_type = self.check_pair_type(new_corr, ipoint, c_select, index)
							pair_type_label = self.pair_type_labels[pair_type]

						if 'ENC' in observable:
							# print('hname is',hname.format(cone_R, observable + str(ipoint) + pair_type_label, jetR, obs_label))
							getattr(self, hname.format(cone_R, observable + str(ipoint) + pair_type_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index])
							getattr(self, hname.format(cone_R, observable + str(ipoint) + pair_type_label + 'Pt', jetR, obs_label, suffix)).Fill(jet_pt, jet_pt*new_corr.correlator(ipoint).rs()[index], new_corr.correlator(ipoint).weights()[index]) # NB: fill pt*RL

						if ipoint==2 and 'EEC_noweight' in observable:
							getattr(self, hname.format(cone_R, observable + pair_type_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index])

						if ipoint==2 and 'EEC_weight2' in observable:
							getattr(self, hname.format(cone_R, observable + pair_type_label, jetR, obs_label, suffix)).Fill(jet_pt, new_corr.correlator(ipoint).rs()[index], pow(new_corr.correlator(ipoint).weights()[index],2))
       

##################################################################
if __name__ == '__main__':
	# Define arguments
	parser = argparse.ArgumentParser(description='Process MC')
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

	# If invalid inputFile is given, exit
	if not os.path.exists(args.inputFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.inputFile))
		sys.exit(0)
	
	# If invalid configFile is given, exit
	if not os.path.exists(args.configFile):
		print('File \"{0}\" does not exist! Exiting!'.format(args.configFile))
		sys.exit(0)

	analysis = ProcessMC_ENC(input_file=args.inputFile, config_file=args.configFile, output_dir=args.outputDir)
	analysis.process_mc()
