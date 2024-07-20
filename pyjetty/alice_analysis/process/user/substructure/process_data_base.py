#!/usr/bin/env python3

"""
Base class to read a ROOT TTree of track information
and do jet-finding, and save basic histograms.
	
To use this class, the following should be done:

	- Implement a user analysis class inheriting from this one, such as in user/james/process_data_XX.py
		You should implement the following functions:
			- initialize_user_output_objects()
			- fill_jet_histograms()
		
	- The histogram of the data should be named h_[obs]_JetPt_R[R]_[subobs]_[grooming setting]
		The grooming part is optional, and should be labeled e.g. zcut01_B0 — from CommonUtils::grooming_label({'sd':[zcut, beta]})
		For example: h_subjet_z_JetPt_R0.4_0.1
		For example: h_subjet_z_JetPt_R0.4_0.1_zcut01_B0

	- You also should modify observable-specific functions at the top of common_utils.py
	
Author: James Mulligan (james.mulligan@berkeley.edu)
"""

from __future__ import print_function

# General
import time

# Data analysis and plotting
import numpy as np
import ROOT
import yaml

# Fastjet via python (from external library heppy)
import fastjet as fj
import fjcontrib

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_base
from pyjetty.mputils.csubtractor import CEventSubtractor

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class ProcessDataBase(process_base.ProcessBase):

	#---------------------------------------------------------------
	# Constructor
	#---------------------------------------------------------------
	def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
	
		# Initialize base class
		super(ProcessDataBase, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)

		# Initialize configuration
		self.initialize_config()
		
	#---------------------------------------------------------------
	# Initialize config file into class members
	#---------------------------------------------------------------
	def initialize_config(self):
		
		# Call base class initialization
		process_base.ProcessBase.initialize_config(self)
		
		# Read config file
		with open(self.config_file, 'r') as stream:
			config = yaml.safe_load(stream)

		self.do_median_subtraction = config['do_median_subtraction']
		if self.do_median_subtraction:
			print("DOING MEDIAN SUBTRACTION")

		self.do_perpendicular_cone = config['do_perpendicular_cone']
		self.randomize_cone = config['randomize_cone']
		if self.do_perpendicular_cone:
			print("DOING PERPENDICULAR CONE")
		
		if self.do_median_subtraction or self.do_perpendicular_cone:
			self.is_pp = False
			self.is_pA = True
		elif self.do_constituent_subtraction:
			self.is_pp = False
			self.is_pA = False
		else:
			self.is_pp = True
			self.is_pA = False
		
		self.mult_threshold = config['mult_threshold']

		if 'ENC_pair_cut' in config:
				self.ENC_pair_cut = config['ENC_pair_cut']
		else:
				self.ENC_pair_cut = False

		self.do_reshuffle = config['do_reshuffle']
		
		# Create dictionaries to store grooming settings and observable settings for each observable
		# Each dictionary entry stores a list of subconfiguration parameters
		#   The observable list stores the observable setting, e.g. subjetR
		#   The grooming list stores a list of grooming settings {'sd': [zcut, beta]} or {'dg': [a]}
		self.observable_list = config['process_observables']
		self.obs_settings = {}
		self.obs_grooming_settings = {}
		for observable in self.observable_list:
		
			obs_config_dict = config[observable]
			obs_config_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
			
			obs_subconfig_list = [name for name in list(obs_config_dict.keys()) if 'config' in name ]
			self.obs_settings[observable] = self.utils.obs_settings(observable, obs_config_dict, obs_subconfig_list)
			self.obs_grooming_settings[observable] = self.utils.grooming_settings(obs_config_dict)
			
		# Construct set of unique grooming settings
		self.grooming_settings = []
		lists_grooming = [self.obs_grooming_settings[obs] for obs in self.observable_list]
		for observable in lists_grooming:
			for setting in observable:
				if setting not in self.grooming_settings and setting != None:
					self.grooming_settings.append(setting)
					
	#---------------------------------------------------------------
	# Main processing function
	#---------------------------------------------------------------
	def process_data(self):
		
		self.start_time = time.time()
		
		# Use IO helper class to convert ROOT TTree into a SeriesGroupBy object of fastjet particles per event
		print('--- {} seconds ---'.format(time.time() - self.start_time))
		io = process_io.ProcessIO(input_file=self.input_file, track_tree_name='tree_Particle',
															is_pp=(self.is_pp or self.is_pA), min_cent=0., max_cent=10., use_ev_id_ext=True)
		self.df_fjparticles, self.df_events = io.load_data(m=self.m)
		self.nEvents = len(self.df_fjparticles.index)
		# print("DF NEVENTS DEEBOOG:",self.nEvents)
		self.nTracks = len(io.track_df.index)
		print('--- {} seconds ---'.format(time.time() - self.start_time))
		
		# Initialize histograms
		self.initialize_output_objects()
		
		# Create constituent subtractor, if configured
		if not self.is_pp and not self.is_pA:
			self.constituent_subtractor = [CEventSubtractor(max_distance=R_max, alpha=self.alpha, max_eta=self.max_eta, bge_rho_grid_size=self.bge_rho_grid_size, max_pt_correct=self.max_pt_correct, ghost_area=self.ghost_area, distance_type=fjcontrib.ConstituentSubtractor.deltaR) for R_max in self.max_distance]
		
		# Create median subtractor, if configured
		if self.do_median_subtraction:
			self.jet_def_medsub = {}
			self.jet_selector_medsub = {}
			self.median_subtractor = {}
			for jetR in self.jetR_list:
				self.jet_def_medsub[jetR] = fj.JetDefinition(fj.kt_algorithm, jetR)
				self.jet_selector_medsub[jetR] = fj.SelectorAbsRapMax(0.9 - jetR) & (~fj.SelectorNHardest(2)) & (~fj.SelectorIsPureGhost())
				self.median_subtractor[jetR] = fj.JetMedianBackgroundEstimator(self.jet_selector_medsub[jetR], self.jet_def_medsub[jetR], fj.AreaDefinition(fj.active_area_explicit_ghosts))


		print(self)

		# Find jets and fill histograms
		print('Analyze events...')
		self.analyze_events()
		
		# Plot histograms
		print('Save histograms...')
		process_base.ProcessBase.save_output_objects(self)

		print('--- {} seconds ---'.format(time.time() - self.start_time))
	
	#---------------------------------------------------------------
	# Initialize histograms
	#---------------------------------------------------------------
	def initialize_output_objects(self):
	
		# Initialize user-specific histograms
		self.initialize_user_output_objects()
		
		# Initialize base histograms
		self.hNevents = ROOT.TH1F('hNevents', 'hNevents', 2, -0.5, 1.5)
		if self.event_number_max < self.nEvents:
			self.hNevents.Fill(1, self.event_number_max)
		else:
			self.hNevents.Fill(1, self.nEvents)
		
		self.hTrackEtaPhi = ROOT.TH2F('hTrackEtaPhi', 'hTrackEtaPhi', 200, -1., 1., 628, 0., 6.28)
		self.hTrackPt = ROOT.TH1F('hTrackPt', 'hTrackPt', 300, 0., 300.)
		
		if not self.is_pp and not self.is_pA:
			self.hRho = ROOT.TH1F('hRho', 'hRho', 1000, 0., 1000.)
		
		self.hMult = ROOT.TH1F('hMult', 'hMult', 1200, 0., 1200.)
				
		for jetR in self.jetR_list:
			
			name = 'hZ_R{}'.format(jetR)
			h = ROOT.TH2F(name, name, 300, 0, 300, 100, 0., 1.)
			setattr(self, name, h)

			if self.is_pA:
				if self.do_constituent_subtraction:
					name = 'hRhoCS_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, 0, 20.)
					setattr(self, name, h)

				if self.do_median_subtraction:
					name = 'hMedRho_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, 0, 20.)
					setattr(self, name, h)

					name = 'hMedRhoCArea_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, 0, 20.)
					setattr(self, name, h)

					name = 'hMedRho_SigJetEvents_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, 0, 20.)
					setattr(self, name, h)

					name = 'hRhoCS_SigJetEvents_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, 0, 20.)
					setattr(self, name, h)

					name = 'hCArea_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 100, 0, 1.)
					setattr(self, name, h)

					name = 'hNJets_MedRho_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 30, 0, 30., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hMedJetpT_Rho_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, 0, 100., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hSigJetpT_MedRho_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 100, 0, 100., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hJetconstituents_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, 0, 200, 30, 0., 30.)
					setattr(self, name, h)

					name = 'hJetconstituent_pT_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, 0, 200, 100, 0., 100.)
					setattr(self, name, h)

				if self.do_perpendicular_cone:
					name = 'hJetArea_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 100, 0, 4.)
					setattr(self, name, h)

					name = 'hPerpConepT_Raw_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 400, 0., 20.)
					setattr(self, name, h)

					name = 'hPerpConepT_MedCorrected_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, 0, 200., 400, -5., 15.)
					setattr(self, name, h)

					name = 'hPerpConepT_MednoC_Corrected_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 400, -5., 15.)
					setattr(self, name, h)

					name = 'hPerpConepT_CSCorrected_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 400, -5., 15.)
					setattr(self, name, h)

					name = 'hPerpConepTRes_MedCorrected_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, -1., 1.)
					setattr(self, name, h)

					name = 'hPerpConepTComp_MedCorrected_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, -5., 15., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hPerpConepTComp_CSCorrected_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, -5., 15., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hPerpRho_R{}'.format(jetR)
					h = ROOT.TH1F(name, name, 200, 0, 20.)
					setattr(self, name, h)

					name = 'hSigJetpt_PerpRho_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 100, 0, 100., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hMedRho_PerpRho_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 200, 0, 20., 200, 0, 20.)
					setattr(self, name, h)

					name = 'hPerpConeMult_Rho_R{}'.format(jetR)
					h = ROOT.TH2F(name, name, 20, 0, 20., 200, 0, 20.)
					setattr(self, name, h)

					for ptlo, pthi in [(20,40),(40,60),(60,80)]:
						name = 'hPerpConeMult_JetMult_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH2F(name, name, 20, 0, 20., 100, 0, 100.)
						setattr(self, name, h)

						name = 'hPerpConeMult_JetMult_150_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH2F(name, name, 20, 0, 20., 100, 0, 100.)
						setattr(self, name, h)

						name = 'hPerpConeMult_JetMult_1GeV_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH2F(name, name, 20, 0, 20., 100, 0, 100.)
						setattr(self, name, h)

						name = 'hPerpConeMult_JetMult_2GeV_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH2F(name, name, 20, 0, 20., 100, 0, 100.)
						setattr(self, name, h)

						name = 'hPerpConeRho_150_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH1F(name, name, 200, 0, 20.)
						setattr(self, name, h)

						name = 'hPerpConeRho_1GeV_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH1F(name, name, 200, 0, 20.)
						setattr(self, name, h)

						name = 'hPerpConeRho_2GeV_{}{}_R{}'.format(ptlo,pthi,jetR)
						h = ROOT.TH1F(name, name, 200, 0, 20.)
						setattr(self, name, h)

	#---------------------------------------------------------------
	# Main function to loop through and analyze events
	#---------------------------------------------------------------
	def analyze_events(self):
		
		# Fill track histograms
		print('--- {} seconds ---'.format(time.time() - self.start_time))
		print('Fill track histograms')
		[[self.fillTrackHistograms(track) for track in fj_particles] for fj_particles in self.df_fjparticles]
		print('--- {} seconds ---'.format(time.time() - self.start_time))
		
		print('Find jets...')
		fj.ClusterSequence.print_banner()
		print()
		self.event_number = 0
	
		# Use list comprehension to do jet-finding and fill histograms
		result = [self.analyze_event(fj_particles) for fj_particles in self.df_fjparticles]
		
		print('--- {} seconds ---'.format(time.time() - self.start_time))
		print('Save thn...')
		process_base.ProcessBase.save_thn_th3_objects(self)
	
	#---------------------------------------------------------------
	# Fill track histograms.
	#---------------------------------------------------------------
	def fillTrackHistograms(self, track):
		
		self.hTrackEtaPhi.Fill(track.eta(), track.phi())
		self.hTrackPt.Fill(track.pt())
	
	#---------------------------------------------------------------
	# Analyze jets of a given event.
	# fj_particles is the list of fastjet pseudojets for a single fixed event.
	#---------------------------------------------------------------
	def analyze_event(self, fj_particles):
	
		self.event_number += 1
		if self.event_number > self.event_number_max:
			return
		if self.debug_level > 1:
			print('-------------------------------------------------')
			print('event {}'.format(self.event_number))
		
		if len(fj_particles) > 1:
			if np.abs(fj_particles[0].pt() - fj_particles[1].pt()) <  1e-10:
				print('WARNING: Duplicate particles may be present')
				print([p.user_index() for p in fj_particles])
				print([p.pt() for p in fj_particles])
	
		# Perform constituent subtraction for each R_max (do this once, for all jetR)
		if not self.is_pp and not self.is_pA:
			fj_particles_subtracted = [self.constituent_subtractor[i].process_event(fj_particles) for i, R_max in enumerate(self.max_distance)]
		
		self.mult = self.df_events['V0Amult'].values[self.event_number-1]
		getattr(self, 'hMult').Fill(self.mult)
		if type(self.mult_threshold) is list:
			self.isHighMult = self.mult > self.mult_threshold[1]
			self.isLowMult = self.mult <= self.mult_threshold[0]
		else:
			self.isHighMult = self.mult > self.mult_threshold
			self.isLowMult = self.mult <= self.mult_threshold

		# Loop through jetR, and process event for each R
		for jetR in self.jetR_list:
		
			# Keep track of whether to fill R-independent histograms
			self.fill_R_indep_hists = (jetR == self.jetR_list[0])

			# Set jet definition and a jet selector
			jet_def = fj.JetDefinition(fj.antikt_algorithm, jetR)
			jet_selector = fj.SelectorPtMin(5.0) & fj.SelectorAbsRapMax(0.9 - jetR)
			if self.debug_level > 2:
				print('jet definition is:', jet_def)
				print('jet selector is:', jet_selector,'\n')

			if self.do_median_subtraction:
				csa_medsub = fj.ClusterSequenceArea(fj_particles, self.jet_def_medsub[jetR], fj.AreaDefinition(fj.active_area_explicit_ghosts))
				self.median_subtractor[jetR].set_cluster_sequence(csa_medsub)
				rho = self.median_subtractor[jetR].rho()
				getattr(self, 'hMedRho_R{}'.format(jetR)).Fill(rho)

				Cjet_selector = fj.SelectorAbsRapMax(0.9 - jetR) & (~fj.SelectorIsPureGhost())
				selected_jets = fj.sorted_by_pt(Cjet_selector(csa_medsub.inclusive_jets()))
				n_activejets = len(selected_jets)
				C_area = 0
				for jet in selected_jets:
					C_area += jet.area()
				C_area = C_area / (2*np.pi*1.8)
				getattr(self, 'hNJets_MedRho_R{}'.format(jetR)).Fill(n_activejets, rho)
				getattr(self, 'hCArea_R{}'.format(jetR)).Fill(C_area)
				getattr(self, 'hMedRhoCArea_R{}'.format(jetR)).Fill(rho*C_area)
				medianjet_pt = np.median([jet.pt() for jet in selected_jets])
				getattr(self, 'hMedJetpT_Rho_R{}'.format(jetR)).Fill(medianjet_pt,rho)

			# Analyze
			if self.is_pp:
			
				# Do jet finding
				cs = fj.ClusterSequenceArea(fj_particles, jet_def, fj.AreaDefinition(fj.VoronoiAreaSpec()))
				jets = fj.sorted_by_pt(cs.inclusive_jets())
				jets_selected = jet_selector(jets)
			
				self.analyze_jets(jets_selected, jetR)

			elif self.is_pA:

				# Perform constituent subtraction
				rhocs = 0
				# if self.do_constituent_subtraction:
				# 	for i, R_max in enumerate(self.max_distance):
				# 		rhocs = self.constituent_subtractor[i].bge_rho.rho()
				# 		getattr(self, 'hRhoCS_R{}'.format(jetR)).Fill(rhocs)

				# Jet finding
				cs = fj.ClusterSequenceArea(fj_particles, jet_def, fj.AreaDefinition(fj.VoronoiAreaSpec()))
				jets = fj.sorted_by_pt(cs.inclusive_jets())
				jets_selected = jet_selector(jets)

				med_rhoC = [0,0]
				if self.do_median_subtraction:
					med_rhoC = [rho,C_area]

					jets_reselected = []
					if rho*C_area > 0:
						for jet in jets_selected:
							if jet.pt()-rho*C_area*jet.area() > 5:
								getattr(self, 'hSigJetpT_MedRho_R{}'.format(jetR)).Fill(jet.pt()-rho*C_area*jet.area(),rho*C_area)
								getattr(self, 'hMedRho_SigJetEvents_R{}'.format(jetR)).Fill(rho)
								jets_reselected.append(jet)

					self.analyze_jets(jets_reselected, jetR, rho_bge=rho*C_area)
					if self.do_perpendicular_cone:
						self.analyze_perpendicular_cone(jets_reselected, jetR, fj_particles, *med_rhoC, rhocs)
						self.analyze_perp_cones(fj_particles, jets_reselected, jetR, rho_bge = rho*C_area)
				
				else:
					self.analyze_jets(jets_selected, jetR)
					if self.do_perpendicular_cone:
						self.analyze_perpendicular_cone(jets_selected, jetR, fj_particles)
						self.analyze_perp_cones(fj_particles, jets_selected, jetR, R_max = R_max)
				
			else:
			
				for i, R_max in enumerate(self.max_distance):
										
					if self.debug_level > 1:
						print('R_max: {}'.format(R_max))
						
					# Keep track of whether to fill R_max-independent histograms
					self.fill_Rmax_indep_hists = (i == 0)
					
					# Perform constituent subtraction
					rho = self.constituent_subtractor[i].bge_rho.rho()
					if self.fill_R_indep_hists and self.fill_Rmax_indep_hists:
						getattr(self, 'hRho').Fill(rho)
					
					# Do jet finding (re-do each time, to make sure matching info gets reset)
					cs = fj.ClusterSequence(fj_particles_subtracted[i], jet_def)
					jets = fj.sorted_by_pt(cs.inclusive_jets())
					jets_selected = jet_selector(jets)
					
					self.analyze_jets(jets_selected, jetR, R_max = R_max)
			
	#---------------------------------------------------------------
	# Analyze jets of a given event.
	#---------------------------------------------------------------
	def analyze_jets(self, jets_selected, jetR, rho_bge=0, R_max = None):
	
		# Set suffix for filling histograms
		if R_max:
			suffix = '_Rmax{}'.format(R_max)
		else:
			suffix = ''
		
		# Loop through jets and call user function to fill histos
		result = [self.analyze_accepted_jet(jet, jetR, suffix, rho_bge=rho_bge) for jet in jets_selected]
	
	#---------------------------------------------------------------
	# Fill histograms
	#---------------------------------------------------------------
	def analyze_accepted_jet(self, jet, jetR, suffix, rho_bge=0):
		
		# Check additional acceptance criteria
		if not self.utils.is_det_jet_accepted(jet):
			return
		
		# Fill base histograms
		jet_pt_subtracted = jet.pt() - rho_bge*jet.area()
		jet_pt_ungroomed = jet.pt()
		if self.is_pp or self.is_pA or self.fill_Rmax_indep_hists:
			if self.do_median_subtraction:
				hZ = getattr(self, 'hZ_R{}'.format(jetR))
				nconstituents_1gev = 0
				for constituent in jet.constituents():
					z = constituent.pt() / jet_pt_subtracted
					hZ.Fill(jet_pt_subtracted, z)
					if constituent.pt() > 1: 
						getattr(self, 'hJetconstituent_pT_R{}'.format(jetR)).Fill(jet_pt_subtracted, constituent.pt())
						nconstituents_1gev += 1

				getattr(self, 'hJetconstituents_R{}'.format(jetR)).Fill(jet_pt_subtracted, nconstituents_1gev)
		
		# Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
		# Note that the subconfigurations are defined by the first observable, if multiple are defined
		observable = self.observable_list[0]
		for i in range(len(self.obs_settings[observable])):
		
			obs_setting = self.obs_settings[observable][i]
			grooming_setting = self.obs_grooming_settings[observable][i]
			obs_label = self.utils.obs_label(obs_setting, grooming_setting)
		
			# Groom jet, if applicable
			if grooming_setting:
				gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
				jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
				if not jet_groomed_lund:
					continue
			else:
				jet_groomed_lund = None

			# Call user function to fill histograms
			self.fill_jet_histograms(jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
															 obs_label, jet_pt_subtracted, suffix)
			
	#---------------------------------------------------------------
	# Calculate rho from perpendicular cone 
	# and call user histograms
	#---------------------------------------------------------------
	def analyze_perpendicular_cone(self, jets_selected, jetR, fj_particles, rho=0, C_area=0, rhocs = 0, R_max = None):
	
		# Set suffix for filling histograms
		if R_max:
			suffix = '_Rmax{}'.format(R_max)
		else:
			suffix = ''

		# Loop through jets and call user function on particles in perp cone
		for jet in jets_selected:
			
			if not self.utils.is_det_jet_accepted(jet):
				continue
			
			jet_area = jet.area()
			getattr(self, 'hJetArea_R{}'.format(jetR)).Fill(jet_area)
			
			jet_Rsq = jet_area/np.pi
			perp_area = np.pi*0.4*0.4
			perp_cone_phis = [(jet.phi()-np.pi/2)%(2*np.pi),(jet.phi()+np.pi/2)%(2*np.pi)]
			perp_cone_eta = jet.eta()
			
			for perp_cone_phi in perp_cone_phis:
				
				perp_cone_particles = fj.vectorPJ()
				perp_cone_pj = fj.PseudoJet()
				for particle in fj_particles:
					if (np.square(particle.eta()-perp_cone_eta)+np.square(particle.phi()-perp_cone_phi) < 0.4*0.4): #jet_Rsq):   --- changed so it's a more precise jet radius
						perp_cone_particles.append(particle)
						perp_cone_pj += particle
			
				perp_cone_pt = perp_cone_pj.pt()
				perp_cone_rho = perp_cone_pt / perp_area #jet_area
				getattr(self, 'hPerpRho_R{}'.format(jetR)).Fill(perp_cone_rho)
				getattr(self, 'hPerpConepT_Raw_R{}'.format(jetR)).Fill(perp_cone_pt)

				nperp_cone_particles = len(perp_cone_particles)
				if self.do_median_subtraction and nperp_cone_particles != 0:
					perp_cone_ptcorr = perp_cone_pt - rho*C_area*perp_area #jet_area
					perp_cone_ptcorr_woC = perp_cone_pt - rho*perp_area #jet_area
					perp_cone_ptcorr_CS = perp_cone_pt - rhocs*perp_area #jet_area
					getattr(self, 'hPerpConepT_MedCorrected_R{}'.format(jetR)).Fill(jet.pt()-rho*jet_area, perp_cone_ptcorr)
					getattr(self, 'hPerpConepT_MednoC_Corrected_R{}'.format(jetR)).Fill(perp_cone_ptcorr_woC)
					getattr(self, 'hPerpConepT_CSCorrected_R{}'.format(jetR)).Fill(perp_cone_ptcorr_CS)
					getattr(self, 'hPerpConepTRes_MedCorrected_R{}'.format(jetR)).Fill(perp_cone_ptcorr/perp_cone_pt)
					getattr(self, 'hPerpConepTComp_MedCorrected_R{}'.format(jetR)).Fill(perp_cone_ptcorr,perp_cone_pt)
					getattr(self, 'hPerpConepTComp_CSCorrected_R{}'.format(jetR)).Fill(perp_cone_ptcorr_CS,perp_cone_pt)
					getattr(self, 'hMedRho_PerpRho_R{}'.format(jetR)).Fill(rho,perp_cone_rho)
			 
				getattr(self, 'hSigJetpt_PerpRho_R{}'.format(jetR)).Fill(jet.pt()-rho*jet_area, perp_cone_rho)
				getattr(self, 'hPerpConeMult_Rho_R{}'.format(jetR)).Fill(nperp_cone_particles, perp_cone_rho)

				nperp_150 = [p for p in perp_cone_particles if p.pt()>0.15]
				nperp_1gev = [p for p in perp_cone_particles if p.pt()>1]
				nperp_2gev = [p for p in perp_cone_particles if p.pt()>2]
				njet_150 = [p for p in jet.constituents() if p.pt()>0.15]
				njet_1gev = [p for p in jet.constituents() if p.pt()>1]
				njet_2gev = [p for p in jet.constituents() if p.pt()>2]
				
				for ptlo, pthi in [(20,40),(40,60),(60,80)]:
					if jet.pt()-rho*jet_area < pthi and jet.pt()-rho*jet_area > ptlo:
						getattr(self, 'hPerpConeMult_JetMult_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(nperp_cone_particles, len(jet.constituents()))
						getattr(self, 'hPerpConeMult_JetMult_150_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(len(nperp_150), len(njet_150))
						getattr(self, 'hPerpConeMult_JetMult_1GeV_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(len(nperp_1gev), len(njet_1gev))
						getattr(self, 'hPerpConeMult_JetMult_2GeV_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(len(nperp_2gev), len(njet_2gev))

						getattr(self, 'hPerpConeRho_150_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(np.sum([p.pt() for p in nperp_150])/perp_area)
						getattr(self, 'hPerpConeRho_1GeV_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(np.sum([p.pt() for p in nperp_1gev])/perp_area)
						getattr(self, 'hPerpConeRho_2GeV_{}{}_R{}'.format(ptlo,pthi,jetR)).Fill(np.sum([p.pt() for p in nperp_2gev])/perp_area)
				
				observable = self.observable_list[0]
				for i in range(len(self.obs_settings[observable])):
					obs_setting = self.obs_settings[observable][i]
					grooming_setting = self.obs_grooming_settings[observable][i]
					obs_label = self.utils.obs_label(obs_setting, grooming_setting)
					self.fill_perpcone_histograms(jet, jetR, perp_cone_particles, obs_setting, obs_label, suffix, rho_bge=rho)

	#---------------------------------------------------------------
	# Helper functions for perpendicular cone 
	# background pair subtraction
	#---------------------------------------------------------------
	def find_parts_around_jet(self, parts, jet, cone_R):
		# select particles around jet axis
		cone_parts = fj.vectorPJ()
		for part in parts:
			if jet.delta_R(part) <= cone_R:
				cone_parts.push_back(part)
		
		return cone_parts

	def rotate_parts(self, parts, rotate_phi):
		# rotate parts in azimuthal direction
		parts_rotated = fj.vectorPJ()
		for part in parts:
			pt_new = part.pt()
			y_new = part.rapidity()
			phi_new = part.phi() + rotate_phi
			m_new = part.m()
			# print('before',part.phi())
			part.reset_PtYPhiM(pt_new, y_new, phi_new, m_new)
			# print('after',part.phi())
			parts_rotated.push_back(part)
		
		return parts_rotated

	#---------------------------------------------------------------
	# Using perpendicular cone to
	# do background pair subtraction
	#---------------------------------------------------------------	
	def analyze_perp_cones(self, parts, jets_selected, jetR, R_max = None, rho_bge = 0):
		# analyze cones perpendicular to jet in the azimuthal plane
		# if R_max and (not self.do_median_subtraction):
		# 	suffix = '_Rmax{}'.format(R_max)
		# else:
		# 	suffix = ''
		suffix = ''


		perpcone_R_list = [jetR] # FIX ME: not sure it's better to just use jet R or calculate from area, For now just use jet R so it's more consistent with the jet cone method
		# if self.do_median_subtraction and rho_bge > 0:
		#   perpcone_R_list = [math.sqrt(jet.area()/np.pi)] # NB: jet area is available only when rho subtraction flag is on

		for jet in jets_selected:

			angle = np.pi/2
			if self.randomize_cone:
				angle = np.random.uniform(low=np.pi/3, high= 2*np.pi/3)
				print(angle)

			# print('jet pt',jet.perp()-rho_bge*jet.area(),'phi',jet.phi(),'eta',jet.eta(),'area',jet.area())
			perp_jet1 = fj.PseudoJet()
			perp_jet1.reset_PtYPhiM(jet.pt(), jet.rapidity(), jet.phi() + angle, jet.m())
			perp_jet2 = fj.PseudoJet()
			perp_jet2.reset_PtYPhiM(jet.pt(), jet.rapidity(), jet.phi() - angle, jet.m())

			for perpcone_R in perpcone_R_list:

				constituents = jet.constituents()
				# FIX ME: current implemetation is to use jet constituents as "signal" for perp cone if cone radius == jetR, else use jet cone as "signal" for perp cone. May want to implement both jet and jet cone later for radius = jet R case
				# if perpcone_R != jetR:
					# constituents = self.find_parts_around_jet(parts, jet, jetcone_R)
				parts_in_perpcone1 = self.find_parts_around_jet(parts, perp_jet1, perpcone_R)
				# for part in parts_in_perpcone1:
				#   print('before rotation',part.phi())
				parts_in_perpcone1 = self.rotate_parts(parts_in_perpcone1, -angle)
				# for part in parts_in_perpcone1:
				#   print('after rotation',part.phi())
				
				parts_in_perpcone2 = self.find_parts_around_jet(parts, perp_jet2, perpcone_R)
				parts_in_perpcone2 = self.rotate_parts(parts_in_perpcone2, +angle)

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

				self.analyze_accepted_cone(True, parts_in_cone1, perpcone_R, jet, jetR, suffix, rho_bge)
				self.analyze_accepted_cone(True, parts_in_cone2, perpcone_R, jet, jetR, suffix, rho_bge)

	def analyze_accepted_cone(self, is_perp, cone_parts, cone_R, jet, jetR, suffix, rho_bge = 0):
		
		# Check additional acceptance criteria
		if not self.utils.is_det_jet_accepted(jet):
			return
					
		# Fill base histograms
		if self.do_median_subtraction and rho_bge > 0:
			jet_pt_ungroomed = jet.pt() - rho_bge*jet.area()
		else:
			jet_pt_ungroomed = jet.pt()
		
		# Loop through each jet subconfiguration (i.e. subobservable / grooming setting)
		# Note that the subconfigurations are defined by the first observable, if multiple are defined
		observable = self.observable_list[0]
		for i in range(len(self.obs_settings[observable])):
		
			obs_setting = self.obs_settings[observable][i]
			grooming_setting = self.obs_grooming_settings[observable][i]
			obs_label = self.utils.obs_label(obs_setting, grooming_setting)
		
			# Groom jet, if applicable
			if grooming_setting:
				gshop = fjcontrib.GroomerShop(jet, jetR, self.reclustering_algorithm)
				jet_groomed_lund = self.utils.groom(gshop, grooming_setting, jetR)
				if not jet_groomed_lund:
					continue
			else:
				jet_groomed_lund = None

			# Call user function to fill histograms
			self.fill_perp_cone_histograms(cone_parts, cone_R, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
															 obs_label, jet_pt_ungroomed, suffix, rho_bge)

	
	#---------------------------------------------------------------
	# This function is called once
	# You must implement this
	#---------------------------------------------------------------
	def initialize_user_output_objects(self):
	
		raise NotImplementedError('You must implement initialize_user_output_objects()!')

	#---------------------------------------------------------------
	# This function is called once for each jet subconfiguration
	# You must implement this
	#---------------------------------------------------------------
	def fill_jet_histograms(self, jet, jet_groomed_lund, jetR, obs_setting, grooming_setting,
													obs_label, jet_pt_ungroomed, suffix):
	
		raise NotImplementedError('You must implement fill_jet_histograms()!')
