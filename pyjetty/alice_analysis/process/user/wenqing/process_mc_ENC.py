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

# Analysis utilities
from pyjetty.alice_analysis.process.base import process_io
from pyjetty.alice_analysis.process.base import process_io_emb
from pyjetty.alice_analysis.process.base import jet_info
from pyjetty.alice_analysis.process.user.substructure import process_mc_base
from pyjetty.alice_analysis.process.base import thermal_generator
from pyjetty.mputils import CEventSubtractor

def linbins(xmin, xmax, nbins):
  lspace = np.linspace(xmin, xmax, nbins+1)
  arr = array.array('f', lspace)
  return arr

def logbins(xmin, xmax, nbins):
  lspace = np.logspace(np.log10(xmin), np.log10(xmax), nbins+1)
  arr = array.array('f', lspace)
  return arr

def calculate_distance(p0, p1):
    
  dphiabs = math.fabs(p0.phi() - p1.phi())
  dphi = dphiabs

  if dphiabs > math.pi:
    dphi = 2*math.pi - dphiabs


  dy = p0.rap() - p1.rap()
  return math.sqrt(dy*dy + dphi*dphi)

class CorrelatorBuilder:
  def __init__(self, particle_list, jetpt):
    self.particle_list = particle_list
    self.mult = len(self.particle_list)
    self.pair_list = np.empty((0,self.mult),float)
    self.scale = jetpt

  def make_pairs(self):
    for i,part_i in enumerate(self.particle_list):
      inner_list = np.array([])
      for j,part_j in enumerate(self.particle_list[i+1:]):
        
        dist = calculate_distance(part_i,part_j)
        # print(' pairing particle i = ', i, ' and j = ', j+i+1, ' with distance ', dist)
        inner_list = np.append(inner_list, dist)

      inner_list = np.pad(inner_list, (i+1, 0), 'constant')
      inner_list.resize(1,self.mult)
      self.pair_list = np.append(self.pair_list,inner_list,axis=0)

      del inner_list

    # print(self.pair_list)

  def construct_EEC(self, hist):
    for ipart1 in range(self.mult):

      for ipart2 in range(ipart1+1,self.mult):
        dist12 = self.pair_list[ipart1][ipart2]
        # print(' EEC combining particle i =', ipart1, 'and j =', ipart2, 'with distance', dist12)
        
        # print(' E(i) =', self.particle_list[ipart1].E(), 'and E(j) =', self.particle_list[ipart2].E(), 'with jet pt =', self.scale)
        eec_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()/math.pow(self.scale,2)
        hist.Fill(self.scale, dist12, eec_weight) # scale is jet pt

  def construct_E3C(self, hist):
    for ipart1 in range(self.mult):

      for ipart2 in range(ipart1+1,self.mult):
        dist12 = self.pair_list[ipart1][ipart2]

        for ipart3 in range(ipart2+1,self.mult):
          dist23 = self.pair_list[ipart2][ipart3]
          dist13 = self.pair_list[ipart1][ipart3]

          dist_list= [dist12, dist23, dist13]
          dist_list_sorted = sorted(dist_list)
          dist_max = dist_list_sorted[len(dist_list)-1]
          # print(' E3C combining particle', ipart1, ipart2, ipart3, 'with distance', dist12, dist23, dist13,'max',dist_max)
        
          e3c_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()*self.particle_list[ipart3].E()/math.pow(self.scale,3)
          hist.Fill(self.scale, dist_max, e3c_weight)

  def construct_E4C(self, hist):
    for ipart1 in range(self.mult):

      for ipart2 in range(ipart1+1,self.mult):
        dist12 = self.pair_list[ipart1][ipart2]

        for ipart3 in range(ipart2+1,self.mult):
          dist23 = self.pair_list[ipart2][ipart3]
          dist13 = self.pair_list[ipart1][ipart3]

          for ipart4 in range(ipart3+1,self.mult):
            dist34 = self.pair_list[ipart3][ipart4]
            dist24 = self.pair_list[ipart2][ipart4]
            dist14 = self.pair_list[ipart1][ipart4]

            dist_list= [dist12, dist23, dist13, dist34, dist24, dist14]
            dist_list_sorted = sorted(dist_list)
            dist_max = dist_list_sorted[len(dist_list)-1]
        
            e4c_weight = self.particle_list[ipart1].E()*self.particle_list[ipart2].E()*self.particle_list[ipart3].E()*self.particle_list[ipart4].E()/math.pow(self.scale,4)
            hist.Fill(self.scale, dist_max, e4c_weight)

################################################################
class ProcessMC_ENC(process_mc_base.ProcessMCBase):

  #---------------------------------------------------------------
  # Constructor
  #---------------------------------------------------------------
  def __init__(self, input_file='', config_file='', output_dir='', debug_level=0, **kwargs):
  
    # Initialize base class
    super(ProcessMC_ENC, self).__init__(input_file, config_file, output_dir, debug_level, **kwargs)
    
    self.observable = self.observable_list[0]

  #---------------------------------------------------------------
  # Initialize histograms
  #---------------------------------------------------------------
  def initialize_user_output_objects_R(self, jetR):

    for observable in self.observable_list:

      for trk_thrd in self.obs_settings[observable]:
        
        obs_label = self.utils.obs_label(trk_thrd, None) 

        for ipoint in range(2, 5):
          name = 'h_{}{}_JetPt_R{}_{}'.format(observable, ipoint, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          RL_bins = logbins(1E-4,1,50)
          h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('R_{L}')
          setattr(self, name, h)

          # Truth histograms
          name = 'h_{}{}_JetPt_Truth_R{}_{}'.format(observable, ipoint, jetR, obs_label)
          pt_bins = linbins(0,200,200)
          RL_bins = logbins(1E-4,1,50)
          h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
          h.GetXaxis().SetTitle('p_{T,ch jet}')
          h.GetYaxis().SetTitle('R_{L}')
          setattr(self, name, h)
          
          if self.thermal_model:
            for R_max in self.max_distance:
              name = 'h_{}{}_JetPt_R{}_{}_Rmax{}'.format(observable, ipoint, jetR, obs_label, R_max)
              pt_bins = linbins(0,200,200)
              RL_bins = logbins(1E-4,1,50)
              h = ROOT.TH2D(name, name, 200, pt_bins, 50, RL_bins)
              h.GetXaxis().SetTitle('p_{T,ch jet}')
              h.GetYaxis().SetTitle('R_{L}')
              setattr(self, name, h)

        # Residuals and responses
        for trk_thrd in self.obs_settings[observable]:
        
          for ipoint in range(2, 5):
            if not self.is_pp:
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

  #---------------------------------------------------------------
  # This function is called once for each jet subconfiguration
  # Fill 2D histogram of (pt, obs)
  #---------------------------------------------------------------
  def fill_observable_histograms(self, hname, jet, jet_groomed_lund, jetR, obs_setting,
                                 grooming_setting, obs_label, jet_pt_ungroomed):
    
    constituents = fj.sorted_by_pt(jet.constituents())
    c_select = []
    trk_thrd = obs_setting
    # print('jet R',jetR,'observable',observable,'thrd',trk_thrd)
    # print('const size',len(constituents))
    for c in constituents:
      if c.pt() < trk_thrd:
        break
      c_select.append(c) # NB: use the break statement since constituents are already sorted

      # print('jet',jet.perp(),'pt',c.pt())
    
    new_corr = CorrelatorBuilder(c_select,jet.perp())
    new_corr.make_pairs()
    for ipoint in range(2, 5):
      if (ipoint==2):
          new_corr.construct_EEC(getattr(self, hname.format(self.observable + str(ipoint),obs_label)))
      if (ipoint==3):
          new_corr.construct_E3C(getattr(self, hname.format(self.observable + str(ipoint),obs_label)))
      if (ipoint==4):
          new_corr.construct_E4C(getattr(self, hname.format(self.observable + str(ipoint),obs_label)))  
        
  #---------------------------------------------------------------
  # Fill matched jet histograms
  #---------------------------------------------------------------
  def fill_matched_jet_histograms(self, jet_det, jet_det_groomed_lund, jet_truth,
                                  jet_truth_groomed_lund, jet_pp_det, jetR,
                                  obs_setting, grooming_setting, obs_label,
                                  jet_pt_det_ungroomed, jet_pt_truth_ungroomed, R_max, suffix, **kwargs):
    pass
       
    # if (jetR - obs_setting) < 1e-3:
    #   return
      
    # # If jetscape, we will need to correct substructure observable for holes (pt is corrected in base class)
    # if self.jetscape:
    #     holes_in_det_jet = kwargs['holes_in_det_jet']
    #     holes_in_truth_jet = kwargs['holes_in_truth_jet']
       
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
