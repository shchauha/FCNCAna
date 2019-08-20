import ROOT as r

def get_first_hist(f1):
    return f1.Get(f1.GetListOfKeys()[0].GetName())

def dump_bins(h2, name, stupid=False, nofabseta=False, default1=True, do_err=False):
    buff = "float %s(float pt, float eta) {\n" % name
    for ix in range(1,h2.GetNbinsX()+1):
        for iy in range(1,h2.GetNbinsY()+1):
            if stupid:
                val = h2.GetBinContent(ix,iy)
                if do_err:
                    val = h2.GetBinError(ix,iy)
                if iy != h2.GetNbinsY():
                    buff += "  if (pt >= %.0f && pt < %.0f && eta >= %.3f && eta < %.3f) return %.4f;\n" \
                            % (h2.GetYaxis().GetBinLowEdge(iy), h2.GetYaxis().GetBinUpEdge(iy), h2.GetXaxis().GetBinLowEdge(ix), h2.GetXaxis().GetBinUpEdge(ix), val)
                else:
                    buff += "  if (pt >= %.0f  && eta >= %.3f && eta < %.3f) return %.4f;\n" \
                            % (h2.GetYaxis().GetBinLowEdge(iy), h2.GetXaxis().GetBinLowEdge(ix), h2.GetXaxis().GetBinUpEdge(ix), val)
            else:
                val = h2.GetBinContent(ix,iy)
                if do_err:
                    val = h2.GetBinError(ix,iy)
                if ix != h2.GetNbinsX():
                    buff += "  if (pt >= %.0f && pt < %.0f && fabs(eta) >= %.3f && fabs(eta) < %.3f) return %.4f;\n" \
                            % (h2.GetXaxis().GetBinLowEdge(ix), h2.GetXaxis().GetBinUpEdge(ix), h2.GetYaxis().GetBinLowEdge(iy), h2.GetYaxis().GetBinUpEdge(iy), val)
                else:
                    buff += "  if (pt >= %.0f  && fabs(eta) >= %.3f && fabs(eta) < %.3f) return %.4f;\n" \
                            % (h2.GetXaxis().GetBinLowEdge(ix), h2.GetYaxis().GetBinLowEdge(iy), h2.GetYaxis().GetBinUpEdge(iy), val)
    if not default1:
        buff += "  return 0.;\n"
    else:
        buff += "  return 1.;\n"
    buff += "}\n"
    if nofabseta:
        buff = buff.replace("fabs(eta)", "eta")
    print buff

# def dump_bins_fastsim(h3, name):
#     print "float %s(float pt, float eta, int nvtx) {" % name
#     print "  if (pt>=200) pt = 199.9;"
#     if "electron" in name:
#         print "  if (fabs(eta)>=2.5) eta = 2.4;"
#     else:
#         print "  if (fabs(eta)>=2.4) eta = 2.3;"
#     print "  if (nvtx>=40) nvtx = 39;"
#     for ix in range(1,h3.GetNbinsX()+1):
#         for iy in range(1,h3.GetNbinsY()+1):
#             for iz in range(1,h3.GetNbinsZ()+1):
#                 val = h3.GetBinContent(ix,iy,iz)
#                 ptlow = h3.GetXaxis().GetBinLowEdge(ix)
#                 pthigh = h3.GetXaxis().GetBinUpEdge(ix)
#                 etalow = h3.GetYaxis().GetBinLowEdge(iy)
#                 etahigh = h3.GetYaxis().GetBinUpEdge(iy)
#                 nvtxlow = h3.GetZaxis().GetBinLowEdge(iz)
#                 nvtxhigh = h3.GetZaxis().GetBinUpEdge(iz)
#                 print "  if (pt >= %.0f && pt < %.0f && fabs(eta) >= %.3f && fabs(eta) < %.3f && nvtx >= %i && nvtx < %i) return %.4f;" \
#                         % (ptlow, pthigh, etalow, etahigh, nvtxlow, nvtxhigh, val)
#     print "  return 0.;"
#     print "}"

# el_sfs.root


if True:
    # FULLSIM ELECTRONS
    electron_file = r.TFile("rootfiles_moriond/scaleFactors.root") # from http://tomc.web.cern.ch/tomc/tagAndProbe/20160726/output/scaleFactors.root

    # tightid2d3d = electron_file.Get("GsfElectronToMVATightTightIP2DSIP3D4")
    tightid2d3d = electron_file.Get("GsfElectronToMVATightIDEmuTightIP2DSIP3D4")
    multi = electron_file.Get("MVATightElectronToMultiIsoT")
    convhit0chg = electron_file.Get("MVATightConvIHit0ElectronToCharge")
    multiemu = electron_file.Get("MVATightElectronToMultiIsoTISOEmu")

    el_noiso = tightid2d3d.Clone("el_noiso")
    el_noiso.Multiply(multi)
    el_noiso.Multiply(convhit0chg)

    el_iso = tightid2d3d.Clone("el_iso")
    el_iso.Multiply(multiemu)
    el_iso.Multiply(convhit0chg)

    f_gsf = r.TFile("rootfiles_moriond/egammaEffi.txt_EGM2D.root")
    h_gsf = f_gsf.Get("EGamma_SF2D")

    # FULLSIM MUONS
    muon_file_dxy = r.TFile("rootfiles_moriond/mu_dxy0p05dz0p1.root") # http://mmarionn.web.cern.ch/mmarionn/SUSY/MuonSFMoriond17/TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root
    muon_file_mediumid = r.TFile("rootfiles_moriond/mu_mediumID.root") # http://mmarionn.web.cern.ch/mmarionn/SUSY/MuonSFMoriond17/TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root
    muon_file_sip3d = r.TFile("rootfiles_moriond/mu_sip3dlt4.root") # http://mmarionn.web.cern.ch/mmarionn/SUSY/MuonSFMoriond17/TnP_NUM_TightIP3D_DENOM_MediumID_VAR_map_pt_eta.root
    muon_file_multiiso = r.TFile("rootfiles_moriond/mu_multiiso.root") # http://mmarionn.web.cern.ch/mmarionn/SUSY/MuonSFMoriond17/TnP_NUM_MultiIsoMedium_DENOM_MediumID_VAR_map_pt_eta.root
    muon_file_dpt = r.TFile("rootfiles_moriond/mu_dpt0p2.root") # http://mmarionn.web.cern.ch/mmarionn/SUSY/MuonSFMoriond17/TnP_NUM_PtError_DENOM_LooseID_VAR_map_pt_eta.root

    dxy0p05dz0p1 = get_first_hist(muon_file_dxy)
    mediumID = get_first_hist(muon_file_mediumid)
    sip3dlt4 = get_first_hist(muon_file_sip3d)
    multiiso = get_first_hist(muon_file_multiiso) 
    dpt = get_first_hist(muon_file_dpt) 

    # muon_tracking = r.TFile("rootfiles_moriond//Tracking_EfficienciesAndSF_BCDEFGH.root")
    # mutrack = muon_tracking.Get("ratio_eff_eta3_dr030e030_corr")

    mu = dxy0p05dz0p1.Clone("mu")
    mu.Multiply(mediumID)
    mu.Multiply(sip3dlt4)
    mu.Multiply(multiiso)
    mu.Multiply(dpt) 

    # FASTSIM ELECTRONS

    el_fs_file_tightid2d3d = r.TFile("rootfiles_moriond/fs/sf_el_tight2d3dIDEmu.root")
    el_fs_file_multi = r.TFile("rootfiles_moriond/fs/sf_el_multi.root")
    el_fs_file_convhit0chg = r.TFile("rootfiles_moriond/fs/sf_el_inhit_eq0_chargeConst.root")
    el_fs_file_multiemu = r.TFile("rootfiles_moriond/fs/sf_el_multiISO.root")


    el_fs_tightid2d3d = get_first_hist(el_fs_file_tightid2d3d)
    el_fs_multi = get_first_hist(el_fs_file_multi)
    el_fs_convhit0chg = get_first_hist(el_fs_file_convhit0chg)
    el_fs_multiemu = get_first_hist(el_fs_file_multiemu)


    el_fs_noiso = el_fs_tightid2d3d.Clone("el_fs_noiso")
    el_fs_noiso.Multiply(el_fs_multi)
    el_fs_noiso.Multiply(el_fs_convhit0chg)

    el_fs_iso = el_fs_tightid2d3d.Clone("el_fs_iso")
    el_fs_iso.Multiply(el_fs_multiemu)
    el_fs_iso.Multiply(el_fs_convhit0chg)


    # FASTSIM MUONS
    muon_fs_file_dxy = r.TFile("rootfiles_moriond/fs/sf_mu_mediumID_tightIP2D.root")
    muon_fs_file_mediumid = r.TFile("rootfiles_moriond/fs/sf_mu_mediumID.root")
    muon_fs_file_sip3d = r.TFile("rootfiles_moriond/fs/sf_mu_mediumID_tightIP3D.root")
    muon_fs_file_multiiso = r.TFile("rootfiles_moriond/fs/sf_mu_mediumID_multiM.root") 
    muon_fs_file_dpt = r.TFile("rootfiles_moriond/fs/sf_mu_deltaPt.root") 

    dxy0p05dz0p1_fs = get_first_hist(muon_fs_file_dxy)
    mediumID_fs = get_first_hist(muon_fs_file_mediumid)
    sip3dlt4_fs = get_first_hist(muon_fs_file_sip3d)
    multiiso_fs = get_first_hist(muon_fs_file_multiiso) 
    dpt_fs = get_first_hist(muon_fs_file_dpt) 

    mu_fastsim = dxy0p05dz0p1_fs.Clone("mu")
    mu_fastsim.Multiply(mediumID_fs)
    mu_fastsim.Multiply(sip3dlt4_fs)
    mu_fastsim.Multiply(multiiso_fs)
    mu_fastsim.Multiply(dpt_fs) 

    ### DO THE DUMPING

    # dump_bins(mutrack, name="trackingSF", nofabseta=True)

    dump_bins(el_noiso, name="electronScaleFactorHighHT", nofabseta=True)
    dump_bins(el_iso, name="electronScaleFactorLowHT", nofabseta=True)
    dump_bins(mu, name="muonScaleFactor")
    dump_bins(h_gsf,"electronGSF", stupid=True)

    dump_bins(el_noiso, name="electronScaleFactorHighHT", nofabseta=True, do_err=True)
    dump_bins(el_iso, name="electronScaleFactorLowHT", nofabseta=True, do_err=True)
    dump_bins(mu, name="muonScaleFactor", do_err=True)
    dump_bins(h_gsf,"electronGSF", stupid=True, do_err=True)

    # dump_bins(el_fs_iso, name="electronScaleFactorFastSimHighHT")
    # dump_bins(el_fs_noiso, name="electronScaleFactorFastSimLowHT")
    # dump_bins(mu_fastsim, name="muonScaleFactorFastSim")
