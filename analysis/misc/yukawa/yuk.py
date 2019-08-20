import os
import numpy as np
import ROOT as r

r.gROOT.SetBatch(1)

def get_kt(sigma, A,B,F, other=False):
    # A: gza
    # B: interference
    # F: higgs
    fact = 1.0 if not other else -1.0
    return np.sqrt(fact*np.sqrt(B*B-4.*A*F+4.*F*sigma)/F  - B/F)/np.sqrt(2.0)

def bifurcated(cent,low,high, nsig=0.):
    # given a central value, a 1-sigma low and 1-sigma high
    # produce a number that differs from central value by
    # nsigma (different for up and down)
    if nsig > 0.:
        return cent + (high-cent)*nsig
    else:
        return cent + (cent-low)*nsig

# use actual likelihood shape from combine
# otherwise, we use a bifurcated shape
use_actual_obspdf = True
if use_actual_obspdf:
    f1 = r.TFile("higgsCombineTest.MultiDimFit.mH120.root")
    tree = f1.Get("limit")
    # tree.Draw("TMath::Exp(-2*deltaNLL):9.2*r>>hobs(25,0,50,100,0,1)","","profx")
    tree.Draw("TMath::Exp(-2*deltaNLL):9.2*r>>hobs(35,0,70,100,0,1)","","profx")
    hobs = r.gDirectory.Get("hobs")

# we measured 16.9 +13.8 -11.4 fb
sigma_obs = 1.8329030 * 9.2
# 95% CL observed UL
sigma_obs_ul = 41.7

gza_8tev = 1.344
int_8tev = -0.224
higgs_8tev = 0.171

gza_13tev = 9.997
int_13tev = -1.547
higgs_13tev = 1.168

# scale variations for 13 tev numbers
gza_13tev_up = 14.104
int_13tev_up = -2.152
higgs_13tev_up = 1.625
gza_13tev_down = 6.378
int_13tev_down = -0.999
higgs_13tev_down = 0.7655

# NLO/LO kfactor from 14TeV
kfactor = 1.27


print "The paper should get kt<{0:.2f} for 13TeV".format(get_kt(10.2*9.2, gza_13tev, int_13tev, higgs_13tev))
print "The paper should get kt<{0:.2f} for 8TeV without kfactor".format(get_kt(23.,gza_8tev,int_8tev,higgs_8tev))
print "The paper should get kt<{0:.2f} for 8TeV with kfactor".format(get_kt(23.,gza_8tev*kfactor,int_8tev*kfactor,higgs_8tev*kfactor))

print "I get kt={0:.2f} using obs xsec of {1:.1f} for 13TeV with kfactor".format(get_kt(sigma_obs, gza_13tev*kfactor, int_13tev*kfactor, higgs_13tev*kfactor), sigma_obs)
print "I get kt<{0:.2f} using obs UL of {1:.1f} for 13TeV with kfactor".format(get_kt(sigma_obs_ul, gza_13tev*kfactor, int_13tev*kfactor, higgs_13tev*kfactor), sigma_obs_ul)

# reset seed for reproducibility
np.random.seed(1)
# N = 10000
N = 10000

# first make N gaussian random numbers to guarantee correlation amongst the 
# 3 terms below
# gausrands = np.random.normal(size=N)


toy_gza = []
toy_h = []
toy_int = []
toy_obs = []
gaus = []
kts = []
weights = []

do_correlated = True
do_uniform = False

for i in range(N):
    if not do_uniform:
        gausrand = np.random.normal()
    else:
        gausrand = np.random.uniform(-1.,1.)
    gausrand2 = np.random.normal()

    if do_correlated:
        gza = bifurcated(gza_13tev*kfactor,gza_13tev_down*kfactor,gza_13tev_up*kfactor,gausrand)
        higgs = bifurcated(higgs_13tev*kfactor,higgs_13tev_down*kfactor,higgs_13tev_up*kfactor,gausrand)
        interference = bifurcated(int_13tev*kfactor,int_13tev_down*kfactor,int_13tev_up*kfactor,gausrand)
    else:
        if not do_uniform:
            rands = [np.random.normal() for _ in range(3)]
        else:
            rands = [np.random.uniform(-1.,1.) for _ in range(3)]
        gza = bifurcated(gza_13tev*kfactor,gza_13tev_down*kfactor,gza_13tev_up*kfactor,rands[0])
        higgs = bifurcated(higgs_13tev*kfactor,higgs_13tev_down*kfactor,higgs_13tev_up*kfactor,rands[1])
        interference = bifurcated(int_13tev*kfactor,int_13tev_down*kfactor,int_13tev_up*kfactor,rands[2])

    if use_actual_obspdf:
        obs = hobs.GetRandom()
    else:
        obs = bifurcated(16.9,16.9-11.4,16.9+13.8,gausrand2)

    gaus.append(gausrand)
    toy_gza.append(gza)
    toy_h.append(higgs)
    toy_int.append(interference)
    toy_obs.append(obs)

    if obs < 0. or gza < 0. or higgs < 0. or interference > 0.:
        # print "bad:",gza, higgs, interference, obs
        continue

    soln1 = get_kt(obs, gza, interference, higgs)
    soln2 = get_kt(obs, gza, interference, higgs, other=True)


    good1 = not np.isnan(soln1)
    good2 = not np.isnan(soln2)

    weight = 1.0
    # if good1 and good2:
    #     weight = 1.0

    if good1:
        kts.append(soln1)
        weights.append(weight)
    if good2:
        kts.append(soln2)
        weights.append(weight)
    # kt = soln1
    # print soln1, soln2

    # kts.append(kt)
    # weights.append(1)


toy_gza = np.array(toy_gza)
toy_h = np.array(toy_h)
toy_int = np.array(toy_int)
toy_obs = np.array(toy_obs)
gaus = np.array(gaus)

# sometimes we get NaN, so remove those
# kts = get_kt(toy_obs, toy_gza, toy_int, toy_h)
# kts = kts[~np.isnan(kts)]

c1 = r.TCanvas()

for name,arr in [
        ("gza",toy_gza),
        ("higgs",toy_h),
        ("int",toy_int),
        ("obs",toy_obs),
        ("gaus",gaus),
        ]:
    htemp = r.TH1F(name,name,100,min(-1,arr.min()-1),max(1,arr.max()+1))
    htemp.Sumw2()
    for val in arr:
        htemp.Fill(val)
    htemp.Draw("pe")
    c1.SaveAs("plots/{0}.pdf".format(name))

h1 = r.TH1F("h1","#kappa_{t}",90,-1.0,3.5)
h1.Sumw2()
for kt,weight in zip(kts,weights):
    h1.Fill(kt,weight)
h1.Scale(1.0/h1.Integral())
h1.SetMinimum(0.)
h1.Draw("pe")

# compute x value for which we collect 95% of the probability
cumulative = 0.
for ib in range(1,h1.GetNbinsX()):
    cumulative += h1.GetBinContent(ib)
    if cumulative > 0.95:
        break
ul95 = h1.GetXaxis().GetBinLowEdge(ib)

latex = r.TLatex()
latex.SetTextSize(0.035)
latex.SetTextColor(r.kBlue)
latex.SetTextFont(42)
latex.DrawLatexNDC(0.65,0.60, "95%% UL: %.2f" % ul95)

line = r.TArrow(ul95,h1.GetMaximum(),ul95,0.0, 0.02, "|>")
line.SetLineWidth(2)
line.SetLineColor(r.kBlue)
line.SetFillColor(r.kBlue)
line.SetAngle(40)
line.Draw()

os.system("mkdir -p plots")
c1.SaveAs("plots/kt.pdf")
os.system("ic plots/kt.pdf")

os.system("niceplots plots yukawa_Aug12")

