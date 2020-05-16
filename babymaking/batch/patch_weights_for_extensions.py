import numpy as np
import os
from array import array
import ROOT as r
from tqdm import tqdm
import sys

# /nfs-7/userdata/namin/tupler_babies/merged/FT/v3.21_fix2017/output/year_2017/TTdilep0jet.root
def patch_files(oldfnames):

    dolds = {}
    for fname in oldfnames:
        dolds[fname] = {}
        f = r.TFile(fname)
        t = f.Get("t")
        hcounts = f.Get("counts")
        dolds[fname]["file"] = f
        dolds[fname]["tree"] = t
        # dolds[fname]["weighthist"] = f.Get("weight")
        dolds[fname]["otherstufftowrite"] = [f.Get(k.GetName()) for k in f.GetListOfKeys() if k.GetName() != "t"]
        dolds[fname]["counts"] = hcounts.GetEntries()
        dolds[fname]["old_scale1fb"] = array("f",[0.])
        t.SetBranchAddress("scale1fb",dolds[fname]["old_scale1fb"])
        t.GetEntry(0)

    print "Input event count and scale1fbs:"
    tot_counts = sum(v["counts"] for v in dolds.values())
    for fname in oldfnames:
        scale1fb = float(dolds[fname]["old_scale1fb"][0])
        counts = dolds[fname]["counts"]
        toscale = counts/tot_counts 
        dolds[fname]["toscale"] = toscale
        print "Scaling {} (N={:.0f}, scale1fb={:.4f}) by {:.4f} to get scale1fb={:.4f}".format(fname, counts, scale1fb, toscale, scale1fb*toscale)

    old_scale1fbs = np.array([float(v["old_scale1fb"][0]) for v in dolds.values()])
    same_oldscale1fb = np.all(np.abs(old_scale1fbs-old_scale1fbs.mean()) < 0.01)
    if same_oldscale1fb:
        print "Hmm, the old scale1fbs are already pretty much the same. You probably already scaled them. Not continuing."
        sys.exit()

    # sys.exit() # FIXME

    for oldfname,dold in dolds.items():
        basename = os.path.basename(oldfname)
        basedir = os.path.dirname(oldfname)

        tmp_basename = "tmp_{}".format(basename)
        fnew = r.TFile(tmp_basename, "recreate")
        oldtree = dold["tree"]
        newtree = oldtree.CloneTree(0)

        [x.Write() for x in dolds[fname]["otherstufftowrite"]] # write these things into the new TFile

        print "Looping to fill {}".format(tmp_basename)
        scale = dold["toscale"]
        old_scale1fb = dold["old_scale1fb"]
        for i in tqdm(xrange(oldtree.GetEntries()),total=oldtree.GetEntries()):
            oldtree.GetEntry(i)
            old_scale1fb[0] = old_scale1fb[0]*scale
            newtree.Fill()
            # if (i > 30000): break # FIXME
        newtree.AutoSave()
        del oldtree
        del newtree
        fnew.Close()
        print "Copying {} over {}".format(tmp_basename,oldfname)
        # FIXME actually do the mv
        os.system("mv {} {}".format(tmp_basename,oldfname))
        print "Done"


if __name__ == "__main__":

    patch_files([
            "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2017/TTdilep0jet.root",
            "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2017/TTdilep0jetext1.root",
            ])
    patch_files([
            "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2017/TTdilep1jet.root",
            "/nfs-7/userdata/namin/tupler_babies/merged/FT/v3.24/output/year_2017/TTdilep1jetext1.root",
            ])
