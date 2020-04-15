import numpy as np
from pytable import Table
from analysis.limits.errors import E

def write_table(data, bgs, outname=None, signal=None, extra_hists=[],precision=2,sep = u"\u00B1".encode("utf-8"), binedge_fmt="{}-{}", fix_negative=True, binlabels=[], show_errors=True, cell_callback=None):
    tab = Table()
    sumbgs = sum(bgs)
    databg = data/sumbgs
    if signal is not None:
        procs = bgs+[sumbgs,data,databg,signal]
        cnames = [bg.get_attr("label") for bg in bgs] + ["Total bkg","Data", "Data/bkg","tttt"]
    else:
        procs = bgs+[sumbgs,data,databg]
        cnames = [bg.get_attr("label") for bg in bgs] + ["Total bkg","Data", "Data/bkg"]
    for eh in extra_hists:
        procs.append(eh)
        cnames.append(eh.get_attr("label"))
    tab.set_column_names(["bin"]+cnames)
    if outname:
        sep = "+-"
    binpairs = zip(data.edges[:-1],data.edges[1:])
    #tab.set_theme_basic()
    tab.set_theme_latex()
    for ibin,binrow in enumerate(binpairs):
        row = [("[%s]"%binedge_fmt).format(binrow[0],binrow[1])]
        if ibin < len(binlabels):
            row = [binlabels[ibin]]
        for iproc,proc in enumerate(procs):
            if fix_negative:
                cent = max(proc.counts[ibin],0.)
            else:
                cent = proc.counts[ibin]
            err = proc.errors[ibin]
            if show_errors:
                tmp = ("{0:5.%if} {1}{2:%i.%if}" % (precision,precision+3,precision)).format(cent,sep,err)
            else:
                tmp = ("{0:5.%if}" % (precision)).format(cent)
            if cell_callback: tmp = cell_callback(tmp)
            row.append(tmp)
        tab.add_row(row)
    tab.add_line()

    row = ["total"]
    for iproc,proc in enumerate(procs):
        if iproc == len(procs)-(1+(signal is not None)+len(extra_hists)):
            totbg = E(sum(sumbgs.counts), np.sum(sumbgs.errors**2.)**0.5)
            totdata = E(sum(data.counts))
            ratio = totdata/totbg
            cent, err = ratio[0], ratio[1]
            precision = max(precision, 2) if precision != 0 else 0
        else:
            cent = sum(proc.counts)
            err = np.sum(proc.errors**2.)**0.5
        if show_errors:
            tmp = ("{0:5.%if} {1}{2:%i.%if}" % (precision,precision+3,precision)).format(cent,sep,err)
        else:
            tmp = ("{0:5.%if}" % (precision)).format(cent)
            if cell_callback: tmp = cell_callback(tmp)
        row.append(tmp)
    tab.add_row(row)

    if outname:
        with open(outname,"w") as fhout:
            # towrite = "".join(tab.get_table_string(show_row_separators=False,show_alternating=False))
            towrite = "".join(tab.get_table_strings(show_row_separators=False,show_alternating=False))
            fhout.write(towrite)
            parts = towrite.split("\n")
            header = parts[:3]
            binparts = parts[3:-4]
            total = parts[-4:-1]
            table_info = { "header":"<br>".join(header), "bins":binparts, "total":"<br>".join(total) }
            return table_info
    return tab

def get_bin_info(bin_coords, table_info):
    info = { "xedges": [], "yedges": [], "table": table_info }
    for ipair,pair in enumerate(zip(bin_coords[::2],bin_coords[1::2])):
        xlow,xhigh = round(pair[0][0],5), round(pair[1][0],5)
        ylow,yhigh = round(pair[0][1],5), round(pair[1][1],5)
        info["xedges"].append([xlow,xhigh])
        info["yedges"].append([ylow,yhigh])
    return info

