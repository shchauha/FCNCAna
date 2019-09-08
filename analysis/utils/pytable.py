#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
import math
import os
import commands

def round_sig(x, sig=2):
    if x < 0.001: return x
    return round(x, sig-int(math.floor(math.log10(x)))-1)

def visible_len(val):
    try: val = val.decode("utf-8")
    except: pass
    lenval = len(val)
    if "\033" in val:
        # lenval -= 12
        lenval -= sum(map(len,re.findall(r"\033\[[0-9\;]+m",val)))
    return lenval


class Table():

    def __init__(self):
        self.matrix = []
        self.colnames = []
        self.colsizes = []
        self.hlines = []
        self.rowcolors = {}
        self.extra_padding = 1
        self.d_style = {}
        self.theme = ""
        self.set_theme_fancy()
        self.use_color = True

    def shorten_string(self, val, length):
        return val[:length//2-1] + "..." + val[-length//2+2:]
    def fmt_string(self, val, length, fill_char=" ", justify="c", bold=False, offcolor=False, color=None):
        ret = ""
        val = str(val)
        lenval = visible_len(val)
            
        if lenval > length: val = self.shorten_string(val, length)
        if justify == "l": 
            nr = (length-lenval-1)
            ret = " " + val + fill_char*nr
            # ret = " "+val.ljust(length-1, fill_char)
        elif justify == "r": ret = val.rjust(length, fill_char)
        elif justify == "c":
            nl = (length-lenval)//2
            nr = (length-lenval)-(length-lenval)//2
            ret = fill_char*nl + val + fill_char*nr
        if bold and self.use_color:
            ret = '\033[1m' + ret + '\033[0m'
        if offcolor and self.use_color:
            ret = '\033[2m' + ret + '\033[0m'
        if self.use_color:
            if color == "green":
                ret = '\033[00;32m' + ret + '\033[0m'
            if color == "blue":
                ret = '\033[00;34m' + ret + '\033[0m'
            if color == "red":
                ret = '\033[00;31m' + ret + '\033[0m'
            if color == "lightblue":
                ret = '\033[38;5;117m' + ret + '\033[0m'
        return ret

    def set_theme_fancy(self):
        self.theme = "fancy"
        self.d_style["INNER_HORIZONTAL"] = '\033(0\x71\033(B'
        self.d_style["INNER_INTERSECT"] = '\033(0\x6e\033(B'
        self.d_style["INNER_VERTICAL"] = '\033(0\x78\033(B'
        self.d_style["OUTER_LEFT_INTERSECT"] = '\033(0\x74\033(B'
        self.d_style["OUTER_LEFT_VERTICAL"] = self.d_style["INNER_VERTICAL"]
        self.d_style["OUTER_RIGHT_INTERSECT"] = '\033(0\x75\033(B'
        self.d_style["OUTER_RIGHT_VERTICAL"] = self.d_style["INNER_VERTICAL"]
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = self.d_style["INNER_HORIZONTAL"]
        self.d_style["OUTER_BOTTOM_INTERSECT"] = '\033(0\x76\033(B'
        self.d_style["OUTER_BOTTOM_LEFT"] = '\033(0\x6d\033(B'
        self.d_style["OUTER_BOTTOM_RIGHT"] = '\033(0\x6a\033(B'
        self.d_style["OUTER_TOP_HORIZONTAL"] = self.d_style["INNER_HORIZONTAL"]
        self.d_style["OUTER_TOP_INTERSECT"] = '\033(0\x77\033(B'
        self.d_style["OUTER_TOP_LEFT"] = '\033(0\x6c\033(B'
        self.d_style["OUTER_TOP_RIGHT"] = '\033(0\x6b\033(B'

    def set_theme_basic(self):
        self.theme = "basic"
        self.use_color = False

        self.d_style["INNER_HORIZONTAL"] = '-'
        self.d_style["INNER_INTERSECT"] = '+'
        self.d_style["INNER_VERTICAL"] = '|'
        self.d_style["OUTER_LEFT_INTERSECT"] = '|'
        self.d_style["OUTER_LEFT_VERTICAL"] = '|'
        self.d_style["OUTER_RIGHT_INTERSECT"] = '+'
        self.d_style["OUTER_RIGHT_VERTICAL"] = '|'
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = '-'
        self.d_style["OUTER_BOTTOM_INTERSECT"] = '+'
        self.d_style["OUTER_BOTTOM_LEFT"] = '+'
        self.d_style["OUTER_BOTTOM_RIGHT"] = '+'
        self.d_style["OUTER_TOP_HORIZONTAL"] = '-'
        self.d_style["OUTER_TOP_INTERSECT"] = '+'
        self.d_style["OUTER_TOP_LEFT"] = '+'
        self.d_style["OUTER_TOP_RIGHT"] = '+'

    def set_theme_csv(self, delimiter=","):
        self.theme = "csv"
        self.use_color = False

        self.d_style["INNER_HORIZONTAL"] = ''
        self.d_style["INNER_INTERSECT"] = delimiter
        self.d_style["INNER_VERTICAL"] = delimiter
        self.d_style["OUTER_LEFT_INTERSECT"] = ''
        self.d_style["OUTER_LEFT_VERTICAL"] = ''
        self.d_style["OUTER_RIGHT_INTERSECT"] = ''
        self.d_style["OUTER_RIGHT_VERTICAL"] = ''
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = ''
        self.d_style["OUTER_BOTTOM_INTERSECT"] = ''
        self.d_style["OUTER_BOTTOM_LEFT"] = ''
        self.d_style["OUTER_BOTTOM_RIGHT"] = ''
        self.d_style["OUTER_TOP_HORIZONTAL"] = ''
        self.d_style["OUTER_TOP_INTERSECT"] = ''
        self.d_style["OUTER_TOP_LEFT"] = ''
        self.d_style["OUTER_TOP_RIGHT"] = ''

    def set_theme_latex(self):
        self.theme = "latex"
        self.use_color = False
        self.d_style["INNER_HORIZONTAL"] = ''
        self.d_style["INNER_INTERSECT"] = ''
        self.d_style["INNER_VERTICAL"] = ' & '
        self.d_style["OUTER_LEFT_INTERSECT"] = ''
        self.d_style["OUTER_LEFT_VERTICAL"] = ''
        self.d_style["OUTER_RIGHT_INTERSECT"] = r'\hline'
        self.d_style["OUTER_RIGHT_VERTICAL"] = r'\\ \hline'
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = ''
        self.d_style["OUTER_BOTTOM_INTERSECT"] = ''
        self.d_style["OUTER_BOTTOM_LEFT"] = ''
        self.d_style["OUTER_BOTTOM_RIGHT"] = ''
        self.d_style["OUTER_TOP_HORIZONTAL"] = ''
        self.d_style["OUTER_TOP_INTERSECT"] = ''
        self.d_style["OUTER_TOP_LEFT"] = ''
        self.d_style["OUTER_TOP_RIGHT"] = r'\hline'

    def compile_latex(self, tstr, name="table", temp_folder_name="tmp"):
        template = r"""
        \documentclass{article}
        \usepackage[table]{xcolor}
        \usepackage{multirow}
        \usepackage{graphicx}
        \usepackage{slashed}
        \newcommand{\met}{\slashed{E}_\mathrm{T}}
        \newcommand{\mt}{m_\mathrm{T}}
        \newcommand{\pt}{p_\mathrm{T}}
        \newcommand{\mtmin}{m_{T}^\mathrm{min}}
        \newcommand{\Ht}{H_\mathrm{T}}
        \renewcommand{\arraystretch}{1.2}
        \usepackage{helvet} \renewcommand{\familydefault}{\sfdefault}
        \definecolor{white}{rgb}{1,1,1}
        \definecolor{black}{rgb}{0,0,0}
        \definecolor{red}{rgb}{0.96,0.10,0.12}
        \definecolor{green}{rgb}{0.25,0.90,0.10}
        \definecolor{blue}{rgb}{0.14,0.10,1.0}
        \definecolor{orange}{rgb}{0.98,0.85,0.22}
        \definecolor{yellow}{rgb}{1.0,0.98,0.28}
        \definecolor{lightred}{rgb}{0.96,0.60,0.52}
        \definecolor{lightgreen}{rgb}{0.75,0.90,0.60}
        \definecolor{lightblue}{rgb}{0.64,0.80,1.0}
        \definecolor{lightorange}{rgb}{0.98,0.85,0.62}
        \definecolor{lightyellow}{rgb}{1.0,0.98,0.68}
        \begin{document}
        \pagenumbering{gobble}
        \begin{table}[h]
        \centering
        \resizebox{\textwidth}{!}{
            \begin{tabular}{%s}
            \rowcolor{white} %% first row has bg so that it's easily selectable in keynote
            %s
            \end{tabular}
        }
        \end{table}
        \end{document}
        """
        colstr = "|"+"c|"*len(self.get_column_names())
        os.system("mkdir -p {}".format(temp_folder_name))
        with open("{tmpname}/{name}.tex".format(tmpname=temp_folder_name,name=name),"w") as fh:
            fh.write(template%(colstr,tstr))
        status,out = commands.getstatusoutput("pdflatex -interaction=nonstopmode -output-directory={tmpname} {tmpname}/{name}.tex".format(tmpname=temp_folder_name,name=name))
        if " Error" in out:
            print "ERROR: Tried to compile, but failed. Last few lines of printout below."
            print "_"*40
            print "\n".join(out.split("\n")[-30:])
        else:
            status,out = commands.getstatusoutput("pdfcrop --margins 4 {tmpname}/{name}.pdf {name}.pdf >& /dev/null".format(tmpname=temp_folder_name,name=name))
            return "{name}.pdf".format(name=name)
        return ""

    def print_pdf(self, **kwargs):
        if self.theme != "latex":
            raise Exception("Can't save themes other than latex to pdf")
        tstr = self.get_table(**kwargs)
        fname = self.compile_latex(tstr)
        os.system("which ic >& /dev/null && ic {fname} || echo Made {fname}".format(fname=fname))

    def set_column_names(self, cnames):
        self.colnames = cnames
        self.update()

    def get_column_names(self):
        return self.colnames

    def add_row(self, row, color=None):
        self.matrix.append(row)
        if color:
            self.rowcolors[len(self.matrix)] = color

    def add_column(self, colname, values):
        # if no matrix to begin with, just add the column, otherwise append to rows
        if len(self.matrix) == 0:
            for val in values:
                self.matrix.append([val])
        else:
            for irow in range(len(self.matrix)):
                if irow < len(values): val = values[irow]
                else: val = ""
                self.matrix[irow].append(val)
        self.colnames.append(colname)


    def add_line(self):
        # draw hlines by making list of the row 
        # indices, which we will check when drawing
        # the matrix
        self.hlines.append(len(self.matrix))

    def update(self):
        if not self.colnames:
            if self.matrix:
                self.colnames = range(1,len(self.matrix[0])+1)
        if self.matrix:
            for ic, cname in enumerate(self.colnames):
                self.colsizes.append( max(
                    max([visible_len(str(r[ic]).decode("utf-8")) for r in self.matrix])+2,
                    len(str(cname))+2
                    ) )

    def sort(self, column=None, descending=True):
        self.update()
        icol = self.colnames.index(column)
        # sort matrix and range of numbers to get sorted indices for later use
        self.matrix, self.sortedidxs = zip(
                *sorted(
                    zip(
                        self.matrix,range(len(self.matrix))
                        ), key=lambda x: x[0][icol], reverse=descending
                    )
                )
        self.matrix = list(self.matrix)

        # now update row colors and hlines to match sorted matrix
        # one based indexing, not zero, so add 1
        oldtonewidx = dict([(self.sortedidxs[i]+1,i+1) for i in range(len(self.sortedidxs))])
        newrowcolors = {}
        for key,val in self.rowcolors.items():
            newrowcolors[oldtonewidx[key]] = val
        self.rowcolors = newrowcolors
        self.hlines = [oldtonewidx[hline] for hline in self.hlines]

    def get_table(self, **kwargs):
        return "".join(self.get_table_strings(**kwargs))

    def print_table(self, **kwargs):
        print self.get_table(**kwargs)

    def get_table_strings(self, bold_title=True, show_row_separators=False, show_alternating=False, ljustall=False, show_colnames=True):
        self.update()
        nrows = len(self.matrix) + 1
        draw_row_separators = True
        draw_hlines = True
        if self.theme in ["csv"]:
            draw_row_separators = False
            draw_hlines = False

        for irow,row in enumerate([self.colnames]+self.matrix):

            # line at very top
            if irow == 0:
                yield self.d_style["OUTER_TOP_LEFT"]
                for icol,col in enumerate(row):
                    yield self.d_style["OUTER_TOP_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                    if icol != len(row)-1: yield self.d_style["OUTER_TOP_INTERSECT"]
                yield self.d_style["OUTER_TOP_RIGHT"]+"\n"

            if not show_colnames and irow == 0: continue

            # lines separating columns
            yield self.d_style["OUTER_LEFT_VERTICAL"]
            oc = False if not show_alternating else (irow%2==1 )
            bold = False if not bold_title else (irow==0)
            color = self.rowcolors.get(irow,None)
            if irow == 0: color = "lightblue"
            for icol,col in enumerate(row):
                j = "l" if icol == 0 else "c"
                if ljustall: j = "l"
                yield self.fmt_string(col, self.colsizes[icol]+self.extra_padding, justify=j, bold=bold,offcolor=oc,color=color)
                if icol != len(row)-1: yield self.d_style["INNER_VERTICAL"]
            yield self.d_style["OUTER_RIGHT_VERTICAL"]+"\n"

            # lines separating rows
            if draw_row_separators and ((show_row_separators and (irow < nrows-1)) or (irow == 0)):
                yield self.d_style["OUTER_LEFT_INTERSECT"]
                for icol,col in enumerate(row):
                    yield self.d_style["INNER_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                    if icol != len(row)-1: yield self.d_style["INNER_INTERSECT"]
                yield self.d_style["OUTER_RIGHT_INTERSECT"]
                yield "\n"

            # line at very bottom
            if irow == nrows-1:
                yield self.d_style["OUTER_BOTTOM_LEFT"]
                for icol,col in enumerate(row):
                    yield self.d_style["OUTER_BOTTOM_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                    if icol != len(row)-1: yield self.d_style["OUTER_BOTTOM_INTERSECT"]
                yield self.d_style["OUTER_BOTTOM_RIGHT"]+"\n"
            else:
                # extra hlines
                if draw_hlines and (irow in self.hlines):
                    yield self.d_style["OUTER_LEFT_INTERSECT"]
                    for icol,col in enumerate(row):
                        yield self.d_style["OUTER_TOP_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                        if icol != len(row)-1: yield self.d_style["INNER_INTERSECT"]
                    yield self.d_style["OUTER_RIGHT_INTERSECT"]+"\n"


if __name__ == "__main__":

    if(sys.stdin.isatty()):

        tab = Table()
        tab.set_column_names(["name", "age", "blahhhhhh"])
        for row in [
                ["Alice", 42, 4293.9923344],
                ["Bob", 1, 0.9999999],
                ["Jim", -3, 4293],
                ["Pam", 4.2, 0.9923344],
                ["David", 4.2, 0.99999923344],
                ["John", 4.2, 0.9923344],
                ["John", 4.2, 0.9923344],
                ["John", 4.2, 0.9923344],
                ]:
            color = "green" if row[0] in ["Bob","Alice"] else None
            tab.add_row(row,color=color)
            if row[0] == "Pam":
                tab.add_line()
        # oh crap, forgot a field. no worries ;)
        tab.add_column("forgot this",range(8))
        tab.sort(column="age", descending=True)
        for w in [
                tab.set_theme_basic,
                tab.set_theme_fancy,
                tab.set_theme_latex,
                tab.set_theme_csv,
                ]:
            w()
            tab.print_table(show_row_separators=False,show_alternating=True)

        tab = Table()
        tab.set_theme_latex()
        tab.set_column_names(["col1","col2","col3","col4","col5","col6","col7"])
        for row in [
                [r"$\slashed{E}_T$", 4.2, r"\color{lightred}{4.2}"   , r"\cellcolor{lightred}{4.2}"   , r"\cellcolor{black}{\color{lightred}{4.2}}"   , r"\cellcolor{white}{\color{red}{4.2}}"   , 4.2],
                [r"$\met$"         , 4.2, r"\color{lightblue}{4.2}"  , r"\cellcolor{lightblue}{4.2}"  , r"\cellcolor{black}{\color{lightblue}{4.2}}"  , r"\cellcolor{white}{\color{blue}{4.2}}"  , 4.2],
                [r"$\mt$"          , 4.2, r"\color{lightgreen}{4.2}" , r"\cellcolor{lightgreen}{4.2}" , r"\cellcolor{black}{\color{lightgreen}{4.2}}" , r"\cellcolor{white}{\color{green}{4.2}}" , 4.2],
                [r"$\pt$"          , 4.2, r"\color{lightorange}{4.2}", r"\cellcolor{lightorange}{4.2}", r"\cellcolor{black}{\color{lightorange}{4.2}}", r"\cellcolor{white}{\color{orange}{4.2}}", 4.2],
                [r"$\Ht$"          , 4.2, r"\color{lightyellow}{4.2}", r"\cellcolor{lightyellow}{4.2}", r"\cellcolor{black}{\color{lightyellow}{4.2}}", r"\cellcolor{white}{\color{orange}{4.2}}", 4.2],
                [r"$\mt$"          , 4.2, r"4.2$\pm$0.1"             , 4.2                            , 4.2                                           , 4.2                                      , 4.2],
                ]:
            tab.add_row(row)
        tab.print_pdf(show_row_separators=False,show_alternating=False)


    else:

        import argparse

        parser = argparse.ArgumentParser()

        parser.add_argument("-f", "--first", help="first line is heading", action="store_true", default=False)
        parser.add_argument("-a", "--alternating", help="don't highlight alternating lines", action="store_true", default=False)
        parser.add_argument("-s", "--separators", help="show row separators", action="store_true", default=False)
        args = parser.parse_args()

        first_row_colnames = args.first
        rows = []
        maxcols = -1
        for row in sys.stdin: 
            row = row.strip()
            rows.append(row)
            maxcols = max(maxcols,len(row.split()))
        rows = [r for r in rows if len(r.split()) == maxcols]
        tab = Table()
        for irow,row in enumerate(rows):
            parts = row.split()
            if first_row_colnames:
                if irow == 0: tab.set_column_names(parts)
                else: tab.add_row(parts)
            else:
                tab.add_row(parts)
        tab.print_table(show_row_separators=args.separators,show_alternating=(not args.alternating))

