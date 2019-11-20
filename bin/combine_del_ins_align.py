#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
combine delAlignment and insAlignment together to form final sequence
"""

__author__ = "LXF"
__version__ = "1.0"


from pathlib import Path
import numpy as np
import pandas as pd
import re
import argparse

parser = argparse.ArgumentParser(description="generate frequency matrix from raw")
parser.add_argument('-i', metavar='infile', type=argparse.FileType('r'), help="input file")
parser.add_argument('-o', metavar='outfile', type=argparse.FileType('w'), help="output file")
args = parser.parse_args()


def main():
    infile_name = args.i
    df = pd.read_csv(infile_name, sep=",", na_filter=False)
    wt_df = df.loc[(df["sample_name"] == "sample") & (df["seqType"] == "WT"), ["gene_name", "finalSeq"]]
    wt_new_df = wt_df.set_index(wt_df["gene_name"])
    # pos from 0, the pos of the space before "|"
    wt_new_df["cut_site"] = wt_new_df["finalSeq"].str.find("|") - 1
    # change NA to wt_seq for wt, not so necessary
    wt_dict = wt_new_df.to_dict("index")
    condition = (df["seqType"] == "WT") & (df["finalSeq"] == "NA")
    df.loc[condition, ["delAlignment", "insAlignment", "finalSeq"]] = \
        df[condition]["gene_name"].map(lambda x: wt_dict[x]["finalSeq"])
    # add cut_site to df
    df["WT_seq_con"] = df["gene_name"].map(lambda x: wt_dict[x]["finalSeq"])
    df["cut_site"] = df["gene_name"].map(lambda x: wt_dict[x]["cut_site"])

    # combine sequence
    # replace target seq
    df["combinedSeq"] = df[["finalSeq", "WT_seq_con", "cut_site", "seqType"]].apply(lambda x: combine_seq(*x), axis=1)

    # add cut_site
    df["combinedSeq"] = df[["combinedSeq", "cut_site"]].apply(lambda x: add_cutsite(*x), axis=1)

    df.to_csv(args.o, header=True, index=False)


# combine ins alt and del alt together
# AT---{CC}AAA -> AT-CCAAA
def combine_seq(seq, original, site, seqType):
    """
        pattern1: insertion before deletion
        pattern1 = re.compile(r'\{[ATCG]+\}[ \| ]*[-]+[ \| ]*[-]*')
        pattern2: deletion before insertion
        pattern2 = re.compile(r'[-]*[ \| ]*[-]+[ \| ]*\{[ATCG]+\}')
        pattern3: insertion between deletion
        pattern4: deletion between insertion
        get replace sequence for INDEL
    """
    match_pat_indel = re.compile(r'\{[ATCG]+\}[ \| ]*[-]+[ \| ]*[-]*|[-]+[ \| ]*[-]*[ \| ]*\{[ATCG]+\}')
    match_pat_del = re.compile(r'[-]*[ \| ]+[-]*')
    match_pat_ins = re.compile(r'\{.+\}')
    if seqType == "DEL":
        match_pat = match_pat_del
    elif seqType == "INDEL":
        match_pat = match_pat_indel
    else:
        match_pat = match_pat_ins

    matchobj = re.search(match_pat, seq)
    if matchobj:
        ms, me = matchobj.span()
        match_seq = matchobj.group(0)
        cut_site_mseq = re.search(' \| ', match_seq)
        if cut_site_mseq:
            cut_s, cut_e = cut_site_mseq.span()
            match_seq_rm = match_seq[:cut_s] + match_seq[cut_e:]
        else:
            match_seq_rm = match_seq
        del_len = match_seq_rm.count("-")
        ins_s = match_seq_rm.find("{")
        ins_e = match_seq_rm.find("}")
        ins_len = ins_e - ins_s - 1
        if ins_len == -1 or del_len == 0:
            alt_seq = match_seq_rm
        else:
            ins_seq = match_seq_rm[ins_s + 1:ins_e]
            diff_len = ins_len - del_len
            # ID
            if match_seq_rm[0] == "{":
                if diff_len > 0:
                    # diff_len = ins_len
                    alt_seq = "{" + ins_seq[:diff_len] + "}" + "[" + ins_seq[diff_len:] + "]"
                else:
                    alt_seq = "[" + ins_seq + "]" + abs(diff_len) * "-"
            # DI
            elif match_seq_rm[0] == "-":
                if diff_len > 0:
                    alt_seq = "[" + ins_seq[:del_len] + "]" + "{" + ins_seq[del_len:] + "}"
                else:
                    alt_seq = abs(diff_len) * "-" + "[" + ins_seq + "]"

        # compare the cut_site of seq with the cut_site of WT seq
        cut_site = re.search(' \| ', seq)
        cut_ss, cut_es = cut_site.span()
        seq_b = seq[:cut_ss]
        ins_s = seq_b.find("{")
        ins_e = seq_b.find("}")
        ins_l = ins_e - ins_s + 1
        if ins_s == -1:
            cut_site_seq = cut_ss
        else:
            cut_site_seq = cut_ss - ins_l
        if cut_site_seq == site:
            new_seq = original[:ms] + alt_seq + original[me:]
        else:
            dis = site - cut_site_seq
            new_seq = original[:ms+dis] + alt_seq + original[me+dis:]
    else:
        new_seq = original

    return new_seq


# {ATCG}-- | -ATG -> [AT] | [C]{G}ATG
def add_cutsite(new_seq, cut_site):
    cut_flag = re.search(' \| ', new_seq)
    if cut_flag:
        final_seq = new_seq
    else:
        snv_s = [i for i, c in enumerate(new_seq) if c == "["]
        snv_e = [j for j, c in enumerate(new_seq) if c == "]"]
        ins_s = [m for m, c in enumerate(new_seq) if c == "{"]
        ins_e = [n for n, c in enumerate(new_seq) if c == "}"]
        flag_snv = 0
        for ss, se in zip(snv_s, snv_e):
            if se <= cut_site + 1:
                cut_site = cut_site + 2
            elif ss < cut_site and se > cut_site + 1:
                cut_site = cut_site + 1
                flag_snv = 1
        for inss, inse in zip(ins_s, ins_e):
            ins_len = inse - inss - 1
            if inse <= cut_site:
                cut_site = cut_site + ins_len + 2
        if flag_snv:
            final_seq = new_seq[:cut_site] + "] | [" + new_seq[cut_site:]
        else:
            final_seq = new_seq[:cut_site] + " | " + new_seq[cut_site:]

    return final_seq


if __name__ == "__main__":
    main()