import re
import glob
import shutil
import os


def parse_chrom(chrs):
    clist = [x.split(":") for x in chrs.split(",")]
    parsed = []
    for chrs in clist:
        if len(chrs) == 2:
            chrs = [str(c) for c in range(int(chrs[0]), int(chrs[1]) + 1)]
        elif len(chrs) != 1:
            raise ValueError("Invalid chromosome list.")
        parsed += chrs
    return parsed


def fix_path(path):
    path = os.path.abspath(path) + '/'
    return path


def build_samp(in_path):
    p = os.path.abspath(in_path)
    p = [directory for directory,y,files in os.walk(p)
            if any("info.gz" in f for f in files)]
    return [os.path.basename(x) for x in p]


def parser(config):
    # Construct chromosome list
    chrom = parse_chrom(config["chroms"])

    # Figure out samples
    in_path = fix_path(config["directory"])
    sample = build_samp(in_path)

    # Make copy file with exclusions, else make empty one
    if "exclude_samp" in config:
        exclude = config["exclude_samp"]
        if exclude and os.path.isfile(exc):
            shutil.copyfile(exclude, "samp.irem")
        elif exclude:
            raise Exception("Invalid exclusion list.")
    open("samp.irem", 'a').close()

    # Copy file with sample keep list and set plink command
    if "include_samp" in config:
        keep = config["include_samp"]
        if keep and os.path.isfile(keep):
            shutil.copyfile(keep, "samp.ikeep")
            keep_command = "--keep samp.ikeep "
        elif keep:
            raise Exception("Invalid inclusion list.")
    if ("include_samp" not in config) or (not keep):
        keep_command = ""

    return chrom, sample, in_path, keep_command
