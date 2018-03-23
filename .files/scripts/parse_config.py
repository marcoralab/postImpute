import re
import glob
import shutil

def parser(config):
    # Construct chromosome list
    if not config["chr"]["chrlist"]:
        chrom_list = []
    else:
        chrom_list = config["chr"]["chrlist"]
    chrom = list(range(config["chr"]["from"],
                 config["chr"]["to"]+1))
    chrom += chrom_list

    # Retrieve input path and make sure it ends in a /
    in_path = config["directory"]
    if in_path[-1] != '/':
        in_path += '/'

    # Get list of all .bed files
    beds = glob.iglob(in_path + "*.bed")

    # Sample is every base name for the .bed files in the in_path directory
    expression = r"(^.*\/)(.*)(?=\.bed)"
    sample =  [re.search(expression, x)[2] for x in  beds]

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
