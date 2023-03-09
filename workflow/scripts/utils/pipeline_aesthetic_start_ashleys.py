class fg:
    BLACK = "\u001b[30m"
    RED = "\u001b[31m"
    GREEN = "\u001b[32m"
    YELLOW = "\u001b[33m"
    BLUE = "\u001b[34m"
    MAGENTA = "\u001b[35m"
    CYAN = "\u001b[36m"
    WHITE = "\u001b[37m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def pipeline_aesthetic_start(config):

    sep = """------------------------------------------------------"""

    smk = """
                     _                        _        
     ___ _ __   __ _| | _____ _ __ ___   __ _| | _____ 
    / __| '_ \ / _` | |/ / _ \ '_ ` _ \ / _` | |/ / _ \\
    \__ \ | | | (_| |   <  __/ | | | | | (_| |   <  __/
    |___/_| |_|\__,_|_|\_\___|_| |_| |_|\__,_|_|\_\___|
    """

    wf_name = """                                                   
              _     _                                                 
     __ _ ___| |__ | | ___ _   _ ___        __ _  ___ 
    / _` / __| '_ \| |/ _ \ | | / __|_____ / _` |/ __|
   | (_| \__ \ | | | |  __/ |_| \__ \_____| (_| | (__
    \__,_|___/_| |_|_|\___|\__, |___/      \__, |\___|
                           |___/              |_|     
    """

    wf_info = "smk-wf-catalog/ashleys-qc-pipeline v{version}".format(version=str(config["version"]))
    print(sep + fg.GREEN + smk)
    print(fg.ENDC)
    print(fg.YELLOW + wf_name)
    print(fg.ENDC)
    print(fg.MAGENTA + wf_info)
    print(fg.ENDC)
    print(sep)

    # Input / Output
    print("\033[1m{}\033[0m".format("Input/Output options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Folder to processed", ": " + str(config["data_location"])),
    ]
    if config["genecore"] is True:
        l.append(
            f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Genecore Folder to processed", ": " + str(config["genecore_date_folder"]))
        )
    [print(e) for e in l]

    print(fg.ENDC)
    # Main options
    print("\033[1m{}\033[0m".format("Main options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("GC analysis", ": " + str(config["GC_analysis"])),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Binning window size", ": " + str(config["window"])),
    ]
    [print(e) for e in l]

    print(fg.ENDC)
    # Behavior options
    print("\033[1m{}\033[0m".format("Behavior options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Genecore mode enabled", ": " + str(config["genecore"])),
    ]
    [print(e) for e in l]

    print(fg.ENDC)
    # Genome & chrom
    chroms = ["chr{e}".format(e=str(e)) for e in range(1, 23)] + ["chrX", "chrY"]
    if config["chromosomes"] == chroms:
        print_chroms = "chr1..22,chrX,chrY"
    else:
        print_chroms = ",".join(config["chromosomes"])
    print("\033[1m{}\033[0m".format("Reference genome & Chromosomes options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("List of chromosomes processed", ": " + print_chroms),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "List of chromosomes to exclude", ": " + ",".join(config["chromosomes_to_exclude"])
        ),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Reference genome selected", ": " + str(config["reference"])),
        # f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Reference FASTA file", ": " + str(config["references_data"][config["reference"]]["reference_file_location"])),
    ]
    [print(e) for e in l]
    print("\n\n")


def argparse_help(config):
    import argparse
    import yaml
    import sys, os

    # config = yaml.safe_load(open("config/config.yaml", "r"))
    # pipeline_aesthetic_start(config)

    sep = """------------------------------------------------------"""

    smk = """
                     _                        _        
     ___ _ __   __ _| | _____ _ __ ___   __ _| | _____ 
    / __| '_ \ / _` | |/ / _ \ '_ ` _ \ / _` | |/ / _ \\
    \__ \ | | | (_| |   <  __/ | | | | | (_| |   <  __/
    |___/_| |_|\__,_|_|\_\___|_| |_| |_|\__,_|_|\_\___|
    """

    wf_name = """                                                   
              _     _                                                 
     __ _ ___| |__ | | ___ _   _ ___        __ _  ___ 
    / _` / __| '_ \| |/ _ \ | | / __|_____ / _` |/ __|
   | (_| \__ \ | | | |  __/ |_| \__ \_____| (_| | (__
    \__,_|___/_| |_|_|\___|\__, |___/      \__, |\___|
                           |___/              |_|     
    """

    wf_info = "smk-wf-catalog/ashleys-qc-pipeline v{version}".format(version=str(config["version"]))
    print(sep + fg.GREEN + smk)
    print(fg.ENDC)
    print(fg.YELLOW + wf_name)
    print(fg.ENDC)
    print(fg.MAGENTA + wf_info)
    print(fg.ENDC)
    print(sep)

    config_metadata = yaml.safe_load(open("config/config_metadata.yaml", "r"))
    # parser = argparse.ArgumentParser(prog="mosaicatcher-pipeline", add_help=False)
    # group = parser.add_argument_group("--config")
    print(f"{fg.BLACK}\n\n\033[1mConfig options available (--config):\033[0m\n\n")
    l1 = [
        f"{fg.BLACK}\033[1m       {{:<20}}{fg.BLACK}       {{:<10}}{fg.BLACK}       {{:<10}}{fg.BLACK}       {{:<20}}{fg.BLACK}{{:<30}}\033[0m".format(
            "Option", "Type", "Required", "Default", "Description"
        )
    ]
    l = [
        f"{fg.BLACK}       {{:<20}}{fg.BLACK}       {{:<10}}{fg.BLACK}       {{:<10}}{fg.BLACK}       {{:<20}}{fg.BLACK}{{:<30}}".format(
            e,
            config_metadata[e]["type"],
            str(config_metadata[e]["required"]),
            str(config_metadata[e]["default"]),
            ": " + config_metadata[e]["desc"],
        )
        for e in config_metadata
    ]
    [print(e) for e in l1 + l]
    print("\n\n")
    # args = parser.parse_args()
    # print(pipeline_aesthetic_start.fg.BOLD)
    # parser.print_help()
    # print(pipeline_aesthetic_start.fg.ENDC)
