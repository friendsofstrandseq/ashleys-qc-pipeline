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

    wf_info = "smk-wf-catalog/ashleys-qc-pipeline v{version}".format(
        version=str(config["version"])
    )
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
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "Folder to processed", ": " + str(config["data_location"])
        ),
    ]
    if config["genecore"] is True:
        l.append(
            f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
                "Genecore Folder to processed",
                ": " + str(config["genecore_date_folder"]),
            )
        )
    [print(e) for e in l]

    print(fg.ENDC)
    # Main options
    print("\033[1m{}\033[0m".format("Main options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "Multistep normalisation module",
            ": " + str(config["multistep_normalisation"]),
        ),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "Binning window size", ": " + str(config["window"])
        ),
    ]
    [print(e) for e in l]

    print(fg.ENDC)
    # Behavior options
    print("\033[1m{}\033[0m".format("Behavior options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "Genecore mode enabled", ": " + str(config["genecore"])
        ),
    ]
    [print(e) for e in l]

    print(fg.ENDC)
    # Genome & chrom
    chroms = (
        ["chr{e}".format(e=str(e)) for e in range(1, 23)] + ["chrX", "chrY"]
        if config["reference"] != "mm10"
        else ["chr{e}".format(e=str(e)) for e in range(1, 20)] + ["chrX", "chrY"]
    )
    if config["chromosomes"] == chroms:
        print_chroms = (
            "chr1..22,chrX,chrY"
            if config["reference"] != "mm10"
            else "chr1..19,chrX,chrY"
        )
    else:
        print_chroms = ",".join(config["chromosomes"])
    print("\033[1m{}\033[0m".format("Reference genome & Chromosomes options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "List of chromosomes processed", ": " + print_chroms
        ),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "List of chromosomes to exclude",
            ": " + ",".join(config["chromosomes_to_exclude"]),
        ),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format(
            "Reference genome selected", ": " + str(config["reference"])
        ),
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
    / __| '_ \ / _` | |/ / _ \ '_ ` _ \ / _` | |/ / _ \
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

    wf_info = "smk-wf-catalog/ashleys-qc-pipeline v{version}".format(
        version=str(config["version"])
    )
    print(sep + "\033[32m" + smk + "\033[0m")
    print("\033[33m" + wf_name + "\033[0m")
    print("\033[35m" + wf_info + "\033[0m")
    print(sep)

    config_metadata = yaml.safe_load(open("config/config_metadata.yaml", "r"))
    config_metadata = dict(
        sorted(config_metadata.items(), key=lambda item: item[1]["category"])
    )
    print("\033[30m\n\n\033[1mConfig options available (--config):\033[0m\n\n")

    header_format = "{:<30} {:<40} {:<10} {:<10} {:<30} {:<30}"
    print(
        "\033[1m"
        + header_format.format(
            "Category", "Option", "Type", "Required", "Default", "Description"
        )
        + "\033[0m"
    )

    for e in config_metadata:
        print(
            header_format.format(
                config_metadata[e]["category"],
                e,
                config_metadata[e]["type"],
                str(config_metadata[e]["required"]),
                str(config_metadata[e]["default"]),
                config_metadata[e]["desc"],
            )
        )

    print("\n\n")
