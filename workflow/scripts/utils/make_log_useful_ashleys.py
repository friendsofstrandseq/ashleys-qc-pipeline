import sys, os


def make_log_useful(log_path, status, config):
    error_buffer = []
    record = 0
    with open(log_path, "r") as logfile:
        for line in logfile:
            if "error" in line.lower() or "exception" in line.lower():
                if record == 0:
                    error_buffer.append("=== Recording ERROR entry")
                record = 10
                error_buffer.append(line.strip())
            elif line.lower().startswith("[w::") or line.lower().startswith("[e::"):
                if record == 0:
                    error_buffer.append("=== Recording library stderr")
                record = 3
                error_buffer.append(line.strip())
            elif record > 0:
                error_buffer.append(line.strip())
                record -= 1
                if record == 0:
                    error_buffer.append("=== STOP recording error entry")
            else:
                continue

    with open(log_path, "w") as logfile:
        _ = logfile.write("\n".join(error_buffer))
        _ = logfile.write("\n\n")

    my_env = dict(os.environ)
    with open(log_path, "a") as logfile:
        _ = logfile.write("=======[{}]=======\n".format(status))
        _ = logfile.write("\n===[{}]===\n".format("Infrastructure information"))
        _ = logfile.write("Host: {}\n".format(my_env.get("HOST", "N/A")))
        _ = logfile.write("Hostname: {}\n".format(my_env.get("HOSTNAME", "N/A")))
        _ = logfile.write("Display: {}\n".format(my_env.get("DISPLAY", "N/A")))
        _ = logfile.write("Shell: {}\n".format(my_env.get("SHELL", "N/A")))
        _ = logfile.write("Terminal: {}\n".format(my_env.get("TERM", "N/A")))
        _ = logfile.write("Screen: {}\n".format(my_env.get("STY", "N/A")))
        _ = logfile.write("Conda ENV: {}\n".format(my_env.get("CONDA_DEFAULT_ENV", "N/A")))

        _ = logfile.write("\n===[{}]===\n".format("Workflow information"))
        _ = logfile.write("smk-wf-catalog/ashleys-qc-pipeline v{version}\n".format(version=str(config["version"])))
        _ = logfile.write("Folder to processed : {}\n".format(str(config["data_location"])))
        _ = logfile.write("Multistep normalisation module : {}\n".format(str(config["multistep_normalisation"])))
        # _ = logfile.write("Multistep normalisation for SV calling : {}\n".format(str(config["multistep_normalisation_for_SV_calling"])))
        # _ = logfile.write("HGSVC based normalisation for SV calling : {}\n".format(str(config["hgsvc_based_normalisation"])))
        _ = logfile.write("Binning window size : {}\n".format(str(config["window"])))
        _ = logfile.write("List of chromosomes processed : {}\n".format(str(config["chromosomes"])))
        _ = logfile.write("Reference genome selected : {}\n".format(config["reference"]))
        _ = logfile.write("MultiQC : {}\n".format(config["MultiQC"]))
        # _ = logfile.write("Split PDF QC plot : {}\n".format(config["split_qc_plot"]))
        _ = logfile.write("\n")

    return


# if __name__ == "__main__":
#     make_log_useful(sys.argv[1], sys.argv[2], sys.argv[3])
