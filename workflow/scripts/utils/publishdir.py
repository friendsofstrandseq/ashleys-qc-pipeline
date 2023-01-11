import subprocess

for file in list(snakemake.input.list_publishdir):
    print(file)
    subprocess.Popen("echo {file}".format(file=file), shell=True, stdout=subprocess.PIPE)
    print(snakemake.config["publishdir"])
    folder_path = snakemake.config["publishdir"] + "/".join(file.replace(snakemake.config["data_location"], "").split("/")[:-1]) + "/"
    # folder_path = snakemake.config["publishdir"] + "/".join(file.replace(snakemake.config["data_location"], ""))
    print(folder_path)

    subprocess.Popen("mkdir -p {folder_path}".format(folder_path=folder_path), shell=True, stdout=subprocess.PIPE)
    print("rsync --ignore-existing -avzh --progress {file} {folder_path}".format(file=file, folder_path=folder_path))
    subprocess.Popen(
        "rsync --ignore-existing -avzh --progress {file} {folder_path}".format(file=file, folder_path=folder_path),
        shell=True,
        stdout=subprocess.PIPE,
    )
