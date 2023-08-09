import time
import os, sys, glob, subprocess, re
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from datetime import datetime
import logging
import json
import pandas as pd
import threading


os.makedirs("watchdog/logs", exist_ok=True)

# Setup the logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("watchdog/logs/watchdog_ashleys.log"),  # File handler to log to a file
        logging.StreamHandler(),  # Stream handler to log to the console
    ],
)


# Set the path you want to watch
path_to_watch = sys.argv[1]

data_location = "/scratch/tweber/DATA/MC_DATA/STOCKS"
publishdir_location = "/g/korbel/WORKFLOW_RESULTS"
genecore_prefix = path_to_watch
profile_slurm = ["--profile", "workflow/snakemake_profiles/HPC/slurm_EMBL/"]
profile_dry_run = ["--profile", "workflow/snakemake_profiles/local/conda/"]
dry_run_options = ["-c", "1", "-n", "-q"]
snakemake_binary = "/g/korbel2/weber/miniconda3/envs/snakemake_latest/bin/snakemake"

# plates_processing_status = pd.read_csv("watchdog/processing_status.json", sep="\t")
# print(plates_processing_status)


# Define the event handler
class MyHandler(FileSystemEventHandler):
    def on_created(self, event):
        if event.is_directory:  # if a directory is created
            logging.info(f"Directory {event.src_path} has been created!")
            self.process_new_directory(event.src_path)

    def check_unprocessed_folder(self):
        unwanted = ["._.DS_Store", ".DS_Store", "config"]
        list_runs_processed = sorted([e for e in os.listdir(data_location) if e not in unwanted])
        total_list_runs = sorted([e for e in os.listdir(path_to_watch) if e not in unwanted])
        unprocessed_plates = set(total_list_runs).difference(list_runs_processed)
        for plate in unprocessed_plates:
            # if plate not in plates_processing_status["plate"].values.tolist():
            # plates_processing_status_plate_dict = collections.defaultdict(dict)
            nb_txt_gz_files = len(glob.glob(f"{path_to_watch}/{plate}/*.txt.gz"))
            if nb_txt_gz_files == 576:
                print(f"PROCESSING {path_to_watch}/{plate}")
                self.process_new_directory(f"{path_to_watch}/{plate}")
            else:
                print(f"Not possible to process {path_to_watch}/{plate}, containing {nb_txt_gz_files} txt.gz files")

    def process_new_directory(self, directory_path):
        """Process the new directory, check for .txt.gz files and execute snakemake command if conditions are met."""

        # Poll the directory until 576 files appear or a timeout is reached
        timeout = 60  # Timeout in seconds
        start_time = time.time()

        while True:
            # Count the number of .txt.gz files in the new directory
            txt_gz_files = glob.glob(directory_path + "/*.txt.gz")
            num_files = len(txt_gz_files)

            # If the desired number of files is found or timeout is reached, break the loop
            if num_files == 576 or time.time() - start_time > timeout:
                break

            # Sleep for a while before the next poll
            time.sleep(5)  # Sleep for 5 seconds

        # Process the found .txt.gz files
        self.process_txt_gz_files(directory_path, txt_gz_files, num_files)

    def process_txt_gz_files(self, directory_path, txt_gz_files, num_files):
        """Process the found .txt.gz files and execute snakemake command if conditions are met."""

        if num_files == 576:
            logging.info(f"The new directory contains exactly 576 .txt.gz files.")
            self.execute_snakemake(directory_path, txt_gz_files)

        else:
            logging.info(f"The new directory contains {str(num_files)} .txt.gz files, not 576.")

    def execute_snakemake(self, directory_path, txt_gz_files):
        """Execute the snakemake command based on the found prefixes."""

        pattern = re.compile(r"(iTRU|PE20)\d{3}")
        prefixes = set()

        for file_path in txt_gz_files:
            match = pattern.search(file_path)
            if match:
                prefix = match.group()[:4]  # Get the first 4 characters, which is the prefix
                prefixes.add(prefix)

        if len(prefixes) > 1:
            logging.info("Multiple different prefixes found: %s", prefixes)
        elif prefixes:
            self.execute_command(directory_path, prefixes.pop())
        else:
            logging.info("No match found in any file.")

    def execute_command(self, directory_path, prefix):
        """Execute the command."""

        # Change directory and run the snakemake command
        date_folder = directory_path.split("/")[-1]

        cmd = [
            f"{snakemake_binary}",
            "--config",
            "genecore=True",
            f"genecore_prefix={genecore_prefix}",
            f"genecore_date_folder={date_folder}",
            f"genecore_regex_element={prefix}",
            "multistep_normalisation=True",
            "MultiQC=True",
            "split_qc_plot=False",
            f"publishdir={publishdir_location}",
            "email=thomas.weber@embl.de",
            f"data_location={data_location}",
            "--nolock",
        ]

        logging.info("Running command: %s", " ".join(cmd + profile_dry_run + dry_run_options))

        process = subprocess.Popen(
            cmd + profile_dry_run + dry_run_options, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True
        )

        # Variable to store the penultimate line
        penultimate_line = ""

        # Read the output line by line in real-time
        for line in iter(process.stdout.readline, ""):
            logging.info(line.strip())  # log line in real-time
            if line.strip():  # If line is not blank
                penultimate_line = line.strip()

        # Wait for the subprocess to finish
        process.wait()
        logging.info("Return code: %s", process.returncode)

        # Check the penultimate line
        if str(process.returncode) == str(0):
            self.run_second_command(cmd, profile_slurm, data_location, date_folder)
        else:
            logging.info("\nThe output is not as expected.")

    def run_second_command(self, cmd, profile_slurm, data_location, date_folder):
        """Run the second command and write the output to a log file."""

        logging.info("\nThe output is as expected.")
        logging.info("Running command: %s", " ".join(cmd + profile_slurm))

        os.makedirs("watchdog/logs/per-run", exist_ok=True)

        # Get the current date and time
        now = datetime.now()

        # Convert it to a string
        current_time = now.strftime("%Y%m%d%H%M%S")

        with open(f"watchdog/logs/per-run/{date_folder}_{current_time}.log", "w") as f:
            process2 = subprocess.Popen(cmd + profile_slurm, stdout=f, stderr=f, universal_newlines=True)
            process2.wait()

            logging.info("Return code: %s", process2.returncode)

        # Change the permissions of the new directory
        subprocess.run(["chmod", "-R", "777", f"{data_location}/{date_folder}"])


def main():
    # Create the event handler
    event_handler = MyHandler()

    # Create an observer
    observer = Observer()

    # Assign the observer to the path and the event handler
    observer.schedule(event_handler, path_to_watch, recursive=False)

    # Start the observer
    observer.start()

    # Start the periodical directory scanning in a separate thread
    def periodic_scan():
        while True:
            event_handler.check_unprocessed_folder()
            time.sleep(3600)  # Scan the directory every hour

    scan_thread = threading.Thread(target=periodic_scan)
    scan_thread.start()

    try:
        while True:
            logging.info("Waiting for new plate ...")
            time.sleep(3600)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()


if __name__ == "__main__":
    main()
