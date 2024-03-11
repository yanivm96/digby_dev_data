import os
import shutil
import json
import pandas as pd
import argparse
import subprocess

import datetime
import time
import sys
sys.path.append('C:/Users/yaniv/Desktop/test')  # The path to the directory containing `digby_backend`.
from digby_backend.make_vdjbase_db import make_vdjbase_db


REQUIRED_FILES = ['haplotype', 'genotype.tsv', 'ogrdb_plots.pdf', 'ogrdb_report.csv']
#SPLIT = '/'
SPLIT= '\\'
# Extracts repertoire, subject, and sample IDs from a JSON file
def get_repertoire_details(file_path):
    
    with open(file_path, 'r') as details:
        data = json.load(details)
        repertoire_id = data['repertoire_id']
        subject_id = data['subject_id']
        sample_id = data['sample_id']
    
    return repertoire_id, subject_id, sample_id

# Merges metadata from various sources into a single JSON file in the destination project
def merge_metadata(metadata_filename, project_dest, tsv_map, pre_processed_map, vdjbase_project_name, repertoire_mapping):
    with open(metadata_filename, 'r') as metadata:
        project_metadata = json.load(metadata)

        for file in tsv_map:
            repertoire_id, subject_id, sample_id = get_repertoire_details(file['repertoire_ids'])
            with open(file['annotation_metadata'], 'r') as annotation_metadata:
                annotation_metadata = json.load(annotation_metadata)
                update_annotated_metadata(project_metadata, repertoire_id, annotation_metadata)
        
        for file in pre_processed_map:
            repertoire_id, subject_id, sample_id = get_repertoire_details(file['repertoire_ids'])
            with open(file['pre_processed_metadata'], 'r') as pre_processed_metadata:
                pre_processed_metadata = json.load(pre_processed_metadata)
                update_pre_processed_metadata(project_metadata, repertoire_id, pre_processed_metadata)

        # Write the updated project_metadata to a new JSON file
        new_metadata_path = os.path.join(project_dest, f'{vdjbase_project_name}.json')
        with open(new_metadata_path, 'w') as new_metadata_file:
            json.dump(project_metadata, new_metadata_file, indent=4)

        copy_required_files(repertoire_mapping, tsv_map, project_dest)

def copy_required_files(repertoire_mapping, tsv_map, project_dest):
    for vdjbase_project in repertoire_mapping:
        vdjbase_project = repertoire_mapping[vdjbase_project]
        vdjbase_project_path = os.path.join(project_dest, 'samples', vdjbase_project['project_number'], vdjbase_project['vdjbase_name'])
        if not os.path.exists(vdjbase_project_path):
            os.makedirs(vdjbase_project_path)
        
        for projcet in tsv_map:
            repertoire_id = None
            with open(projcet['repertoire_ids'], 'r') as repertoire:
                repertoire_data = json.load(repertoire)
                repertoire_id = repertoire_data['repertoire_id']

            if repertoire_id == vdjbase_project['airr_repertoire_id']:
                for file in projcet['required_files']:
                    file_name = file.split(SPLIT)[-1]
                    new_file_name = change_file_name_to_vdjbase(vdjbase_project['vdjbase_name'], file_name)
                    destination_path = os.path.join(vdjbase_project_path, new_file_name)

                    copy_file(source_path=file, destination_path= destination_path)
        

def change_file_name_to_vdjbase(vdjbase_project_name ,file_name):
    for file in REQUIRED_FILES:
        if file in file_name:
            if file == 'haplotype':
                gene = file_name.split('Finale_')[1]
                return vdjbase_project_name + '_' + gene
            
            return vdjbase_project_name + '_' + file

# Updates project metadata with annotated metadata for a specific repertoire
def update_annotated_metadata(project_metadata, repertoire_id, annotation_metadata):
    new_data = annotation_metadata['sample']['data_processing']
    for repertoire in project_metadata['Repertoire']:
        if repertoire['repertoire_id'] == repertoire_id:
            original_data = repertoire['data_processing'][0]
            repertoire['data_processing'][0] = merge_json_data_recursive(original_data, new_data)

# Updates project metadata with pre-processed metadata for a specific repertoire
def update_pre_processed_metadata(project_metadata, repertoire_id, pre_processed_metadata):
    new_data = pre_processed_metadata['sample']['data_processing']
    for repertoire in project_metadata['Repertoire']:
        if repertoire['repertoire_id'] == repertoire_id:
            original_data = repertoire['data_processing'][0]
            repertoire['data_processing'][0] = merge_json_data_recursive(original_data, new_data)


def merge_json_data_recursive(original_data, new_data):
    """
    Recursively merges new_data into original_data. If a key in new_data already exists in original_data
    and both values are dictionaries, it merges them recursively. If both are lists, it appends the items
    from the new list to the old list. Otherwise, the value in original_data is updated with the value from new_data.
    """
    for key, value in new_data.items():
        if key in original_data:
            if isinstance(original_data[key], dict) and isinstance(value, dict):
                merge_json_data_recursive(original_data[key], value)
            elif isinstance(original_data[key], list) and isinstance(value, list):
                original_data[key].extend(value)
            else:
                original_data[key] = value
        else:
            original_data[key] = value

    return original_data


# Copies content from a source directory to a destination directory and merges metadata
def copy_folder_content(source_folder, target_repo_path, vdjbase_project_name,project_number, metadata_filename, repertoire_mapping):
    # Create the destination directory if it does not exist
    if not os.path.exists(target_repo_path):
        os.makedirs(target_repo_path)
    
    tsv_files_paths, pre_processed_files = find_project_tsv_files(source_folder)

    merge_metadata(metadata_filename, target_repo_path, tsv_files_paths, pre_processed_files, vdjbase_project_name, repertoire_mapping)


def copy_file(source_path, destination_path):
    """
    Copy a file from source_path to destination_path.

    """
    # Check if the source file exists
    if not os.path.isfile(source_path):
        raise FileNotFoundError(f"The source file does not exist: {source_path}")
    
    # Ensure the destination directory exists, if not, create it
    os.makedirs(os.path.dirname(destination_path), exist_ok=True)
    
    # Copy the file
    shutil.copy2(source_path, destination_path)
    #print(f"File copied from {source_path} to {destination_path}")


# Finds TSV files and pre-processed files within a project directory
def find_project_tsv_files(project_path):
    pre_processed_folders = []
    try:
        annotated_folder_path = os.path.join(project_path, 'annotated')
        annotated_folders = os.listdir(annotated_folder_path)
        pre_processed_folder_path = os.path.join(project_path, 'pre_processed')
        if os.path.exists(pre_processed_folder_path):
            pre_processed_folders = os.listdir(pre_processed_folder_path)
            
        tsv_files = start_scan(annotated_folder_path,annotated_folders , False)
        pre_processed_files = start_scan(pre_processed_folder_path, pre_processed_folders , True)

    except Exception as e:
        print(e)

    return tsv_files, pre_processed_files

# Initiates scanning of folders for TSV and metadata files
def start_scan(folder_path,folders, pre_processed):
    files = []
    for subject in folders:
        subject_path = os.path.join(folder_path, subject)
        res = scan_subject_folder(subject_path, pre_processed)
        for file in res:
            files.append(file)
    
    return files

# Scans a subject folder for TSV and metadata files
def scan_subject_folder(subject_path, pre_processed):
    tsv_files = []
    samples = os.listdir(subject_path)
    for sample in samples:
        sample_path = os.path.join(subject_path, sample)
        files = scan_run_folder(sample_path, pre_processed)
        for file in files:
            tsv_files.append(file)
    
    return tsv_files

# Scans a run folder for TSV and metadata files
def scan_run_folder(sample_path, pre_processed):
    tsv_files = []
    repertoires = os.listdir(sample_path)
    for rep in repertoires:
        rep_path = os.path.join(sample_path, rep)
        if not pre_processed:
            file = find_tsv_and_metadata_for_annotated(rep_path)
        else:
            file = find_metadata_for_pre_processed(rep_path)

        if file != None:
            tsv_files.append(file[0])
    
    return tsv_files

# Finds metadata for pre-processed results
def find_metadata_for_pre_processed(result_path):
    res_list = []
    res = {
        'repertoire_ids': None,
        'pre_processed_metadata': None
    }
    result_folders = os.listdir(result_path)
    for folder in result_folders:
        folder_path = os.path.join(result_path, folder)
        folder_files = os.listdir(folder_path)
        
        if 'meta_data' in folder:
            if 'pre_processed_metadata.json' in folder_files:
                res['pre_processed_metadata'] = os.path.join(folder_path, 'pre_processed_metadata.json')
            
            if 'repertoire_id.json' in folder_files:
                res['repertoire_ids'] = os.path.join(folder_path, 'repertoire_id.json')

    check_result_fileds(res, result_path)
    if all(value is not None for value in res.values()):
        res_list.append(res)
        return res_list
    
    return None     
                

# Finds TSV files and their corresponding metadata for annotated results
def find_tsv_and_metadata_for_annotated(result_path):
    res_list = []
    res = {
            'file_path': None,
            'file_name': None,
            'repertoire_ids': None,
            'annotation_metadata': None,
            'required_files' : []
        }
    
    result_folders = os.listdir(result_path)
    for folder in result_folders:
        folder_path = os.path.join(result_path, folder)
        folder_files = os.listdir(folder_path)
        for file in folder_files:
            if 'Finale' in file:
                res['file_path'] = os.path.join(folder_path, file)
                res['file_name'] = file
            
            if file == 'repertoire_id.json':
                res['repertoire_ids'] = os.path.join(folder_path, file)
                
            for required_file in REQUIRED_FILES:
                if required_file in file:
                    res['required_files'].append(os.path.join(folder_path, file))


        if 'meta_data' in folder:
            if 'annotation_metadata.json' in folder_files:
                res['annotation_metadata'] = os.path.join(folder_path, 'annotation_metadata.json')

    check_result_fileds(res, result_path)
    if all(value is not None for value in res.values()):
        res_list.append(res)
        return res_list
    
    return None 

# Checks if all required fields in a result are present
def check_result_fileds(result, folder):
    for key, value in result.items():
        if value == None:
            print(f"{key} was not found in the {folder}")



def verify_directory_exists(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Directory does not exist: {path}")

def verify_and_clear_project_directory(target_repo_path, project_name):
    airr_correspondence_path = os.path.join(target_repo_path, 'airr_correspondence.csv')
    if not os.path.isfile(airr_correspondence_path):
        raise FileNotFoundError(f"airr_correspondence.csv not found in the target repo at {target_repo_path}")

    # Read the CSV file into a DataFrame
    correspondence_df = pd.read_csv(airr_correspondence_path)

    # Check if the project name is in the airr_file column and extract the first matching vdjbase_name
    vdjbase_project_name = None
    for airr_file in correspondence_df['airr_file']:
        if project_name in airr_file:
            vdjbase_project_name = correspondence_df[correspondence_df['airr_file'] == airr_file]['vdjbase_name'].str.split('_').str[0].iloc[0]
            break

    if vdjbase_project_name is None:
        raise ValueError(f"No matching project name {project_name} found in airr_correspondence.csv")

    # Path to the project directory within the target repo
    project_dir_path = os.path.join(target_repo_path, 'samples' ,vdjbase_project_name)
    
    if not os.path.isdir(project_dir_path):
        raise FileNotFoundError(f"Project directory {project_dir_path} does not exist")

    # Clear the project directory
    for filename in os.listdir(project_dir_path):
        file_path = os.path.join(project_dir_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

    print(f"Project directory {project_dir_path} has been cleared.")
    return airr_correspondence_path


def derive_vdjbase_project_mapping(airr_correspondence_path, project_name):
    correspondence_df = pd.read_csv(airr_correspondence_path)
    
    # Filter the DataFrame for rows where the vdjbase_name contains the specific project name
    project_specific_df = correspondence_df[correspondence_df['airr_file'].str.contains(f"{project_name}", regex=True)]

    mapping = {}
    for index, row in project_specific_df.iterrows():
        vdjbase_name = row['vdjbase_name']
        airr_repertoire_id = row['airr_repertoire_id']
        
        # Extract individual and sample from vdjbase_name
        project_number, individual, sample = vdjbase_name.split('_')

        mapping[vdjbase_name] = {
            'project_name' : project_name,
            'vdjbase_name' : vdjbase_name, 
            'project_number': project_number,
            'individual': individual,
            'sample': sample,
            'airr_repertoire_id': airr_repertoire_id
        }

    return mapping


def verify_annotations_exist(sequence_data_store_path, airr_correspondence_path, project_name):
    # Read the airr_correspondence.csv to get the expected repertoire IDs and corresponding file patterns
    correspondence_df = pd.read_csv(airr_correspondence_path)
    filtered_df = correspondence_df[correspondence_df['airr_file'].str.contains(project_name)]
    expected_repertoires = filtered_df['airr_repertoire_id'].tolist()
    #expected_repertoires = correspondence_df['airr_repertoire_id'].to_list()
    
    # Dictionary to hold the found 'Final' files
    found_files = []

    # Walk through the annotated directory to find 'Final' files
    annotated_path = os.path.join(sequence_data_store_path, 'annotated')
    if not os.path.isdir(annotated_path):
        raise FileNotFoundError(r"there is no annotated file for {sequence_data_store_path}")
        

    for root, dirs, files in os.walk(annotated_path):
        for file in files:
            if "Final" in file:
                found_files.append(file)

    # Verify that each expected file pattern has at least one matching 'Final' file
    missing_annotations = []
    for repertoire in expected_repertoires:
        if not any(repertoire in f for f in found_files):
            missing_annotations.append(repertoire)

    if missing_annotations:
        raise FileNotFoundError(f"Missing 'Final' annotation files for the following patterns: {', '.join(missing_annotations)}")
    
    return True  # Return True if all checks pass


def consolidate_metadata(repertoire_metadata_file, additional_metadata_paths):
    with open(repertoire_metadata_file, 'r') as f:
        metadata = json.load(f)
    
    for metadata_path in additional_metadata_paths:
        if os.path.isfile(metadata_path):
            with open(metadata_path, 'r') as f:
                additional_metadata = json.load(f)
            metadata.update(additional_metadata)  # Assumes that additional metadata should overwrite
    
    return metadata

def check_file(path):
    # Check if the file exists and is not empty
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return True
    else:
        return False

def last_modified_time(path):
    # Get the last modified time of the file
    return time.ctime(os.path.getmtime(path))

def check_files_updated():
    db_path = os.path.join(target_repo_path, "db.sqlite3")
    samples_path = os.path.join(target_repo_path, "samples.zip")

    # Check db.sqlite3
    if not check_file(db_path):
        raise Exception("db.sqlite3 does not exist or is empty.")

    # Check samples.zip
    if not check_file(samples_path):
        raise Exception("samples.zip does not exist or is empty.")

def run_git_command(command):
    try:
        result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True, universal_newlines=True)
        return result.strip()
    except subprocess.CalledProcessError as e:
        print("Error executing Git command:", e.output)
        sys.exit(1)

def is_repo_up_to_date(repo_path):
    # Change to the repository directory
    original_dir = subprocess.os.getcwd()
    subprocess.os.chdir(repo_path)

    # Fetch the latest updates from remote without merging
    run_git_command("git fetch")

    # Get the latest commit hash of the local branch
    local_commit = run_git_command("git rev-parse HEAD")

    # Get the latest commit hash of the remote branch
    remote_commit = run_git_command("git rev-parse @{u}")

    # Change back to the original directory
    subprocess.os.chdir(original_dir)

    return local_commit == remote_commit



def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[((n//10%10!=1)*(n%10<4)*n%10)::4])

def update_description_file(target_repo_path):
    description_path = os.path.join(target_repo_path, "db_description.txt")
    res = target_repo_path.split('\\') #need to change!!!
    with open(description_path, 'w') as file:
        date = datetime.datetime.now()
        formatted_date = date.strftime(f"{ordinal(date.day)} %B %Y")
        file.write(f'Analysis of {res[-2]} {res[-1]} datsets, compiled {formatted_date}')

def main(project_name, source_folder, metadata_filename, target_repo_path):
    try:
        repo_path = r'C:\Users\yaniv\Desktop\test\digby_data'
        #Verify that the source folder and target repo path exist
        if is_repo_up_to_date(repo_path):
            print("The repository is up-to-date with the remote GitHub repository.")
        else:
            raise Exception("The repository is not up-to-date with the remote GitHub repository.")

        verify_directory_exists(source_folder)
        verify_directory_exists(target_repo_path)

        # Verify airr_correspondence.csv file exists and get the mapping
        airr_correspondence_path = verify_and_clear_project_directory(target_repo_path, project_name) #need to add the check of contains references to a file matching the project, in the airr_file column.
        repertoire_mapping = derive_vdjbase_project_mapping(airr_correspondence_path, project_name)
        
        # Verify annotations exist
        verify_annotations_exist(source_folder, airr_correspondence_path, project_name)
        first_key = next(iter(repertoire_mapping.keys()))# Access the first key
        project_number = repertoire_mapping[first_key]['project_number']
        vdjbase_project_name = project_number + ('_' + project_name)
        copy_folder_content(source_folder, target_repo_path, vdjbase_project_name, project_number, metadata_filename, repertoire_mapping)

        print("Data copy completed successfully.")
        os.chdir(target_repo_path)
        bat_file = os.path.join(target_repo_path, 'make.bat')
        # subprocess.run(bat_file, shell=True)
        with open("make.log", "w") as log_file:
            original_stdout = sys.stdout
            sys.stdout = log_file
            make_vdjbase_db("Human", "IGH")
            sys.stdout = original_stdout

        # Change directory to 'samples'
        os.chdir("samples")

        # Execute the 7z command to zip the contents
        subprocess.run([r"C:\Program Files\7-Zip\7z.exe", "a", "-tzip", r"..\samples.zip", "*"])

        # Change directory back to the original
        os.chdir("..")
        # Change directory back to the original
        check_files_updated()
        update_description_file(target_repo_path)

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    # Create the parser
    # parser = argparse.ArgumentParser(description='Process some inputs.')

    # # Add arguments
    # parser.add_argument('project_name', type=str, help='Name of the project')
    # parser.add_argument('source_folder', type=str, help='Path to the source folder')
    # parser.add_argument('metadata_filename', type=str, help='Path to the metadata file')
    # parser.add_argument('target_repo_path', type=str, help='Path to the target repository')

    # # Parse the arguments
    # args = parser.parse_args()
    #main(args.project_name, args.source_folder, args.metadata_filename, args.target_repo_path)
   
   # Hardcoded for demonstration purposes
    project_name = r"PRJNA248475"
    source_folder = r"C:\Users\yaniv\Desktop\PRJNA248475\runs\current"
    metadata_filename = r"C:\Users\yaniv\Desktop\PRJNA248475\project_metadata\metadata.json"
    target_repo_path = r"C:\Users\yaniv\Desktop\test\digby_data\AIRR-seq\Human\IGH"
    main(project_name, source_folder, metadata_filename, target_repo_path)
#python your_script.py "PRJNA248411" "/home/bcrlab/malachy7/sequence_data_store_test/PRJNA248411/runs/current/" "/home/bcrlab/malachy7/sequence_data_store_test/PRJNA248411/project_metadata/metadata.json" "/home/bcrlab/malachy7/digby_data/AIRR-seq/Human/IGH/"
