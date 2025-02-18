#!/bin/bash

# run by typing this into the terminal on linux:
# cd Documents/code/IlseNeuro/trace_c4
# bash process_all_mice.sh

source /home/no1/Documents/code/IlseNeuro/trace_c4/.venv/bin/activate  

# Define the Dropbox path
DROPBOX_PATH="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder"

# List of mouse folders to process
MICE=("Venice" "Tallinn" "Rotterdam" "Queens")
# "Reno" "Porto" "Kyiv" "Istanbul" "Amsterdam" "Copenhagen" "Missouri"

# Function to synchronize (force Dropbox to download the folder)
sync_folder() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"
    
    echo "Ensuring $folder is synchronized..."
    dropbox exclude remove "$folder_path" 
    
    sleep 10
    
    # Find all subfolders inside the mouse folder
    find "$folder_path" -mindepth 1 -maxdepth 1 -type d | while read subfolder; do
        
        # Now, find all folders inside each subfolder
        find "$subfolder" -mindepth 1 -maxdepth 1 -type d | while read inner_folder; do
            folder_name=$(basename "$inner_folder")

            # List of folders we want to KEEP (not exclude)
            specific_folder_names=("Extraction2Bin" "Data" "c4")

            # If folder is NOT in the allowed list, exclude it
            if [[ ! " ${specific_folder_names[@]} " =~ " $folder_name " ]]; then
                echo "Excluding folder: $inner_folder"
                dropbox exclude add "$inner_folder"
            fi
        done  
        
        # Look inside each subfolder for Extraction2Bin
        extraction_folder="$subfolder/Extraction2Bin"
        
        if [[ -d "$extraction_folder" ]]; then
            echo "Checking $extraction_folder for files to exclude..."

            # Files to exclude
            excluded_files=("Data4KS2.bin" "temp_wh.dat")            
            
            for file in "${excluded_files[@]}"; do
                file_path="$extraction_folder/$file"
                
                echo "Excluding $file_path from Dropbox sync..."
                dropbox exclude add "$file_path"                
            done           
            
        fi
    done
}


# Function to check if Dropbox is done syncing
wait_for_sync() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"
    
    sleep 10

    echo "Checking sync status for: $folder"

    # Collect all subfolder paths that have a "Data" subdirectory
    mapfile -t data_folders < <(find "$folder_path" -mindepth 1 -maxdepth 1 -type d -exec test -d "{}/Data" \; -print)

    if [[ ${#data_folders[@]} -eq 0 ]]; then
        echo "No valid subfolders found in $folder. Exiting."
        return 1
    fi

    while true; do
        all_synced=true  # Assume everything is synced unless proven otherwise
        echo "Waiting for data folders to be synchronized"
        for data_folder in "${data_folders[@]}"; do
            data_path="$data_folder/Data"
            status=$(dropbox filestatus "$data_path")

            if [[ "$status" != *"up to date"* ]]; then
                all_synced=false
            fi
        done

        if $all_synced; then
            echo "All subfolders in $folder are fully synced!"
            break
        fi

        sleep 300  # Wait before checking again
    done
}


# Function to process data using Python scripts
process_data() {
    local folder="$1"
    echo "Processing data for $folder..."
    
    /home/no1/Documents/code/IlseNeuro/trace_c4/.venv/bin/python3.10 OpenEphys_wrapper_IK.py "$folder"
    
    if [[ $? -ne 0 ]]; then
        echo "Error: Data processing failed for $folder. Exiting."
        exit 1
    fi

    echo "Processing complete for $folder."
}

desync_no_longer_needed_folders() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"

      
    # Find all subfolders inside the mouse folder
    find "$folder_path" -mindepth 1 -maxdepth 1 -type d | while read subfolder; do
        
        # Now, find all folders inside each subfolder
        find "$subfolder" -mindepth 1 -maxdepth 1 -type d | while read inner_folder; do
            folder_name=$(basename "$inner_folder")

            # List of folders we want to KEEP (not exclude)
            specific_folder_names=("c4")

            # If folder is NOT in the allowed list, exclude it
            if [[ ! " ${specific_folder_names[@]} " =~ " $folder_name " ]]; then
                echo "Excluding folder: $inner_folder"
                dropbox exclude add "$inner_folder"
            fi
        done  
    done
}


# Function to process data using Python scripts
run_c4() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"

    echo "Running c4 for $folder..."
    
    /home/no1/Documents/code/IlseNeuro/trace_c4/.venv/bin/python3.10 run_run_cell_types_classifier.py "$folder"

    if [[ $? -ne 0 ]]; then
        echo "Error: c4 processing failed for $folder. Exiting."
        exit 1
    fi

    echo "c4 processing complete for $folder."
}

desync_no_longer_needed_folders_inside_c4() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"

    # Ensure c4 exists before proceeding
    local c4_path="$folder_path/c4"
    if [[ ! -d "$c4_path" ]]; then
        echo "Folder 'c4' not found inside $folder_path. Exiting."
        return
    fi

    # Iterate over all subfolders inside c4
    find "$c4_path" -mindepth 1 -maxdepth 1 -type d | while read inner_folder; do
        folder_name=$(basename "$inner_folder")
        
        keep_folder_names=("cell_type_classification")

        # If the folder is NOT "cell_type_classification", exclude it
        if [[ ! " ${keep_folder_names[@]} " =~ " $folder_name " ]]; then
            echo "Excluding folder: $inner_folder"
            dropbox exclude add "$inner_folder"
        fi
    done
}

# Function to desynchronize (remove local copy but keep cloud data)
desync_folder() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"

    #echo "Desynchronizing $folder..."
    #dropbox exclude add "$folder_path"
}

# Main loop to iterate through all mice
for MOUSE in "${MICE[@]}"; do
    echo "Starting process for $MOUSE..."

    # Ensure folder is set to sync
    sync_folder "$MOUSE"
    
    # Wait until Dropbox finishes syncing
    wait_for_sync "$MOUSE"

    # Process the data using Python
    process_data "$MOUSE"

    desync_no_longer_needed_folders "$MOUSE"    
    
    run_c4 "$MOUSE"

    #desync_no_longer_needed_folders_inside_c4 "$MOUSE"
    
    # Desynchronize after processing
    desync_folder "$MOUSE"
    
    sleep 1800

    echo "Completed processing for $MOUSE."
    echo "----------------------------------"
done

echo "All mice processed successfully!"




