#!/bin/bash

# run by typing this into the terminal on linux:
# cd Documents/code/IlseNeuro/trace_c4
# bash process_all_mice.sh


case "$(hostname)" in
      "sphinx")
          source /home/no1/Documents/code/IlseNeuro/trace_c4/.venv/bin/activate
          # Add commands specific to Sphinx here
          ;;
      "hydra")
          source /home/devika/Documents/code/IlseNeuro/trace_c4/.venv/bin/activate
          # Add commands specific to Hydra here
          ;;
      *)
          echo "Unknown host: $(hostname)"
          # Add default commands here
          ;;
esac

DROPBOX_PATH=$(python3.10 -c "from get_dropbox_path import get_dropbox_path; print(get_dropbox_path() + '/ExperimentOutput/Ephys4Trace1/MainFolder')")

# Define the Dropbox path
#DROPBOX_PATH="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/MainFolder"

# List of mouse folders to process
MICE=("Ana1" "Ana3" "Uppsala" "Tallinn" "Porto" "York" "Xanthi" "Willemstad")
#MICE=("Quimper" "Madrid" "Dallas" "Reno" "Ana2" "Ana4" "Ana5" "Willemstad" "Zurich" "Uppsala" "Tallinn" "York" "Xanthi" "Iowa" "Ana1" "Ana3" "Zachary" "Kyiv" "Istanbul" "Copenhagen" "Rotterdam")
# "Reno" "Amsterdam"

# Function to synchronize (force Dropbox to download the folder)
sync_folder() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"
    mouse_skip_list=("Ana1" "Seattle")
    #mouse_skip_list=("Dallas" "Iowa" "Quimper" "Pittsburg" "Houston" "Yosemite" "Venice" "Seattle" "Newark" "Lisbon" "Jackson")
    # Check if folder is in the skip list
    if [[ " ${mouse_skip_list[@]} " =~ " $folder " ]]; then
        echo "Skipping synchronization for $folder"
        return
    fi

    echo "Ensuring $folder is synchronized..."
    dropbox exclude add "$folder_path"
    sleep 5
    dropbox exclude remove "$folder_path"
    # Save exclude list to a temporary file
    dropbox exclude list > exclude_list_tmp.txt

    # Filter for matching folders and save to another temporary file
    grep "$folder" exclude_list_tmp.txt | grep "$folder/$folder" > filtered_exclude_list.txt

    # Debugging: Check if anything was found
    cat filtered_exclude_list.txt  # Show what was found before proceeding

    # Remove each matching file from the exclude list
    while IFS= read -r file; do
      echo "Folder found: $file"
      dropbox exclude remove "$file"
    done < filtered_exclude_list.txt

    sleep 15
    
    # Find all subfolders inside the mouse folder
    find "$folder_path" -mindepth 1 -maxdepth 1 -type d | while read subfolder; do

        # Now, find all folders inside each subfolder
        find "$subfolder" -mindepth 1 -maxdepth 1 -type d | while read inner_folder; do
            folder_name=$(basename "$inner_folder")

            # List of folders we want to KEEP (not exclude)
            specific_folder_names=("c4") #specific_folder_names=("Extraction2Bin" "Data" "c4")

            # If folder is NOT in the allowed list, exclude it
            if [[ ! " ${specific_folder_names[@]} " =~ " $folder_name " ]]; then
                echo "Excluding folder: $inner_folder"
                dropbox exclude add "$inner_folder"
            else
                echo "Including folder: $inner_folder"
                dropbox exclude remove "$inner_folder"
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
    #echo "Folder path: $folder_path"
    folder_to_check="c4"

    echo "Checking sync status for: $folder"
    sleep 5

    # Collect all subfolder paths that have a "Data" subdirectory
    mapfile -t data_folders < <(find "$folder_path" -mindepth 1 -maxdepth 1 -type d -exec test -d "{}/$folder_to_check" \; -print)

    if [[ ${#data_folders[@]} -eq 0 ]]; then
        echo "No valid subfolders found in $folder. Exiting."
        return
    fi

    while true; do
        all_synced=true  # Assume everything is synced unless proven otherwise
        echo "Waiting for data folders to be synchronized"
        for data_folder in "${data_folders[@]}"; do
            data_path="$data_folder/$folder_to_check"
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

    case "$(hostname)" in
      "sphinx")
          /home/no1/Documents/code/IlseNeuro/trace_c4/.venv/bin/python3.10 run_script_3.py "$folder"
          # Add commands specific to Sphinx here
          ;;
      "hydra")
          /home/devika/Documents/code/IlseNeuro/trace_c4/.venv/bin/python3.10 run_script_3.py "$folder"
          # Add commands specific to Hydra here
          ;;
      *)
          echo "Unknown host: $(hostname)"
          # Add default commands here
          ;;
    esac

    
    if [[ $? -ne 0 ]]; then
        echo "Error: Data processing failed for $folder. Exiting."
        return 1
    fi

    echo "Processing complete for $folder."
}

# Function to desynchronize (remove local copy but keep cloud data)
desync_folder() {
    local folder="$1"
    local folder_path="$DROPBOX_PATH/$folder"
    mouse_skip_list=("Ana1" "Ana3" "Yosemite" "Venice" "Seattle" "Newark" "Lisbon" "Lincoln" "Jackson")
    # Check if folder is in the skip list
    if [[ " ${mouse_skip_list[@]} " =~ " $folder " ]]; then
        echo "Skipping desynchronization for $folder"
        return
    fi
    #echo "Desynchronizing $folder..."
    dropbox exclude add "$folder_path"
}

# Main loop to iterate through all mice
# Main loop to iterate through all mice
for MOUSE in "${MICE[@]}"; do
    echo "Starting process for $MOUSE..."

    # Ensure folder is set to sync
    sync_folder "$MOUSE" || { echo "Failed to sync folder for $MOUSE. Skipping..."; continue; }

    # Wait until Dropbox finishes syncing
    wait_for_sync "$MOUSE" || { echo "Failed wait for sync to finish for $MOUSE. Skipping..."; continue; }

    # Process the data using Python
    process_data "$MOUSE" || { echo "Processing failed for $MOUSE. Skipping..."; continue; }

    # Desynchronize after processing
    desync_folder "$MOUSE" || { echo "Failed to desync folder for $MOUSE. Continuing..."; }

    sleep 1800

    echo "Completed processing for $MOUSE."
    echo "----------------------------------"
done

DROPBOX_PATH=$(python3.10 -c "from get_dropbox_path import get_dropbox_path; print(get_dropbox_path() + '/ExperimentOutput/Ephys4Trace1/ReserveFolder')")

#DROPBOX_PATH="/home/no1/Lucas Bayones/BayesLab Dropbox/Lucas Bayones/TraceExperiments/ExperimentOutput/Ephys4Trace1/ReserveFolder"
RESERVE_MICE=("ReserveMouse3")
# Main loop to iterate through all mice
# Main loop to iterate through all mice
for MOUSE in "${RESERVE_MICE[@]}"; do
    echo "Starting process for $MOUSE..."

    # Ensure folder is set to sync
    sync_folder "$MOUSE" || { echo "Failed to sync folder for $MOUSE. Skipping..."; continue; }

    # Wait until Dropbox finishes syncing
    wait_for_sync "$MOUSE" || { echo "Failed wait for sync to finish for $MOUSE. Skipping..."; continue; }

    # Process the data using Python
    process_data "$MOUSE" || { echo "Processing failed for $MOUSE. Skipping..."; continue; }

    # Desynchronize after processing
    desync_folder "$MOUSE" || { echo "Failed to desync folder for $MOUSE. Continuing..."; }

    sleep 1800

    echo "Completed processing for $MOUSE."
    echo "----------------------------------"
done


echo "All mice processed successfully!"




