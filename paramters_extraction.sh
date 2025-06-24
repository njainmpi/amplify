#!/bin/sh

#Following script has been made by Naman Jain with following features included in the
#different version upgrades


#Read All the parameters in the scan and compiles them in a text file for easy reading


## Main Script Starts from here
File_with_Dataset_Names="/Users/njain/Desktop/data.txt"

indices=(1) #enter the index number of the file name that you would like to analyse

for datasets in "${indices[@]}"; do
    
    DatasetName=$(awk -F "\"*,\"*" -v var="$datasets" 'NR == var {print $1}' $File_with_Dataset_Names)
    echo "Dataset Currently Being Analysed is": $DatasetName

    #Locate the source of Raw Data on the server, this needs to be changed by the user based on the paths defined in their system#
    Raw_Data_Path="/Users/njain/Desktop/$DatasetName"
    cd ${Raw_Data_Path}

    for run in *; do
        regx='^[0-9]+$'           # Regular Expression to check for numerics
        if [[ $run =~ $regx ]]; then  # Check if folder name is numeric
            echo "Extracting Parameters from Run Number: E$run"
            Raw_Data_Path_Run=$Raw_Data_Path/$run
        

            ## Extracting Relevant Parameters

            Sequence="`grep -A1 'ACQ_protocol_name=( 64 )' $Raw_Data_Path_Run/acqp | grep -v -e 'ACQ_protocol_name=( 64 )' -e '--'`"
            #here the grep command will read the no of slices from the Bruker raw file acquired during acquisiton
            SequenceName=$(echo "$Sequence" | sed 's/[<>]//g')
            

            Field="`grep -A1 '##$PVM_Fov' $Raw_Data_Path_Run/method | grep -v -e '##$PVM_Fov' -e '--' | head -n 1`"
            FieldOfView=$(echo "$Field" | sed 's/[<>]//g')

            Matrix="`grep -A1 '##$PVM_Matrix' $Raw_Data_Path_Run/method | grep -v -e '##$PVM_Matrix' -e '--'`"
            MatrixSize=$(echo "$Matrix" | sed 's/[<>]//g')

            Orient="`grep -A1 '##$PVM_SPackArrSliceOrient' $Raw_Data_Path_Run/method | grep -v -e '##$PVM_SPackArrSliceOrient' -e '--'`"
            Orientation=$(echo "$Orient" | sed 's/[<>]//g')

              
            echo "E$run\t$SequenceName\t$FieldOfView\t$MatrixSize\t$Orientation" >> Parameters.txt
        else    
            echo "It is a Bruker Scan related file"
        fi
        
    done
done
