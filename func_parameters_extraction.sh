#!/bin/sh


#07th November 2024: $$Naman Jain$$ This function is to extract parameters from the methods file

FUNC_PARAM_EXTARCT () {

    Sequence="`grep -A1 'ACQ_protocol_name=( 64 )' $1/acqp | grep -v -e 'ACQ_protocol_name=( 64 )' -e '--'`"
    #here the grep command will read the no of slices from the Bruker raw file acquired during acquisiton
    SequenceName=$(echo "$Sequence" | sed 's/[<>]//g')
    
    # NoOfRepetitions=$(awk '/PVM_NRepetitions=/ {print substr($0,21,4)}' $1/method)
    NoOfRepetitions=$(awk -F'=' '/##\$PVM_NRepetitions=/{print $2}' $1/method)
    TotalScanTime=$(awk -F'=' '/##\$PVM_ScanTime=/{print $2}' $1/method)

    # TotalScanTime=$(awk '/PVM_ScanTime=/ {print substr($0,17,6)}' $1/method)
    #here the awk will look at the number of slices acquired using the information located in the methods file    
                    
    # 07.08.2024 Estimating Volume TR
    VolTR_msec=$(echo "scale=0; $TotalScanTime/$NoOfRepetitions" | bc)
    VolTR=$(echo "scale=0; $VolTR_msec/1000" | bc)
        
    Baseline_TRs=$(awk '/PreBaselineNum=/ {print substr($0,19,3)}' $1/method)
    StimOn_TRs=$(awk '/StimNum=/ {print substr($0,12,2); exit}' $1/method)  
    StimOff_TRs=$(awk '/InterStimNum=/ {print substr($0,17)}' $1/method)  
    NoOfEpochs=$(awk '/NEpochs=/ {print substr($0,12)}' $1/method)  

       
    # TaskDuration=$(echo "$Baseline_TRs + ($StimOn_TRs + $StimOff_TRs) * $NoOfEpochs" | bc)
    
    # BlockLength=$(($StimOn_TRs + $StimOff_TRs))
    MiddleVolume=$(($NoOfRepetitions / 2))

    
    export SequenceName NoOfRepetitions TotalScanTime VolTR_msec VolTR
    export Baseline_TRs StimOn_TRs StimOff_TRs NoOfEpochs TaskDuration
    export BlockLength MiddleVolume

}                    