#!/bin/sh


#06th June 2024: $$Naman Jain$$ This function is to check if a folder exists or not. If it doesn't exist
#                               then it creates a new and stores file in it.


#Function 1
CHECK_FILE_EXISTENCE () {
        if [ -d $1 ]; then
            echo "Analysed Data folder exists, Proceeding to Analyse the data"
            return 1
        else 
            mkdir $1
            return 0
        fi

}
