#!/bin/bash

path=${1}
subfolder=${2}

cd ${path}

foldersToCreate=(Results Results/Opening Results/Pressure)

sectionCreate(){

    folderList=$(find . -type d | sed -e "s/\.\//"\$INSTDIR"\/${subfolder}//g" | tr '/' '\' | awk '{print "  SetOutPath \"" $1 "\"\n  File \"" $1 "\\*\""}' | sed -e "s/File \"\$INSTDIR\\${subfolder}/File \"Builds\\${subfolder}/g")

    numberOfFolders=$(($(echo "${folderList}" | wc -l) - 2))

    for (( i = 0; i < ${#foldersToCreate[@]}; ++i )); do
        echo "  CreateDirectory \"\$INSTDIR\\${subfolder}\\${foldersToCreate[i]}\"" | tr '/' '\'
    done

    echo ""

    echo "  SetOutPath \"\$INSTDIR\\${subfolder}\""
    echo "  File \"Builds\\${subfolder}\*"
    echo "${folderList}" | tail -n ${numberOfFolders}
}

sectionDelete(){
    folderList=$(find . -type d | sed -e "s/\.\//"\$INSTDIR"\//g" | tr '/' '\' | awk '{print "  Delete \"" $1 "\\*\""}')
    numberOfFolders=$(($(echo "${folderList}" | wc -l) - 1))

    for (( i = ${#foldersToCreate[@]} - 1; i >= 0; --i )); do
        echo "  Delete \"\$INSTDIR\\${foldersToCreate[i]}\*\"" | tr '/' '\'
    done

    echo "${folderList}" | tail -r | head -n ${numberOfFolders}
    echo "  Delete \"\$INSTDIR\\*\""

    for (( i = ${#foldersToCreate[@]} - 1; i >= 0; --i )); do
        echo "  RMDir /r \"\$INSTDIR\\${foldersToCreate[i]}\"" | tr '/' '\'
    done

    echo "${folderList}" | sed -e "s/Delete/RMDir \/r/g"  | rev | cut -c 4- | rev | tail -r | head -n ${numberOfFolders} | awk '{print "  " $1 " " $2 " " $3 "\""}'
    echo "  RMDir /r \"\$INSTDIR\""
}

sectionCreate
