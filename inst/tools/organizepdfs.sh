#!/bin/bash
# organizepdfs.sh
# Organizes a folder of SPADE output PDFs into subfolders for each parameter
#
# Instructions:  Copy this script to the "pdf" folder in a SPADE output folder, and run it.
#
# Version 0.1.0
# Erin Simonds
# December 14, 2010

(
IFS=$'\n'
filenames=(`ls -1 *.pdf`)


echo "=== Unique parameters: ==="
uniqueparams=(`printf "%s\n" "${filenames[@]}" | sed 's/.*\.gml\.median//g' | sed 's/.*\.gml\.fold//g' | sed 's/\_clust//g' | sed 's/.*\.gml\.//g' | sed 's/\.pdf//g' |sort -u`)
printf "%s\n" "${uniqueparams[@]}"

echo "=== Making folder for each parameter: ==="
mkdir ../pdf_organized_by_parameter
COUNT="0"
for i in ${uniqueparams[*]}; do
    mkdir -v ../pdf_organized_by_parameter/$i
    let COUNT++
done

echo "=== Copying each PDF to the appropriate parameter folder ==="
COUNT="0"
for i in ${uniqueparams[*]}; do
    cp -v *$i*\.pdf ../pdf_organized_by_parameter/$i
    let COUNT++
done

echo "=== Unique FCS files: ==="
uniquefcs=(`printf "%s\n" "${filenames[@]}" | sed 's/\.fcs.*//g' | sort -u`)
printf "%s\n" "${uniquefcs[@]}"

echo "=== Making folder for each FCS file: ==="
mkdir ../pdf_organized_by_FCS
COUNT="0"
for i in ${uniquefcs[*]}; do
    mkdir -v ../pdf_organized_by_FCS/$i
    let COUNT++
done

echo "=== Copying each PDF to the appropriate FCS file folder ==="
COUNT="0"
for i in ${uniquefcs[*]}; do
    cp -v $i*\.pdf ../pdf_organized_by_FCS/$i
    let COUNT++
done

echo "=== Finished ==="
)