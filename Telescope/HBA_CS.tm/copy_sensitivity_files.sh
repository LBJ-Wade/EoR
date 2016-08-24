#!/bin/bash

# Copy into folders
for i in $(ls -d CS*/); do
    #echo "cp sensitivity_CS.txt ${i%%/}/.";
    cp sensitivity_CS.txt ${i%%/}/sensitivity.txt
    #echo "rm -f ${i%%/}/sensitivity.txt.";
    #rm -f ${i%%/}/sensitivity_CS.txt
    #rm -f ${i%%/}/sensitivity.txt.
done