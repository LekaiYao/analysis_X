#!/bin/bash

mkdir -p FOM
> FOM/FOM_summary.txt  

TEMPLATE="optimization_template.C"

cat input_BDT.txt | while read rows
do
	var=$(echo $rows | awk 'BEGIN{FS="/"} {print $1}')
	min=$(echo $rows | awk 'BEGIN{FS="/"} {print $2}')
	max=$(echo $rows | awk 'BEGIN{FS="/"} {print $3}')
	step=$(echo $rows | awk 'BEGIN{FS="/"} {print $4}')
    script_name="optimization_${var}.C"
    cp $TEMPLATE $script_name

    sed -i -e "s/@VAR@/${var}/g" \
        -e "s/@MIN@/${min}/g" \
        -e "s/@MAX@/${max}/g" \
        -e "s/@STEP@/${step}/g" $script_name
    echo "Running $script_name..."
    root -l -b -q "$script_name"
    rm "$script_name"
done

