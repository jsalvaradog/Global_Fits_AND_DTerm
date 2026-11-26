
#!/usr/bin/env bash

set -euo pipefail

TEMPLATE="computeManyKinematicsForDVCSComptonFormFactor.xml"
DIR="CFF_NN_Predictions"
XML_DIR="CFF_NN_Predictions/XML"
LOG_DIR="CFF_NN_Predictions/Logs"
OUT_DIR="CFF_NN_Predictions"
NREPS=101

mkdir -p "$XML_DIR" "$LOG_DIR" "$OUT_DIR"
mkdir -p "$DIR/Replicas"

for i in $(seq 0 100); do
    echo "Processing replica $i..."

    XML_FILE="${XML_DIR}/computeManyKinematicsForDVCSComptonFormFactor_${i}.xml"
    LOG_FILE="${LOG_DIR}/Rep_${i}.log"

    # Create the XML file for this replica
    sed "s/RepID/${i}/g" "$TEMPLATE" > "$XML_FILE"

    # Run PARTONS
    ./bin/PARTONS_example "$XML_FILE" > "$LOG_FILE" 2>&1

    # Create eight output files for this replica
    for cff in H E Ht Et; do
        for part in Re Im; do
            OUT_FILE="${OUT_DIR}/Replicas/${part}_${cff}_${i}.txt"

            awk -v cff="$cff" -v part="$part" '
                $1 == "CFF" && $2 == cff":" {
                    for (j = 1; j <= NF; j++) {
                        if ($j == part ":") {
                            value = $(j+1)
                            gsub(/,/, "", value)
                            print value
                        }
                    }
                }
            ' "$LOG_FILE" > "$OUT_FILE"
        done
    done
done

echo "All replicas processed."



for cff in H E Ht Et; do
    for part in Re Im; do

        output_mean="${DIR}/${part}${cff}_Mean.txt"
        output_std="${DIR}/${part}${cff}_StdDev.txt"

        awk -v nreps="$NREPS" '
            {
                x = $1
                sum[FNR] += x
                sumsq[FNR] += x * x
                count[FNR]++
            }
            END {
                for (i = 1; i <= FNR; i++) {
                    if (count[i] != nreps) {
                        printf "Error: line %d has %d values; expected %d\n", \
                               i, count[i], nreps > "/dev/stderr"
                        exit 1
                    }

                    mean = sum[i] / nreps
                    variance = sumsq[i] / nreps - mean * mean

                    # Protect against tiny negative roundoff
                    if (variance < 0 && variance > -1e-12)
                        variance = 0

                    if (variance < 0) {
                        printf "Error: negative variance at line %d\n", \
                               i > "/dev/stderr"
                        exit 1
                    }

                    printf "%.15g\n", mean > meanfile
                    printf "%.15g\n", sqrt(variance) > stdfile
                }
            }
        ' meanfile="$output_mean" stdfile="$output_std" \
          "${DIR}/Replicas/${part}_${cff}_"{0..100}".txt"

        echo "Created $output_mean and $output_std"
    done
done

echo "Means and standard deviations computed."


