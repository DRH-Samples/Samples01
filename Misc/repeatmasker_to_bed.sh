#!/bin/env bash
awk '{ OFS="\t"; sub(/C/, "-", $9); sub(/\(/, "", $10); sub(/\)/, "_", $10); print $5, $6-1, $7, $11"|"$10, $2*10, $9 }' $1 | tail -n +4
