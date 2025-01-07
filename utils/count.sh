#!/usr/bin/env bash
# @AUTHOR: Chun-Jie Liu
# @CONTACT: chunjie.sam.liu.at.gmail.com
# @DATE: 2024-07-18 16:50:51
# @DESCRIPTION:

# Number of input parameters
param=$#
inputfile=$1

awk '
BEGIN { OFS = "\t" }
{
    pos = $2
    bases = $5
    a_count = g_count = c_count = t_count = 0

    for (i = 1; i <= length(bases); i++) {
        base = substr(bases, i, 1)
        if (base == "A" || base == "a") a_count++
        else if (base == "G" || base == "g") g_count++
        else if (base == "C" || base == "c") c_count++
        else if (base == "T" || base == "t") t_count++
    }

    print pos, a_count, g_count, c_count, t_count
}
' $inputfile
