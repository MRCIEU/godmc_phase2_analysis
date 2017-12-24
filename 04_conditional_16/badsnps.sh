#!/bin/bash

cat bad/* | cut -f 1 | sort -u > badsnps.txt

