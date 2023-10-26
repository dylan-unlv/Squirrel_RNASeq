#!/bin/bash

fastqc raw_data/*/*.fq.gz -o data/fastqc/ -t 25
