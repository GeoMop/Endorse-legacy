#!/bin/bash
root=$1
find $root -name "FINISHED" -exec rm -rf {} \;
