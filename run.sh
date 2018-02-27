#!/bin/sh

make
echo rm -rf data.xyz
rm -rf data.xyz
./Quake
