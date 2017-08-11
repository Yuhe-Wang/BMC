#!/bin/sh

cd RunBMC
rm *.so*
rm *.app
rm *.encrypt
cd .. # to BMC
cd ..
rm -r Build
