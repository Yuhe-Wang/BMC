#!/bin/sh

cd RunBMC
rm *.so*
rm *.app
rm *.encrypt
rm -r Materials
cd .. # to BMC
cd ..
rm -r Build
cd BMC
