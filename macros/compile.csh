#!/bin/tcsh

# builds libraries from scratch
#
# Currently builds:
# 1) StEfficiencyAssessor

echo "[i] loading embedding library"
starver pro

echo "[i] Remove any existing libs"
rm -v libStEfficiencyAssessor.so
rm -v libs/libStEfficiencyAssessor.so
rm -v sandbox/libStEfficiencyAssessor.so
rm -v libStRefMultCorr.so
rm -v libs/libStRefMultCorr.so
rm -v sandbox/libStRefMultCorr.so

setenv CXXFLAGSNEW "-pipe -fPIC -Wall -Woverloaded-virtual -ansi -Wno-long-long -pthread -m32 -std=c++11"
#LDFLAGS       += -m32

echo "[i] Changing CXXFLAGS to: "${CXXFLAGSNEW}

echo "[i] Running cons for StEfficiencyAssessor"

cons CXXFLAGS="${CXXFLAGSNEW}" +StEfficiencyAssessor +StRefMultCorr


# places copies of the libraries into the local lib directory
# as well as into sandbox/
echo "[i] Copying libraries to the lib & sandbox"
find .sl*/lib -name "libStEfficiencyAssessor.so" -exec cp -v {} ./libs/ \;
find .sl*/lib -name "libStEfficiencyAssessor.so" -exec cp -v {} ./sandbox/ \;
find .sl*/lib -name "libStRefMultCorr.so" -exec cp -v {} ./libs/ \;
find .sl*/lib -name "libStRefMultCorr.so" -exec cp -v {} ./sandbox/ \;


