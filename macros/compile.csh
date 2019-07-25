#!/bin/tcsh

# builds libraries from scratch
#
# Currently builds:
# 1) StEfficiencyAssessor

echo "[i] loading embedding library"
starver SL18h

echo "[i] Remove any existing libs"
rm -v libStEfficiencyAssessor.so
rm -v StEfficiencyAssessor.so
rm -v libs/libStEfficiencyAssessor.so
rm -v libs/StEfficiencyAssessor.so
rm -v sandbox/libStEfficiencyAssessor.so
rm -v sandbox/StEfficiencyAssessor.so

setenv CXXFLAGSNEW "-pipe -fPIC -Wall -Woverloaded-virtual -ansi -Wno-long-long -pthread -m32 -std=c++11"
#LDFLAGS       += -m32

echo "[i] Changing CXXFLAGS to: "${CXXFLAGSNEW}

echo "[i] Running cons for StEfficiencyAssessor"

cons CXXFLAGS="${CXXFLAGSNEW}" +StEfficiencyAssessor


# places copies of the libraries into the local lib directory
# as well as into sandbox/
echo "[i] Copying libraries to the lib & sandbox"
find .sl*/lib -name "libStEfficiencyAssessor.so" -exec cp -v {} ./libs/ \;
find .sl*/lib -name "libStEfficiencyAssessor.so" -exec cp -v {} ./sandbox/ \;


