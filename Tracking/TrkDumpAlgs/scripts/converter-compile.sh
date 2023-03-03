source athena/Tracking/TrkDumpAlgs/scripts/setup_cvmfs_lcg97.sh
cmake -Hathena/Tracking/TrkDumpAlgs/scripts/ -Bbuild -Sathena/Tracking/TrkDumpAlgs/scripts -DCMAKE_CXX_STANDARD=17
cmake --build build -- -j3
cp build/libDumper_rdict.pcm bin/libDumper_rdict.pcm
