/*
  Copyright (C) 2002-2020 CERN for the benefit of the ATLAS collaboration
*/

/********************************************************
 Class def for MuonGeoModel DblQ00/ASZT
 *******************************************************/

 //  author: S Spagnolo
 // entered: 07/28/04
 // comment: MUON STATION ELEMENT

#ifndef DBLQ00_ASZT_H
#define DBLQ00_ASZT_H



class IRDBAccessSvc;
#include <string>

class AmdcDb;

namespace MuonGM {
class DblQ00Aszt {
public:
    DblQ00Aszt();
    ~DblQ00Aszt();
    DblQ00Aszt(IRDBAccessSvc *pAccessSvc, const std::string & GeoTag="", const std::string & GeoNode="");
    DblQ00Aszt(const std::string& asciiFileName);
    DblQ00Aszt(AmdcDb* aszt);
    
    void WriteAsztToAsciiFile(const std::string& filename);

    // data members for DblQ00/ASZT fields
    struct ASZT {
        int version; // VERSION
        int line; // LINE NUMBER
        char type[8]; // STATION TYPE
        int jff; // PHI POSITION
        int jzz; // Z POSITION
        int job; // JOB POSITION
        float tras; // S TRANSLATION [MM]
        float traz; // Z TRANSLATION
        float trat; // T TRANSLATION
        float rots; // S ROTATION [RAD]
        float rotz; // Z ROTATION
        float rott; // T ROTATION
        int i; // STATION AMDB INDEX
    };

    const ASZT* data() const { return m_d; };
    unsigned int size() const { return m_nObj; };
    const char* getName() const { return "ASZT"; };
    const char* getDirName() const { return "DblQ00"; };
    const char* getObjName() const { return "ASZT"; };

private:
    ASZT* m_d;
    unsigned int m_nObj; // > 1 if array; 0 if error in retrieve.
    DblQ00Aszt & operator=(const DblQ00Aszt &right);
    DblQ00Aszt(const DblQ00Aszt&);
};


} // end of MuonGM namespace

#endif // DBLQ00_ASZT_H

