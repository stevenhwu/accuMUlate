#pragma once
#ifndef parsers_H
#define parsers_H

//#include "utils/bamtools_pileup_engine.h"
#include "local_bamtools/bamtools_pileup_engine.h"
#include "local_bamtools/bamtools_fasta.h"

#include <unordered_map>
#include <iostream>
#include <string>

//using namespace std;
//typedef vector< string > SampleNames;
typedef std::unordered_map<std::string, uint16_t> SampleMap;

struct FastaReferenceData{
    std::string name;
    uint32_t length; 
    uint64_t end; //endpoint relative to entire length of whole ref
};

struct BedInterval{
    std::string chr;
    uint64_t start;
    uint64_t end;
};

typedef std::vector<FastaReferenceData> FastaReferenceVector;

class FastaReference{
        //string ref_file_name;
    public:
        FastaReference(std::string ref_file_name);
        FastaReferenceVector chromosomes;
        void get_ref_id(std::string name, int& chr_id);
    private:
    std::ifstream ref_file;
};

class BedFile{
//        string bed_file_name;
    public:
        BedFile(std::string bed_file_name);
        int get_interval(BedInterval& current_interval);
    std::ifstream bed_file;
        
};
//Helper functions

bool include_site(LocalBamToolsUtils::PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut);
/* include_site vs include_site_3(4) and with lots of cache (V2)
Total time V1: 8.884579 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 6.720386 Count:100000 G_count:6927746 Base_count:16143
Total time V1: 8.730865 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 6.683539 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 6.702135 Count:100000 G_count:6927746 Base_count:16143
Total time V1: 8.788570 Count:100000 G_count:6927746 Base_count:16143
Total time V1: 8.855847 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 7.087987 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 6.983055 Count:100000 G_count:6927746 Base_count:16143
Total time V1: 8.986538 Count:100000 G_count:6927746 Base_count:16143
Total time V1: 9.053064 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 6.796710 Count:100000 G_count:6927746 Base_count:16143
Total time V2: 6.986318 Count:100000 G_count:6927746 Base_count:16143
Total time V1: 9.067891 Count:100000 G_count:6927746 Base_count:16143

 */

bool include_site_2(const LocalBamToolsUtils::PileupAlignment & pileup, uint16_t map_cut, uint16_t qual_cut);
bool include_site_3(const LocalBamToolsUtils::PileupAlignment & pileup, uint16_t map_cut, char qual_cut);
//using namespace BamTools;
bool include_site_4(const BamTools::BamAlignment & alignment, const int &pos, const uint16_t &map_cut, const char &qual_cut);
uint16_t base_index(char b);
std::string get_sample(std::string& tag);

extern int base_index2[128];

//uint32_t find_sample_index(string, SampleNames);

#endif



