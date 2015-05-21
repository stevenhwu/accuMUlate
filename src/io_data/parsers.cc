#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>

#include "parsers.h"



using namespace BamTools;

FastaReference::FastaReference(std::string ref_file_name){
    std::ifstream ref_file (ref_file_name);
    std::string L;
    uint64_t cummulative_len = 0;
    while( getline(ref_file, L)){
        size_t i = 0;
        std::string chrom;
        for(; L[i] != '\t'; i++){
            chrom.push_back(L[i]);
        }
        i += 1;
        std::string s_len;
        for (; L[i] !='\t'; i++){
            s_len.push_back(L[i]);
        }
        uint32_t chrom_len = std::stoul(s_len);
        cummulative_len += chrom_len;
        chromosomes.push_back(FastaReferenceData{ chrom, 
                                                  chrom_len, 
                                                  cummulative_len
                                                 });
    }
}
    
void FastaReference::get_ref_id(std::string search_name, int& chr_id){
    auto it = find_if(chromosomes.begin(), chromosomes.end(),
                     [&search_name](const FastaReferenceData& chrom){
                        return chrom.name == search_name;
                      });
    int idx  = distance(chromosomes.begin(), it);
    chr_id = idx;
}


BedFile::BedFile(std::string bed_file_name){
     bed_file.open(bed_file_name);
}

int BedFile::get_interval(BedInterval& current_interval){
    std::string L;
    if(getline(bed_file, L)){
        std::stringstream Lstream(L);
        std::string chrom;
        std::string start_s;
        std::string end_s;
        getline(Lstream, chrom, '\t');
        getline(Lstream, start_s, '\t');
        getline(Lstream, end_s, '\t');
        current_interval = BedInterval{ chrom,
                std::stoul(start_s),
                std::stoul(end_s)};
        return 0;
    }
   else {return 1;}
    
}

//
//Helper functions for parsing data out of BAMs

uint16_t base_index(char b){
    switch(b){        
        case 'A':
        case 'a':    
            return 0;
         case 'T':
         case 't':
             return 3;
         case 'C':
         case 'c':
             return 1;
         case 'G':
         case 'g':
             return 2;
         case '-':
         case 'N':
             return -1 ;
         default: // Unknown base, alert in debug mode?
             return -1;
    }
}


//uint32_t find_sample_index(string s, SampleNames sv){
//    for (size_t i=0; i < sv.size(); i++){
//        if(s.compare(sv[i])==0){
//            return i;
//        }
//    }
//    return(-1); //TODO refactor this to  update sample in place
//}

bool include_site(LocalBamToolsUtils::PileupAlignment pileup, uint16_t map_cut, uint16_t qual_cut){
    const BamAlignment *ali = &pileup.Alignment;
    if(ali->MapQuality > map_cut){
        uint16_t bqual = static_cast<short>(ali->Qualities[pileup.PositionInAlignment]) - 33;
        if(bqual > qual_cut){
            return(not (ali->IsDuplicate()) && not(ali->IsFailedQC()) && ali->IsPrimaryAlignment());
        }
    }
    return false;
}


bool include_site_2(const LocalBamToolsUtils::PileupAlignment & pileup, uint16_t map_cut, uint16_t qual_cut){
    const BamAlignment *ali = &(pileup.Alignment);
    if(ali->MapQuality > map_cut){
        uint16_t bqual = static_cast<short>(ali->Qualities[pileup.PositionInAlignment]) - 33;
        if(bqual > qual_cut){
            return(not (ali->IsDuplicate()) && not(ali->IsFailedQC()) && ali->IsPrimaryAlignment());
        }
    }

    return false;
}



bool include_site_3(const LocalBamToolsUtils::PileupAlignment & pileup, uint16_t map_cut, char qual_cut){
//    const BamAlignment *ali = &(pileup.Alignment);
    if(pileup.Alignment.MapQuality > map_cut){
        char reference = pileup.Alignment.Qualities[pileup.PositionInAlignment];

//        uint16_t bqual = static_cast<short>(ali->Qualities[pileup.PositionInAlignment]) - 33;
//
//        char rc = (char)(qual_cut+33);
//        if(  (bqual > qual_cut) != (reference > rc)    ){
//            cout << reference << "\t" <<
//                    (static_cast<short>(ali->Qualities[pileup.PositionInAlignment]) - 33) << "\t" << qual_cut << "\t" <<
//                    ((char)(qual_cut+33)) <<  (bqual > qual_cut) << "\t" <<  (reference > rc) << endl;
//        }
        if(reference > qual_cut){
            return(not (pileup.Alignment.IsDuplicate()) && not(pileup.Alignment.IsFailedQC()) && pileup.Alignment.IsPrimaryAlignment());
        }
    }

    return false;
}


bool include_site_4(const BamAlignment & alignment, const int &pos, const uint16_t &map_cut, const char &qual_cut){
//    const BamAlignment *ali = &(pileup.Alignment);
    if(alignment.MapQuality > map_cut){
        char reference = alignment.Qualities[pos];

        if(reference > qual_cut){
            return(not (alignment.IsDuplicate()) && not(alignment.IsFailedQC()) && alignment.IsPrimaryAlignment());
        }
    }

    return false;
}






//int main() {
//    FastaReference reference_g ("test/test.fai");
//    int chr_idx;
//    reference_g.get_ref_id("scf_8254727", chr_idx);
//    cout << "chr index = " << chr_idx << " (should be 9)" << endl;
//
//    BedFile bed ("test/test.bed");
//    BedInterval current_line;
//    while(bed.get_interval(current_line) == 0){
//         cout << current_line.chr << endl;
//    }
//    return 0;
//}
//







int base_index2[128] ={
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 0-15
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 16-31
////                                          -
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,	// 32-47
////                                                ?
//    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,	// 48-63
////	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
//    17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 64-79
////	 p  q  R  S  T  U  V  W  x  Y  z
//    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17,	// 80-95
////	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
//    17, 0,11, 1,12,16,16, 2,13,16,16,10,16, 7,15,16,	// 96-111
////	 p  q  R  S  T  U  V  W  x  Y  z
//    16,16, 5, 9, 3, 3,14, 8,16, 6,16,17,17,17,17,17		// 112-127

        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 0-15
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 16-31
//                                          -
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 32-47
//                                                ?
        -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 48-63
//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,    // 64-79
//	 p  q  R  S  T  U  V  W  x  Y  z
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,	// 80-95
//	    A  B  C  D  e  f  G  H  i  j  K  l  M  N  o
        -1, 0,-1, 1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,   	// 96-111
//	 p  q  R  S  T  U  V  W  x  Y  z
        -1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1		// 112-127
};


