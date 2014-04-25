#ifndef parsers_H
#define parsers_H

using namespace std;

typedef vector< string > SampleNames;

struct FastaReferenceData{
    string name;
    uint32_t length; 
    uint64_t end; //endpoint relative to entire length of whole ref
};

struct BedInterval{ // make clas?
    string chr;
    uint64_t start;
    uint64_t end;
    bool contains(int pos);
    void set_open();
};

typedef vector<FastaReferenceData> FastaReferenceVector;

class FastaReference{
        //string ref_file_name;
    public:
        FastaReference(string ref_file_name);
        FastaReferenceVector chromosomes;
        void get_ref_id(string name, int& chr_id);
    private:
        ifstream ref_file;
};

class BedFile{
//        string bed_file_name;
    public:
        BedFile(string bed_file_name);
        int get_interval(BedInterval& current_interval);
        ifstream bed_file;
        
};
//Helper functions

//operations on genomic intervals

uint32_t genomic_distance(const BedInterval& a, const BedInterval& b){
    if(a.chr != b.chr){
       return -1;
    }
    else{
        //onyl for zero base indices!!
        int d =   max(a.start, b.start) - min(a.end , b.end);
        if(d < 0){
            return -1;
        }
        return d;        
    }
}

uint32_t genomic_overlap(const BedInterval& a, const BedInterval& b){
    if(a.chr != b.chr){ 
        return 0;
    }
    else{
       int overlap = min(a.end - b.end) - max(a.start - b.start);
       if (overlap < 0){
           return 0;
       }
       return overlap;
    }
}

inline bool operator==(const BedInterval& a, const BedInterval& b){ 
    if (a.chr == b.chr){
        if (a.start == a.start){
            return a.end == b.end;
        }
    }
    return false;
}












        


uint16_t base_index(char b);
uint32_t find_sample_index(string, SampleNames);

#endif



