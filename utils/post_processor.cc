#include <iostream>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <vector>
#include <numeric> 

#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "boost/program_options.hpp"

#include "model.h"
#include "parsers.h"

using namespace std;
using namespace BamTools;


class SiteData{

    public:
        string name;
        vector<uint16_t> BQs;
        vector<uint16_t> MQs;
        uint32_t fwd_reads;
        uint32_t rev_reads;
        ReadData base_calls;

        SiteData(string sname){
          name = sname;
          base_calls.key = 0;
        }
                  
        int get_genotype(){
           //These have already been called for mutation-ness, and are haploid
           //so, to make a start, we are just calling the most common base
           return distance(base_calls.reads, max_element(base_calls.reads, base_calls.reads + 4 ));
        }
        void import_alignment(const BamAlignment& al, const int& pos, const int& bindex){
            uint16_t b_index = base_index(al.AlignedBases[pos]);
            base_calls.reads[bindex] += 1;
            BQs.push_back(al.Qualities[pos]);
            MQs.push_back(al.MapQuality);
            if(al.IsReverseStrand()){ 
                rev_reads += 1; 
            } 
            else{ 
                fwd_reads += 1; 
            }
        }
};


 

typedef vector<SiteData> SiteDataVector;

class FilterVisitor: public PileupVisitor{
    public: 
        FilterVisitor(BamAlignment& ali, 
                      ostream *out_stream,
                      SampleNames samples, 
                      int ref_pos,
                      string initial_data):

            PileupVisitor(), m_samples(samples), m_ostream(out_stream), 
                             m_ref_pos(ref_pos), m_initial_data(initial_data){
                nsamp = m_samples.size(); 
            } 
        ~FilterVisitor(void) { }


    public:
        void Visit(const PileupPosition& pileupData){
            SiteDataVector sample_data;
            for (size_t i = 0; i < nsamp; i++){
                sample_data.push_back(SiteData(m_samples[i]));
            }
            for (auto it =  pileupData.PileupAlignments.begin();
                      it != pileupData.PileupAlignments.end();
                      it++){
                if(pileupData.Position == m_ref_pos){
                    int const *pos = &it->PositionInAlignment;
                    if(it->Alignment.Qualities[*pos] > 46){//TODO user-defined qual cut to match first call?
                        uint16_t b_index = base_index(it->Alignment.AlignedBases[*pos]);
                        if (b_index < 4){
                            string tag_id;
                            it->Alignment.GetTag("RG", tag_id);
                            uint32_t sindex = find_sample_index(get_sample(tag_id), m_samples);
                            sample_data[sindex].import_alignment(it->Alignment, *pos, b_index);                        
                        }
                    }
                }
            }

            vector<uint16_t> genotypes;
            int gfreqs[4];
            for(size_t i=0; i < nsamp; i++){
                uint16_t g = sample_data[i].get_genotype();
                gfreqs[g] += 1;
                genotypes.push_back(g);                        
            }
            if(count(gfreqs, gfreqs+4, 1) == 1){//should only be one mutant
                    // OK let's collect all that data....
                    auto it = find_if(genotypes.begin(), genotypes.end(),
                        [](int v) {return v==1;});
                    unsigned mutant = *it;
                    uint16_t mutant_base = genotypes[mutant];
                    int mutant_alleles;
                    int mutant_allele_denom;
                    int wt_MQs;
                    int wt_MQ_denom;
                    int mutant_MQs =  accumulate(sample_data[mutant].MQs.begin(), sample_data[mutant].MQs.end(), 0) /
                                      sample_data[mutant].MQs.size();
                    for (size_t i=0; i < nsamp; i++){
                        if(i != mutant){
                            mutant_alleles += sample_data[i].base_calls.reads[mutant_base];
                            mutant_allele_denom += sample_data[i].fwd_reads + sample_data[i].rev_reads;
                            wt_MQs += accumulate(sample_data[i].MQs.begin(), sample_data[i].MQs.end(), 0); 
                            wt_MQ_denom += sample_data[i].MQs.size();     
                        }
                        
                    }
                    double xbar_MQs = double(wt_MQs/wt_MQ_denom);
                    double mutant_freq = (double)mutant_alleles/mutant_allele_denom;
                    *m_ostream << m_initial_data << '\t'
                               << mutant_base << '\t'
                               << sample_data[mutant].name << '\t'
                               << mutant_freq << '\t' 
                               << sample_data[mutant].fwd_reads << '\t'
                               << sample_data[mutant].rev_reads << '\t'
                               << mutant_MQs << '\t'
                               << xbar_MQs  << endl;

            }
        }

    private:
        SampleNames m_samples;
        int m_ref_pos;
        ostream* m_ostream;
        int nsamp;
        string m_initial_data;
    
};
                                 

int main(int argc, char* argv[]){
    string bam_path;
    string input_path;
    string ref_path;
    namespace po = boost::program_options;
    po::options_description cmd("Command line args");
    cmd.add_options()
        ("help,h", "Print a help message")
        ("bam,b", po::value<string>(&bam_path)->required(), "Path to BAM file")
        ("reference,r", po::value<string>(&ref_path)->required(), "Path to reference genome")
        ("input,i", po::value<string>(&input_path)->required(), "Path to results file")
        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags")
        ("config,c", po::value<string>(), "Path to config file")
        ("out,o", po::value<string>()->default_value("filtered_result.tsv"),
                    "Out file name");

 
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);

    if (vm.count("help")){
        cout << cmd <<endl;
        return 0;
    }

    if (vm.count("config")){
        ifstream config_stream(vm["config"].as<string>());
        po::store(po::parse_config_file(config_stream, cmd, false), vm);
    }

    vm.notify();

    ofstream outfile (vm["out"].as<string>());
    ifstream putations(input_path);
    BamReader experiment;
    experiment.Open(bam_path);
    //FastaReference ref_genome (ref_path + ".fai");
    string L;
    PileupEngine pileup;
    BamAlignment ali;

    while(getline(putations, L)){    

        // chr \t pos \t ref \t prob
        size_t i = 0;
        string chr;
        for(; L[i] != '\t'; ++i){
            chr.push_back(L[i]);
        }
        string pos_s;
        i += 1;
        for(; L[i] !='\t'; i++){
            pos_s.push_back(L[i]);
        }
        int pos = stoul(pos_s);
        int ref_id = experiment.GetReferenceID(chr);
        experiment.SetRegion(ref_id, pos, ref_id, pos+1);
        FilterVisitor *f = new FilterVisitor(ali, 
                                             &outfile,
                                             vm["sample-name"].as<vector< string> >(), 
                                             pos, L);
        pileup.AddVisitor(f);
        while( experiment.GetNextAlignment(ali) ) {
            pileup.AddAlignment(ali);
        }
    }
    pileup.Flush();
    return 0;
}

        


            






