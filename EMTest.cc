#include <iostream>
#include <stdint.h>
#include <map>
#include <vector>
#include <string>

#include <boost/program_options.hpp>
#include "api/BamReader.h"
#include "utils/bamtools_pileup_engine.h"
#include "utils/bamtools_fasta.h"

#include "model.h"
#include "parsers.h"
#include "SequenceProb.h"
#include "MutationProb.h"

#include <algorithm>

using namespace std;
using namespace BamTools;


class VariantVisitor : public PileupVisitor{
public:
    VariantVisitor(const RefVector& bam_references,
            const SamHeader& header,
            const Fasta& idx_ref,
            GenomeData& all_the_data,
//                      const ModelParams& p,  
            const SampleMap& samples,
            BamAlignment& ali,
            int qual_cut,
            int mapping_cut,
            double prob_cut):

            PileupVisitor(), m_idx_ref(idx_ref), m_bam_ref(bam_references),
            m_header(header), m_samples(samples),
            m_qual_cut(qual_cut), m_ali(ali),
            m_all_the_data(all_the_data), m_prob_cut(prob_cut),
            m_mapping_cut(mapping_cut)
    { }
    ~VariantVisitor(void) { }
public:
    void Visit(const PileupPosition& pileupData) {
        string chr = m_bam_ref[pileupData.RefId].RefName;
        uint64_t pos  = pileupData.Position;
        m_idx_ref.GetBase(pileupData.RefId, pos, current_base);
        ReadDataVector bcalls (m_samples.size(), ReadData{{ 0,0,0,0 }});
        string tag_id;
        for(auto it = begin(pileupData.PileupAlignments);
            it !=  end(pileupData.PileupAlignments);
            ++it){
            if( include_site(*it, m_mapping_cut, m_qual_cut) ){
                it->Alignment.GetTag("RG", tag_id);
                string sm =  m_header.ReadGroups[tag_id].Sample;
                uint32_t sindex = m_samples[sm]; //TODO check samples existed!
                uint16_t bindex  = base_index(it->Alignment.QueryBases[it->PositionInAlignment]);
                if (bindex < 4 ){
                    bcalls[sindex].reads[bindex] += 1;
                }
            }
        }
        uint16_t ref_base_idx = base_index(current_base);
        if (ref_base_idx < 4  ){ //TODO Model for bases at which reference is 'N'
//                ModelInput d = {ref_base_idx, bcalls};
            m_all_the_data.push_back(ModelInput{ ref_base_idx, bcalls });
//                std::cout << m_all_the_data.size() << endl;
//                ModelInput d = {ref_base_idx, bcalls};
//                double prob_one = TetMAProbOneMutation(m_params,d);
//                double prob = TetMAProbability(m_params, d);
//                if(prob >= m_prob_cut){
//                     *m_ostream << chr << '\t' 
//                                << pos << '\t' 
//                                << current_base << '\t' 
//                                << prob << '\t' 
//                                << prob_one << '\t' 
//                               << endl;          
//                }
        }
    }
private:
    RefVector m_bam_ref;
    SamHeader m_header;
    Fasta m_idx_ref;
    GenomeData& m_all_the_data;
    SampleMap m_samples;
    BamAlignment& m_ali;
//        ModelParams m_params;
    int m_qual_cut;
    int m_mapping_cut;
    double m_prob_cut;
    char current_base;
    uint64_t chr_index;
};

int main(int argc, char** argv){

//	exit(-1);
    namespace po = boost::program_options;
    string ref_file;
    string config_path;
    po::options_description cmd("Command line options");
    cmd.add_options()
            ("help,h", "Print a help message")
            ("bam,b", po::value<string>()->required(), "Path to BAM file")
            ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
            ("reference,r", po::value<string>(&ref_file)->required(),  "Path to reference genome")
//       ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
//        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags")
            ("qual,q", po::value<int>()->default_value(13),
                    "Base quality cuttoff")

            ("mapping-qual,m", po::value<int>()->default_value(13),
                    "Mapping quality cuttoff")

            ("prob,p", po::value<double>()->default_value(0.1),
                    "Mutaton probability cut-off")
            ("out,o", po::value<string>()->default_value("acuMUlate_result.tsv"),
                    "Out file name")
            ("intervals,i", po::value<string>(), "Path to bed file")
            ("config,c", po::value<string>(), "Path to config file")
            ("theta", po::value<double>()->required(), "theta")
            ("nfreqs", po::value<vector<double> >()->multitoken(), "")
            ("mu", po::value<double>()->required(), "")
            ("seq-error", po::value<double>()->required(), "")
            ("phi-haploid",     po::value<double>()->required(), "")
            ("phi-diploid",     po::value<double>()->required(), "");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cmd), vm);

    if (vm.count("help")){
        cout << cmd << endl;
        return 0;
    }

    if (vm.count("config")){
        ifstream config_stream (vm["config"].as<string>());
        po::store(po::parse_config_file(config_stream, cmd, false), vm);
    }

    vm.notify();
//    ModelParams params = {
//        vm["theta"].as<double>(),
//        vm["nfreqs"].as<vector< double> >(),
//        vm["mu"].as<double>(),
//        vm["seq-error"].as<double>(), 
//        vm["phi-haploid"].as<double>(), 
//        vm["phi-diploid"].as<double>(),
//    };
    string bam_path = vm["bam"].as<string>();
    string index_path = vm["bam-index"].as<string>();
    if(index_path == ""){
        index_path = bam_path + ".bai";
    }

//    ofstream result_stream (vm["out"].as<string>());
    //TODO: check sucsess of all these opens/reads:
    BamReader experiment;
    experiment.Open(bam_path);
    experiment.OpenIndex(index_path);
    RefVector references = experiment.GetReferenceData();
    SamHeader header = experiment.GetHeader();
    Fasta reference_genome; // BamTools::Fasef_file);
    reference_genome.Open(ref_file, ref_file+ ".fai");
//    reference_genome.CreateIndex(ref_file + ".fai");
    PileupEngine pileup;
    BamAlignment ali;

    uint32_t total_len = 0;
    for_each(references.begin(),references.end(),[&](RefData chrom){
        total_len += chrom.RefLength;
    });

    GenomeData base_counts;
    base_counts.reserve(total_len);

    SampleMap samples;
    uint16_t sindex = 0;
    for(auto it = header.ReadGroups.Begin(); it!= header.ReadGroups.End(); it++){
        if(it->HasSample()){
            auto s  = samples.find(it->Sample);
            if( s == samples.end()){ // not in there yet
                samples[it->Sample] = sindex;
                sindex += 1;
            }
        }
    }
    for (auto t : samples) {
        cout << t.first << "\t" << t.second << endl;
    }

    VariantVisitor *v = new VariantVisitor(
            references,
            header,
            reference_genome,
            base_counts,
//            &result_stream,
            samples,
//            params, 
            ali,
            vm["qual"].as<int>(),
            vm["mapping-qual"].as<int>(),
            vm["prob"].as<double>()
    );
    pileup.AddVisitor(v);
//  TODO: Only allocate interval-sized memory vector   
//  if intervals are set
    if (vm.count("intervals")){
        BedFile bed (vm["intervals"].as<string>());
        BedInterval region;
        while(bed.get_interval(region) == 0){
            int ref_id = experiment.GetReferenceID(region.chr);
            experiment.SetRegion(ref_id, region.start, ref_id, region.end);
            while( experiment.GetNextAlignment(ali) ){
                pileup.AddAlignment(ali);
            }
        }
    }
    else{
        cout << "here"<< endl;
        BamAlignment ali;
        while( experiment.GetNextAlignment(ali)){
            pileup.AddAlignment(ali);
        };

        pileup.Flush();
        cout << base_counts.size() << endl;
        cout << "here3"<< endl;

        ModelParams params = {
                vm["theta"].as<double>(),
                vm["nfreqs"].as<vector< double> >(),
                vm["mu"].as<double>(),
                vm["seq-error"].as<double>(),
                vm["phi-haploid"].as<double>(),
                vm["phi-diploid"].as<double>(),
        };
        int count = 0;
        for (auto base : base_counts) {

//		ModelInput d = {ref_base_idx, bcalls};
//		double prob_one = TetMAProbOneMutation(params,d);
            for (auto t : base.all_reads) {
                cout << t.reads[0] <<t.reads[1]<<t.reads[2]<<t.reads[3] << endl;
            }
            cout <<""<< endl;
//    	exit(-1);
            double prob = TetMAProbability(params, base);
            if(prob >= 0.1){
                cout.precision(10);
                cout <<"=================" <<endl;
                cout << prob << "\t" << (1-prob) << endl;
                cout << base.reference  << "\tSIZE:\t" << base.all_reads.size() << endl;
                for (auto t : base.all_reads) {
                    cout << t.reads[0] <<t.reads[1]<<t.reads[2]<<t.reads[3] << endl;
                }
                cout<<count<<endl;
                MutationProb muProb = MutationProb(params);
                SequenceProb sp (params, base, muProb);
                sp.UpdateLikelihood();
                double likelihood = sp.GetLikelihood();
                cout << likelihood << (1-likelihood) << endl;

                int size = 20;
                double muArray[size];
                muArray[0] = 0.1;
                for (int i = 1; i < size; ++i) {
                    muArray[i] = muArray[i-1]*0.1;
                }
                for (int i =0; i < size; ++i) {
                    muProb.UpdateMu(muArray[i]);
                    sp.UpdateMuProb(muProb);
                    likelihood = sp.GetLikelihood();

                    printf("I:%d %e %e %e\n", i, muArray[i], likelihood, (1-likelihood));
                }



//			sp.UpdateMu(1e-7);
//			 likelihood = sp.GetLikelihood();
//			cout << likelihood << (1-likelihood) << endl;
//
//			sp.UpdateMu(1e-8);
//			 likelihood = sp.GetLikelihood();
//			cout << likelihood << (1-likelihood) << endl;
//
//			sp.UpdateMu(1e-60);
//			 likelihood = sp.GetLikelihood();
//			cout << likelihood << (1-likelihood) << endl;


                cout <<"=================" <<endl;
////					std::cout << base.reference << base.all_reads << endl;
////			 *m_ostream << chr << '\t'
////						<< pos << '\t'
////						<< current_base << '\t'
////						<< prob << '\t'
////						<< prob_one << '\t'
////					   << endl;
            }

            count++;
            if(count == 600){ //598 4195 5391 6589
                break;
            }
        }


        cout << count << endl;
        MutationProb muProb = MutationProb(params);
        int xx[3];
        xx[0] = 1;
//		int sp_count = 10;
        SequenceProb sp[5];

        for (int i = 0; i < 5; ++i) {
            sp[i] = SequenceProb(params, base_counts[500+i], muProb);
        }
        sp[0] = SequenceProb(params, base_counts[598], muProb);
//		sp[1] = SequenceProb(params, base_counts[4195], muProb);
        sp[1] = SequenceProb(params, base_counts[5391], muProb);
//		sp[2] = SequenceProb(params, base_counts[5391], muProb);
//		sp[3] = SequenceProb(params, base_counts[6589], muProb);

        ModelInput base_custom =  base_counts[598];
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 4; ++j) {
                base_custom.all_reads[i].reads[j] = (uint16_t) i*2+j*3;
            }
        }
//		base_custom.all_reads[0].reads[0] = 12;
//		base_custom.all_reads[1].reads = {42,65,12,35};
//		base_custom.all_reads[2].reads = {2,6,84,5};
//		base_custom.all_reads[3].reads = {7,48,51,21};
//		base_custom.all_reads[4].reads = {24,24,23,29};

        sp[2] = SequenceProb(params, base_custom, muProb);
//		double d = sp[0].GetLikelihood();
//		cout << "D:\t" << d << endl;
//		for (auto t : sp) {
//			t.UpdateLikelihood();
//			double likelihood = t.GetLikelihood();
//			cout << likelihood << (1 - likelihood) << endl;
//		}

        int size = 20;
        double muArray[size];
        double likelihood[size];
        muArray[0] = 1;
        for (int i = 1; i < size; ++i) {
            muArray[i] = muArray[i - 1] * 0.2;
        }
        for (int i = 0; i < size; ++i) {
            muProb.UpdateMu(muArray[i]);
            for (auto t : sp) {
                t.UpdateMuProb(muProb);
                cout << log( t.GetLikelihood() )<< " ";

                likelihood[i] += log( t.GetLikelihood() );
            }
            printf("\tI:%d\t%.1e %e %f \n", i, muArray[i], likelihood[i], likelihood[i]
            );
        }
        double mma = *std::max_element( likelihood, likelihood+size);
        double mmi = *std::min_element( likelihood, likelihood+size);
        cout << mma <<"\t"<< mmi << endl;


        const int cat = 2;
        muArray[cat];
        muArray[0] = 1e-2;
        muArray[0] = 1e-5;

        likelihood[cat];
        double proportion[cat];
        double weight[2];
        for (auto t : sp) {
            for (int i = 0; i < cat; ++i) {
                muProb.UpdateMu(muArray[i]);

                t.UpdateMuProb(muProb);
//				cout << log( t.GetLikelihood() )<< " ";
                likelihood[i] = t.GetLikelihood();
//				likelihood[i] += log( t.GetLikelihood() );
            }
            double sum = likelihood[0] + likelihood[1];

            proportion[0] = likelihood[0]/sum;
            proportion[1] = likelihood[1]/sum;



            cout << proportion[0] << "\t" << proportion[1] << endl;

            ModelInput m = t.GetData();
//			m.all_reads

            ModelInput mA, mB;
            ReadData vA, vB;
            //Might have to redo this part, with double[4] TODO:
            mA.all_reads = vector<ReadData>(7);
            for (int b = 0; b < m.all_reads.size(); ++b) {
                ReadData v = m.all_reads[b];

                cout << v.reads[0] << "\t" << v.reads[1] << "\t"<< v.reads[2] << "\t"<< v.reads[3] << "\t" <<endl;
                for (int bb = 0; bb < 3; ++bb) {
                    vA.reads[bb] = v.reads[bb] * proportion[0];
                    vB.reads[bb] = v.reads[bb] * proportion[1];
                }
                cout << vA.reads[0] << "\t" << vA.reads[1] << "\t"<< vA.reads[2] << "\t"<< vA.reads[3] << "\t" <<endl;
                cout << "AA" << endl;
                mA.all_reads[b] = vA;
                cout << "BB" << endl;
                mB.all_reads.push_back( vB );
                cout << "??" << endl;
            }
        }
        cout << muArray[0] << "\t" << muArray[2] << endl;
        return 0;


    }
}

