/*
 * boost_input_utils.cc
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */



#include "boost_input_utils.h"

namespace BoostUtils {
    using namespace std;


    static bool file_exists_test3(const std::string &name) {
        struct stat buffer;
        return (stat(name.c_str(), &buffer) == 0);
    }

    static bool file_exists_test1(const std::string &name) {
        ifstream ifile(name.c_str());
        return ifile.good();

//    ifstream f(name.c_str());
//    if (f.good()) {
//        f.close();
//        return true;
//    } else {
//        f.close();
//        return false;
//    }
    }


    void ParseCommandLinkeInput(int argc, char **argv, boost::program_options::variables_map &vm) {
//    boost::program_options::variables_map vm;
        namespace po = boost::program_options;
        std::string ref_file;
        std::string config_path;
        po::options_description cmd("Command line options");
        cmd.add_options()
                ("help,h", "Print a help message")
                ("bam,b", po::value<string>()->required(), "Path to BAM file")
                ("bam-index,x", po::value<string>()->default_value(""), "Path to BAM index, (defalult is <bam_path>.bai")
                ("reference,r", po::value<string>(&ref_file)->required(), "Path to reference genome")
//       ("ancestor,a", po::value<string>(&anc_tag), "Ancestor RG sample ID")
//        ("sample-name,s", po::value<vector <string> >()->required(), "Sample tags")
                ("qual,q", po::value<int>()->default_value(13), "Base quality cuttoff")

                ("mapping-qual,m", po::value<int>()->default_value(13), "Mapping quality cuttoff")

                ("prob,p", po::value<double>()->default_value(0.1), "Mutaton probability cut-off")
                ("out,o", po::value<string>()->default_value("acuMUlate_result.tsv"),
                        "Out file name")
                ("intervals,i", po::value<string>(), "Path to bed file")
                ("config,c", po::value<string>(), "Path to config file")
                ("theta", po::value<double>()->required(), "theta")
                ("nfreqs", po::value<vector<double> >()->multitoken(), "")
                ("mu", po::value<double>()->required(), "")
                ("seq-error", po::value<double>()->required(), "")
                ("phi-haploid", po::value<double>()->required(), "")
                ("phi-diploid", po::value<double>()->required(), "")
                ("output_binary_file,", po::value<string>()->default_value("default_binary_output.bin"),
                        "Output the parsed BAM file to binary file")
                ("outfile", po::value<string>()->required(), "outfile prefix")
                ("thread,t", po::value<int>()->default_value(1), "thread count")
    ;

        po::store(po::parse_command_line(argc, argv, cmd), vm);

        if (vm.count("help")) {
            cout << cmd << endl;
            exit(-1);
        }

        if (vm.count("config")) {
            ifstream config_stream(vm["config"].as<string>());
            po::store(po::parse_config_file(config_stream, cmd, false), vm);
        }

        vm.notify();
//    string bam_path = vm["bam"].as<string>();
//    string index_path = vm["bam-index"].as<string>();
//    string ref_file2 = vm["reference"].as<string>();
//    cout << bam_path << "\t" << index_path << ref_file << "\t" << ref_file2 << endl;

    }


    void ExtractInputVariables(boost::program_options::variables_map &vm, GenomeData &genome_data,
            BamTools::BamReader &experiment, BamTools::RefVector &references,
            BamTools::SamHeader &header, LocalBamToolsUtils::Fasta &reference_genome) {


        std::string ref_file = vm["reference"].as<std::string>();
        std::string bam_path = vm["bam"].as<std::string>();
        std::string index_path = vm["bam-index"].as<std::string>();
        if (index_path == "") {
            index_path = bam_path + ".bai";
        }
        experiment.Open(bam_path);
        experiment.OpenIndex(index_path);
        references = experiment.GetReferenceData();
        header = experiment.GetHeader();

//    ofstream result_stream (vm["out"].as<std::string>());
        //TODO: check sucsess of all these opens/reads:experiment.Open(bam_path);
//        experiment.OpenIndex(index_path);// BamTools::Fasef_file);


        if (!file_exists_test1(ref_file + ".fai")) {
            reference_genome.Open(ref_file);
            reference_genome.CreateIndex(ref_file + ".fai");
        }
        else {
            reference_genome.Open(ref_file, ref_file + ".fai");
        }

//  Assign some memory for the big list
        uint32_t total_len = 0;
        if (vm.count("intervals")) {
            BedFile bed(vm["intervals"].as<std::string>());
            BedInterval region;
            while (bed.get_interval(region) == 0) {
                total_len += (region.end - region.start);
            }
        }
        else {

            for_each(references.begin(), references.end(), [&](BamTools::RefData chrom) {
                total_len += chrom.RefLength;
            });
        }
        genome_data.clear();
        genome_data.reserve(total_len);


    }


    ModelParams CreateModelParams(boost::program_options::variables_map variables_map) {

        ModelParams params = {
                variables_map["theta"].as<double>(),
                variables_map["nfreqs"].as<vector<double> >(),
                variables_map["mu"].as<double>(),
                variables_map["seq-error"].as<double>(),
                variables_map["phi-haploid"].as<double>(),
                variables_map["phi-diploid"].as<double>(),
        };
        return params;

    }
}
