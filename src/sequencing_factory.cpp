//
// Created by steven on 3/31/15.
//

#include "sequencing_factory.h"


SequencingFactory::SequencingFactory(ModelParams const &model_params2) :model_params(model_params2) {

    phi_haploid = model_params.phi_haploid;
    phi_diploid = model_params.phi_diploid;
    error_prob = model_params.error_prob;
    theta = model_params.theta;
	frequency_prior = model_params.nuc_freq;

//    size_t i=0;
//    DiploidProbs pop_genotypes = DiploidPopulation(site_data.reference);
////    ancestor_data = site_data.all_reads[i];
////	ancestor_genotypes = DiploidSequencing(ancestor_data);
////	ancestor_genotypes *= pop_genotypes;
//
//	descendant_count = site_data.all_reads.size() - 1;
//	descendant_index.assign(descendant_count, -1);
//
//	all_descendant_data.reserve(descendant_count);
//	all_descendant_genotypes.reserve(descendant_count);
//
//    for (i = 1; i < (descendant_count+1); ++i) {
//        ReadData data = site_data.all_reads[i];
//        all_descendant_data.push_back(data);
////		all_descendant_data[i-1] =data ;
//
////        HaploidProbs p = HaploidSequencing(data);
//		all_descendant_genotypes.emplace_back(HaploidSequencing(data));
////        all_descendant_genotypes.push_back(p);
////		all_descendant_genotypes[i-1]=p;
//
//    }


    double alphas_total = (1.0 - phi_haploid) / phi_haploid;
	for (int i : { 0, 1, 2, 3 }) {
        std::fill(haploid_alphas[i], haploid_alphas[i] + 4, error_prob / 3.0 * alphas_total);
        for (int k : {0, 1, 2, 3}) {
            if (k == i) {
                haploid_alphas[i][k] = (1.0 - error_prob) * alphas_total;
                break;
            }
        }

    }


}



SequenceProb & SequencingFactory::InitSequenceProb(SequenceProb &seq_prob, ModelInput &site_data) {

    seq_prob.AddModel(model_params);
    seq_prob.SetupDiploid(site_data);

//    CalculateHaploidProb(seq_prob, site_data);

    return seq_prob;
}

void SequencingFactory::CalculateHaploidProb(SequenceProb &seq_prob, ModelInput &site_data) {

    std::vector<HaploidProbs> all_descendant_genotypes;
    all_descendant_genotypes.reserve(seq_prob.GetDescendantCount());

//    ReadDataVector GetDescendantReadData();
//    ReadDataVector const & GetDescendantReadData2();
//    ReadDataVector const * GetDescendantReadData3();



    for (size_t i = 0; i < seq_prob.GetDescendantCount(); ++i) {
//        ReadData data = seq_prob.GetDescendantReadData(i);
//        HaploidProbs p = HaploidSequencing(data);

        int index = seq_prob.GetDescendantIndex(i);
//        convert_index_key_to_haploid.push_back(HaploidSequencing(rd));

        all_descendant_genotypes.emplace_back(convert_index_key_to_haploid[index]);
//        all_descendant_genotypes.push_back(p);
//		all_descendant_genotypes[i-1]=p;

    }

    seq_prob.SetupHaploid(all_descendant_genotypes);
}



HaploidProbs SequencingFactory::HaploidSequencing(ReadData const &data) {
    HaploidProbs result;

//    double alphas_total = (1.0 - phi_haploid) / phi_haploid;
//    map_readdata_haploidprobs.
    for (int i : { 0, 1, 2, 3 }) {
//        double alphas[4];
//        std::fill(alphas, alphas+4, error_prob / 3.0 * alphas_total);
//        for (int k : { 0, 1, 2, 3 }) {
//            if (k == i) {
//                alphas[k] = (1.0 - error_prob) * alphas_total;
//                break;
//            }
////			else[k] = error_prob / 3.0 * alphas_total;
//        }
        result[i] = DirichletMultinomialLogProbability(haploid_alphas[i], data);
    }

    double scale = result.maxCoeff();
    return (result - scale).exp();
//    result -= scale;
//    result = ScaleLogArray(result);
//    return result.exp();
}

//void temp() {
//    int index = 0;
//    for (size_t i = 0; i < all_sequence_prob.size(); ++i) {
//        SequenceProb &item = all_sequence_prob[i];
//
//        for (int j = 0; j < item.GetDescendantCount(); ++j) {
//            ReadData rd = item.GetDescendantReadData(j);
//            auto rd_key = rd.key;
//            count[rd_key]++;
//            auto find_key = map_rd_to_index.find(rd_key);
//
//            if (find_key == map_rd_to_index.end()) {
//
//                map_rd_key_to_haploid[rd_key] = item.GetDescendantGenotypes(j);
//
////                std::array<std::array<double, 2>, 10> temp;
//                std::array<std::pair<double, double>, 10> temp;
//
//
////                cache_read_data_to_all_index.push_back( temp);
//                cache_read_data_to_all_index.push_back(temp);
//                map_index_key_to_haploid.push_back(item.GetDescendantGenotypes(j));
//                map_rd_to_index[rd_key] = index;
//                index++;
////                std::cout << map_rd_key_to_haploid.size() << "\t" << cache_read_data_to_all.size() << std::endl;
////                std::cout << "\t" << rd.reads[0] << "\t" << rd.reads[1] << "\t" <<rd.reads[2] << "\t" <<rd.reads[3] << "\t" << rd_key << std::endl;
//            }
//            item.SetDescendantIndex(j, map_rd_to_index[rd_key]);
//        }
//
//    }
//    std::cout << map_rd_key_to_haploid.size() << "\t" << cache_read_data_to_all.size() << "\t" << index << std::endl;
//}


void SequencingFactory::CreateSequenceProbsVector(std::vector<SequenceProb> &sp, GenomeData &genome_data) {


    sp.reserve(genome_data.size());
    for (size_t i = 0; i < genome_data.size(); ++i){

//        std::cout << i << "\t" << std::endl;
        SequenceProb ss;
        InitSequenceProb(ss, genome_data[i]);
        sp.emplace_back(ss);//copy, how to move?? std::move?
    }

    std::unordered_map<uint64_t, int> map_rd_to_index;

    int index = 0;
    for (size_t i = 0; i < sp.size(); ++i) {
        SequenceProb &item = sp[i];

        for (int j = 0; j < item.GetDescendantCount(); ++j) {
            ReadData rd = item.GetDescendantReadData(j);
            auto rd_key = rd.key;
            auto find_key = map_rd_to_index.find(rd_key);

            if (find_key == map_rd_to_index.end()) {

//                map_rd_key_to_haploid[rd_key] = item.GetDescendantGenotypes(j);
//                std::cout << map_rd_key_to_haploid.size() << "\t" << index << std::endl;
                convert_index_key_to_haploid.push_back(HaploidSequencing(rd));

                map_rd_to_index[rd_key] = index;
                index++;
//                std::cout << map_rd_key_to_haploid.size() << "\t" << index << std::endl;
//                std::cout << "\t" << rd.reads[0] << "\t" << rd.reads[1] << "\t" <<rd.reads[2] << "\t" <<rd.reads[3] << "\t" << rd_key << std::endl;
            }
            item.SetDescendantIndex(j, map_rd_to_index[rd_key]);
        }
        CalculateHaploidProb(item, genome_data[i]);

    }
    std::cout << map_rd_key_to_haploid.size() << "\t" << index << "\t" << convert_index_key_to_haploid.size() << "\t" << map_rd_to_index.size() <<std::endl;


//    return SequenceProb();
}

void SequencingFactory::CreateSequenceProbV1(std::vector<SequenceProbV1> &sp, GenomeData &genome_data) {
    sp.reserve(genome_data.size());
    for (size_t i = 0; i < genome_data.size(); ++i) {
//        std::cout << i << "\t" << std::endl;
//        SequenceProbV1 ss = SequenceProbV1(genome_data[i], model_params);
//        sp.emplace_back(ss); //copy
//        sp.push_back(ss); //copy
//        sp.emplace_back(SequenceProbV1(genome_data[i], model_params)); //move
        sp.emplace_back(genome_data[i], model_params);

    }
}