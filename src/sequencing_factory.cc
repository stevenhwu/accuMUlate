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

    for (int i : { 0, 1, 2, 3 }) {
        ref_diploid_probs[i] = CreateRefDiploidProbs(i);
    }
    CalculateAncestorPrior();

    index = 0;
    index_ancestor = 0;

}


DiploidProbs SequencingFactory::CreateRefDiploidProbs(int ref_allele) {
    ReadData d;
    DiploidProbs result;

    double alphas[4];
    for(int i : {0,1,2,3}){
        alphas[i] = theta * frequency_prior[i];
    }

    for(int i : {0,1,2,3}) {
        for(int j=0;j<i;++j) {
            d.key = 0;
            d.reads[ref_allele] = 1;
            d.reads[i] += 1;
            d.reads[j] += 1;
//			printf("  %d %d: %u %u %u %u\n", i, j, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
            result[i+j*4] = DirichletMultinomialLogProbability(alphas, d);
            result[j+i*4] = result[i+j*4];
        }
        d.key = 0;
        d.reads[ref_allele] = 1;
        d.reads[i] += 2;
//		printf("D:%d %d: %u %u %u %u\n", i, i, d.reads[0], d.reads[1], d.reads[2], d.reads[3]);
        result[i+i*4] = DirichletMultinomialLogProbability(alphas, d);
    }
//	std::cout << result << std::endl;
    return result.exp();
}


DiploidProbs SequencingFactory::DiploidSequencing(ReadData const &data) {
    DiploidProbs result;
    double alphas_total = (1.0 - phi_diploid) / phi_diploid;
    for (int i : { 0, 1, 2, 3 }) {
        for (int j = 0; j < i; ++j) {
            double alphas[4];
            for (int k : { 0, 1, 2, 3 }) {
                if (k == i || k == j)
                    alphas[k] = (0.5 - error_prob / 3.0) * alphas_total;
                else
                    alphas[k] = (error_prob / 3.0) * alphas_total;
            }
            result[i * 4 + j] = DirichletMultinomialLogProbability(alphas, data);
            result[j * 4 + i] = result[i * 4 + j];
        }
        double alphas[4];
        for (int k : { 0, 1, 2, 3 }) {
            if (k == i)
                alphas[k] = (1.0 - error_prob) * alphas_total;
            else
                alphas[k] = error_prob / 3.0 * alphas_total;
        }
        result[i * 4 + i] = DirichletMultinomialLogProbability(alphas, data);
    }


    double scale = result.maxCoeff();
    return result;
//    return (result - scale).exp();

}

HaploidProbs SequencingFactory::HaploidSequencing(ReadData const &data) {
    HaploidProbs result;
    for (int i : { 0, 1, 2, 3 }) {
        result[i] = DirichletMultinomialLogProbability(haploid_alphas[i], data);
    }
    double scale = result.maxCoeff();
    return result;
//    return (result - scale).exp();
}


//
//void SequencingFactory::CreateSequenceProbsVector(std::vector<SiteGenotypes> &sp, GenomeData &genome_data) {
//
//    //    ReadDataVector GetDescendantReadData();
////    ReadDataVector const & GetDescendantReadData2();
////    ReadDataVector const * GetDescendantReadData3();
//
//    sp.reserve(genome_data.size());
//    for (size_t i = 0; i < genome_data.size(); ++i){
//        sp.emplace_back(genome_data[i]);
//    }
//
//    int index = 0;
//    int index_ancestor = 0;
//
//
//    sp.reserve(genome_data.size());
//
//    for (size_t i = 0; i < genome_data.size(); ++i) {
//
////        uint16_t ref = genome_data[i].reference;
////        int descendant_conut = genome_data[i].all_reads.size()-1;
//
//        ModelInput &data = genome_data[i];
//        int descendant_conut = data.all_reads.size() - 1;
//
////        sp.emplace_back(data.reference, descendant_conut);
//        SiteGenotypes &item = sp[i];
//
//        for (int j = 0; j < item.GetDescendantCount(); ++j) {
////            std::cout << i <<  "\t" << j << std::endl;
////            ReadData rd = item.GetDescendantReadData(j);
////            ReadData &rd = genome_data[i].all_reads[j+1];
//            ReadData &rd = data.all_reads[j+1];
//            auto rd_key = rd.key;
//
//            auto find_key = map_rd_to_index.find(rd_key);
//
//            if (find_key == map_rd_to_index.end()) {
//                convert_index_key_to_haploid.push_back(HaploidSequencing(rd));
//                map_rd_to_index[rd_key] = index;
//                index++;
//            }
//            item.SetDescendantIndex(j, map_rd_to_index[rd_key]);
//        }
//
//
////        ReadData rd = item.GetAncestorReadData();
////        ReadData &rd = genome_data[i].all_reads[0];
//        ReadData &rd = data.all_reads[0];
//        auto rd_key = rd.key;
//        auto find_key = map_ancestor_to_index[data.reference].find(rd_key);
//        if (find_key == map_ancestor_to_index[data.reference].end()) {
//            convert_index_key_to_diploid.push_back(DiploidSequencing(rd));
//            map_ancestor_to_index[data.reference][rd_key] = index_ancestor;
//            index_ancestor++;
//        }
//        item.SetAncestorIndex(map_ancestor_to_index[data.reference][rd_key]);
//
//
//        CalculateDescendantGenotypes(item);
//        CalculateAncestorGenotype(item);
//
//    }
//    std::cout << index << "\t" << convert_index_key_to_haploid.size() << "\t" << map_rd_to_index.size() <<std::endl;
//    std::cout << index_ancestor << "\t" << convert_index_key_to_diploid.size() << "\t" << map_ancestor_to_index.size() <<std::endl;
//
//
////    return SiteGenotypes();
//}



//void SequencingFactory::CreateSequenceProb(SequenceProb &sp, ModelInput const &data, ModelParams const params){
//
//}std::vector<HaploidProbs> & SequencingFactory::RemoveConvertIndexKeyToHaploid(){
//return <#initializer#>;
//}

const std::vector<HaploidProbs> SequencingFactory::RemoveConvertIndexKeyToHaploid() {
    return std::move(convert_index_key_to_haploid);
}

const std::vector<DiploidProbsIndex10> SequencingFactory::RemoveConvertIndexKeyToDiploidIndex10() {
    return std::move(convert_index_key_to_diploid_10);
}

std::vector<double> && SequencingFactory::RemoveConvertIndexKeyToHaploidScaler() {
    return std::move(convert_index_key_to_haploid_scaler);
}

std::vector<double> && SequencingFactory::RemoveConvertIndexKeyToDiploidIndex10Scaler() {
    return std::move(convert_index_key_to_diploid_10_scaler);
}

const std::vector<HaploidProbs> SequencingFactory::RemoveConvertIndexKeyToHaploidUnnormalised() {

    return std::move(convert_index_key_to_haploid_unnormalised);
}

const std::vector<DiploidProbsIndex10> SequencingFactory::RemoveConvertIndexKeyToDiploidIndex10Unnormalised() {

    return std::move(convert_index_key_to_diploid_10_unnormalised);
}


void SequencingFactory::CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sgi, GenomeData &genome_data) {

//
//    std::unordered_map<uint64_t, uint> map_rd_to_index;
//    std::array<std::unordered_map<uint64_t, uint>, 4> map_ancestor_to_index;


    sgi.reserve(genome_data.size());
//    for (size_t i = 0; i < genome_data.size(); ++i){
//        sgi.emplace_back(genome_data[i]);
//    }

    uint index = 0;
    uint index_ancestor = 0;


//    sgi.reserve(genome_data.size());

    for (size_t i = 0; i < genome_data.size(); ++i) {

//        uint16_t ref = genome_data[i].reference;
//        int descendant_conut = genome_data[i].all_reads.size()-1;

        ModelInput data = std::move(genome_data[i]);//FIXME:: This can reduce the memeory even further
//        ModelInput data = genome_data[i];
        int descendant_conut = data.all_reads.size() - 1;

//        sgi.emplace_back(descendant_conut);
//        SiteGenotypesIndex &item = sgi[i];

        SiteGenotypesIndex item (descendant_conut);

        for (int j = 0; j < item.GetDescendantCount(); ++j) {
//            std::cout << i <<  "\t" << j << std::endl;
//            ReadData rd = item.GetDescendantReadData(j);
//            ReadData &rd = genome_data[i].all_reads[j+1];
            ReadData &rd = data.all_reads[j+1];
            auto rd_key = rd.key;

            auto find_key = map_rd_to_index.find(rd_key);

            if (find_key == map_rd_to_index.end()) {

                HaploidProbs prob = HaploidSequencing(rd);
                {
                    HaploidProbs prob_exp = prob.exp();
                    convert_index_key_to_haploid_unnormalised.push_back(prob_exp);
                }

                double scale = prob.maxCoeff();
                convert_index_key_to_haploid_scaler.push_back(scale);
                prob = (prob - scale).exp();
                convert_index_key_to_haploid.push_back(prob);

                map_rd_to_index[rd_key] = index;
                index++;

            }
            item.SetDescendantIndex(j, map_rd_to_index[rd_key]);
        }


//        ReadData rd = item.GetAncestorReadData();
//        ReadData &rd = genome_data[i].all_reads[0];
        ReadData &rd = data.all_reads[0];

        auto rd_key = rd.key;
        auto find_key = map_ancestor_to_index[data.reference].find(rd_key);
        if (find_key == map_ancestor_to_index[data.reference].end()) {
            DiploidProbs temp_prob = DiploidSequencing(rd);
            {
                DiploidProbsIndex10 temp_dip_10 = ConvertDiploid16ToDiploid10(temp_prob.exp(), data.reference);
                convert_index_key_to_diploid_10_unnormalised.push_back(temp_dip_10);
            }

            double scale = temp_prob.maxCoeff();
            convert_index_key_to_diploid_10_scaler.push_back(scale);
            temp_prob = (temp_prob - scale).exp();
            DiploidProbsIndex10 temp_dip_10 = ConvertDiploid16ToDiploid10(temp_prob, data.reference);
            convert_index_key_to_diploid_10.push_back(temp_dip_10);



            map_ancestor_to_index[data.reference][rd_key] = index_ancestor;
            index_ancestor++;
        }
        item.SetAncestorIndex(map_ancestor_to_index[data.reference][rd_key]);

        sgi.push_back(std::move(item));

//        if(i==5){
//            break;
//        }
    }
    std::cout << "Cache descendant size: " << index << "\t" << convert_index_key_to_haploid.size() << "\t" << map_rd_to_index.size() <<std::endl;
    std::cout << "Cache ancestor size: " << index_ancestor << "\t" << convert_index_key_to_diploid_10.size() << "\t" << map_ancestor_to_index.size() <<std::endl;

//    return SiteGenotypes();
}

void SequencingFactory::CalculateAncestorPrior() {
    for (int i = 0; i < 4; ++i) {
        for (int j = i; j < 4; ++j) {
            int index10 = LookupTable::index_converter_4_4_to_10[i][j];
            ancestor_prior[index10] = frequency_prior[i] * frequency_prior[j];
            if(i != j){
                ancestor_prior[index10] *= 2; //Count both AC and CA
            }
        }
    }
}


//
//void SequencingFactory::CalculateAncestorGenotypeIndex(SiteGenotypesIndex &seq_prob) {
//
//    auto ancestor_genotypes = convert_index_key_to_diploid[seq_prob.GetAncestorIndex()];
//    ancestor_genotypes *= ref_diploid_probs[seq_prob.GetReference()];
////    seq_prob.SetAncestorGenotypes(ancestor_genotypes);
//
////    seq_prob.AddModel(model_params);
////    seq_prob.SetupDiploid(site_data);
////        DiploidProbs pop_genotypes = CreateRefDiploidProbs(site_data.reference);
////    auto ancestor_genotypes = convert_index_key_to_diploid[seq_prob.GetAncestorIndex()];
////    ancestor_genotypes *= ref_diploid_probs[site_data.reference];
////    seq_prob.SetAncestorGenotypes(ancestor_genotypes);
////    seq_prob.SetAncestorGenotypes(pop_genotypes);
////    CalculateDescendantGenotypes(seq_prob, site_data);
//
//
//}
//
//void SequencingFactory::CalculateDescendantGenotypesIndex(SiteGenotypesIndex &seq_prob) {
//
//    std::vector<HaploidProbs> all_descendant_genotypes;
//    all_descendant_genotypes.reserve(seq_prob.GetDescendantCount());
//
//    for (size_t i = 0; i < seq_prob.GetDescendantCount(); ++i) {
//        int index = seq_prob.GetDescendantIndex(i);
//        all_descendant_genotypes.emplace_back(convert_index_key_to_haploid[index]);
//    }
//    seq_prob.SetDescendantGenotypes(all_descendant_genotypes);
//}


DiploidProbsIndex10 SequencingFactory::ConvertDiploid16ToDiploid10(DiploidProbs diploid_16, int reference) {
    diploid_16 *= ref_diploid_probs[reference];
    DiploidProbsIndex10 temp_diploid_10;
    for (int index10 = 0; index10 < ANCESTOR_COUNT; ++index10) {
        int index16 = LookupTable::index_converter_10_to_16[index10];
        temp_diploid_10[index10] = diploid_16[index16] * ancestor_prior[index10];
    }

    return temp_diploid_10;
}

std::vector<SiteGenotypesIndex> SequencingFactory::CreateSequenceProbsVector(GenomeData &data) {
    std::vector<SiteGenotypesIndex> sgi;
    CreateSequenceProbsVector(sgi, data);
    return sgi;
}

void SequencingFactory::CreateSequenceProbsVector(std::vector<SiteGenotypesIndex> &sgi, ModelInput &data) {



    int descendant_conut = data.all_reads.size() - 1;
    SiteGenotypesIndex item (descendant_conut);

    for (int j = 0; j < item.GetDescendantCount(); ++j) {
        ReadData &rd = data.all_reads[j+1];
        auto rd_key = rd.key;

        auto find_key = map_rd_to_index.find(rd_key);

        if (find_key == map_rd_to_index.end()) {

            HaploidProbs prob = HaploidSequencing(rd);
            {
                HaploidProbs prob_exp = prob.exp();
                convert_index_key_to_haploid_unnormalised.push_back(prob_exp);
            }

            double scale = prob.maxCoeff();
            prob = (prob - scale).exp();
            convert_index_key_to_haploid.push_back(prob);

            map_rd_to_index[rd_key] = index;
            index++;
//                std::cout << map_rd_to_index[rd_key] << std::endl;
        }
        item.SetDescendantIndex(j, map_rd_to_index[rd_key]);
    }


//        ReadData rd = item.GetAncestorReadData();
//        ReadData &rd = genome_data[i].all_reads[0];
    ReadData &rd = data.all_reads[0];

    auto rd_key = rd.key;
    auto find_key = map_ancestor_to_index[data.reference].find(rd_key);
    if (find_key == map_ancestor_to_index[data.reference].end()) {
        DiploidProbs temp_prob = DiploidSequencing(rd);
        {
            DiploidProbsIndex10 temp_dip_10 = ConvertDiploid16ToDiploid10(temp_prob.exp(), data.reference);
            convert_index_key_to_diploid_10_unnormalised.push_back(temp_dip_10);
        }

        double scale = temp_prob.maxCoeff();
        temp_prob = (temp_prob - scale).exp();
        DiploidProbsIndex10 temp_dip_10 = ConvertDiploid16ToDiploid10(temp_prob, data.reference);
        convert_index_key_to_diploid_10.push_back(temp_dip_10);



        map_ancestor_to_index[data.reference][rd_key] = index_ancestor;
        index_ancestor++;
    }
    item.SetAncestorIndex(map_ancestor_to_index[data.reference][rd_key]);

//        std::cout << "TEST: " << i << std::endl;
//        CalculateDescendantGenotypesIndex(item);
//        CalculateAncestorGenotypeIndex(item);
//        SiteGenotypesIndex tt = std::move(item);
    sgi.push_back(std::move(item));
//        sgi.push_back(item);
//        sgi.emplace_back(item);
//        std::cout << sgi[i].GetDescendantIndex(0) << "\t" << item.GetDescendantIndex(0) << std::endl;
//        std::cout << sgi[i].GetAncestorIndex() << "\t" << sgi[i].GetDescendantIndex().size() << "\t" << std::endl;
//        std::cout << item.GetAncestorIndex() << "\t" << item.GetDescendantIndex().size() << "\t" << std::endl;

//        if(i==5){
//            break;
//        }
}
