/*
 * genome_data_stream.h
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */


#pragma once
#ifndef GENOME_DATA_STREAM_H_
#define GENOME_DATA_STREAM_H_

#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include <cstdio>

#include <cstring>
#include <cmath> // for log2f

#include <mutations/model.h>
#include <errno.h>
#include <stddef.h>
#include <stdint.h>


#define BUFFER_SIZE 16384 // same as kseq.h

#define WRITE_BUFFER  32768// 16384  //800000
#define TAILLE_NOM 1024

off_t fsize(const char *filename);

class GenomeDataStream {

protected:
    char filename[TAILLE_NOM];
    FILE *binary_read_file;
    int sizeElement;
    void *buffer;
    int cpt_buffer;
    int buffer_size_nelem;


protected:
    void ReadModelInput(ModelInput &model_input);

    void ReadReference(uint16_t &ref);

    void ReadReadDataKey(uint64_t &key);

    void write(void const *element, int size);

    size_t read(void *element, int size);

    void write(void *element, int size);

public:
    GenomeDataStream() {
    }

    GenomeDataStream(std::string filename, bool write);
//    GenomeDataStream(char *filename, int sizeElement, bool write);
//    void write_element(void *element);
//    void write_element(void *element, int size);
//    size_t read_element(void *element);
//    size_t read_element_buffered(void *element);

//    void write_element_buffered(void *element);
//    off_t nb_elements();
    void rewind_all();
//int get_sizeElement(){ return sizeElement; };

    void close();

    void open(bool write);

    off_t file_size();

    size_t read_kmer(void *element);

    size_t read_colour(void *element);

    size_t read_kmer_colour(void *element, void *element2);

    size_t read_kmer_skip_colour(void *element);


    ~GenomeDataStream();

    void WriteHeader(uint64_t const &sequence_count);

    void WriteHeader(uint64_t const &base_count, uint64_t const &sequence_count);

    void WriteReference(const uint16_t &ref);

    void WriteReadDataKey(uint64_t const &key);

    void WriteModelInput(ModelInput &input);

    void ReadHeader(uint64_t &site_count0, uint64_t &read_data_count0);

    void ReadGenomeData(GenomeData &data);

    uint64_t read_data_count;
    uint64_t site_count;

    void End();


    ModelInput ReadOne();
};


#endif //GENOME_DATA_STREAM_H_
