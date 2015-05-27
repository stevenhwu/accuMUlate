/*
 * genome_data_stream.cc
 *
 *  Created on: 12/12/14
 *      Author: Steven Wu
 */


#include <stddef.h>
#include <stdint.h>


#include "genome_data_stream.h"


off_t fsize(const char *filename) {
    struct stat st;

    if (stat(filename, &st) == 0)
        return st.st_size;

    return -1;
}
GenomeDataStream::GenomeDataStream(std::string given_filename, bool write) : sizeElement(0)
{
    strcpy(filename,given_filename.c_str());

    open(write);
//    buffer_size_nelem= (WRITE_BUFFER/given_sizeElement);
//    buffer = (void *) malloc(given_sizeElement * buffer_size_nelem);
//    buffer_size_nelem= (WRITE_BUFFER/given_sizeElement);
    buffer = (void *) malloc(WRITE_BUFFER);
    cpt_buffer=0;
}
//GenomeDataStream::GenomeDataStream(char *given_filename, int given_sizeElement, bool write) : sizeElement(given_sizeElement)
//{
//    strcpy(filename,given_filename);
//    open(write);
//    buffer_size_nelem= (WRITE_BUFFER/given_sizeElement);
//    buffer = (void *) malloc(given_sizeElement * buffer_size_nelem);
//    cpt_buffer=0;
//}
//
//
//void GenomeDataStream::write_element( void *element)
//{
//    //  flockfile(binary_read_file);ca
//    // fprintf(stderr,"write elem %lli \n",*(int64_t *)element);
//    if (!fwrite(element, sizeElement, 1, binary_read_file))
//    {
//        // funlockfile(binary_read_file);
//        printf("error: can't fwrite (disk full?)\n");
//        exit(1);
//    }
//    //  funlockfile(binary_read_file);
//}
//
//void GenomeDataStream::write_element( void *element, int size)
//{
//
//    if (!fwrite(element, size, 1, binary_read_file))
//    {
//        printf("error: can't fwrite (disk full?)\n");
//        exit(1);
//    }
//
//}
//
//void GenomeDataStream::write_element_buffered( void *element)
//{
//
//    if(cpt_buffer==buffer_size_nelem)
//    {
//        if (!fwrite(buffer, sizeElement, buffer_size_nelem, binary_read_file))
//        {
//            printf("error: can't fwrite (disk full?)\n");
//            exit(1);
//        }
//        cpt_buffer=0;
//    }
//
//    //((kmer_type *)buffer)[cpt_buffer]= *((kmer_type *)element);
//    memcpy((unsigned char *)buffer + (cpt_buffer * sizeElement), element, sizeElement);
//    cpt_buffer++;
//
//
//
//}
//
//
//
//size_t GenomeDataStream::read_element( void *element)
//{
//    return fread(element, sizeElement,1, binary_read_file);
//}
//
//size_t GenomeDataStream::read_element_buffered( void *element)
//{
//    if(cpt_buffer==0)
//    {
//        cpt_buffer=fread(buffer, sizeElement,buffer_size_nelem, binary_read_file);
//        if (cpt_buffer==0) return 0;
//    }
//    memcpy(element, (unsigned char *)buffer + (cpt_buffer-1) * sizeElement, sizeElement);
//    cpt_buffer --;
//    return cpt_buffer+1; // nb remaining before read
//}

// used to read/write raw information to the binary file (e.g. kmer count)

void GenomeDataStream::write(const void *element, int size)
{
    if (!fwrite(element, size, 1, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
}

void GenomeDataStream::write( void *element, int size)
{
    if (!fwrite(element, size, 1, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
}

size_t GenomeDataStream::read( void *element, int size)
{

    return fread(element, size,1, binary_read_file);
}


void GenomeDataStream::rewind_all()
{
    rewind(binary_read_file);
}

void GenomeDataStream::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    if(cpt_buffer)
    {
        if (!fwrite(buffer, sizeElement, cpt_buffer, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
    }
    cpt_buffer=0;

    fclose(binary_read_file);
}

void GenomeDataStream::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",buffer,write,filename);
        free(buffer);
        exit(1);
    }

}
void GenomeDataStream::End(){
    long i = ftell(binary_read_file);
    std::cout << "CUR: " << i << "\t" ;

    fseek(binary_read_file, 0, SEEK_END);
    i = ftell(binary_read_file);
    std::cout << "CUR: " << i;
}
off_t GenomeDataStream::file_size()
{
    return fsize(filename);
}

//off_t GenomeDataStream::nb_elements()
//{
//    return fsize(filename)/sizeElement;
//}

//
//size_t GenomeDataStream::read_kmer(void *element)
//{
//    return fread(element, kSizeOfKmerType, 1, binary_read_file);
//}
//
//size_t GenomeDataStream::read_kmer_colour(void *element_kmer, void* element_colour)
//{
//    size_t size = fread(element_kmer, kSizeOfKmerType, 1, binary_read_file);
//    size += fread(element_colour, kSizeOfKmerColour, 1, binary_read_file);
//    return size;
//}
//
//size_t GenomeDataStream::read_kmer_skip_colour( void *element)
//{
//    size_t size = fread(element, kSizeOfKmerType, 1, binary_read_file);
//    fseek(binary_read_file, kSizeOfKmerColour, SEEK_CUR);
//    return size;
//}
//
//size_t GenomeDataStream::read_colour( void *element)
//{
//    return fread(element, kSizeOfKmerColour, 1, binary_read_file);
//}

GenomeDataStream::~GenomeDataStream()
{
    if(buffer!=NULL)
    {
        free (buffer); //buffer =NULL;
    }


}

void GenomeDataStream::WriteHeader(uint64_t const &sequence_count) {
    uint64_t zero = 0;
    write(&zero, sizeof(uint64_t));
    write(&sequence_count, sizeof(uint64_t));
}


void GenomeDataStream::WriteHeader(uint64_t const &base_count, uint64_t const &sequence_count) {
    write(&base_count, sizeof(uint64_t));
    write(&sequence_count, sizeof(uint64_t));
}


void GenomeDataStream::WriteReference(const uint16_t &ref) {
    write(&ref, sizeof(uint16_t));
}

void GenomeDataStream::WriteReadDataKey(uint64_t const &key) {
    write(&key, sizeof(uint64_t));
}

void GenomeDataStream::ReadHeader(uint64_t &site_count0, uint64_t &read_data_count0) {
    read(&site_count, sizeof(uint64_t));
    read(&read_data_count, sizeof(uint64_t));

    if(read_data_count == 0){
        printf("Error: Something wrong with the header, read_data_count==0\n");
        exit(109);
    }
    if(site_count == 0){
        off_t total = file_size();
        total -= sizeof(uint64_t) * 2; //2 refs
        size_t each_model_input = sizeof(uint16_t) + sizeof(uint64_t) * read_data_count; //ref base + n*key

        site_count = total / each_model_input;
    }
    site_count0 = site_count;
    read_data_count0 = read_data_count;
}

void GenomeDataStream::ReadGenomeData(GenomeData &data) {

    if(site_count == 0 | read_data_count == 0){
        ReadHeader(site_count, read_data_count);
        
    }
    data.clear();
    data.reserve(site_count);

    for (size_t s = 0; s < site_count; ++s) {
        ModelInput m(s, read_data_count);
        ReadModelInput(m);
        data.push_back(std::move(m));
    }

}

void GenomeDataStream::ReadModelInput(ModelInput &model_input) {
    uint16_t ref;
    uint64_t key;
    ReadReference(ref);

    model_input.reference = ref;

    for (size_t i = 0; i < read_data_count; ++i) {
        ReadReadDataKey(key);
        model_input.all_reads[i].key = key;
//        model_input.all_reads.emplace_back(key);
    }

}

void GenomeDataStream::ReadReference(uint16_t &ref) {
    read(&ref, sizeof(uint16_t));
}

void GenomeDataStream::ReadReadDataKey(uint64_t &key) {
    read(&key, sizeof(uint64_t));
}

void GenomeDataStream::WriteModelInput(ModelInput &input) {

//    std::cout << input.reference << " : ";
    WriteReference(input.reference);
    for (auto &item : input.all_reads) {
//        std::cout << item.key << " : ";
        WriteReadDataKey(item.key);
    }
//    std::cout << "\n";

}

ModelInput GenomeDataStream::ReadOne() {

    uint16_t ref;
    uint64_t key;
    ReadReference(ref);
    ModelInput model_input(-1, read_data_count);
    model_input.reference = ref;

    for (size_t i = 0; i < read_data_count; ++i) {
        ReadReadDataKey(key);
        model_input.all_reads[i].key = key;
//        model_input.all_reads.emplace_back(key);
    }

    return model_input;
}
