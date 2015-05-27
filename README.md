#Calling mutations from MA lines


Travis: [![Build Status](https://travis-ci.org/stevenhwu/accuMUlate.svg?branch=EM_develop)](https://travis-ci.org/stevenhwu/accuMUlate)

Jenkins: [![Build Status](http://jenkins.scit.us/job/accuMUlate/badge/icon)](http://jenkins.scit.us/job/accuMUlate/)

<!---
your comment goes here
and here
-->
[//]: # (To improve platform compatibility (and to save one keystroke) it is also possible to use # (which is a legitimate hyperlink target) instead of <>:)

##Required libs

Needs to be compiled/on path:
* [bamtools](https://github.com/pezmaster31/bamtools)
* [eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [Boost](http://www.boost.org/)
  * [Boost::program_options](http://www.boost.org/doc/libs/1_58_0/doc/html/program_options.html)
  * [Boost::thread](http://www.boost.org/doc/libs/1_58_0/doc/html/thread.html)
  * [Boost::system](http://www.boost.org/doc/libs/1_58_0/libs/system/doc/index.html)

Provided in this repo
* bamtools headers  (for bamtools-utils, not in the normal include headers)

##Compilation

At present `cmake` handles a basic build. Presuming bamtools and eigein are on
your path you can do this

```sh
mkdir -p build
cd build
cmake ..
make
./accuMUlate -c test/test_params.ini -b test/test.bam -r test/test.fasta -o test/test.out
./pp -c test/pp_params.ini -b test/test.bam -i test/test.out 

```

The last line will  run the caller on a test dateset with 6 000 bases, and
show find mutations in each gene.

##TODO


* ~~have the main loop run the analysis from~~
    * ~~BAM file~~
    * ~~kist of sample IDs (one ancestor, rest descendants)~~
    * ~~Model Params:~~
        *~~theta~~
        *~~mu~~
        *~~nfreqs~~
        *~~phi-haploid~~
        *~~phi-diploid~~
    * ~~[BED file of coords to include]~~
    * ~~[probability threshold]~~
    * ~~[quality cut-off]~~
* ~~Boost PO/YAML/some other format to handle run variables~~
* ~~Manage build across platforms w/ CMake~~

* Optimize model code
* Clean up sample mapping/base calling
* Post-processor for putative mutations
* Create `.vcf` and annotations files from snps
* Make real test suite

