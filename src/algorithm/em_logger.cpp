/*
 * em_logger.cc
 *
 *  Created on: 4/9/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/9/15.
//

#include <iostream>
#include "em_logger.h"

EmLogger::EmLogger(){
}

void EmLogger::Log(std::string log){
    if(can_log) {
        outfile << log;
    }
}

void EmLogger::LogLine(std::string log){
    if(can_log) {
        outfile << log << "\n";
        outfile.flush();
    }
    else{
        std::cout << "Warning: No Header" << std::endl;
    }
}

void EmLogger::Stop() {
    outfile.close();
}

void EmLogger::Header(std::string header) {
    can_log = true;
    LogLine(header);
}

void EmLogger::SetOutifle(std::string const &outfile_prefix0) {
    outfile_prefix = outfile_prefix0;
        outfile.open(outfile_prefix);
}

void EmLogger::LogLine(int i, std::string line) {
    outfile << i << "\t";
    LogLine(line);

}
