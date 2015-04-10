/*
 * em_logger.cc
 *
 *  Created on: 4/9/15
 *      Author: Steven Wu
 */

//
// Created by steven on 4/9/15.
//

#pragma once
#ifndef _EM_LOGGER_H_
#define _EM_LOGGER_H_


#include <string>
#include <fstream>

class EmLogger {


public:
    EmLogger();

    void Stop();

    void LogLine(std::string log);

    void Header(std::string header);

    void SetOutifle(std::string const &outfile_prefix0);

    void LogLine(int i, std::string string);

private:
    std::ofstream outfile;
    std::string outfile_prefix;

    void Log(std::string log);
    bool can_log;
};


#endif //_EM_LOGGER_H_
