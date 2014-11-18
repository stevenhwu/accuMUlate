#pragma once
#ifndef __EvolutionModel_H_
#define __EvolutionModel_H_


class EvolutionModel {

public:
    virtual ~EvolutionModel();

    virtual EvolutionModel() ;

    virtual void UpdateMu();
};


#endif //__EvolutionModel_H_
