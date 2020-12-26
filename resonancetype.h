#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H

#include <iostream>

class ResonanceType : public ParticleType
{
public:
   ResonanceType(std::string name, double mass, int charge, double width); 
   
   double GetWidth() const override;
   void Print() const override; 

private: 
    const double fWidth_;
};

#endif
