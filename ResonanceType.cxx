#include <iostream>
#include <cmath>
#include <string>
#include <array>
#include "ResonanceType.h"

ResonanceType::ResonanceType (std::string Name, double Mass, int Charge, double Width)
              : ParticleType(Name, Mass, Charge), fWidth_{Width} {};
              
double ResonanceType::GetWidth() const { return fWidth_; }
 
void ResonanceType::Print() const {

    std::cout << "----------" << '\n';
    ParticleType::Print();
    std::cout << "Width " << fWidth_ << '\n';
}

