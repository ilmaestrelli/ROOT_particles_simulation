#include "particletype.h"
#include "resonancetype.h"

#include <iostream>

ResonanceType::ResonanceType(std::string name, double mass, int charge, double width) : ParticleType{name, mass, charge}, fWidth_{width} {};

void ResonanceType::Print() const
{
    ParticleType::Print();
    std::cout << "The width of the particle is " << fWidth_ << '\n';
}

double ResonanceType::GetWidth() const { return fWidth_; }

