#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <cmath>
#include "particletype.h"
#include "resonancetype.h"

class Particle
{
public:
    Particle(std::string name, double fPx, double fPy_, double fPz_);
    Particle() = default;

    static void AddParticleType(std::string name, double mass, int charge, double width = 0);

    int Decay2body(Particle &dau1, Particle &dau2) const;

    int GetIndex() const;
    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    double GetImpulse() const;

    double GetMass() const;
    double GetEnergy() const;
    double GetInvMass(Particle &p);

    void SetIndex(int index);
    void SetParticleID(std::string name);
    void SetImpulse(double Px, double Py, double Pz);
    
    void PrintParticle() const;
    void static PrintArray();

private:
    static int FindParticle(std::string name);

    void Boost(double bx, double by, double bz);
    
    
    int fIndex_ = 0;
    double fPx_ = 0;     
    double fPy_ = 0;
    double fPz_ = 0;

    static const int fMaxNumParticleType_ = 10;                
    static ParticleType *fParticleType_[fMaxNumParticleType_]; 
    static int fNParticleType_;                                

};

#endif
