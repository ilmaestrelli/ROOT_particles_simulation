#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

#include <iostream>

class ParticleType
{
public:
   ParticleType(std::string name, double mass, int charge); 
    
   std::string GetName() const; 
   double GetMass() const;
   int GetCharge() const;
   virtual double GetWidth() const;

   virtual void Print() const; 

private: 
   std::string const fName_;
   const double fMass_;
   const int fCharge_;
};

#endif
