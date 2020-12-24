#ifndef PARTICLE_H
#define PARTICLE_H 

#include<array>
#include<string>
#include "ResonanceType.h"

static const int fMaxNumParticleType=10;

class Particle
{
 public:
  Particle(std::string name, double Px, double Py, double Pz);
  
  Particle();
  
  int GetIndex() const;
  void SetIndex(int Index);
  void SetParticleID(std::string Name);
  
  static void AddParticleType(std::string name, double Mass, int Charge, double Width);
  
  void ParticleIdentifier() const;

  void PrintParticleID() const;
  static void PrintParticle();

  double GetPx() const;
  double GetPy() const;
  double GetPz() const;
  
  std::string GetName() const;
  double GetMass() const;
  int GetCharge() const;
  
  double GetEnergy() const;
  double InvMass(Particle &p);
  
  void SetP(double Px, double Py, double Pz);

  ~Particle();
  
  int Decay2body(Particle &dau1,Particle &dau2);

private:
  static std::array<ParticleType*, fMaxNumParticleType> fParticleType_;
  static int fNParticleType_;
  int fIndex_;
  double  fPx_, fPy_, fPz_;

  static int FindParticle(int fIndex);
  static int FindParticle(std::string Name);
  
  void Boost(double bx, double by, double bz);
  
  
};

#endif
