#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H 

#include "ParticleType.h"

class ResonanceType : public ParticleType
{
public:
  ResonanceType(std::string Name, double Mass, int Charge, double Width);
  double GetWidth() const;
  void Print() const;

private:
  const double fWidth_;
};

#endif

