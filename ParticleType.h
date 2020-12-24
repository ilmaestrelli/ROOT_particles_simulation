#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H 

class ParticleType
{
public:
  ParticleType(const char* Name, const double Mass, const int Charge);
  const char* GetName() const;
  double GetMass() const;
  int GetCharge() const;
  virtual double GetWidth() const;
  void Print() const;

private:
  const char* fName_;
  const double fMass_;
  const int fCharge_;
  const double fWidth_ = 0;
};

#endif 
