#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H 

class ParticleType
{
public:
  ParticleType(std::string const Name, const double Mass, const int Charge);
  std::string GetName() const;
  double GetMass() const;
  int GetCharge() const;
  virtual double GetWidth() const;
  virtual void Print() const;

private:
  std::string const fName_;
  const double fMass_;
  const int fCharge_;
  const double fWidth_ = 0;
};

#endif 
