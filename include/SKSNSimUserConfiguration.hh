/************************************
 * File: SKSNSimUserConfiguration.hh
 ***********************************/

#ifndef __SKSNSIMUSERCONFIGURATION_H_INCLUDED__
#define __SKSNSIMUSERCONFIGURATION_H_INCLUDED__

class SKSNSimUserConfiguration{
  private:
    double m_energy_min;
    double m_energy_max;
    double m_time_min;
    double m_time_max;

  public:
    SKSNSimUserConfiguration(){}
    ~SKSNSimUserConfiguration(){}

    SKSNSimUserConfiguration &SetFluxEnergyMin(double e){ m_energy_min = e; return *this;}
    SKSNSimUserConfiguration &SetFluxEnergyMax(double e){ m_energy_max = e; return *this;}
    double GetFluxEnergyMin() const { return m_energy_min;}
    double GetFluxEnergyMax() const { return m_energy_max;}

};


#endif
