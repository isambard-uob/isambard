#include <math.h>
#include <string>

#include "dataClasses.hpp"

namespace bude {

    AtomDataFF::AtomDataFF(std::string AtomName, std::string ElectrostaticType, double Radius, double HydrophobicPotential,
                           double Hardness, double DistNpNp, double DistNpP, double RadiusScale,
                           double ElectrostaticPotential) {
         atomName = AtomName;
         electrostaticType = electrostaticType;
         radius = Radius;
         hydrophobicPotential = HydrophobicPotential;
         hardness = Hardness;
         distNpNp = DistNpNp;
         distNpP = DistNpP;
         radiusScale = RadiusScale;
         electrostaticPotential = ElectrostaticPotential;
    }

    AtomDataFF::~AtomDataFF(){ }

    std::vector<double> calculatePairEnergy(const AtomDataFF * atom1FF, const AtomDataFF * atom2FF,
                                  const double atomsDistance) {

        double stericEnergy, desolvationEnergy, chargeEnergy;
        double atomsRadius, electrostaticDistance;
        std::vector<double> splitEnergies;

        const double compZero = 0.0;
        const double columbicConstant = 22.5;
        const std::string dipoleChargeDC = "D";
        const std::string formalChargeDC = "F";
        const std::string donorAcceptorDC = "E";
        const std::string noRelevantDC = "-";

        stericEnergy = desolvationEnergy = chargeEnergy = 0.0;
        atomsRadius = electrostaticDistance = 0.0;

        atomsRadius = atom1FF->radius + atom2FF->radius;

        // Steric Energy
        if (atomsDistance < atomsRadius)
            stericEnergy = ((atom1FF->hardness + atom2FF->hardness) / 2) * (1 - atomsDistance / atomsRadius);

        // Electrostatic Energy
        if (atom1FF->electrostaticType != noRelevantDC && atom2FF->electrostaticType != noRelevantDC) {
            if (atom1FF->electrostaticType == formalChargeDC
                && atom2FF->electrostaticType == formalChargeDC)
                electrostaticDistance = 4.0;
            else
                electrostaticDistance = 2.0;
            if (atomsDistance < atomsRadius)
                chargeEnergy = columbicConstant * atom1FF->electrostaticPotential
                * atom2FF->electrostaticPotential;
            else if (atomsDistance < (atomsRadius + electrostaticDistance)) {
                chargeEnergy = columbicConstant * atom1FF->electrostaticPotential
                * atom2FF->electrostaticPotential;
                chargeEnergy -= chargeEnergy * (atomsDistance - atomsRadius)/ electrostaticDistance;
            }
            if ((atom1FF->electrostaticType == donorAcceptorDC || atom2FF->electrostaticType == donorAcceptorDC)
                && chargeEnergy > compZero)
                chargeEnergy = -chargeEnergy;
            }

        // Desolvation Energy
        // Non Polar - Polar
        if (atom1FF->hydrophobicPotential > compZero
            && atom2FF->hydrophobicPotential < compZero) {
            if (atomsDistance < atomsRadius)
                desolvationEnergy = (atom1FF->hydrophobicPotential - atom2FF->hydrophobicPotential) / 2;
            else if (atomsDistance < (atomsRadius + atom1FF->distNpP)) {
                desolvationEnergy = (atom1FF->hydrophobicPotential - atom2FF->hydrophobicPotential) / 2;
                desolvationEnergy -= desolvationEnergy * (atomsDistance - atomsRadius) / atom1FF->distNpP;
            }
            } else if (atom1FF->hydrophobicPotential < compZero
                       && atom2FF->hydrophobicPotential > compZero) {
        if (atomsDistance < atomsRadius)
            desolvationEnergy = (-atom1FF->hydrophobicPotential + atom2FF->hydrophobicPotential) / 2;
        else if (atomsDistance < (atomsRadius + atom2FF->distNpP)) {
            desolvationEnergy = (-atom1FF->hydrophobicPotential + atom2FF->hydrophobicPotential) / 2;
            desolvationEnergy -= desolvationEnergy * (atomsDistance - atomsRadius) / atom2FF->distNpP;
        }
        // Non Polar - Non Polar
        } else if (atom1FF->hydrophobicPotential < compZero && atom2FF->hydrophobicPotential < compZero) {
            if (atomsDistance < atomsRadius) {
                desolvationEnergy = (atom1FF->hydrophobicPotential + atom2FF->hydrophobicPotential) / 2;
            } else if (atomsDistance < (atomsRadius + atom1FF->distNpNp)) {
                desolvationEnergy = (atom1FF->hydrophobicPotential + atom2FF->hydrophobicPotential) / 2;
                desolvationEnergy -= desolvationEnergy * (atomsDistance - atomsRadius) / atom1FF->distNpNp;
                }
        }

        splitEnergies.push_back(stericEnergy);
        splitEnergies.push_back(desolvationEnergy);
        splitEnergies.push_back(chargeEnergy);

        return splitEnergies;
    }
}
