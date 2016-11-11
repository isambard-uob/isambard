#ifndef DATACLASSES_HPP
#define DATACLASSES_HPP

#include <string>

namespace bude {

    class AtomDataFF {
    public:
        std::string atomName;
        std::string electrostaticType;
        double radius;
        double hydrophobicPotential;
        double hardness;
        double distNpNp;
        double distNpP;
        double radiusScale;
        double electrostaticPotential;

        AtomDataFF(std::string atomName, std::string electrostaticType, double radius, double hydrophobicPotential,
                   double hardness, double distNpNp, double distNpP, double radiusScale,
                   double electrostaticPotential);
        ~AtomDataFF();
    };
}
#endif
