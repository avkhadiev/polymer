// 2018 Artur Avkhadiev
/*! \file geodesic_path_analyzer.h
*/
#ifndef GEODESIC_PATH_ANALYZER_H
#define GEODESIC_PATH_ANALYZER_H
namespace geodesic{
    /***************************************************************************
    *                          GEODESIC PATH ANALYZER
    ***************************************************************************/
    class PathAnalyzer{
    public:
        PathAnalyzer(double landscape_energy);
        ~PathAnalyzer();
        double landscape_energy() const {return E_l;};
        double set_landscape_energy(double energy);
    private:
        double E_l;
    };
} // namespace geodesic

#endif
