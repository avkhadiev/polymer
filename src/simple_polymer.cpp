// 2017 Artur Avkhadiev
/*! \file simple_polymer.cpp
*/
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include "../include/vector.h"
#include "../include/simple_atom.h"
#include "../include/simple_bond.h"
#include "../include/default_macros.h"
#include "../include/simple_polymer.h"
namespace simple {
    /***************************************************************************
    *                               BASE POLYMER
    ***************************************************************************/
    int BasePolymer::_nb = DEFAULT_NB;
    double BasePolymer::_m = DEFAULT_M;
    double BasePolymer::_d = DEFAULT_D;
    BasePolymer::BasePolymer(Vector new_rcm, Vector new_vcm) :
        _rcm (new_rcm),
        _vcm (new_vcm){}
    BasePolymer::~BasePolymer(){}
    bool BasePolymer::operator==(const BasePolymer &other) const{
        return (_rcm == other.rcm()
            && _vcm == other.vcm());
    }
    bool BasePolymer::operator!=(const BasePolymer &other) const{
        return !(*this == other);
    }
    std::string BasePolymer::header_str(bool verbose) const{
        std::string m_str_rcm_header = "RCM = ";
        std::string m_str_rcm = vector_to_string(rcm());
        std::string m_str_vcm_header = "VCM = ";
        std::string m_str_vcm = vector_to_string(vcm());
        std::string m_str;
        if (verbose){
            m_str = m_str_rcm_header + m_str_rcm + "\n"
                + m_str_vcm_header + m_str_vcm;
         }
        else{
            m_str = m_str_rcm + " " + m_str_vcm;
        }
        return m_str;
    }
    /***************************************************************************
    *                               BOND POLYMER
    ***************************************************************************/
    BondPolymer::~BondPolymer(){}
    void BondPolymer::subtract_vel_projections(){
        for (int i = 0; i < nb(); ++i){
            // \dot{Omega} = \dot{Omega} - (\dot{Omega} \cdot \Omega) \Omega
            bonds.at(i).velocity
                -= multiply(bonds.at(i).position,
                    dot(bonds.at(i).velocity, bonds.at(i).position));
        }
    }
    BondPolymer::BondPolymer(const std::vector<Bond>& new_bonds,
        Vector rcm,
        Vector vcm) :
        BasePolymer(rcm, vcm),
        bonds (new_bonds)
        {
            // default value is empty
            if(bonds.empty()){
                bonds.resize(nb());
            }
            // all polymers have to have the same number of atoms and bonds
            else if(new_bonds.size() != nb()){
                fprintf(stderr, "%s (%d); %s\n",
                    "BondPolymer: given vector of bonds does not match the required size",
                    nb(),
                    "the vector will be resized.");
                bonds.resize(nb());
            }
            subtract_vel_projections();
        }
    bool BondPolymer::operator==(const BondPolymer &other) const{
        bool cm = (rcm() == other.rcm() && vcm() == other.vcm());
        bool molecules = true;
        for(int i = 0; i < nb(); ++i){
            molecules = molecules && (bonds.at(i) == other.bonds.at(i));
        }
        return cm && molecules;
    }
    bool BondPolymer::operator!=(const BondPolymer &other) const{
        return !(*this == other);
    }
    /***************************************************************************
    *                               ATOM POLYMER
    ***************************************************************************/
    AtomPolymer::~AtomPolymer(){}
    void AtomPolymer::subtract_average_vel(){
        Vector excess = vector(0.0, 0.0, 0.0);
        // calculate average velocity
        for (int i = 0; i < nb() + 1; ++i){
            // \dot{Omega} = \dot{Omega} - (\dot{Omega} \cdot \Omega) \Omega
            excess += atoms.at(i).velocity;
        }
        excess = divide(excess, nb() + 1);
        // subtract average velocity
        for (int i = 0; i < nb() + 1; ++i){
            // \dot{Omega} = \dot{Omega} - (\dot{Omega} \cdot \Omega) \Omega
            atoms.at(i).velocity -= excess;
        }
    }
    AtomPolymer::AtomPolymer(const std::vector<Atom>& new_atoms,
        Vector rcm,
        Vector vcm) :
        BasePolymer(rcm, vcm),
        atoms (new_atoms)
        {
            // default value is empty
            if(atoms.empty()){
                atoms.resize(nb() + 1);
            }
            // all polymers have to have the same number of atoms and bonds
            else if(new_atoms.size() - 1 != nb()){
                fprintf(stderr, "%s (%d + 1); %s\n",
                    "AtomPolymer: given vector of atoms does not match the required size",
                    nb(),
                    "the vector will be resized.");
                atoms.resize(nb() + 1);
            }
            subtract_average_vel();
        }
    bool AtomPolymer::operator==(const AtomPolymer &other) const{
        bool cm = (rcm() == other.rcm() && vcm() == other.vcm());
        bool molecules = true;
        for(int i = 0; i < nb() + 1; ++i){
            molecules = molecules && (atoms.at(i) == other.atoms.at(i));
        }
        return cm && molecules;
    }
    bool AtomPolymer::operator!=(const AtomPolymer &other) const{
        return !(*this == other);
    }
    /**************************************************************************/
    BondPolymer::BondPolymer(const AtomPolymer& polymer) :
        BasePolymer(polymer.rcm(), polymer.vcm()),
        bonds (polymer.nb())
        {
            // convert atoms to bonds
            for (int i = 0; i < nb(); ++i){
                bonds.at(i) = simple::Bond(polymer.atoms.at(i),
                    polymer.atoms.at(i+1));
            }
        }
    AtomPolymer::AtomPolymer(const BondPolymer& polymer) :
        BasePolymer(polymer.rcm(), polymer.vcm()),
        atoms (nb() + 1)
        {
            // r_0^{CM} = - \frac{d}{N+1}\sum_{i=1}^{N}(N+1-i)\Omega_i
            // v_0^{CM} = - \frac{d}{N+1}\sum_{i=1}^{N}(N+1-i)\dot{\Omega}_i
            Vector sum_pos = vector(0.0, 0.0, 0.0);
            Vector sum_vel = vector(0.0, 0.0, 0.0);
            for(int i = 0; i < nb(); ++i){
                sum_pos += multiply(polymer.bonds.at(i).position, nb() - i);
                sum_vel += multiply(polymer.bonds.at(i).velocity, nb() - i);
            }
            double frac = - 1.0 * d()/((double)(nb() + 1));
            atoms.at(0).position = multiply(sum_pos, frac);
            atoms.at(0).velocity = multiply(sum_vel, frac);
            // r_{i+1} = r_i + d\Omega_i
            // v_{i+1} = v_i + d\dot{\Omega}_i
            for(int i = 0; i < nb(); ++i){
                atoms.at(i+1).position = add(atoms.at(i).position,
                    multiply(polymer.bonds.at(i).position, d()));
                atoms.at(i+1).velocity = add(atoms.at(i).velocity,
                    multiply(polymer.bonds.at(i).velocity, d()));
            }
            // velocities are automatically in CM, no need to subtract average
            // velocity
        }
    std::string BondPolymer::to_string(bool verbose) const{
        std::string m_str_name = "SIMPLE BOND POLYMER";
        std::string m_str_bonds = "";
        for (int i = 0; i < nb(); ++i) {
            // verbose output
            if (verbose) {
                m_str_bonds += std::to_string(i)
                    + " "
                    +  bonds.at(i).to_string(verbose)
                    + "\n";
            }
            else {
            // non-verbose output
                m_str_bonds += bonds.at(i).to_string(verbose) + "\n";
            }
        }
        std::string m_str;
        if (verbose) {
            // verbose output
            m_str = m_str_name + "\n"
                + header_str(verbose) + "\n"
                + m_str_bonds;
        }
        else {
            // non-verbose output
            m_str = header_str(verbose) + "\n" + m_str_bonds;
        }
        return m_str;
    }
    std::string AtomPolymer::to_string(bool verbose) const{
        std::string m_str_name = "SIMPLE ATOM POLYMER:";
        std::string m_str_atoms = "";
        for (int i = 0; i < nb() + 1; ++i) {
            // verbose output
            if (verbose) {
                m_str_atoms += std::to_string(i)
                    + " "
                    +  atoms.at(i).to_string(verbose)
                    + "\n";
            }
            else {
            // non-verbose output
                m_str_atoms += atoms.at(i).to_string(verbose) + "\n";
            }
        }
        std::string m_str;
        if (verbose) {
            // verbose output
            m_str = m_str_name + "\n"
                + header_str(verbose) + "\n"
                + m_str_atoms;
        }
        else {
            // non-verbose output
            m_str = header_str(verbose) + "\n" + m_str_atoms;
        }
        return m_str;
    }
    ::std::ostream& operator<<(::std::ostream& os,
        const simple::BondPolymer& polymer){
        return os << polymer.to_string().c_str();
    }
    ::std::ostream& operator<<(::std::ostream& os,
        const simple::AtomPolymer& polymer){
        return os << polymer.to_string().c_str();
    }
    void BasePolymer::read_header(std::string non_verbose_header){
        // in non-verbose lines, the information about an atom is one line
        std::istringstream ss(non_verbose_header.c_str());
        std::istream_iterator<std::string> begin(ss);
        std::istream_iterator<std::string> end;
        std::vector<std::string> words(begin, end);
        // convert strings to data
        double rx = atof(words.at(0).c_str());
        double ry = atof(words.at(1).c_str());
        double rz = atof(words.at(2).c_str());
        double vx = atof(words.at(3).c_str());
        double vy = atof(words.at(4).c_str());
        double vz = atof(words.at(5).c_str());
        _rcm = vector(rx, ry, rz);
        _vcm = vector(vx, vy, vz);
    }
    BondPolymer string_to_bond_polymer(std::ifstream& input_stream){
        BondPolymer polymer;
        // data format in non-verbose lines :
        // ... header line ...
        // ... bond vector ... (nb lines)
        std::string line;
        // get header line
        std::getline(input_stream, line);
        polymer.read_header(line);
        // get all bonds
        for(int nlines = 0; nlines < BondPolymer::nb(); ++nlines) {
            std::getline(input_stream, line);
            polymer.bonds.at(nlines) = string_to_bond(line);
        }
        return polymer;
    }
    AtomPolymer string_to_atom_polymer(std::ifstream& input_stream){
        AtomPolymer polymer;
        // data format in non-verbose lines :
        // ... header line ...
        // ... atom vector ... (nb lines)
        std::string line;
        // get header line
        std::getline(input_stream, line);
        polymer.read_header(line);
        for(int nlines = 0; nlines < BondPolymer::nb() + 1; ++nlines) {
            std::getline(input_stream, line);
            polymer.atoms.at(nlines) = string_to_atom(line);
        }
        return polymer;
    }
}   // namespace simple
