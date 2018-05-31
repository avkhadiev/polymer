// 2018 Artur Avkhadiev
/*! \file geodesic_path.h
*/
#ifndef GEODESIC_PATH_H
#define GEODESIC_PATH_H
#include <string>
#include <list>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include "default_macros.h"
#include "geodesic_record.h"
#include "geodesic_observables.h"
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC PATH
    ***************************************************************************/
    class Path {
    protected:
        Record _initial;                         /**> head                  */
        Record _final;                           /**> tail                  */
        std::list<Record> _path;                 /**> omits boundary values */
    public:
        /* path length */
        /* TODO add Euclidian endpoint distance calculation upon construction */
        //class Length {
        //     friend class Path;
        // private:
        //bool _calculate_diag_length;
        //// TODO
        ///** euclidian distance, including and excluding cross terms      */
        //double _euclidian_full; /**> full euclidian distance             */
        //double _euclidian_diag; /**> diagonal-term euclidian distance    */
        ///** full kinematic length, including and excluding cross terms   */
        //double _ell_full;       /**> full length                         */
        //double _ell_diag;       /**> diagonal-term (SLERP) length        */
        ///** kinematic length of solvents and polymers separately         */
        //double _g_pol_diag;     /**> polymer component, diagonal terms   */
        //double _g_pol_full;     /**> polymer component, all terms        */
        //    // further break down of polymer into
        //    // translational/rotational
        //double _g_pol_trans;    /**> translational component (polymer)   */
        //double _g_pol_rot_full; /**> rotational component, all terms     */
        //double _g_pol_rot_diag; /**> rotational component, diagonal terms*/
        //double _g_sol;          /**> solvent component                   */
        ///** corresponding square accumulators                            */
        //double _g_pol_trans_sq;
        //double _g_pol_rot_full_sq;
        //double _g_pol_rot_diag_sq;
        //double _g_sol_sq;
        ////** auxilliary quanities                                        */
        //double _I0_over_M;      /**> (md^2)/m (polymer)                  */
        public:
        //  bool calculate_diag_length() const {return _calculate_diag_length;};
         //  /** euclidian distance, including and excluding cross terms       */
        //  double euclidan() const {return _euclidian_full;};
         //  double euclidian_diag() const {return _euclidian_diag;}
        //  /** full kinematic length, including and excluding cross terms    */
         //  double ell() const {return _ell_full;};
        //  double ell_diag() const {return _ell_diag;};
         //  /** kinematic length of solvents and polymers separately          */
        //  double g_polymer() const {return _g_pol_full;};
         //  double g_polymer_diag() const {return _g_pol_diag;};
        //  double g_solvent() const {return _g_sol;};
         //  /** further break down of polymer into translational/rotational   */
        //  double g_polymer_trans() const {return _g_pol_trans;};
         //  double g_polymer_rot() const {return _g_pol_rot_full;};
        //  double g_polymer_rot_diag() const {return _g_pol_rot_diag;};
            /** FIXME TODO */
            /** g_pol_tran g_pol_rot g_pol_rot_diag g_sol_trans */
        //    std::string to_string() const;
        private:
            /** g_pol_tran g_pol_rot g_pol_rot_diag g_sol_trans */
           // void _update_from_string(std::string str);
           //void _zero();            /**> zero accumulators & their squares  */
           //void _zero_sq_increments(); /**> zero squares of accumulators    */
          // void _increment(         /**> increment length at new step   */
        //        simple::BondState& last_state,
        //        simple::BondState& next_state);
           //void _calculate_square_increments(  /**> at new timestep        */
           //    simple::BondState& last_state,
           //    simple::BondState& next_state
           //);
           //void _calculate_g_pol_trans_sq_increment(
           //    simple::BondState& last_state,
           //    simple::BondState& next_state);
           //void _calculate_g_pol_rot_sq_increment(
           //    simple::BondState& last_state,
           //    simple::BondState& next_states);
           //void _calculate_g_sol_sq_increment(
           //    simple::BondState& last_state,
           //    simple::BondState& next_state);
            /** returns (r_next - r_last)^2 */
            //double _trans_sq_increment(
            //    const Vector& r_last, const Vector& r_next);
            /**
            * returns del Psi_i del Psi_j p_ij, where
            *   del Psi = arccos(omega_next * omega_last)
            *   p_ij = del_Omega_i * del_Omega_j, where
            *       del Omega = omega_next - omega_last
            */
            // double _rot_sq_increment(
            //    const Vector& omega_last_i, const Vector& omega_next_i,
            //    const Vector& omega_last_j, const Vector& omega_next_j);
        public:
            // TODO
        //    Length();
        //    Length( simple::BondState& first,
        //            simple::BondState& last,
        //            bool calculate_diag_length = true);
        //    Length( simple::BondState& first,
        //            simple::BondState& last,
        //            double g_pol_trans,
        //            double g_pol_rot,
        //            double g_pol_rot_diag,
        //            double g_sol_trans,
        //            bool calculate_diag_length = true);
        //    ~Length();
        //};
         Length length;
        Record initial() const {return _initial;};  /**> boundary value 1   */
        Record final() const {return _final;};      /**> boundary value 2   */
    public:
        Path& operator= (const Path &fraction);
        size_t nrecords() const;                    /** number of states    */
        // NO TIME IN CONFIGURATION SPACE
        //double time() const;                      /**> abs(t_f - t_i)     */
        /* operations on a path */
        /* current tail is the record preceding the final configuration */
        Record current_tail() const;
        /* returns a copy of the entire path, including boundary values     */
        std::list<Record> get_path() const;
        Length* get_length() {return &length;};
        void append(const Record& new_record);
        void reverse();                             /**> A->B to B->A       */
        /* TODO */
        void recompute_length();                    /**> may take a while   */
        /* TODO */
        void merge(const Path& new_path);           /**> append new path    */
        /* if overwrite = true, truncates and outputs header string */
        void write(std::string file, bool overwrite);
    protected:
        // tracking geodesic length
        /* TODO */
        void _increment_length( const Record& last_record,
                                const Record& next_record);
        // NO TIME IN CONFIGURATION SPACE
        /**
        * if t_f - t_i > 0, then preceding.t < next.t < final.t
        * if t_f - t_i < 0, then preceding.t > next.t > final.t
        * if t_f - t_i = 0, there can be no intermediate values
        */
        // bool _check_time_flow(  const Record& last_record,
        //                        const Record& next_record) const;
        /* performs _check_time_flow(record) along the entire path */
        // bool _check_time_flow() const;
        /* state header string                                              */
        virtual std::string _header_str() const;
        virtual void _read_header(std::ifstream& readout);
    public:
        Path();
        Path(Record initial, Record final);         /**> minimal path       */
        Path(std::list<Record>& records);           /**> general path       */
        Path(std::string file);                     /**> read from file     */
        ~Path();
    };
} // namespace geodesic
#endif
