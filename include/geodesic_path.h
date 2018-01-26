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
namespace geodesic{
    /***************************************************************************
    *                               GEODESIC PATH
    ***************************************************************************/
    class Path {
    public:
        Path();
        Path(Record initial, Record final);         /**> minimal path       */
        Path(std::list<Record>& records);           /**> general path       */
        Path(std::string file);                     /**> read from file     */
        ~Path();
        Record initial() const {return _initial;};  /**> boundary value 1   */
        Record final() const {return _final;};      /**> boundary value 2   */
        /* path length */
        double kinematic_length() const {return _ell;};
        size_t nrecords() const;                    /** number of states    */
        double time() const;                        /**> abs(t_f - t_i)     */
        /* operations on a path */
        /* current tail is the record preceding the final configuration */
        const Record current_tail() const;
        void append(const Record& new_record);
        void reverse();                             /**> A->B to B->A       */
        void recompute_length();                    /**> may take a while   */
        Path& operator= (const Path &fraction);
        void reset_length() {_ell = 0.0;};
        /* TODO */
        void merge(const Path& new_path);           /**> append new path    */
        /* if overwrite = true, truncates and outputs header string */
        void write(std::string file, bool overwrite);
    protected:
        Record _initial;                            /**> head               */
        Record _final;                              /**> tail               */
        double _ell;                                /**> kinematic length   */
        virtual void _increment_length( const Record& preceding_record,
                                        const Record& next_record);
        std::list<Record> _path;               /**> omits boundary values */
        /* TODO */
        /**
        * if t_f - t_i > 0, then preceding.t < next.t < final.t
        * if t_f - t_i < 0, then preceding.t > next.t > final.t
        * if t_f - t_i = 0, there can be no intermediate values
        */
        bool _check_time_flow(  const Record& preceding_record,
                                const Record& next_record) const;
        /* performs _check_time_flow(record) along the entire path */
        bool _check_time_flow() const;
        /* (1st line) kinematic_length nrecords time         */
        /* (2nd line) state header string                    */
        virtual std::string _header_str() const;
        virtual void _read_header(std::ifstream& readout);
    };
} // namespace geodesic
#endif
