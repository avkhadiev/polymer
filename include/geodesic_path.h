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
        Length length;
        Record initial() const {return _initial;};  /**> boundary value 1   */
        Record final() const {return _final;};      /**> boundary value 2   */
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
        void recompute_length();                    /**> may take a while   */
        void merge(const Path& new_path);           /**> append new path    */
        /* if overwrite = true, truncates and outputs header string */
        void write(std::string file, bool overwrite);
        /* allows to slice the path in n parts and get the mth slice        */
        /* n should be less than or equal to one less the number of records */
        /* m should be less than or equal to n                              */
        Path get_slice(size_t n_slices, size_t mth_slice); 
    protected:
        void _increment_length( const Record& last_record,
                                const Record& next_record );
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
