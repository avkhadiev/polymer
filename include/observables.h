#ifndef POLYMER_INCLUDE_OBSERVABLES_H
#define POLYMER_INCLUDE_OBSERVABLES_H
#include <basics.h>
namespace observables {
inline namespace abstract {
    class Observable {
    private:
        basics::Units *_units;
        basics::Time *_time;
    public:
        virtual ~Observable();
        const basics::Units get_units();
        basics::Time get_time();
};
    template<class NumType> class ScalarObservable: public Observable {
    };
    template<class NumType> class VectorObservable: public Observable {
    };
}   // namespace abstract
}   // namespace observables
#endif
