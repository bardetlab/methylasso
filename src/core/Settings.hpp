#ifndef SETTINGS_HPP
#define SETTINGS_HPP

/*! \file Settings.hpp
 * @brief TODO
 */

namespace MethyLasso {
namespace base {

//class to store all settings used in the c++ side 
//to add settings for a new class, forward-declare the class
//and then specialize the Settings template on it
template<typename T> struct Settings {};

//settings for GFLLibraryImpl
class GFLLibraryImpl;
template<> struct Settings<GFLLibraryImpl> {
    //initial adamts step size
    static const double get_alpha() { return 5.; }
    //inflation factor for adamts step update
    static const double get_inflate() { return 2.; }
    //maximum number of adamts steps
    static const int get_ninner_max() { return 30000; }
};


//settings for FusedLassoGaussianEstimator
template<typename Library> class FusedLassoGaussianEstimator;
template<typename Library> struct Settings<FusedLassoGaussianEstimator<Library> > {
    //at which value to clamp computed beta. Negative to turn off
    static const double get_clamp() { return 35; } // 1-inv_logit(x) != 0 when x < -36.7
};

}
}

#endif

