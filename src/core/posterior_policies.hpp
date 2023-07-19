#ifndef POSTERIOR_POLICIES_HPP
#define POSTERIOR_POLICIES_HPP

namespace MethyLasso {

/// @brief The Basic policy never orders damping of coefficients
/*!
 * Note: Use BasicIRLSPolicy or BasicHyperPolicy to pass it to the Posterior
 */
class BasicPolicy {
protected:
  bool condition_check(unsigned irls_step, unsigned damp_step, double mlogp, double old_mlogp) const {
    return false;
  }
  double get_new_param_weight() const { return 1; }
};

class BasicIRLSPolicy : public BasicPolicy {};

class BasicPriorPolicy : public BasicPolicy {};
                                            
class BasicHyperPolicy : public BasicPolicy {};
            
                                            
// @brief The fixed damping policy orders damping of coefficients after a number of IRLS steps have passed
/*! Damping results to setting the new parameter to \f$(f*p_{old}+(1-f)*p_{new})/10\f$, where \f$p_{new}\f$ 
 *  (resp. \f$p_{old}\f$) is the new (resp. old) parameter estimate.
 *  @tparam damp_after the number of steps after which damping occurs
 *  @tparam damp_factor the damping strength \f$f\f$, between 1 and 9.
 *  
 * Note: Use FixedDampingIRLSPolicy or FixedDampingHyperPolicy to pass it to the Posterior
 */
template<unsigned damp_after = 5, unsigned damp_factor = 5>
class FixedDampingPolicy {
protected:
  bool condition_check(unsigned irls_step, unsigned damp_step, double mlogp, double old_mlogp) const {
    if (damp_step>=1) return false;
    return irls_step > damp_after;
  }
   
   // the higher the factor, the closer we stay to the previous value
  double get_new_param_weight() const { return 1-damp_factor/10.; }
};

template<unsigned damp_after = 5, unsigned damp_factor = 5>
struct FixedDampingIRLSPolicy : public FixedDampingPolicy<damp_after,damp_factor> {};

template<unsigned damp_after = 5, unsigned damp_factor = 5>
struct FixedDampingPriorPolicy : public FixedDampingPolicy<damp_after,damp_factor> {};

template<unsigned damp_after = 5, unsigned damp_factor = 5>
struct FixedDampingHyperPolicy : public FixedDampingPolicy<damp_after,damp_factor> {};



/// @brief The adaptive damping policy orders damping of coefficients when the posterior gets worse
template<unsigned damp_after = 5, unsigned damp_factor = 5, unsigned max_reductions = 10>
class AdaptiveDampingPolicy {
protected:
  bool condition_check(unsigned irls_step, unsigned damp_step, double mlogp, double old_mlogp) const {
    if (damp_step>=max_reductions) return false;
    return mlogp > old_mlogp;
  }
  
  // the higher the factor, the closer we stay to the previous value
  double get_new_param_weight() const { return 1-damp_factor/10.; }

};

template<unsigned damp_after = 5, unsigned damp_factor = 5>
class AdaptiveDampingIRLSPolicy : public AdaptiveDampingPolicy<damp_after,damp_factor> {};

template<unsigned damp_after = 5, unsigned damp_factor = 5>
class AdaptiveDampingPriorPolicy : public AdaptiveDampingPolicy<damp_after,damp_factor> {};

template<unsigned damp_after = 5, unsigned damp_factor = 5>
class AdaptiveDampingHyperPolicy : public AdaptiveDampingPolicy<damp_after,damp_factor> {};


}

#endif
