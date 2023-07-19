#include "predictor_helpers.hpp"

namespace MethyLasso {
namespace predictor {

LinearResponse get_linear_response(const binned::Estimate& est, const binned::Binner& B) {
  return LinearResponse{B.project_back(est.get_phi())};
}

double get_max_absolute_deviation(const MeanVector& mu1, const MeanVector& mu2) {
  return (mu1.get_mu()-mu2.get_mu()).abs().maxCoeff();
}

}
}
