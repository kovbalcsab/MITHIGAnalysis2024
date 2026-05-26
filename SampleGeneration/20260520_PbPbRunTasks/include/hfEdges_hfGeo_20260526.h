#pragma once

#include <vector>

namespace HFPFBinEdges {
inline double PhiShift18() { return 0.1745329252; }
inline double PhiStep18() { return 0.3490658504; }
inline double PhiStep36() { return 0.1745329252; }

// Eta edges derived from /afs/cern.ch/user/w/wangj/public/forBalazs/hfGeo.root
// using depth==1 eta_corner values. Adjacent HF rings overlap in the full
// corner envelope, so the 1D eta partition is made consistent by placing each
// shared boundary at the midpoint of the adjacent-ring overlap.
inline std::vector<double> EtaEdges18Negative() { return {-5.3221402, -4.9515372, -4.7781482}; }

inline std::vector<double> EtaEdges18Positive() { return {4.7781482, 4.9515372, 5.3221402}; }

inline std::vector<double> EtaEdges36Negative() {
  return {-4.7781482, -4.6008436, -4.4255354, -4.2526436, -4.0753205, -3.9017789,
          -3.7258617, -3.5513379, -3.3757032, -3.2006227, -3.0258648, -2.8456032};
}

inline std::vector<double> EtaEdges36Positive() {
  return {2.8456032, 3.0258648, 3.2006227, 3.3757032, 3.5513379, 3.7258617,
          3.9017789, 4.0753205, 4.2526436, 4.4255354, 4.6008436, 4.7781482};
}

// Phi edges already agree with the HF condition file at the 1e-7 rad level,
// so the existing 20260520 values are retained here.
inline std::vector<double> PhiEdges36() {
  return {0 - M_PI,           0.1745329252 - M_PI, 0.3490658504 - M_PI, 0.5235987756 - M_PI, 0.6981317008 - M_PI,
          0.872664626 - M_PI, 1.047197551 - M_PI,  1.221730476 - M_PI,  1.396263402 - M_PI,  1.570796327 - M_PI,
          1.745329252 - M_PI, 1.919862177 - M_PI,  2.094395102 - M_PI,  2.268928028 - M_PI,  2.443460953 - M_PI,
          2.617993878 - M_PI, 2.792526803 - M_PI,  2.967059728 - M_PI,  3.141592654 - M_PI,  3.316125579 - M_PI,
          3.490658504 - M_PI, 3.665191429 - M_PI,  3.839724354 - M_PI,  4.01425728 - M_PI,   4.188790205 - M_PI,
          4.36332313 - M_PI,  4.537856055 - M_PI,  4.71238898 - M_PI,   4.886921906 - M_PI,  5.061454831 - M_PI,
          5.235987756 - M_PI, 5.410520681 - M_PI,  5.585053606 - M_PI,  5.759586532 - M_PI,  5.934119457 - M_PI,
          6.108652382 - M_PI, 6.283185307 - M_PI};
}

inline std::vector<double> PhiEdges18Logical() {
  return {0.1745329252 - M_PI, 0.5235987804 - M_PI, 0.8726646308 - M_PI, 1.221730481 - M_PI, 1.570796332 - M_PI,
          1.919862182 - M_PI,  2.268928032 - M_PI,  2.617993883 - M_PI,  2.967059733 - M_PI, 3.316125584 - M_PI,
          3.665191434 - M_PI,  4.014257284 - M_PI,  4.363323135 - M_PI,  4.712388985 - M_PI, 5.061454836 - M_PI,
          5.410520686 - M_PI,  5.759586536 - M_PI,  6.108652387 - M_PI,  6.457718237 - M_PI};
}

inline std::vector<double> PhiEdges18Main() {
  return {0.1745329252 - M_PI, 0.5235987804 - M_PI, 0.8726646308 - M_PI, 1.221730481 - M_PI, 1.570796332 - M_PI,
          1.919862182 - M_PI,  2.268928032 - M_PI,  2.617993883 - M_PI,  2.967059733 - M_PI, 3.316125584 - M_PI,
          3.665191434 - M_PI,  4.014257284 - M_PI,  4.363323135 - M_PI,  4.712388985 - M_PI, 5.061454836 - M_PI,
          5.410520686 - M_PI,  5.759586536 - M_PI,  6.108652387 - M_PI};
}

inline std::vector<double> PhiEdges18WrapLow() { return {0 - M_PI, 0.1745329252 - M_PI}; }

inline std::vector<double> PhiEdges18WrapHigh() { return {6.108652377 - M_PI, 6.283185307 - M_PI}; }

} // namespace HFPFBinEdges
