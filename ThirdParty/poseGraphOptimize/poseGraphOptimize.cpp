#include "poseGraphOptimize.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include "coder_bounded_array.h"
#include "coder_posix_time.h"
#include "cs.h"
#include "makeCXSparseMatrix.h"
#include "rt_defines.h"
#include "solve_from_lu.h"
#include "solve_from_qr.h"
#include <algorithm>
#include <cmath>
#include <cstring>

namespace poseGraphOptimize {
namespace coder {
class sparse;

}
} // namespace poseGraphOptimize

namespace poseGraphOptimize {
namespace coder {
namespace robotics {
namespace core {
namespace internal {
enum class NLPSolverExitFlags : int
{
  LocalMinimumFound = 1,
  IterationLimitExceeded,
  TimeLimitExceeded,
  StepSizeBelowMinimum,
  ChangeInErrorBelowMinimum,
  SearchDirectionInvalid,
  HessianNotPositiveSemidefinite,
  TrustRegionRadiusBelowMinimum
};

}
} // namespace core
} // namespace robotics
} // namespace coder
struct struct_T {
  int xstart;
  int xend;
  int depth;
};

struct b_struct_T {
  ::coder::array<double, 2U> edgeNodePairs;
  ::coder::array<double, 2U> edgeMeasurements;
  ::coder::array<double, 2U> edgeInfoMats;
  double tformSize[2];
  double infoMatSize[2];
  double poseDeltaLength;
  ::coder::array<double, 1U> nodeMap;
  ::coder::array<double, 1U> nodeDims;
  ::coder::array<boolean_T, 1U> IsLandmarkNode;
};

struct c_struct_T {
  ::coder::array<int, 1U> a;
  ::coder::array<int, 1U> b;
};

namespace coder {
namespace robotics {
namespace core {
namespace internal {
class BlockMatrix {
public:
  void replaceBlock(double i, const double blockij[9]);
  void extractBlock(double i, ::coder::array<double, 2U> &B) const;
  void replaceBlock();
  ::coder::array<double, 2U> Matrix;
  double NumRowBlocks;
  double NumColBlocks;
  double BlockSize[2];
};

} // namespace internal
} // namespace core
} // namespace robotics
class poseGraph {
public:
  void addRelativePose(const double varargin_1[3], double varargin_3,
                       double varargin_4);
  void findEdgeID(const double nodePair[2],
                  ::coder::array<double, 1U> &edgeID) const;
  double NumNodes;
  double NumEdges;
  double NumLoopClosureEdges;
  robotics::core::internal::BlockMatrix *NodeEstimates;
  ::coder::array<double, 1U> NodeMap;
  ::coder::array<double, 1U> NodeDims;
  ::coder::array<boolean_T, 1U> IsLandmarkNode;
  ::coder::array<double, 2U> EdgeNodePairs;
  ::coder::array<double, 2U> LoopClosureEdgeNodePairs;
  robotics::core::internal::BlockMatrix *EdgeMeasurements;
  robotics::core::internal::BlockMatrix *EdgeInfoMatrices;
  ::coder::array<double, 2U> LoopClosureEdgeIDsInternal;
  double MaxNumEdges;
  boolean_T __MaxNumEdges_AssignmentSentinel;
  double MaxNumNodes;
  boolean_T __MaxNumNodes_AssignmentSentinel;
};

namespace robotics {
namespace core {
namespace internal {
class SystemTimeProvider {
public:
  coderTimespec StartTime;
};

class TrustRegionIndefiniteDogLegSE2 {
public:
  double solve(const ::coder::array<double, 2U> &seed, BlockMatrix &iobj_0,
               BlockMatrix **xSol, sparse &hess, double &solutionInfo_Error,
               double &solutionInfo_ExitFlag);
  boolean_T computeBasicSteps(const ::coder::array<double, 1U> &grad,
                              const sparse &B,
                              ::coder::array<double, 1U> &stepSD,
                              ::coder::array<double, 1U> &stepGN) const;

protected:
  void incrementX(const ::coder::array<double, 2U> &x,
                  const ::coder::array<double, 1U> &epsilons,
                  ::coder::array<double, 2U> &xNew) const;

public:
  b_struct_T ExtraArgs;
  double MaxNumIteration;
  double MaxTime;
  ::coder::array<double, 2U> SeedInternal;
  double MaxTimeInternal;
  double MaxNumIterationInternal;
  double StepTolerance;
  SystemTimeProvider TimeObj;
  double GradientTolerance;
  double FunctionTolerance;
  double InitialTrustRegionRadius;
  double TrustRegionRadiusTolerance;
};

} // namespace internal
} // namespace core
} // namespace robotics
class sparse {
public:
  void ctranspose(sparse &y) const;
  void fillIn();
  ::coder::array<double, 1U> d;
  ::coder::array<int, 1U> colidx;
  ::coder::array<int, 1U> rowidx;
  int m;
  int n;
};

namespace nav {
namespace algs {
namespace internal {
class BlockInserter2 {
public:
  void insertGradientBlock(double i, const double blocki_data[],
                           int blocki_size);
  ::coder::array<double, 1U> Gradient;
  ::coder::array<double, 1U> NodeDims;
  ::coder::array<double, 1U> NodeMap;
  ::coder::array<double, 2U> HessianCSC;
  double HessianCSCCount;
};

} // namespace internal
} // namespace algs
} // namespace nav
class anonymous_function {
public:
  c_struct_T workspace;
};

namespace robotics {
namespace core {
namespace internal {
class b_BlockMatrix {
public:
  ::coder::array<double, 2U> Matrix;
  double BlockSize[2];
};

} // namespace internal
} // namespace core
} // namespace robotics
namespace internal {
class stack {
public:
  ::coder::bounded_array<struct_T, 120U, 1U> d;
  int n;
};

} // namespace internal
namespace nav {
namespace algs {
namespace internal {
class PoseGraphOptimizer {
public:
  static double parseOptimizePoseGraphInputs(
      double &paramStruct_MaxTime, double &paramStruct_FunctionTolerance,
      boolean_T &paramStruct_IsVerbose, double &paramStruct_GradientTolerance,
      double &paramStruct_StepTolerance,
      double &paramStruct_InitialTrustRegionRadius,
      double paramStruct_FirstNodePose[3],
      double &paramStruct_TrustRegionRadiusTolerance,
      double &paramStruct_SolverID);
};

class PoseGraphHelpers {
public:
  static double
  poseGraphCost(const ::coder::array<double, 2U> &posesMat,
                const ::coder::array<double, 2U> &args_edgeNodePairs,
                const ::coder::array<double, 2U> &args_edgeMeasurements,
                const ::coder::array<double, 2U> &args_edgeInfoMats,
                const double args_tformSize[2],
                const double args_infoMatSize[2], double args_poseDeltaLength,
                const ::coder::array<double, 1U> &args_nodeMap,
                const ::coder::array<double, 1U> &args_nodeDims,
                const ::coder::array<boolean_T, 1U> &args_IsLandmarkNode,
                ::coder::array<double, 1U> &gradient, sparse &hessian);
  static double costBetweenTwoNodes(
      const ::coder::array<double, 2U> &Toi,
      const ::coder::array<double, 2U> &Toj,
      const ::coder::array<double, 2U> &measurement,
      const ::coder::array<double, 2U> &Omega, boolean_T nodejIsLandmark,
      double gradi_data[], int &gradi_size, double gradj_data[],
      int &gradj_size, double hessii_data[], int hessii_size[2],
      double hessij_data[], int hessij_size[2], double hessji_data[],
      int hessji_size[2], double hessjj_data[], int hessjj_size[2]);
};

} // namespace internal
} // namespace algs
} // namespace nav
namespace robotics {
namespace core {
namespace internal {
class SEHelpers {
public:
  static void veelogmSE3(const double T[16], double vec[6]);
  static void expSE3hat(const double e[6], double T[16]);
};

class Sim3Helpers {
public:
  static void multiplyLogSim3(const double S1[16], const double S2[16],
                              const double S3[16], double e[7]);
  static void sim3ToSform(const double minVecSim3[7], double S[16]);
};

} // namespace internal
} // namespace core
} // namespace robotics
namespace internal {
class CXSparseAPI {
public:
  static void iteratedQR(const sparse &A, const ::coder::array<double, 1U> &b,
                         int n, ::coder::array<double, 1U> &out);
};

} // namespace internal
} // namespace coder
} // namespace poseGraphOptimize

namespace poseGraphOptimize {
static double freq;

static boolean_T freq_not_empty;

static const signed char iv[9]{1, 0, 0, 0, 1, 0, 0, 0, 1};

static boolean_T isInitialized_poseGraphOptimize{false};

} // namespace poseGraphOptimize

namespace poseGraphOptimize {
static void binary_expand_op(coder::nav::algs::internal::BlockInserter2 *in1,
                             int in2, int in4, int in5, const double in6_data[],
                             const int &in6_size);

static double binary_expand_op(double in1, double in2,
                               const ::coder::array<double, 1U> &in3,
                               const ::coder::array<double, 1U> &in4,
                               const ::coder::array<double, 1U> &in5);

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in2, double in3);

static void binary_expand_op(double in1[3],
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 1U> &in3, int in4,
                             int in5);

static void binary_expand_op(::coder::array<double, 2U> &in1,
                             const ::coder::array<double, 1U> &in2, int in3,
                             int in4);

namespace coder {
static double b_norm(const ::coder::array<double, 1U> &x);

static void b_sparse(const ::coder::array<double, 1U> &varargin_1,
                     const ::coder::array<double, 1U> &varargin_2,
                     const ::coder::array<double, 1U> &varargin_3, sparse &y);

namespace internal {
static void b_heapsort(::coder::array<int, 1U> &x, int xstart, int xend,
                       const anonymous_function &cmp);

namespace blas {
static void mtimes(const double A_data[], const int A_size[2],
                   const double B_data[], const int B_size[2], double C_data[],
                   int C_size[2]);

static void mtimes(const double A_data[], const int A_size[2],
                   const ::coder::array<double, 2U> &B, double C_data[],
                   int C_size[2]);

static void xaxpy(double a, const double x[3], double y[9], int iy0);

static void xaxpy(int n, double a, int ix0, double y[9], int iy0);

static void xaxpy(double a, const double x[9], int ix0, double y[3]);

static double xdotc(int n, const double x[9], int ix0, const double y[9],
                    int iy0);

static double xnrm2(const double x[3]);

static double xnrm2(int n, const double x[9], int ix0);

static void xrot(double x[9], int ix0, int iy0, double c, double s);

static double xrotg(double &a, double &b, double &s);

static void xswap(double x[9], int ix0, int iy0);

} // namespace blas
static void heapify(::coder::array<int, 1U> &x, int idx, int xstart, int xend,
                    const anonymous_function &cmp);

static void insertionsort(::coder::array<int, 1U> &x, int xstart, int xend,
                          const anonymous_function &cmp);

static void introsort(::coder::array<int, 1U> &x, int xend,
                      const anonymous_function &cmp);

namespace scalar {
static void b_sqrt(creal_T &x);

}
static void svd(const double A[9], double U[9], double s[3], double V[9]);

} // namespace internal
static poseGraph *
optimizePoseGraph(poseGraph &b_poseGraph,
                  robotics::core::internal::BlockMatrix &iobj_0,
                  poseGraph &iobj_1);

static double tic(double &tstart_tv_nsec);

static double toc(double tstart_tv_sec, double tstart_tv_nsec);

} // namespace coder
static int div_s32(int numerator, int denominator);

static void minus(::coder::array<double, 1U> &in1,
                  const ::coder::array<double, 1U> &in2,
                  const ::coder::array<double, 1U> &in3);

static double rt_atan2d_snf(double u0, double u1);

static double rt_hypotd_snf(double u0, double u1);

} // namespace poseGraphOptimize

namespace poseGraphOptimize {
namespace coder {
void poseGraph::addRelativePose(const double varargin_1[3], double varargin_3,
                                double varargin_4)
{
  robotics::core::internal::BlockMatrix *obj;
  ::coder::array<double, 2U> T;
  ::coder::array<double, 2U> a;
  ::coder::array<double, 2U> s;
  ::coder::array<double, 1U> id;
  double Omega[9];
  double RR[9];
  double Trel[9];
  double nodePair[2];
  double T_tmp;
  double d;
  double edgeId;
  double fromNodeId;
  double toNodeId;
  int coffset;
  int i;
  int iRowS;
  int inner;
  int j;
  int k;
  boolean_T constraintNeedsInversion;
  boolean_T exitg1;
  boolean_T isLoopClosure;
  boolean_T needNewPoseNode;
  boolean_T tf;
  constraintNeedsInversion = false;
  if ((!(varargin_3 == varargin_4)) &&
      ((!(varargin_3 <= NumNodes)) ||
       (!IsLandmarkNode[static_cast<int>(varargin_3) - 1])) &&
      ((!(varargin_4 <= NumNodes)) ||
       (!IsLandmarkNode[static_cast<int>(varargin_4) - 1]))) {
    if ((varargin_3 <= NumNodes) && (varargin_4 <= NumNodes)) {
      boolean_T guard1{false};
      nodePair[0] = varargin_3;
      nodePair[1] = varargin_4;
      findEdgeID(nodePair, id);
      guard1 = false;
      if (id.size(0) == 0) {
        guard1 = true;
      } else {
        nodePair[0] = varargin_3;
        nodePair[1] = varargin_4;
        tf = false;
        iRowS = 0;
        exitg1 = false;
        while ((!exitg1) && (iRowS <= LoopClosureEdgeNodePairs.size(0) - 1)) {
          boolean_T exitg2;
          needNewPoseNode = true;
          j = 0;
          exitg2 = false;
          while ((!exitg2) && (j < 2)) {
            if (nodePair[j] !=
                LoopClosureEdgeNodePairs[iRowS +
                                         LoopClosureEdgeNodePairs.size(0) *
                                             j]) {
              needNewPoseNode = false;
              exitg2 = true;
            } else {
              j++;
            }
          }
          if (needNewPoseNode) {
            tf = true;
            exitg1 = true;
          } else {
            iRowS++;
          }
        }
        if (tf) {
          guard1 = true;
        } else {
          fromNodeId = std::fmin(varargin_3, varargin_4);
          toNodeId = std::fmax(varargin_3, varargin_4);
          needNewPoseNode = false;
          constraintNeedsInversion = (fromNodeId != varargin_3);
        }
      }
      if (guard1) {
        fromNodeId = varargin_3;
        toNodeId = varargin_4;
        needNewPoseNode = false;
      }
    } else {
      edgeId = std::fmin(varargin_3, varargin_4);
      if ((edgeId <= NumNodes) &&
          (std::fmax(varargin_3, varargin_4) - NumNodes == 1.0)) {
        fromNodeId = edgeId;
        toNodeId = NumNodes + 1.0;
        needNewPoseNode = true;
        constraintNeedsInversion = (edgeId != varargin_3);
      }
    }
  }
  nodePair[0] = fromNodeId;
  nodePair[1] = toNodeId;
  findEdgeID(nodePair, id);
  isLoopClosure = false;
  if (!needNewPoseNode) {
    if (id.size(0) == 0) {
      isLoopClosure = true;
    } else {
      s.set_size(1, LoopClosureEdgeIDsInternal.size(1));
      iRowS = LoopClosureEdgeIDsInternal.size(1);
      for (i = 0; i < iRowS; i++) {
        s[i] = LoopClosureEdgeIDsInternal[i];
      }
      tf = false;
      k = 0;
      exitg1 = false;
      while ((!exitg1) && (k <= s.size(1) - 1)) {
        if (id[0] == s[k]) {
          tf = true;
          exitg1 = true;
        } else {
          k++;
        }
      }
      if (tf) {
        isLoopClosure = true;
      }
    }
  }
  edgeId = std::sin(varargin_1[2]);
  T_tmp = std::cos(varargin_1[2]);
  Trel[0] = T_tmp;
  Trel[3] = -edgeId;
  Trel[6] = varargin_1[0];
  Trel[1] = edgeId;
  Trel[4] = T_tmp;
  Trel[7] = varargin_1[1];
  Trel[2] = 0.0;
  Trel[5] = 0.0;
  Trel[8] = 1.0;
  for (i = 0; i < 9; i++) {
    Omega[i] = iv[i];
  }
  iRowS = static_cast<int>(NumEdges + 1.0);
  EdgeNodePairs[iRowS - 1] = fromNodeId;
  EdgeNodePairs[(iRowS + EdgeNodePairs.size(0)) - 1] = toNodeId;
  if (constraintNeedsInversion) {
    double b_RR[9];
    Trel[0] = T_tmp;
    Trel[1] = -edgeId;
    Trel[6] = T_tmp * -varargin_1[0] + edgeId * -varargin_1[1];
    Trel[3] = edgeId;
    Trel[4] = T_tmp;
    Trel[7] = -edgeId * -varargin_1[0] + T_tmp * -varargin_1[1];
    Trel[2] = 0.0;
    Trel[5] = 0.0;
    Trel[8] = 1.0;
    std::memset(&RR[0], 0, 9U * sizeof(double));
    RR[0] = T_tmp;
    RR[1] = edgeId;
    RR[3] = -edgeId;
    RR[4] = T_tmp;
    RR[8] = 1.0;
    for (i = 0; i < 3; i++) {
      d = RR[i];
      edgeId = RR[i + 3];
      inner = static_cast<int>(RR[i + 6]);
      for (coffset = 0; coffset < 3; coffset++) {
        b_RR[i + 3 * coffset] =
            (d * static_cast<double>(iv[3 * coffset]) +
             edgeId * static_cast<double>(iv[3 * coffset + 1])) +
            static_cast<double>(inner * iv[3 * coffset + 2]);
      }
      d = b_RR[i];
      edgeId = b_RR[i + 3];
      T_tmp = b_RR[i + 6];
      for (inner = 0; inner < 3; inner++) {
        Omega[i + 3 * inner] =
            (d * RR[inner] + edgeId * RR[inner + 3]) + T_tmp * RR[inner + 6];
      }
    }
  }
  EdgeMeasurements->replaceBlock(NumEdges + 1.0, Trel);
  EdgeInfoMatrices->replaceBlock(NumEdges + 1.0, Omega);
  NumEdges++;
  edgeId = NumEdges;
  if (needNewPoseNode) {
    NodeEstimates->extractBlock(fromNodeId, a);
    iRowS = a.size(0);
    inner = a.size(1);
    T.set_size(a.size(0), 3);
    for (j = 0; j < 3; j++) {
      int boffset;
      coffset = j * iRowS;
      boffset = j * 3;
      for (int b_i{0}; b_i < iRowS; b_i++) {
        T[coffset + b_i] = 0.0;
      }
      for (k = 0; k < inner; k++) {
        int aoffset;
        aoffset = k * a.size(0);
        edgeId = Trel[boffset + k];
        for (int b_i{0}; b_i < iRowS; b_i++) {
          i = coffset + b_i;
          T[i] = T[i] + a[aoffset + b_i] * edgeId;
        }
      }
    }
    obj = NodeEstimates;
    edgeId = NumNodes + 1.0;
    edgeId = obj->BlockSize[0] * (edgeId - 1.0) + 1.0;
    d = (edgeId + obj->BlockSize[0]) - 1.0;
    if (edgeId > d) {
      i = 1;
    } else {
      i = static_cast<int>(edgeId);
    }
    iRowS = T.size(0);
    for (inner = 0; inner < 3; inner++) {
      for (coffset = 0; coffset < iRowS; coffset++) {
        obj->Matrix[((i + coffset) + obj->Matrix.size(0) * inner) - 1] =
            T[coffset + T.size(0) * inner];
      }
    }
    NumNodes++;
    NodeDims[static_cast<int>(NumNodes) - 1] = 3.0;
    NodeMap[static_cast<int>(NumNodes) - 1] =
        NodeMap[static_cast<int>(NumNodes - 1.0) - 1] +
        NodeDims[static_cast<int>(NumNodes - 1.0) - 1];
    IsLandmarkNode[static_cast<int>(NumNodes) - 1] = false;
  } else if (isLoopClosure) {
    iRowS = static_cast<int>(NumLoopClosureEdges + 1.0);
    LoopClosureEdgeNodePairs[iRowS - 1] = fromNodeId;
    LoopClosureEdgeNodePairs[(iRowS + LoopClosureEdgeNodePairs.size(0)) - 1] =
        toNodeId;
    LoopClosureEdgeIDsInternal[static_cast<int>(NumLoopClosureEdges + 1.0) -
                               1] = edgeId;
    NumLoopClosureEdges++;
  }
}

void sparse::ctranspose(sparse &y) const
{
  ::coder::array<int, 1U> counts;
  int i;
  int loop_ub;
  int nl;
  int numalloc;
  int outridx;
  nl = n;
  y.m = n;
  y.n = m;
  outridx = colidx[colidx.size(0) - 1];
  if (outridx - 1 >= 1) {
    numalloc = outridx - 2;
  } else {
    numalloc = 0;
  }
  y.d.set_size(numalloc + 1);
  i = m + 1;
  y.colidx.set_size(i);
  y.colidx[0] = 1;
  y.rowidx.set_size(numalloc + 1);
  for (int c{0}; c <= numalloc; c++) {
    y.d[c] = 0.0;
    y.rowidx[c] = 0;
  }
  loop_ub = m;
  for (int c{0}; c < loop_ub; c++) {
    y.colidx[c + 1] = 1;
  }
  y.fillIn();
  if ((m != 0) && (n != 0)) {
    numalloc = y.colidx.size(0);
    for (int c{0}; c < numalloc; c++) {
      y.colidx[c] = 0;
    }
    for (numalloc = 0; numalloc <= outridx - 2; numalloc++) {
      y.colidx[rowidx[numalloc]] = y.colidx[rowidx[numalloc]] + 1;
    }
    y.colidx[0] = 1;
    for (numalloc = 2; numalloc <= i; numalloc++) {
      y.colidx[numalloc - 1] = y.colidx[numalloc - 1] + y.colidx[numalloc - 2];
    }
    counts.set_size(m);
    for (i = 0; i < loop_ub; i++) {
      counts[i] = 0;
    }
    for (int c{0}; c < nl; c++) {
      for (numalloc = colidx[c] - 1; numalloc + 1 < colidx[c + 1]; numalloc++) {
        loop_ub = counts[rowidx[numalloc] - 1];
        outridx = (loop_ub + y.colidx[rowidx[numalloc] - 1]) - 1;
        y.d[outridx] = d[numalloc];
        y.rowidx[outridx] = c + 1;
        counts[rowidx[numalloc] - 1] = loop_ub + 1;
      }
    }
  }
}

void sparse::fillIn()
{
  int i;
  int idx;
  idx = 1;
  i = colidx.size(0);
  for (int c{0}; c <= i - 2; c++) {
    int ridx;
    ridx = colidx[c];
    colidx[c] = idx;
    int exitg1;
    int i1;
    do {
      exitg1 = 0;
      i1 = colidx[c + 1];
      if (ridx < i1) {
        double val;
        int currRowIdx;
        val = 0.0;
        currRowIdx = rowidx[ridx - 1];
        while ((ridx < i1) && (rowidx[ridx - 1] == currRowIdx)) {
          val += d[ridx - 1];
          ridx++;
        }
        if (val != 0.0) {
          d[idx - 1] = val;
          rowidx[idx - 1] = currRowIdx;
          idx++;
        }
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  colidx[colidx.size(0) - 1] = idx;
}

void poseGraph::findEdgeID(const double nodePair[2],
                           ::coder::array<double, 1U> &edgeID) const
{
  ::coder::array<int, 1U> ii;
  ::coder::array<boolean_T, 1U> tf;
  int iRowA;
  int idx;
  int nx;
  boolean_T exitg1;
  tf.set_size(EdgeNodePairs.size(0));
  idx = EdgeNodePairs.size(0);
  for (iRowA = 0; iRowA < idx; iRowA++) {
    boolean_T p;
    tf[iRowA] = false;
    p = true;
    nx = 0;
    exitg1 = false;
    while ((!exitg1) && (nx < 2)) {
      if (EdgeNodePairs[iRowA + EdgeNodePairs.size(0) * nx] != nodePair[nx]) {
        p = false;
        exitg1 = true;
      } else {
        nx++;
      }
    }
    if (p) {
      tf[iRowA] = true;
    }
  }
  nx = tf.size(0);
  idx = 0;
  ii.set_size(tf.size(0));
  iRowA = 0;
  exitg1 = false;
  while ((!exitg1) && (iRowA <= nx - 1)) {
    if (tf[iRowA]) {
      idx++;
      ii[idx - 1] = iRowA + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        iRowA++;
      }
    } else {
      iRowA++;
    }
  }
  if (tf.size(0) == 1) {
    if (idx == 0) {
      ii.set_size(0);
    }
  } else {
    if (idx < 1) {
      idx = 0;
    }
    ii.set_size(idx);
  }
  edgeID.set_size(ii.size(0));
  idx = ii.size(0);
  for (nx = 0; nx < idx; nx++) {
    edgeID[nx] = ii[nx];
  }
}

namespace internal {
void CXSparseAPI::iteratedQR(const sparse &A,
                             const ::coder::array<double, 1U> &b, int n,
                             ::coder::array<double, 1U> &out)
{
  cs_di *cxA;
  cs_din *N;
  cs_dis *S;
  sparse in;
  ::coder::array<double, 1U> outBuff;
  double tol;
  int loop_ub;
  if (A.m < A.n) {
    A.ctranspose(in);
    cxA = makeCXSparseMatrix(in.colidx[in.colidx.size(0) - 1] - 1, in.n, in.m,
                             &(in.colidx.data())[0], &(in.rowidx.data())[0],
                             &(in.d.data())[0]);
  } else {
    cxA =
        makeCXSparseMatrix(A.colidx[A.colidx.size(0) - 1] - 1, A.n, A.m,
                           &(((::coder::array<int, 1U> *)&A.colidx)->data())[0],
                           &(((::coder::array<int, 1U> *)&A.rowidx)->data())[0],
                           &(((::coder::array<double, 1U> *)&A.d)->data())[0]);
  }
  S = cs_di_sqr(2, cxA, 1);
  N = cs_di_qr(cxA, S);
  cs_di_spfree(cxA);
  qr_rank_di(N, &tol);
  out.set_size(n);
  if (b.size(0) < n) {
    outBuff.set_size(n);
  } else {
    outBuff.set_size(b.size(0));
  }
  loop_ub = b.size(0);
  for (int i{0}; i < loop_ub; i++) {
    outBuff[i] = b[i];
  }
  solve_from_qr_di(N, S, (double *)&(outBuff.data())[0], b.size(0), n);
  if (n < 1) {
    loop_ub = 0;
  } else {
    loop_ub = n;
  }
  for (int i{0}; i < loop_ub; i++) {
    out[i] = outBuff[i];
  }
  cs_di_sfree(S);
  cs_di_nfree(N);
}

} // namespace internal
namespace nav {
namespace algs {
namespace internal {
double PoseGraphHelpers::costBetweenTwoNodes(
    const ::coder::array<double, 2U> &Toi,
    const ::coder::array<double, 2U> &Toj,
    const ::coder::array<double, 2U> &measurement,
    const ::coder::array<double, 2U> &Omega, boolean_T nodejIsLandmark,
    double gradi_data[], int &gradi_size, double gradj_data[], int &gradj_size,
    double hessii_data[], int hessii_size[2], double hessij_data[],
    int hessij_size[2], double hessji_data[], int hessji_size[2],
    double hessjj_data[], int hessjj_size[2])
{
  static const signed char b_N[36]{1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                   0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1};
  double Jaci_data[49];
  double Jacj_data[49];
  double Jaci[36];
  double Jacj[36];
  double b_y_tmp_data[21];
  double y_tmp_data[21];
  double dv[18];
  double f_R[9];
  double e_data[7];
  double e_R[3];
  double cost;
  int Jaci_size[2];
  int Jacj_size[2];
  int b_y_tmp_size[2];
  int y_tmp_size[2];
  int aoffset;
  int b_i;
  int boffset;
  int coffset;
  int i;
  int j;
  if (Omega.size(0) == 2) {
    double R[9];
    double N[6];
    double cosTheta;
    double sinTheta;
    for (i = 0; i < 3; i++) {
      R[3 * i] = measurement[measurement.size(0) * i];
      R[3 * i + 1] = measurement[measurement.size(0) * i + 1];
      R[3 * i + 2] = measurement[measurement.size(0) * i + 2];
    }
    double b_measurement;
    double d;
    double q;
    double s;
    double sinThetaDy;
    d = rt_atan2d_snf(Toi[1], Toi[0]);
    sinTheta = std::sin(d);
    cosTheta = std::cos(d);
    q = Toj[Toj.size(0) * 2 + 1] - Toi[Toi.size(0) * 2 + 1];
    sinThetaDy = q * sinTheta;
    b_measurement = Toj[Toj.size(0) * 2] - Toi[Toi.size(0) * 2];
    s = b_measurement * cosTheta;
    b_i = 2;
    e_data[0] = (sinThetaDy + s) - R[6];
    q = q * cosTheta - b_measurement * sinTheta;
    e_data[1] = q - R[7];
    N[0] = -cosTheta;
    N[2] = -sinTheta;
    N[4] = q;
    N[1] = sinTheta;
    N[3] = -cosTheta;
    N[5] = -sinThetaDy - s;
    Jaci_size[0] = 2;
    Jaci_size[1] = 3;
    for (i = 0; i < 6; i++) {
      Jaci_data[i] = N[i];
    }
    Jacj_size[0] = 2;
    Jacj_size[1] = 2;
    Jacj_data[0] = cosTheta;
    Jacj_data[1] = -sinTheta;
    Jacj_data[2] = sinTheta;
    Jacj_data[3] = cosTheta;
  } else if (Omega.size(0) == 3) {
    if (nodejIsLandmark) {
      double Sio[16];
      double Tio[16];
      double deltatform[16];
      double R[9];
      double b_R[9];
      double t[3];
      double d;
      double d1;
      double d2;
      for (i = 0; i < 4; i++) {
        coffset = i << 2;
        Sio[coffset] = Toi[Toi.size(0) * i];
        deltatform[coffset] = measurement[measurement.size(0) * i];
        Sio[coffset + 1] = Toi[Toi.size(0) * i + 1];
        deltatform[coffset + 1] = measurement[measurement.size(0) * i + 1];
        Sio[coffset + 2] = Toi[Toi.size(0) * i + 2];
        deltatform[coffset + 2] = measurement[measurement.size(0) * i + 2];
        Sio[coffset + 3] = Toi[Toi.size(0) * i + 3];
        deltatform[coffset + 3] = measurement[measurement.size(0) * i + 3];
      }
      for (i = 0; i < 3; i++) {
        R[3 * i] = Sio[i];
        R[3 * i + 1] = Sio[i + 4];
        R[3 * i + 2] = Sio[i + 8];
      }
      for (i = 0; i < 9; i++) {
        b_R[i] = -R[i];
      }
      d = Sio[12];
      d1 = Sio[13];
      d2 = Sio[14];
      for (i = 0; i < 3; i++) {
        coffset = i << 2;
        Tio[coffset] = R[3 * i];
        Tio[coffset + 1] = R[3 * i + 1];
        Tio[coffset + 2] = R[3 * i + 2];
        Tio[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
      }
      Tio[3] = 0.0;
      Tio[7] = 0.0;
      Tio[11] = 0.0;
      Tio[15] = 1.0;
      for (i = 0; i < 3; i++) {
        f_R[3 * i] = -deltatform[i];
        f_R[3 * i + 1] = -deltatform[i + 4];
        f_R[3 * i + 2] = -deltatform[i + 8];
        t[i] = ((Tio[i] * Toj[Toj.size(0) * 3] +
                 Tio[i + 4] * Toj[Toj.size(0) * 3 + 1]) +
                Tio[i + 8] * Toj[Toj.size(0) * 3 + 2]) +
               Tio[i + 12];
      }
      b_i = 3;
      d = deltatform[12];
      d1 = deltatform[13];
      d2 = deltatform[14];
      for (i = 0; i < 3; i++) {
        dv[3 * i] = iv[3 * i];
        j = 3 * i + 1;
        dv[j] = iv[j];
        j = 3 * i + 2;
        dv[j] = iv[j];
        e_data[i] = t[i] + ((f_R[i] * d + f_R[i + 3] * d1) + f_R[i + 6] * d2);
      }
      dv[9] = -0.0;
      dv[12] = t[2];
      dv[15] = -t[1];
      dv[10] = -t[2];
      dv[13] = -0.0;
      dv[16] = t[0];
      dv[11] = t[1];
      dv[14] = -t[0];
      dv[17] = -0.0;
      Jaci_size[0] = 3;
      Jaci_size[1] = 6;
      std::copy(&dv[0], &dv[18], &Jaci_data[0]);
      Jacj_size[0] = 3;
      Jacj_size[1] = 3;
      for (i = 0; i < 3; i++) {
        coffset = i << 2;
        Jacj_data[3 * i] = Tio[coffset];
        Jacj_data[3 * i + 1] = Tio[coffset + 1];
        Jacj_data[3 * i + 2] = Tio[coffset + 2];
      }
    } else {
      double R[9];
      double b_R[9];
      double c_R[9];
      double dRoidtheta[4];
      double d_R[4];
      double y[4];
      double y_tmp[4];
      double b_y_tmp[2];
      double e_tmp[2];
      double b_measurement;
      double cosTheta;
      double d;
      double d1;
      double d2;
      double d3;
      double d4;
      double q;
      double sinTheta;
      for (i = 0; i < 3; i++) {
        R[3 * i] = Toi[Toi.size(0) * i];
        b_R[3 * i] = Toj[Toj.size(0) * i];
        c_R[3 * i] = measurement[measurement.size(0) * i];
        boffset = 3 * i + 1;
        R[boffset] = Toi[Toi.size(0) * i + 1];
        b_R[boffset] = Toj[Toj.size(0) * i + 1];
        c_R[boffset] = measurement[measurement.size(0) * i + 1];
        boffset = 3 * i + 2;
        R[boffset] = Toi[Toi.size(0) * i + 2];
        b_R[boffset] = Toj[Toj.size(0) * i + 2];
        c_R[boffset] = measurement[measurement.size(0) * i + 2];
      }
      for (i = 0; i < 2; i++) {
        d = c_R[3 * i];
        d1 = c_R[3 * i + 1];
        d2 = R[0] * d + R[3] * d1;
        d = R[1] * d + R[4] * d1;
        y[i] = d2 * b_R[0] + d * b_R[1];
        y[i + 2] = d2 * b_R[3] + d * b_R[4];
        e_tmp[i] = b_R[i + 6] - R[i + 6];
      }
      b_measurement = rt_atan2d_snf(y[1], y[0]);
      d_R[0] = R[0];
      d_R[1] = R[3];
      d_R[2] = R[1];
      d_R[3] = R[4];
      if (std::isnan(b_measurement + 3.1415926535897931) ||
          std::isinf(b_measurement + 3.1415926535897931)) {
        cosTheta = rtNaN;
      } else if (b_measurement + 3.1415926535897931 == 0.0) {
        cosTheta = 0.0;
      } else {
        boolean_T rEQ0;
        cosTheta =
            std::fmod(b_measurement + 3.1415926535897931, 6.2831853071795862);
        rEQ0 = (cosTheta == 0.0);
        if (!rEQ0) {
          q = std::abs((b_measurement + 3.1415926535897931) /
                       6.2831853071795862);
          rEQ0 =
              !(std::abs(q - std::floor(q + 0.5)) > 2.2204460492503131E-16 * q);
        }
        if (rEQ0) {
          cosTheta = 0.0;
        } else if (b_measurement + 3.1415926535897931 < 0.0) {
          cosTheta += 6.2831853071795862;
        }
      }
      dRoidtheta[0] = R[3];
      dRoidtheta[2] = -R[0];
      dRoidtheta[1] = R[0];
      dRoidtheta[3] = R[3];
      d = c_R[0];
      d1 = c_R[1];
      d2 = c_R[3];
      d3 = c_R[4];
      for (j = 0; j < 2; j++) {
        coffset = j << 1;
        b_measurement = dRoidtheta[j + 2];
        d4 = dRoidtheta[j];
        y[coffset] = d * d4 + d1 * b_measurement;
        y[coffset + 1] = d2 * d4 + d3 * b_measurement;
      }
      d = e_tmp[0];
      d1 = e_tmp[1];
      for (j = 0; j < 2; j++) {
        coffset = j << 1;
        b_measurement = R[j % 2];
        q = R[(j + 2) % 2 + 3];
        dRoidtheta[coffset] = c_R[0] * b_measurement + c_R[1] * q;
        dRoidtheta[coffset + 1] = c_R[3] * b_measurement + c_R[4] * q;
        b_y_tmp[j] = (d_R[j] * d + d_R[j + 2] * d1) - c_R[j + 6];
      }
      b_i = 3;
      e_data[0] = c_R[0] * b_y_tmp[0] + b_y_tmp[1] * c_R[1];
      e_data[1] = b_y_tmp[0] * c_R[3] + b_y_tmp[1] * c_R[4];
      e_data[2] = cosTheta - 3.1415926535897931;
      y_tmp[0] = -c_R[0];
      y_tmp[1] = -c_R[3];
      y_tmp[2] = -c_R[1];
      y_tmp[3] = -c_R[4];
      d = R[0];
      d1 = R[3];
      d2 = R[1];
      d3 = R[4];
      d4 = e_tmp[0];
      sinTheta = e_tmp[1];
      for (i = 0; i < 2; i++) {
        q = y_tmp[i + 2];
        b_measurement = y_tmp[i];
        d_R[i] = b_measurement * d + q * d1;
        d_R[i + 2] = b_measurement * d2 + q * d3;
        b_y_tmp[i] = y[i] * d4 + y[i + 2] * sinTheta;
      }
      f_R[0] = d_R[0];
      f_R[1] = d_R[1];
      f_R[6] = b_y_tmp[0];
      f_R[3] = d_R[2];
      f_R[4] = d_R[3];
      f_R[7] = b_y_tmp[1];
      f_R[2] = 0.0;
      f_R[5] = 0.0;
      f_R[8] = -1.0;
      Jaci_size[0] = 3;
      Jaci_size[1] = 3;
      std::copy(&f_R[0], &f_R[9], &Jaci_data[0]);
      f_R[0] = dRoidtheta[0];
      f_R[1] = dRoidtheta[1];
      f_R[6] = 0.0;
      f_R[3] = dRoidtheta[2];
      f_R[4] = dRoidtheta[3];
      f_R[7] = 0.0;
      f_R[2] = 0.0;
      f_R[5] = 0.0;
      f_R[8] = 1.0;
      Jacj_size[0] = 3;
      Jacj_size[1] = 3;
      std::copy(&f_R[0], &f_R[9], &Jacj_data[0]);
    }
  } else if (Omega.size(0) == 6) {
    double Sio[16];
    double Sji[16];
    double Tio[16];
    double Tji[16];
    double Tjo[16];
    double Tjop[16];
    double deltatform[16];
    double R[9];
    double b_R[9];
    double dv1[6];
    double d;
    double d1;
    double d2;
    double d3;
    for (i = 0; i < 4; i++) {
      coffset = i << 2;
      Sio[coffset] = Toi[Toi.size(0) * i];
      Sji[coffset] = Toj[Toj.size(0) * i];
      deltatform[coffset] = measurement[measurement.size(0) * i];
      Sio[coffset + 1] = Toi[Toi.size(0) * i + 1];
      Sji[coffset + 1] = Toj[Toj.size(0) * i + 1];
      deltatform[coffset + 1] = measurement[measurement.size(0) * i + 1];
      Sio[coffset + 2] = Toi[Toi.size(0) * i + 2];
      Sji[coffset + 2] = Toj[Toj.size(0) * i + 2];
      deltatform[coffset + 2] = measurement[measurement.size(0) * i + 2];
      Sio[coffset + 3] = Toi[Toi.size(0) * i + 3];
      Sji[coffset + 3] = Toj[Toj.size(0) * i + 3];
      deltatform[coffset + 3] = measurement[measurement.size(0) * i + 3];
    }
    for (i = 0; i < 3; i++) {
      R[3 * i] = Sio[i];
      R[3 * i + 1] = Sio[i + 4];
      R[3 * i + 2] = Sio[i + 8];
    }
    for (i = 0; i < 9; i++) {
      b_R[i] = -R[i];
    }
    d = Sio[12];
    d1 = Sio[13];
    d2 = Sio[14];
    for (i = 0; i < 3; i++) {
      coffset = i << 2;
      Tio[coffset] = R[3 * i];
      aoffset = 3 * i + 1;
      Tio[coffset + 1] = R[aoffset];
      boffset = 3 * i + 2;
      Tio[coffset + 2] = R[boffset];
      Tio[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
      R[3 * i] = deltatform[i];
      R[aoffset] = deltatform[i + 4];
      R[boffset] = deltatform[i + 8];
    }
    Tio[3] = 0.0;
    Tio[7] = 0.0;
    Tio[11] = 0.0;
    Tio[15] = 1.0;
    for (i = 0; i < 9; i++) {
      b_R[i] = -R[i];
    }
    d = deltatform[12];
    d1 = deltatform[13];
    d2 = deltatform[14];
    for (i = 0; i < 3; i++) {
      coffset = i << 2;
      Tji[coffset] = R[3 * i];
      Tji[coffset + 1] = R[3 * i + 1];
      Tji[coffset + 2] = R[3 * i + 2];
      Tji[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
    }
    Tji[3] = 0.0;
    Tji[7] = 0.0;
    Tji[11] = 0.0;
    Tji[15] = 1.0;
    for (i = 0; i < 4; i++) {
      d = Tji[i];
      d1 = Tji[i + 4];
      d2 = Tji[i + 8];
      d3 = Tji[i + 12];
      for (j = 0; j < 4; j++) {
        coffset = j << 2;
        deltatform[i + coffset] = ((d * Tio[coffset] + d1 * Tio[coffset + 1]) +
                                   d2 * Tio[coffset + 2]) +
                                  d3 * Tio[coffset + 3];
      }
    }
    for (i = 0; i < 3; i++) {
      R[3 * i] = Sji[i];
      R[3 * i + 1] = Sji[i + 4];
      R[3 * i + 2] = Sji[i + 8];
    }
    for (i = 0; i < 9; i++) {
      b_R[i] = -R[i];
    }
    d = Sji[12];
    d1 = Sji[13];
    d2 = Sji[14];
    for (i = 0; i < 3; i++) {
      coffset = i << 2;
      Tjo[coffset] = R[3 * i];
      Tjo[coffset + 1] = R[3 * i + 1];
      Tjo[coffset + 2] = R[3 * i + 2];
      Tjo[i + 12] = (b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2;
    }
    Tjo[3] = 0.0;
    Tjo[7] = 0.0;
    Tjo[11] = 0.0;
    Tjo[15] = 1.0;
    for (b_i = 0; b_i < 6; b_i++) {
      double Tm[16];
      double g_R[16];
      double c_R[9];
      double N[6];
      for (i = 0; i < 6; i++) {
        N[i] = static_cast<double>(b_N[i + 6 * b_i]) * 1.0E-5;
      }
      robotics::core::internal::SEHelpers::expSE3hat(N, Sio);
      for (i = 0; i < 6; i++) {
        N[i] = -static_cast<double>(b_N[i + 6 * b_i]) * 1.0E-5;
      }
      robotics::core::internal::SEHelpers::expSE3hat(N, Tm);
      for (i = 0; i < 4; i++) {
        d = Sio[i];
        d1 = Sio[i + 4];
        d2 = Sio[i + 8];
        d3 = Sio[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * Tio[coffset] + d1 * Tio[coffset + 1]) +
                               d2 * Tio[coffset + 2]) +
                              d3 * Tio[coffset + 3];
        }
      }
      for (i = 0; i < 4; i++) {
        d = Tji[i];
        d1 = Tji[i + 4];
        d2 = Tji[i + 8];
        d3 = Tji[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          g_R[i + coffset] = ((d * Tjop[coffset] + d1 * Tjop[coffset + 1]) +
                              d2 * Tjop[coffset + 2]) +
                             d3 * Tjop[coffset + 3];
        }
      }
      for (i = 0; i < 4; i++) {
        d = g_R[i];
        d1 = g_R[i + 4];
        d2 = g_R[i + 8];
        d3 = g_R[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * Sji[coffset] + d1 * Sji[coffset + 1]) +
                               d2 * Sji[coffset + 2]) +
                              d3 * Sji[coffset + 3];
        }
      }
      robotics::core::internal::SEHelpers::veelogmSE3(Tjop, dv1);
      for (i = 0; i < 4; i++) {
        d = Tm[i];
        d1 = Tm[i + 4];
        d2 = Tm[i + 8];
        d3 = Tm[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * Tio[coffset] + d1 * Tio[coffset + 1]) +
                               d2 * Tio[coffset + 2]) +
                              d3 * Tio[coffset + 3];
        }
      }
      for (i = 0; i < 4; i++) {
        d = Tji[i];
        d1 = Tji[i + 4];
        d2 = Tji[i + 8];
        d3 = Tji[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          g_R[i + coffset] = ((d * Tjop[coffset] + d1 * Tjop[coffset + 1]) +
                              d2 * Tjop[coffset + 2]) +
                             d3 * Tjop[coffset + 3];
        }
      }
      for (i = 0; i < 4; i++) {
        d = g_R[i];
        d1 = g_R[i + 4];
        d2 = g_R[i + 8];
        d3 = g_R[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * Sji[coffset] + d1 * Sji[coffset + 1]) +
                               d2 * Sji[coffset + 2]) +
                              d3 * Sji[coffset + 3];
        }
      }
      robotics::core::internal::SEHelpers::veelogmSE3(Tjop, N);
      for (i = 0; i < 6; i++) {
        Jaci[i + 6 * b_i] = (dv1[i] - N[i]) / 2.0E-5;
      }
      for (i = 0; i < 4; i++) {
        d = Sio[i];
        d1 = Sio[i + 4];
        d2 = Sio[i + 8];
        d3 = Sio[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * Tjo[coffset] + d1 * Tjo[coffset + 1]) +
                               d2 * Tjo[coffset + 2]) +
                              d3 * Tjo[coffset + 3];
        }
        d = Tm[i];
        d1 = Tm[i + 4];
        d2 = Tm[i + 8];
        d3 = Tm[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Sio[i + coffset] = ((d * Tjo[coffset] + d1 * Tjo[coffset + 1]) +
                              d2 * Tjo[coffset + 2]) +
                             d3 * Tjo[coffset + 3];
        }
      }
      for (i = 0; i < 3; i++) {
        R[3 * i] = Tjop[i];
        b_R[3 * i] = Sio[i];
        boffset = 3 * i + 1;
        R[boffset] = Tjop[i + 4];
        b_R[boffset] = Sio[i + 4];
        boffset = 3 * i + 2;
        R[boffset] = Tjop[i + 8];
        b_R[boffset] = Sio[i + 8];
      }
      for (i = 0; i < 9; i++) {
        c_R[i] = -R[i];
      }
      d = Tjop[12];
      d1 = Tjop[13];
      d2 = Tjop[14];
      for (i = 0; i < 3; i++) {
        boffset = i << 2;
        g_R[boffset] = R[3 * i];
        g_R[boffset + 1] = R[3 * i + 1];
        g_R[boffset + 2] = R[3 * i + 2];
        g_R[i + 12] = (c_R[i] * d + c_R[i + 3] * d1) + c_R[i + 6] * d2;
      }
      g_R[3] = 0.0;
      g_R[7] = 0.0;
      g_R[11] = 0.0;
      g_R[15] = 1.0;
      for (i = 0; i < 4; i++) {
        d = deltatform[i];
        d1 = deltatform[i + 4];
        d2 = deltatform[i + 8];
        d3 = deltatform[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * g_R[coffset] + d1 * g_R[coffset + 1]) +
                               d2 * g_R[coffset + 2]) +
                              d3 * g_R[coffset + 3];
        }
      }
      robotics::core::internal::SEHelpers::veelogmSE3(Tjop, dv1);
      for (i = 0; i < 9; i++) {
        R[i] = -b_R[i];
      }
      d = Sio[12];
      d1 = Sio[13];
      d2 = Sio[14];
      for (i = 0; i < 3; i++) {
        boffset = i << 2;
        g_R[boffset] = b_R[3 * i];
        g_R[boffset + 1] = b_R[3 * i + 1];
        g_R[boffset + 2] = b_R[3 * i + 2];
        g_R[i + 12] = (R[i] * d + R[i + 3] * d1) + R[i + 6] * d2;
      }
      g_R[3] = 0.0;
      g_R[7] = 0.0;
      g_R[11] = 0.0;
      g_R[15] = 1.0;
      for (i = 0; i < 4; i++) {
        d = deltatform[i];
        d1 = deltatform[i + 4];
        d2 = deltatform[i + 8];
        d3 = deltatform[i + 12];
        for (j = 0; j < 4; j++) {
          coffset = j << 2;
          Tjop[i + coffset] = ((d * g_R[coffset] + d1 * g_R[coffset + 1]) +
                               d2 * g_R[coffset + 2]) +
                              d3 * g_R[coffset + 3];
        }
      }
      robotics::core::internal::SEHelpers::veelogmSE3(Tjop, N);
      for (i = 0; i < 6; i++) {
        Jacj[i + 6 * b_i] = (dv1[i] - N[i]) / 2.0E-5;
      }
    }
    for (i = 0; i < 4; i++) {
      d = deltatform[i];
      d1 = deltatform[i + 4];
      d2 = deltatform[i + 8];
      d3 = deltatform[i + 12];
      for (j = 0; j < 4; j++) {
        coffset = j << 2;
        Tjop[i + coffset] = ((d * Sji[coffset] + d1 * Sji[coffset + 1]) +
                             d2 * Sji[coffset + 2]) +
                            d3 * Sji[coffset + 3];
      }
    }
    robotics::core::internal::SEHelpers::veelogmSE3(Tjop, dv1);
    b_i = 6;
    for (i = 0; i < 6; i++) {
      e_data[i] = dv1[i];
    }
    Jaci_size[0] = 6;
    Jaci_size[1] = 6;
    Jacj_size[0] = 6;
    Jacj_size[1] = 6;
    std::copy(&Jaci[0], &Jaci[36], &Jaci_data[0]);
    std::copy(&Jacj[0], &Jacj[36], &Jacj_data[0]);
  } else {
    double Sio[16];
    double Sji[16];
    double Tji[16];
    double Tjop[16];
    double g_R[16];
    double R[9];
    double b_R[9];
    double deltavec[7];
    double b_measurement;
    double cosTheta;
    double d;
    double d1;
    double d2;
    double q;
    double sinThetaDy;
    for (i = 0; i < 4; i++) {
      coffset = i << 2;
      Tjop[coffset] = Toi[Toi.size(0) * i];
      Tji[coffset] = Toj[Toj.size(0) * i];
      Sio[coffset] = measurement[measurement.size(0) * i];
      Tjop[coffset + 1] = Toi[Toi.size(0) * i + 1];
      Tji[coffset + 1] = Toj[Toj.size(0) * i + 1];
      Sio[coffset + 1] = measurement[measurement.size(0) * i + 1];
      Tjop[coffset + 2] = Toi[Toi.size(0) * i + 2];
      Tji[coffset + 2] = Toj[Toj.size(0) * i + 2];
      Sio[coffset + 2] = measurement[measurement.size(0) * i + 2];
      Tjop[coffset + 3] = Toi[Toi.size(0) * i + 3];
      Tji[coffset + 3] = Toj[Toj.size(0) * i + 3];
      Sio[coffset + 3] = measurement[measurement.size(0) * i + 3];
    }
    for (i = 0; i < 3; i++) {
      R[3 * i] = Tjop[i];
      b_R[3 * i] = Sio[i];
      boffset = 3 * i + 1;
      R[boffset] = Tjop[i + 4];
      b_R[boffset] = Sio[i + 4];
      boffset = 3 * i + 2;
      R[boffset] = Tjop[i + 8];
      b_R[boffset] = Sio[i + 8];
    }
    b_measurement = measurement[measurement.size(0) * 3 + 3];
    d = Sio[12];
    d1 = Sio[13];
    d2 = Sio[14];
    for (i = 0; i < 3; i++) {
      coffset = i << 2;
      Sji[coffset] = b_R[3 * i];
      Sji[coffset + 1] = b_R[3 * i + 1];
      Sji[coffset + 2] = b_R[3 * i + 2];
      Sji[i + 12] =
          -((b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2) / b_measurement;
    }
    Sji[3] = 0.0;
    Sji[7] = 0.0;
    Sji[11] = 0.0;
    Sji[15] = 1.0 / measurement[measurement.size(0) * 3 + 3];
    for (i = 0; i < 3; i++) {
      b_R[3 * i] = Tjop[i];
      b_R[3 * i + 1] = Tjop[i + 4];
      b_R[3 * i + 2] = Tjop[i + 8];
    }
    q = Toi[Toi.size(0) * 3 + 3];
    d = Tjop[12];
    d1 = Tjop[13];
    d2 = Tjop[14];
    for (i = 0; i < 3; i++) {
      coffset = i << 2;
      Sio[coffset] = b_R[3 * i];
      Sio[coffset + 1] = b_R[3 * i + 1];
      Sio[coffset + 2] = b_R[3 * i + 2];
      Sio[i + 12] = -((b_R[i] * d + b_R[i + 3] * d1) + b_R[i + 6] * d2) / q;
    }
    Sio[3] = 0.0;
    Sio[7] = 0.0;
    Sio[11] = 0.0;
    Sio[15] = 1.0 / Toi[Toi.size(0) * 3 + 3];
    q = Toi[Toi.size(0) * 3 + 3];
    b_measurement = Toj[Toj.size(0) * 3 + 3];
    cosTheta = Toi[Toi.size(0) * 3 + 3];
    sinThetaDy = Toj[Toj.size(0) * 3 + 3];
    d = Toi[Toi.size(0) * 3 + 3];
    d1 = Toj[Toj.size(0) * 3 + 3];
    for (int k{0}; k < 7; k++) {
      double Tjo[16];
      double Tm[16];
      double deltatform[16];
      double c_R[9];
      double t[3];
      double d3;
      double d4;
      double s;
      double sinTheta;
      for (b_i = 0; b_i < 7; b_i++) {
        deltavec[b_i] = 0.0;
      }
      deltavec[k] = 1.0E-9;
      robotics::core::internal::Sim3Helpers::sim3ToSform(deltavec, Tjo);
      for (i = 0; i < 3; i++) {
        b_R[3 * i] = Tjo[i];
        b_R[3 * i + 1] = Tjo[i + 4];
        b_R[3 * i + 2] = Tjo[i + 8];
      }
      d2 = Tjo[12];
      d3 = Tjo[13];
      d4 = Tjo[14];
      sinTheta = Tjo[15];
      for (i = 0; i < 3; i++) {
        coffset = i << 2;
        deltatform[coffset] = b_R[3 * i];
        deltatform[coffset + 1] = b_R[3 * i + 1];
        deltatform[coffset + 2] = b_R[3 * i + 2];
        deltatform[i + 12] =
            -((b_R[i] * d2 + b_R[i + 3] * d3) + b_R[i + 6] * d4) / sinTheta;
      }
      deltatform[3] = 0.0;
      deltatform[7] = 0.0;
      deltatform[11] = 0.0;
      deltatform[15] = 1.0 / Tjo[15];
      for (i = 0; i < 3; i++) {
        d2 = 0.0;
        for (j = 0; j < 3; j++) {
          coffset = j << 2;
          f_R[i + 3 * j] = (Tjop[i] * deltatform[coffset] +
                            Tjop[i + 4] * deltatform[coffset + 1]) +
                           Tjop[i + 8] * deltatform[coffset + 2];
          d2 += q * Tjop[i + coffset] * deltatform[j + 12];
        }
        e_R[i] = d2 + Tjop[i + 12];
      }
      for (i = 0; i < 3; i++) {
        coffset = i << 2;
        Tjo[coffset] = f_R[3 * i];
        Tjo[coffset + 1] = f_R[3 * i + 1];
        Tjo[coffset + 2] = f_R[3 * i + 2];
        Tjo[i + 12] = e_R[i];
      }
      Tjo[3] = 0.0;
      Tjo[7] = 0.0;
      Tjo[11] = 0.0;
      Tjo[15] = d * deltatform[15];
      s = d1 * deltatform[15];
      deltavec[k] = -1.0E-9;
      robotics::core::internal::Sim3Helpers::sim3ToSform(deltavec, Tm);
      for (i = 0; i < 3; i++) {
        d2 = 0.0;
        for (j = 0; j < 3; j++) {
          boffset = j << 2;
          coffset = i + boffset;
          aoffset = j + 3 * i;
          b_R[aoffset] = Tjo[coffset];
          c_R[i + 3 * j] = (Tji[i] * deltatform[boffset] +
                            Tji[i + 4] * deltatform[boffset + 1]) +
                           Tji[i + 8] * deltatform[boffset + 2];
          d2 += b_measurement * Tji[coffset] * deltatform[j + 12];
          f_R[aoffset] = Tm[coffset];
        }
        t[i] = d2 + Tji[i + 12];
      }
      d2 = Tm[12];
      d3 = Tm[13];
      d4 = Tm[14];
      sinTheta = Tm[15];
      for (i = 0; i < 3; i++) {
        coffset = i << 2;
        deltatform[coffset] = f_R[3 * i];
        deltatform[coffset + 1] = f_R[3 * i + 1];
        deltatform[coffset + 2] = f_R[3 * i + 2];
        deltatform[i + 12] =
            -((f_R[i] * d2 + f_R[i + 3] * d3) + f_R[i + 6] * d4) / sinTheta;
      }
      deltatform[3] = 0.0;
      deltatform[7] = 0.0;
      deltatform[11] = 0.0;
      deltatform[15] = 1.0 / Tm[15];
      for (i = 0; i < 3; i++) {
        d2 = 0.0;
        for (j = 0; j < 3; j++) {
          coffset = j << 2;
          f_R[i + 3 * j] = (Tjop[i] * deltatform[coffset] +
                            Tjop[i + 4] * deltatform[coffset + 1]) +
                           Tjop[i + 8] * deltatform[coffset + 2];
          d2 += cosTheta * Tjop[i + coffset] * deltatform[j + 12];
        }
        e_R[i] = d2 + Tjop[i + 12];
      }
      for (i = 0; i < 3; i++) {
        coffset = i << 2;
        Tm[coffset] = f_R[3 * i];
        Tm[coffset + 1] = f_R[3 * i + 1];
        Tm[coffset + 2] = f_R[3 * i + 2];
        Tm[i + 12] = e_R[i];
      }
      Tm[3] = 0.0;
      Tm[7] = 0.0;
      Tm[11] = 0.0;
      Tm[15] = d * deltatform[15];
      d2 = Tjo[12];
      d3 = Tjo[13];
      d4 = Tjo[14];
      sinTheta = Tjo[15];
      for (i = 0; i < 3; i++) {
        f_R[3 * i] = Tm[i];
        boffset = i << 2;
        g_R[boffset] = b_R[3 * i];
        coffset = 3 * i + 1;
        f_R[coffset] = Tm[i + 4];
        g_R[boffset + 1] = b_R[coffset];
        coffset = 3 * i + 2;
        f_R[coffset] = Tm[i + 8];
        g_R[boffset + 2] = b_R[coffset];
        g_R[i + 12] =
            -((b_R[i] * d2 + b_R[i + 3] * d3) + b_R[i + 6] * d4) / sinTheta;
      }
      g_R[3] = 0.0;
      g_R[7] = 0.0;
      g_R[11] = 0.0;
      g_R[15] = 1.0 / Tjo[15];
      robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, g_R, Tji,
                                                             deltavec);
      d2 = Tm[12];
      d3 = Tm[13];
      d4 = Tm[14];
      sinTheta = Tm[15];
      for (i = 0; i < 3; i++) {
        boffset = i << 2;
        g_R[boffset] = f_R[3 * i];
        g_R[boffset + 1] = f_R[3 * i + 1];
        g_R[boffset + 2] = f_R[3 * i + 2];
        g_R[i + 12] =
            -((f_R[i] * d2 + f_R[i + 3] * d3) + f_R[i + 6] * d4) / sinTheta;
      }
      g_R[3] = 0.0;
      g_R[7] = 0.0;
      g_R[11] = 0.0;
      g_R[15] = 1.0 / Tm[15];
      robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, g_R, Tji,
                                                             e_data);
      for (i = 0; i < 7; i++) {
        Jaci_data[i + 7 * k] =
            (deltavec[i] - e_data[i]) * 4.9999999999999994E+8;
      }
      for (i = 0; i < 3; i++) {
        boffset = i << 2;
        g_R[boffset] = c_R[3 * i];
        g_R[boffset + 1] = c_R[3 * i + 1];
        g_R[boffset + 2] = c_R[3 * i + 2];
        g_R[i + 12] = t[i];
      }
      g_R[3] = 0.0;
      g_R[7] = 0.0;
      g_R[11] = 0.0;
      g_R[15] = s;
      robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, Sio, g_R,
                                                             deltavec);
      for (i = 0; i < 3; i++) {
        d2 = 0.0;
        for (j = 0; j < 3; j++) {
          coffset = j << 2;
          f_R[i + 3 * j] = (Tji[i] * deltatform[coffset] +
                            Tji[i + 4] * deltatform[coffset + 1]) +
                           Tji[i + 8] * deltatform[coffset + 2];
          d2 += sinThetaDy * Tji[i + coffset] * deltatform[j + 12];
        }
        e_R[i] = d2 + Tji[i + 12];
      }
      for (i = 0; i < 3; i++) {
        boffset = i << 2;
        g_R[boffset] = f_R[3 * i];
        g_R[boffset + 1] = f_R[3 * i + 1];
        g_R[boffset + 2] = f_R[3 * i + 2];
        g_R[i + 12] = e_R[i];
      }
      g_R[3] = 0.0;
      g_R[7] = 0.0;
      g_R[11] = 0.0;
      g_R[15] = d1 * deltatform[15];
      robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, Sio, g_R,
                                                             e_data);
      for (i = 0; i < 7; i++) {
        Jacj_data[i + 7 * k] =
            (deltavec[i] - e_data[i]) * 4.9999999999999994E+8;
      }
    }
    q = Toi[Toi.size(0) * 3 + 3];
    d = Tjop[12];
    d1 = Tjop[13];
    d2 = Tjop[14];
    for (i = 0; i < 3; i++) {
      boffset = i << 2;
      g_R[boffset] = R[3 * i];
      g_R[boffset + 1] = R[3 * i + 1];
      g_R[boffset + 2] = R[3 * i + 2];
      g_R[i + 12] = -((R[i] * d + R[i + 3] * d1) + R[i + 6] * d2) / q;
    }
    g_R[3] = 0.0;
    g_R[7] = 0.0;
    g_R[11] = 0.0;
    g_R[15] = 1.0 / Toi[Toi.size(0) * 3 + 3];
    robotics::core::internal::Sim3Helpers::multiplyLogSim3(Sji, g_R, Tji,
                                                           deltavec);
    b_i = 7;
    for (i = 0; i < 7; i++) {
      e_data[i] = deltavec[i];
    }
    Jaci_size[0] = 7;
    Jaci_size[1] = 7;
    Jacj_size[0] = 7;
    Jacj_size[1] = 7;
  }
  coffset = Omega.size(1);
  aoffset = Omega.size(1);
  for (j = 0; j < coffset; j++) {
    boffset = j * Omega.size(0);
    e_R[j] = 0.0;
    for (int k{0}; k < b_i; k++) {
      e_R[j] += e_data[k] * Omega[boffset + k];
    }
  }
  cost = 0.0;
  for (i = 0; i < aoffset; i++) {
    cost += e_R[i] * e_data[i];
  }
  ::poseGraphOptimize::coder::internal::blas::mtimes(
      Jaci_data, Jaci_size, Omega, y_tmp_data, y_tmp_size);
  coffset = y_tmp_size[0] - 1;
  i = y_tmp_size[1];
  gradi_size = y_tmp_size[0];
  std::memset(&gradi_data[0], 0,
              static_cast<unsigned int>(coffset + 1) * sizeof(double));
  for (int k{0}; k < i; k++) {
    aoffset = k * y_tmp_size[0];
    for (b_i = 0; b_i <= coffset; b_i++) {
      gradi_data[b_i] += y_tmp_data[aoffset + b_i] * e_data[k];
    }
  }
  ::poseGraphOptimize::coder::internal::blas::mtimes(
      Jacj_data, Jacj_size, Omega, b_y_tmp_data, b_y_tmp_size);
  coffset = b_y_tmp_size[0] - 1;
  i = b_y_tmp_size[1];
  gradj_size = b_y_tmp_size[0];
  std::memset(&gradj_data[0], 0,
              static_cast<unsigned int>(coffset + 1) * sizeof(double));
  for (int k{0}; k < i; k++) {
    aoffset = k * b_y_tmp_size[0];
    for (b_i = 0; b_i <= coffset; b_i++) {
      gradj_data[b_i] += b_y_tmp_data[aoffset + b_i] * e_data[k];
    }
  }
  ::poseGraphOptimize::coder::internal::blas::mtimes(
      y_tmp_data, y_tmp_size, Jaci_data, Jaci_size, hessii_data, hessii_size);
  ::poseGraphOptimize::coder::internal::blas::mtimes(
      y_tmp_data, y_tmp_size, Jacj_data, Jacj_size, hessij_data, hessij_size);
  ::poseGraphOptimize::coder::internal::blas::mtimes(b_y_tmp_data, b_y_tmp_size,
                                                     Jaci_data, Jaci_size,
                                                     hessji_data, hessji_size);
  ::poseGraphOptimize::coder::internal::blas::mtimes(b_y_tmp_data, b_y_tmp_size,
                                                     Jacj_data, Jacj_size,
                                                     hessjj_data, hessjj_size);
  return cost;
}

void BlockInserter2::insertGradientBlock(double i, const double blocki_data[],
                                         int blocki_size)
{
  ::coder::array<double, 1U> obj;
  double d;
  double rowStart;
  int b_i;
  int i1;
  int i2;
  int loop_ub;
  rowStart = NodeMap[static_cast<int>(i) - 1];
  d = (rowStart + NodeDims[static_cast<int>(i) - 1]) - 1.0;
  if (rowStart > d) {
    b_i = 0;
    i1 = 0;
    i2 = 0;
  } else {
    b_i = static_cast<int>(rowStart) - 1;
    i1 = static_cast<int>(d);
    i2 = static_cast<int>(rowStart) - 1;
  }
  loop_ub = i1 - b_i;
  if (loop_ub == blocki_size) {
    obj.set_size(loop_ub);
    for (i1 = 0; i1 < loop_ub; i1++) {
      obj[i1] = Gradient[b_i + i1] + blocki_data[i1];
    }
    loop_ub = obj.size(0);
    for (b_i = 0; b_i < loop_ub; b_i++) {
      Gradient[i2 + b_i] = obj[b_i];
    }
  } else {
    binary_expand_op(this, i2, b_i, i1 - 1, blocki_data, blocki_size);
  }
}

double PoseGraphOptimizer::parseOptimizePoseGraphInputs(
    double &paramStruct_MaxTime, double &paramStruct_FunctionTolerance,
    boolean_T &paramStruct_IsVerbose, double &paramStruct_GradientTolerance,
    double &paramStruct_StepTolerance,
    double &paramStruct_InitialTrustRegionRadius,
    double paramStruct_FirstNodePose[3],
    double &paramStruct_TrustRegionRadiusTolerance,
    double &paramStruct_SolverID)
{
  static const char cv[128]{
      '\x00', '\x01', '\x02', '\x03', '\x04', '\x05', '\x06', '\a',   '\b',
      '\t',   '\n',   '\v',   '\f',   '\r',   '\x0e', '\x0f', '\x10', '\x11',
      '\x12', '\x13', '\x14', '\x15', '\x16', '\x17', '\x18', '\x19', '\x1a',
      '\x1b', '\x1c', '\x1d', '\x1e', '\x1f', ' ',    '!',    '\"',   '#',
      '$',    '%',    '&',    '\'',   '(',    ')',    '*',    '+',    ',',
      '-',    '.',    '/',    '0',    '1',    '2',    '3',    '4',    '5',
      '6',    '7',    '8',    '9',    ':',    ';',    '<',    '=',    '>',
      '?',    '@',    'a',    'b',    'c',    'd',    'e',    'f',    'g',
      'h',    'i',    'j',    'k',    'l',    'm',    'n',    'o',    'p',
      'q',    'r',    's',    't',    'u',    'v',    'w',    'x',    'y',
      'z',    '[',    '\\',   ']',    '^',    '_',    '`',    'a',    'b',
      'c',    'd',    'e',    'f',    'g',    'h',    'i',    'j',    'k',
      'l',    'm',    'n',    'o',    'p',    'q',    'r',    's',    't',
      'u',    'v',    'w',    'x',    'y',    'z',    '{',    '|',    '}',
      '~',    '\x7f'};
  static const char cv1[3]{'o', 'f', 'f'};
  static const char cv2[2]{'o', 'n'};
  ::coder::array<char, 2U> str;
  double paramStruct_MaxNumIteration;
  int exitg1;
  int kstr;
  char out_data[3];
  boolean_T b_bool;
  paramStruct_MaxNumIteration = 300.0;
  paramStruct_MaxTime = 500.0;
  paramStruct_FunctionTolerance = 1.0E-8;
  str.set_size(1, 3);
  str[0] = 'o';
  str[1] = 'f';
  str[2] = 'f';
  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr <= 2) {
      if (cv[static_cast<int>(str[kstr])] != cv[static_cast<int>(cv1[kstr])]) {
        exitg1 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  if (b_bool) {
    char partial_match_data[3];
    partial_match_data[0] = 'o';
    partial_match_data[1] = 'f';
    partial_match_data[2] = 'f';
    kstr = 3;
    for (int i{0}; i < 3; i++) {
      out_data[i] = partial_match_data[i];
    }
  } else {
    kstr = 2;
    out_data[0] = ' ';
    out_data[1] = ' ';
  }
  paramStruct_IsVerbose = false;
  if (kstr == 2) {
    kstr = 0;
    do {
      exitg1 = 0;
      if (kstr < 2) {
        if (out_data[kstr] != cv2[kstr]) {
          exitg1 = 1;
        } else {
          kstr++;
        }
      } else {
        paramStruct_IsVerbose = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  paramStruct_GradientTolerance = 5.0E-9;
  paramStruct_StepTolerance = 1.0E-12;
  paramStruct_InitialTrustRegionRadius = 100.0;
  paramStruct_FirstNodePose[0] = 0.0;
  paramStruct_FirstNodePose[1] = 0.0;
  paramStruct_FirstNodePose[2] = 0.0;
  paramStruct_TrustRegionRadiusTolerance = 1.0E-10;
  paramStruct_SolverID = -1.0;
  return paramStruct_MaxNumIteration;
}

double PoseGraphHelpers::poseGraphCost(
    const ::coder::array<double, 2U> &posesMat,
    const ::coder::array<double, 2U> &args_edgeNodePairs,
    const ::coder::array<double, 2U> &args_edgeMeasurements,
    const ::coder::array<double, 2U> &args_edgeInfoMats,
    const double args_tformSize[2], const double args_infoMatSize[2],
    double args_poseDeltaLength, const ::coder::array<double, 1U> &args_nodeMap,
    const ::coder::array<double, 1U> &args_nodeDims,
    const ::coder::array<boolean_T, 1U> &args_IsLandmarkNode,
    ::coder::array<double, 1U> &gradient, sparse &hessian)
{
  BlockInserter2 bi;
  robotics::core::internal::BlockMatrix edgeInfoMats;
  robotics::core::internal::BlockMatrix edgeMeasurements;
  robotics::core::internal::BlockMatrix poses;
  ::coder::array<double, 2U> OmegaIn;
  ::coder::array<double, 2U> Tij;
  ::coder::array<double, 2U> Toi;
  ::coder::array<double, 2U> Toj;
  ::coder::array<double, 2U> varargin_1;
  ::coder::array<double, 2U> y;
  ::coder::array<double, 1U> b_bi;
  ::coder::array<double, 1U> c_bi;
  ::coder::array<double, 1U> d_bi;
  ::coder::array<int, 2U> b_m1;
  ::coder::array<int, 2U> m2;
  ::coder::array<int, 2U> m3;
  ::coder::array<int, 2U> m4;
  ::coder::array<int, 2U> n1;
  ::coder::array<int, 2U> n2;
  ::coder::array<int, 2U> n3;
  ::coder::array<int, 2U> v1;
  ::coder::array<int, 2U> vk;
  ::coder::array<signed char, 2U> b_I;
  double H_data[196];
  double hessii_data[49];
  double hessij_data[49];
  double hessji_data[49];
  double hessjj_data[49];
  double gradi_data[7];
  double gradj_data[7];
  double cost;
  double maxNodeDim;
  double numEntries1;
  int i;
  int i1;
  int i2;
  int idx;
  int k;
  int last;
  int loop_ub;
  edgeMeasurements.Matrix.set_size(args_edgeMeasurements.size(0), 3);
  loop_ub = args_edgeMeasurements.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    edgeMeasurements.Matrix[i] = args_edgeMeasurements[i];
  }
  edgeMeasurements.BlockSize[0] = args_tformSize[0];
  edgeMeasurements.BlockSize[1] = args_tformSize[1];
  edgeMeasurements.NumRowBlocks =
      static_cast<double>(args_edgeMeasurements.size(0)) / args_tformSize[0];
  edgeMeasurements.NumColBlocks = 3.0 / args_tformSize[1];
  edgeInfoMats.Matrix.set_size(args_edgeInfoMats.size(0), 3);
  loop_ub = args_edgeInfoMats.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    edgeInfoMats.Matrix[i] = args_edgeInfoMats[i];
  }
  edgeInfoMats.BlockSize[0] = args_infoMatSize[0];
  edgeInfoMats.BlockSize[1] = args_infoMatSize[1];
  edgeInfoMats.NumRowBlocks =
      static_cast<double>(args_edgeInfoMats.size(0)) / args_infoMatSize[0];
  edgeInfoMats.NumColBlocks = 3.0 / args_infoMatSize[1];
  cost = 0.0;
  poses.Matrix.set_size(posesMat.size(0), 3);
  loop_ub = posesMat.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    poses.Matrix[i] = posesMat[i];
  }
  poses.BlockSize[0] = args_tformSize[0];
  poses.BlockSize[1] = args_tformSize[1];
  poses.NumRowBlocks =
      static_cast<double>(posesMat.size(0)) / args_tformSize[0];
  poses.NumColBlocks = 3.0 / args_tformSize[1];
  last = static_cast<int>(poses.NumRowBlocks * args_poseDeltaLength);
  bi.Gradient.set_size(last);
  for (i = 0; i < last; i++) {
    bi.Gradient[i] = 0.0;
  }
  bi.NodeDims.set_size(args_nodeDims.size(0));
  loop_ub = args_nodeDims.size(0);
  for (i = 0; i < loop_ub; i++) {
    bi.NodeDims[i] = args_nodeDims[i];
  }
  bi.NodeMap.set_size(args_nodeMap.size(0));
  loop_ub = args_nodeMap.size(0);
  for (i = 0; i < loop_ub; i++) {
    bi.NodeMap[i] = args_nodeMap[i];
  }
  last = args_nodeDims.size(0);
  if (args_nodeDims.size(0) <= 2) {
    if (args_nodeDims.size(0) == 1) {
      maxNodeDim = args_nodeDims[0];
    } else {
      maxNodeDim = args_nodeDims[args_nodeDims.size(0) - 1];
      if ((!(args_nodeDims[0] < maxNodeDim)) &&
          ((!std::isnan(args_nodeDims[0])) || std::isnan(maxNodeDim))) {
        maxNodeDim = args_nodeDims[0];
      }
    }
  } else {
    if (!std::isnan(args_nodeDims[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(args_nodeDims[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      maxNodeDim = args_nodeDims[0];
    } else {
      maxNodeDim = args_nodeDims[idx - 1];
      i = idx + 1;
      for (k = i; k <= last; k++) {
        numEntries1 = args_nodeDims[k - 1];
        if (maxNodeDim < numEntries1) {
          maxNodeDim = numEntries1;
        }
      }
    }
  }
  i = static_cast<int>(
      (4.0 * static_cast<double>(args_edgeNodePairs.size(0)) + 1.0) *
      maxNodeDim * maxNodeDim);
  bi.HessianCSC.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    bi.HessianCSC[i] = 0.0;
  }
  bi.HessianCSCCount = 1.0;
  i = args_edgeNodePairs.size(0);
  for (k = 0; k < i; k++) {
    double b_i;
    double j;
    double numEntries2;
    double numEntries2_tmp;
    int hessii_size[2];
    int hessij_size[2];
    int hessji_size[2];
    int hessjj_size[2];
    int b_loop_ub;
    int c_loop_ub;
    int m1;
    int m2_idx_0;
    int m3_idx_0;
    int varargin_1_tmp;
    signed char i3;
    boolean_T b;
    b_i = args_edgeNodePairs[k];
    j = args_edgeNodePairs[k + args_edgeNodePairs.size(0)];
    edgeMeasurements.extractBlock(static_cast<double>(k) + 1.0, Tij);
    edgeInfoMats.extractBlock(static_cast<double>(k) + 1.0, OmegaIn);
    b = args_IsLandmarkNode[static_cast<int>(j) - 1];
    if (b) {
      maxNodeDim = args_nodeDims[static_cast<int>(j) - 1];
      if (maxNodeDim < 1.0) {
        last = 0;
      } else {
        last = static_cast<int>(maxNodeDim);
      }
      for (i1 = 0; i1 < last; i1++) {
        for (i2 = 0; i2 < last; i2++) {
          OmegaIn[i2 + last * i1] = OmegaIn[i2 + OmegaIn.size(0) * i1];
        }
      }
      OmegaIn.set_size(last, last);
    }
    poses.extractBlock(b_i, Toi);
    poses.extractBlock(j, Toj);
    maxNodeDim = PoseGraphHelpers::costBetweenTwoNodes(
        Toi, Toj, Tij, OmegaIn, b, gradi_data, idx, gradj_data, last,
        hessii_data, hessii_size, hessij_data, hessij_size, hessji_data,
        hessji_size, hessjj_data, hessjj_size);
    cost += maxNodeDim;
    bi.insertGradientBlock(b_i, gradi_data, idx);
    bi.insertGradientBlock(j, gradj_data, last);
    maxNodeDim = bi.NodeDims[static_cast<int>(b_i) - 1];
    numEntries1 = maxNodeDim * maxNodeDim;
    if (std::isnan(numEntries1)) {
      y.set_size(1, 1);
      y[0] = rtNaN;
    } else if (numEntries1 < 1.0) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, static_cast<int>(numEntries1 - 1.0) + 1);
      loop_ub = static_cast<int>(numEntries1 - 1.0);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        y[i1] = static_cast<double>(i1) + 1.0;
      }
    }
    v1.set_size(1, y.size(1));
    loop_ub = y.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      v1[i1] = static_cast<int>(y[i1]) - 1;
    }
    vk.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      vk[i1] = div_s32(v1[i1], static_cast<int>(maxNodeDim));
    }
    v1.set_size(1, v1.size(1));
    loop_ub = v1.size(1) - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      v1[i1] = v1[i1] - vk[i1] * static_cast<int>(maxNodeDim);
    }
    b_m1.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    n1.set_size(1, vk.size(1));
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_m1[i1] = v1[i1] + 1;
      n1[i1] = vk[i1] + 1;
    }
    numEntries2_tmp = bi.NodeDims[static_cast<int>(j) - 1];
    numEntries2 = maxNodeDim * numEntries2_tmp;
    b = std::isnan(numEntries2);
    if (b) {
      y.set_size(1, 1);
      y[0] = rtNaN;
    } else if (numEntries2 < 1.0) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, static_cast<int>(numEntries2 - 1.0) + 1);
      loop_ub = static_cast<int>(numEntries2 - 1.0);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        y[i1] = static_cast<double>(i1) + 1.0;
      }
    }
    v1.set_size(1, y.size(1));
    loop_ub = y.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      v1[i1] = static_cast<int>(y[i1]) - 1;
    }
    vk.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      vk[i1] = div_s32(v1[i1], static_cast<int>(maxNodeDim));
    }
    v1.set_size(1, v1.size(1));
    loop_ub = v1.size(1) - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      v1[i1] = v1[i1] - vk[i1] * static_cast<int>(maxNodeDim);
    }
    m2.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    n2.set_size(1, vk.size(1));
    for (i1 = 0; i1 < loop_ub; i1++) {
      m2[i1] = v1[i1] + 1;
      n2[i1] = vk[i1] + 1;
    }
    if (b) {
      y.set_size(1, 1);
      y[0] = rtNaN;
    } else if (numEntries2 < 1.0) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, static_cast<int>(numEntries2 - 1.0) + 1);
      loop_ub = static_cast<int>(numEntries2 - 1.0);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        y[i1] = static_cast<double>(i1) + 1.0;
      }
    }
    v1.set_size(1, y.size(1));
    loop_ub = y.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      v1[i1] = static_cast<int>(y[i1]) - 1;
    }
    vk.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      vk[i1] = div_s32(v1[i1], static_cast<int>(numEntries2_tmp));
    }
    v1.set_size(1, v1.size(1));
    loop_ub = v1.size(1) - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      v1[i1] = v1[i1] - vk[i1] * static_cast<int>(numEntries2_tmp);
    }
    m3.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    n3.set_size(1, vk.size(1));
    for (i1 = 0; i1 < loop_ub; i1++) {
      m3[i1] = v1[i1] + 1;
      n3[i1] = vk[i1] + 1;
    }
    maxNodeDim = numEntries2_tmp * numEntries2_tmp;
    if (std::isnan(maxNodeDim)) {
      y.set_size(1, 1);
      y[0] = rtNaN;
    } else if (maxNodeDim < 1.0) {
      y.set_size(1, 0);
    } else {
      y.set_size(1, static_cast<int>(maxNodeDim - 1.0) + 1);
      loop_ub = static_cast<int>(maxNodeDim - 1.0);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        y[i1] = static_cast<double>(i1) + 1.0;
      }
    }
    v1.set_size(1, y.size(1));
    loop_ub = y.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      v1[i1] = static_cast<int>(y[i1]) - 1;
    }
    vk.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      vk[i1] = div_s32(v1[i1], static_cast<int>(numEntries2_tmp));
    }
    v1.set_size(1, v1.size(1));
    loop_ub = v1.size(1) - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      v1[i1] = v1[i1] - vk[i1] * static_cast<int>(numEntries2_tmp);
    }
    m4.set_size(1, v1.size(1));
    loop_ub = v1.size(1);
    y.set_size(1, vk.size(1));
    for (i1 = 0; i1 < loop_ub; i1++) {
      m4[i1] = v1[i1] + 1;
      y[i1] = vk[i1] + 1;
    }
    loop_ub = hessii_size[0] * hessii_size[1];
    idx = hessij_size[0] * hessij_size[1];
    last = hessji_size[0] * hessji_size[1];
    b_loop_ub = hessjj_size[0] * hessjj_size[1];
    if (loop_ub - 1 >= 0) {
      std::copy(&hessii_data[0], &hessii_data[loop_ub], &H_data[0]);
    }
    for (i1 = 0; i1 < idx; i1++) {
      H_data[i1 + loop_ub] = hessij_data[i1];
    }
    for (i1 = 0; i1 < last; i1++) {
      H_data[(i1 + loop_ub) + idx] = hessji_data[i1];
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      H_data[((i1 + loop_ub) + idx) + last] = hessjj_data[i1];
    }
    numEntries1 =
        bi.HessianCSCCount + ((numEntries1 + 2.0 * numEntries2) + maxNodeDim);
    if (bi.HessianCSCCount > numEntries1 - 1.0) {
      i1 = 0;
    } else {
      i1 = static_cast<int>(bi.HessianCSCCount) - 1;
    }
    m1 = b_m1.size(1);
    m2_idx_0 = m2.size(1);
    m3_idx_0 = m3.size(1);
    varargin_1.set_size(((b_m1.size(1) + m2.size(1)) + m3.size(1)) + m4.size(1),
                        2);
    c_loop_ub = b_m1.size(1);
    for (i2 = 0; i2 < c_loop_ub; i2++) {
      maxNodeDim = bi.NodeMap[static_cast<int>(b_i) - 1] - 1.0;
      varargin_1[i2] = maxNodeDim + static_cast<double>(b_m1[i2]);
      varargin_1[i2 + varargin_1.size(0)] =
          maxNodeDim + static_cast<double>(n1[i2]);
    }
    c_loop_ub = m2.size(1);
    for (i2 = 0; i2 < c_loop_ub; i2++) {
      varargin_1_tmp = i2 + m1;
      varargin_1[varargin_1_tmp] =
          (bi.NodeMap[static_cast<int>(b_i) - 1] - 1.0) +
          static_cast<double>(m2[i2]);
      varargin_1[varargin_1_tmp + varargin_1.size(0)] =
          (bi.NodeMap[static_cast<int>(j) - 1] - 1.0) +
          static_cast<double>(n2[i2]);
    }
    c_loop_ub = m3.size(1);
    for (i2 = 0; i2 < c_loop_ub; i2++) {
      varargin_1_tmp = (i2 + m1) + m2_idx_0;
      varargin_1[varargin_1_tmp] = (bi.NodeMap[static_cast<int>(j) - 1] - 1.0) +
                                   static_cast<double>(m3[i2]);
      varargin_1[varargin_1_tmp + varargin_1.size(0)] =
          (bi.NodeMap[static_cast<int>(b_i) - 1] - 1.0) +
          static_cast<double>(n3[i2]);
    }
    c_loop_ub = m4.size(1);
    for (i2 = 0; i2 < c_loop_ub; i2++) {
      maxNodeDim = bi.NodeMap[static_cast<int>(j) - 1] - 1.0;
      varargin_1_tmp = ((i2 + m1) + m2_idx_0) + m3_idx_0;
      varargin_1[varargin_1_tmp] = maxNodeDim + static_cast<double>(m4[i2]);
      varargin_1[varargin_1_tmp + varargin_1.size(0)] = maxNodeDim + y[i2];
    }
    if (varargin_1.size(0) != 0) {
      idx = varargin_1.size(0);
    } else {
      idx = ((loop_ub + idx) + last) + b_loop_ub;
    }
    if ((idx == 0) || (varargin_1.size(0) != 0)) {
      i3 = 2;
    } else {
      i3 = 0;
    }
    loop_ub = i3;
    for (i2 = 0; i2 < loop_ub; i2++) {
      for (last = 0; last < idx; last++) {
        bi.HessianCSC[(i1 + last) + bi.HessianCSC.size(0) * i2] =
            varargin_1[last + idx * i2];
      }
    }
    for (i2 = 0; i2 < idx; i2++) {
      bi.HessianCSC[(i1 + i2) + bi.HessianCSC.size(0) * i3] = H_data[i2];
    }
    bi.HessianCSCCount = numEntries1;
  }
  if (args_poseDeltaLength < 0.0) {
    maxNodeDim = 0.0;
    idx = 0;
  } else {
    maxNodeDim = args_poseDeltaLength;
    idx = static_cast<int>(args_poseDeltaLength);
  }
  b_I.set_size(static_cast<int>(maxNodeDim), static_cast<int>(maxNodeDim));
  last = static_cast<int>(maxNodeDim) * static_cast<int>(maxNodeDim);
  for (i = 0; i < last; i++) {
    b_I[i] = 0;
  }
  if (static_cast<int>(maxNodeDim) > 0) {
    for (k = 0; k < idx; k++) {
      b_I[k + b_I.size(0) * k] = 1;
    }
  }
  if (last < 1) {
    y.set_size(1, 0);
  } else {
    y.set_size(1, last);
    loop_ub = last - 1;
    for (i = 0; i <= loop_ub; i++) {
      y[i] = static_cast<double>(i) + 1.0;
    }
  }
  idx = b_I.size(0);
  v1.set_size(1, y.size(1));
  loop_ub = y.size(1);
  for (i = 0; i < loop_ub; i++) {
    v1[i] = static_cast<int>(y[i]) - 1;
  }
  vk.set_size(1, v1.size(1));
  loop_ub = v1.size(1);
  for (i = 0; i < loop_ub; i++) {
    vk[i] = div_s32(v1[i], idx);
  }
  v1.set_size(1, v1.size(1));
  loop_ub = v1.size(1) - 1;
  for (i = 0; i <= loop_ub; i++) {
    v1[i] = v1[i] - vk[i] * idx;
  }
  y.set_size(1, v1.size(1));
  loop_ub = v1.size(1);
  m4.set_size(1, vk.size(1));
  for (i = 0; i < loop_ub; i++) {
    y[i] = v1[i] + 1;
    m4[i] = vk[i] + 1;
  }
  numEntries1 = bi.HessianCSCCount + static_cast<double>(last);
  if (bi.HessianCSCCount > numEntries1 - 1.0) {
    i = 0;
  } else {
    i = static_cast<int>(bi.HessianCSCCount) - 1;
  }
  loop_ub = y.size(1);
  for (i1 = 0; i1 < loop_ub; i1++) {
    i2 = i + i1;
    bi.HessianCSC[i2] = (bi.NodeMap[0] + y[i1]) - 1.0;
    bi.HessianCSC[i2 + bi.HessianCSC.size(0)] =
        (bi.NodeMap[0] + static_cast<double>(m4[i1])) - 1.0;
  }
  for (i1 = 0; i1 < last; i1++) {
    bi.HessianCSC[(i + i1) + bi.HessianCSC.size(0) * 2] = b_I[i1];
  }
  maxNodeDim = (args_nodeMap[static_cast<int>(poses.NumRowBlocks) - 1] +
                args_nodeDims[static_cast<int>(poses.NumRowBlocks) - 1]) -
               1.0;
  if (maxNodeDim < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(maxNodeDim);
  }
  gradient.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    gradient[i] = bi.Gradient[i];
  }
  if (numEntries1 - 1.0 < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(numEntries1 - 1.0);
  }
  b_bi.set_size(loop_ub);
  c_bi.set_size(loop_ub);
  d_bi.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    b_bi[i] = bi.HessianCSC[i];
    c_bi[i] = bi.HessianCSC[i + bi.HessianCSC.size(0)];
    d_bi[i] = bi.HessianCSC[i + bi.HessianCSC.size(0) * 2];
  }
  b_sparse(b_bi, c_bi, d_bi, hessian);
  return cost;
}

} // namespace internal
} // namespace algs
} // namespace nav
namespace robotics {
namespace core {
namespace internal {
boolean_T TrustRegionIndefiniteDogLegSE2::computeBasicSteps(
    const ::coder::array<double, 1U> &grad, const sparse &B,
    ::coder::array<double, 1U> &stepSD,
    ::coder::array<double, 1U> &stepGN) const
{
  sparse in;
  ::coder::array<double, 2U> a;
  ::coder::array<double, 1U> b_stepGN;
  double b_grad;
  double cd;
  int apend;
  int i;
  int loop_ub_tmp;
  a.set_size(1, B.n);
  loop_ub_tmp = B.n;
  for (i = 0; i < loop_ub_tmp; i++) {
    a[i] = 0.0;
  }
  if ((grad.size(0) != 0) && (B.n != 0) &&
      (B.colidx[B.colidx.size(0) - 1] - 1 != 0)) {
    for (int k{0}; k < loop_ub_tmp; k++) {
      cd = 0.0;
      apend = B.colidx[k + 1] - 1;
      i = B.colidx[k];
      for (int ap{i}; ap <= apend; ap++) {
        cd += B.d[ap - 1] * grad[B.rowidx[ap - 1] - 1];
      }
      a[k] = cd;
    }
  }
  b_grad = 0.0;
  apend = grad.size(0);
  for (i = 0; i < apend; i++) {
    cd = grad[i];
    b_grad += cd * cd;
  }
  cd = 0.0;
  apend = a.size(1);
  for (i = 0; i < apend; i++) {
    cd += a[i] * grad[i];
  }
  cd = b_grad / cd;
  stepGN.set_size(grad.size(0));
  apend = grad.size(0);
  for (i = 0; i < apend; i++) {
    stepGN[i] = -grad[i];
  }
  stepSD.set_size(stepGN.size(0));
  apend = stepGN.size(0);
  for (i = 0; i < apend; i++) {
    stepSD[i] = cd * stepGN[i];
  }
  if ((B.m == 0) || (B.n == 0) || (stepGN.size(0) == 0)) {
    stepGN.set_size(B.n);
    for (i = 0; i < loop_ub_tmp; i++) {
      stepGN[i] = 0.0;
    }
  } else if (stepGN.size(0) == B.n) {
    cs_di *cxA;
    cs_din *N;
    cs_dis *S;
    if (B.m < B.n) {
      B.ctranspose(in);
      cxA = makeCXSparseMatrix(in.colidx[in.colidx.size(0) - 1] - 1, in.n, in.m,
                               &(in.colidx.data())[0], &(in.rowidx.data())[0],
                               &(in.d.data())[0]);
    } else {
      cxA = makeCXSparseMatrix(
          B.colidx[B.colidx.size(0) - 1] - 1, B.n, B.m,
          &(((::coder::array<int, 1U> *)&B.colidx)->data())[0],
          &(((::coder::array<int, 1U> *)&B.rowidx)->data())[0],
          &(((::coder::array<double, 1U> *)&B.d)->data())[0]);
    }
    S = cs_di_sqr(2, cxA, 0);
    N = cs_di_lu(cxA, S, 1);
    cs_di_spfree(cxA);
    if (N == nullptr) {
      cs_di_sfree(S);
      cs_di_nfree(N);
      b_stepGN.set_size(stepGN.size(0));
      apend = stepGN.size(0) - 1;
      for (i = 0; i <= apend; i++) {
        b_stepGN[i] = stepGN[i];
      }
      ::poseGraphOptimize::coder::internal::CXSparseAPI::iteratedQR(
          B, b_stepGN, B.n, stepGN);
    } else {
      solve_from_lu_di(N, S, (double *)&(stepGN.data())[0], stepGN.size(0));
      cs_di_sfree(S);
      cs_di_nfree(N);
    }
  } else {
    b_stepGN.set_size(stepGN.size(0));
    apend = stepGN.size(0) - 1;
    for (i = 0; i <= apend; i++) {
      b_stepGN[i] = stepGN[i];
    }
    ::poseGraphOptimize::coder::internal::CXSparseAPI::iteratedQR(B, b_stepGN,
                                                                  B.n, stepGN);
  }
  return b_norm(grad) < GradientTolerance;
}

void SEHelpers::expSE3hat(const double e[6], double T[16])
{
  double Sphi[9];
  double V[9];
  double absxk;
  double c;
  double c_tmp;
  double d;
  double d1;
  double d2;
  double scale;
  double t;
  double theta;
  int T_tmp;
  Sphi[0] = 0.0;
  Sphi[3] = -e[5];
  Sphi[6] = e[4];
  Sphi[1] = e[5];
  Sphi[4] = 0.0;
  Sphi[7] = -e[3];
  Sphi[2] = -e[4];
  Sphi[5] = e[3];
  Sphi[8] = 0.0;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(e[3]);
  if (absxk > 3.3121686421112381E-170) {
    theta = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    theta = t * t;
  }
  absxk = std::abs(e[4]);
  if (absxk > scale) {
    t = scale / absxk;
    theta = theta * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    theta += t * t;
  }
  absxk = std::abs(e[5]);
  if (absxk > scale) {
    t = scale / absxk;
    theta = theta * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    theta += t * t;
  }
  theta = scale * std::sqrt(theta);
  scale = theta * theta;
  absxk = 1.0 - std::cos(theta);
  t = absxk / scale;
  c_tmp = std::sin(theta);
  c = (theta - c_tmp) / (scale * theta);
  if ((t < 2.2204460492503131E-16) || (std::isinf(t) || std::isnan(t))) {
    std::memset(&V[0], 0, 9U * sizeof(double));
    V[0] = 1.0;
    V[4] = 1.0;
    V[8] = 1.0;
    std::copy(&V[0], &V[9], &Sphi[0]);
  } else {
    double b_Sphi[9];
    signed char V_tmp[9];
    for (int i{0}; i < 9; i++) {
      V_tmp[i] = 0;
    }
    V_tmp[0] = 1;
    V_tmp[4] = 1;
    V_tmp[8] = 1;
    for (int i{0}; i < 3; i++) {
      d = Sphi[i + 3];
      d1 = Sphi[i + 6];
      for (T_tmp = 0; T_tmp < 3; T_tmp++) {
        double d3;
        double d4;
        int b_V_tmp;
        d2 = Sphi[3 * T_tmp];
        d3 = Sphi[3 * T_tmp + 1];
        d4 = Sphi[3 * T_tmp + 2];
        b_V_tmp = i + 3 * T_tmp;
        V[b_V_tmp] = (static_cast<double>(V_tmp[b_V_tmp]) + t * Sphi[b_V_tmp]) +
                     ((c * Sphi[i] * d2 + c * d * d3) + c * d1 * d4);
        b_Sphi[b_V_tmp] = ((Sphi[i] * d2 + d * d3) + d1 * d4) / scale;
      }
    }
    for (int i{0}; i < 9; i++) {
      Sphi[i] = (static_cast<double>(V_tmp[i]) + Sphi[i] * c_tmp / theta) +
                b_Sphi[i] * absxk;
    }
  }
  d = e[0];
  d1 = e[1];
  d2 = e[2];
  for (int i{0}; i < 3; i++) {
    T_tmp = i << 2;
    T[T_tmp] = Sphi[3 * i];
    T[T_tmp + 1] = Sphi[3 * i + 1];
    T[T_tmp + 2] = Sphi[3 * i + 2];
    T[i + 12] = (V[i] * d + V[i + 3] * d1) + V[i + 6] * d2;
  }
  T[3] = 0.0;
  T[7] = 0.0;
  T[11] = 0.0;
  T[15] = 1.0;
}

void BlockMatrix::extractBlock(double i, ::coder::array<double, 2U> &B) const
{
  double colStart;
  double d;
  double rowStart;
  int b_i;
  int b_loop_ub;
  int i1;
  int loop_ub;
  rowStart = BlockSize[0] * (i - 1.0) + 1.0;
  colStart = BlockSize[1] * 0.0 + 1.0;
  d = (rowStart + BlockSize[0]) - 1.0;
  if (rowStart > d) {
    b_i = 0;
    i1 = 0;
  } else {
    b_i = static_cast<int>(rowStart) - 1;
    i1 = static_cast<int>(d);
  }
  d = (colStart + BlockSize[1]) - 1.0;
  if (colStart > d) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(d);
  }
  b_loop_ub = i1 - b_i;
  B.set_size(b_loop_ub, loop_ub);
  for (i1 = 0; i1 < loop_ub; i1++) {
    for (int i2{0}; i2 < b_loop_ub; i2++) {
      B[i2 + B.size(0) * i1] = Matrix[(b_i + i2) + Matrix.size(0) * i1];
    }
  }
}

void TrustRegionIndefiniteDogLegSE2::incrementX(
    const ::coder::array<double, 2U> &x,
    const ::coder::array<double, 1U> &epsilons,
    ::coder::array<double, 2U> &xNew) const
{
  BlockMatrix xBlk;
  b_BlockMatrix xBlkNew;
  ::coder::array<double, 2U> T;
  double b_T[3];
  int i;
  int i1;
  int loop_ub;
  xBlk.Matrix.set_size(x.size(0), 3);
  loop_ub = x.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    xBlk.Matrix[i] = x[i];
  }
  xBlk.BlockSize[0] = 3.0;
  xBlk.BlockSize[1] = 3.0;
  xBlk.NumRowBlocks = static_cast<double>(x.size(0)) / 3.0;
  xBlk.NumColBlocks = 1.0;
  i = static_cast<int>(xBlk.NumRowBlocks * 3.0);
  xBlkNew.Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    xBlkNew.Matrix[i] = 0.0;
  }
  i = static_cast<int>(xBlk.NumRowBlocks);
  for (int b_i{0}; b_i < i; b_i++) {
    double idx;
    double len;
    xBlk.extractBlock(static_cast<double>(b_i) + 1.0, T);
    idx = ExtraArgs.nodeMap[b_i];
    len = ExtraArgs.nodeDims[b_i];
    if (len == 3.0) {
      double T_tmp[9];
      int i2;
      int rowStart;
      if (idx > (idx + 3.0) - 1.0) {
        i1 = 0;
        i2 = -2;
      } else {
        i1 = static_cast<int>(idx) - 1;
        i2 = static_cast<int>(static_cast<unsigned int>(idx));
      }
      if ((i2 - i1) + 2 == 3) {
        b_T[0] = T[T.size(0) * 2];
        b_T[1] = T[T.size(0) * 2 + 1];
        b_T[2] = rt_atan2d_snf(T[1], T[0]);
        b_T[0] += epsilons[i1];
        b_T[1] += epsilons[i1 + 1];
        b_T[2] += epsilons[i1 + 2];
      } else {
        binary_expand_op(b_T, T, epsilons, i1, i2 + 1);
      }
      len = std::sin(b_T[2]);
      idx = std::cos(b_T[2]);
      rowStart = 3 * b_i + 1;
      if (static_cast<unsigned int>(rowStart) >
          static_cast<unsigned int>(rowStart) + 2U) {
        i1 = 0;
        i2 = 0;
      } else {
        i1 = rowStart - 1;
        i2 = static_cast<int>(static_cast<unsigned int>(rowStart) + 2U);
      }
      T_tmp[0] = idx;
      T_tmp[3] = -len;
      T_tmp[6] = b_T[0];
      T_tmp[1] = len;
      T_tmp[4] = idx;
      T_tmp[7] = b_T[1];
      T_tmp[2] = 0.0;
      T_tmp[5] = 0.0;
      T_tmp[8] = 1.0;
      loop_ub = i2 - i1;
      for (i2 = 0; i2 < 3; i2++) {
        for (rowStart = 0; rowStart < loop_ub; rowStart++) {
          xBlkNew.Matrix[(i1 + rowStart) + xBlkNew.Matrix.size(0) * i2] =
              T_tmp[rowStart + loop_ub * i2];
        }
      }
    } else {
      int i2;
      int rowStart;
      len = (idx + len) - 1.0;
      if (idx > len) {
        i1 = 0;
        i2 = 0;
      } else {
        i1 = static_cast<int>(idx) - 1;
        i2 = static_cast<int>(len);
      }
      if (i2 - i1 == 2) {
        T[T.size(0) * 2] = T[T.size(0) * 2] + epsilons[i1];
        T[T.size(0) * 2 + 1] = T[T.size(0) * 2 + 1] + epsilons[i1 + 1];
      } else {
        binary_expand_op(T, epsilons, i1, i2 - 1);
      }
      rowStart = 3 * b_i;
      if (static_cast<unsigned int>(rowStart + 1) >
          static_cast<unsigned int>(rowStart) + 3U) {
        rowStart = 0;
      }
      loop_ub = T.size(1);
      for (i1 = 0; i1 < loop_ub; i1++) {
        int b_loop_ub;
        b_loop_ub = T.size(0);
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          xBlkNew.Matrix[(rowStart + i2) + xBlkNew.Matrix.size(0) * i1] =
              T[i2 + T.size(0) * i1];
        }
      }
    }
  }
  xNew.set_size(xBlkNew.Matrix.size(0), 3);
  loop_ub = xBlkNew.Matrix.size(0);
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      xNew[i1 + xNew.size(0) * i] =
          xBlkNew.Matrix[i1 + xBlkNew.Matrix.size(0) * i];
    }
  }
}

void Sim3Helpers::multiplyLogSim3(const double S1[16], const double S2[16],
                                  const double S3[16], double e[7])
{
  double S12[16];
  double S123[16];
  double omega2[9];
  double w[9];
  double omega1[3];
  double a;
  double a1;
  double b1;
  double c;
  double d;
  double sigma;
  int b_omega2_tmp;
  int omega2_tmp;
  int r1;
  int r2;
  int r3;
  int rtemp;
  for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
    a1 = 0.0;
    for (int k{0}; k < 3; k++) {
      rtemp = k << 2;
      omega2[omega2_tmp + 3 * k] =
          (S1[omega2_tmp] * S2[rtemp] + S1[omega2_tmp + 4] * S2[rtemp + 1]) +
          S1[omega2_tmp + 8] * S2[rtemp + 2];
      a1 += S1[15] * S1[omega2_tmp + rtemp] * S2[k + 12];
    }
    omega1[omega2_tmp] = a1 + S1[omega2_tmp + 12];
  }
  for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
    rtemp = omega2_tmp << 2;
    S12[rtemp] = omega2[3 * omega2_tmp];
    S12[rtemp + 1] = omega2[3 * omega2_tmp + 1];
    S12[rtemp + 2] = omega2[3 * omega2_tmp + 2];
    S12[omega2_tmp + 12] = omega1[omega2_tmp];
  }
  S12[3] = 0.0;
  S12[7] = 0.0;
  S12[11] = 0.0;
  S12[15] = S1[15] * S2[15];
  for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
    a1 = 0.0;
    for (int k{0}; k < 3; k++) {
      rtemp = k << 2;
      omega2[omega2_tmp + 3 * k] =
          (S12[omega2_tmp] * S3[rtemp] + S12[omega2_tmp + 4] * S3[rtemp + 1]) +
          S12[omega2_tmp + 8] * S3[rtemp + 2];
      a1 += S12[15] * S12[omega2_tmp + rtemp] * S3[k + 12];
    }
    omega1[omega2_tmp] = a1 + S12[omega2_tmp + 12];
  }
  for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
    rtemp = omega2_tmp << 2;
    S123[rtemp] = omega2[3 * omega2_tmp];
    S123[rtemp + 1] = omega2[3 * omega2_tmp + 1];
    S123[rtemp + 2] = omega2[3 * omega2_tmp + 2];
    S123[omega2_tmp + 12] = omega1[omega2_tmp];
  }
  S123[3] = 0.0;
  S123[7] = 0.0;
  S123[11] = 0.0;
  S123[15] = S12[15] * S3[15];
  sigma = std::log(S123[15]);
  d = 0.5 * (((S123[0] + S123[5]) + S123[10]) - 1.0);
  if (std::abs(sigma) < 1.0E-5) {
    c = 1.0 - sigma / 2.0;
    if (d > 0.99999) {
      for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        b_omega2_tmp = omega2_tmp << 2;
        omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
        omega2[3 * omega2_tmp + 1] =
            S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
        omega2[3 * omega2_tmp + 2] =
            S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
      }
      omega1[0] = 0.5 * omega2[5];
      omega1[1] = 0.5 * omega2[6];
      omega1[2] = 0.5 * omega2[1];
      omega2[0] = 0.0;
      omega2[3] = -omega1[2];
      omega2[6] = omega1[1];
      omega2[1] = omega1[2];
      omega2[4] = 0.0;
      omega2[7] = -omega1[0];
      omega2[2] = -omega1[1];
      omega2[5] = omega1[0];
      omega2[8] = 0.0;
      a = 0.5;
      a1 = 0.16666666666666666;
    } else {
      double th;
      double thSquare;
      th = std::acos(d);
      thSquare = th * th;
      a = th / (2.0 * std::sqrt(1.0 - d * d));
      for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        b_omega2_tmp = omega2_tmp << 2;
        omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
        omega2[3 * omega2_tmp + 1] =
            S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
        omega2[3 * omega2_tmp + 2] =
            S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
      }
      omega1[0] = a * omega2[5];
      omega1[1] = a * omega2[6];
      omega1[2] = a * omega2[1];
      omega2[0] = 0.0;
      omega2[3] = -omega1[2];
      omega2[6] = omega1[1];
      omega2[1] = omega1[2];
      omega2[4] = 0.0;
      omega2[7] = -omega1[0];
      omega2[2] = -omega1[1];
      omega2[5] = omega1[0];
      omega2[8] = 0.0;
      a = (1.0 - std::cos(th)) / thSquare;
      a1 = (th - std::sin(th)) / (thSquare * th);
    }
  } else {
    c = (S123[15] - 1.0) / sigma;
    if (d > 0.99999) {
      a1 = sigma * sigma;
      for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        b_omega2_tmp = omega2_tmp << 2;
        omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
        omega2[3 * omega2_tmp + 1] =
            S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
        omega2[3 * omega2_tmp + 2] =
            S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
      }
      omega1[0] = 0.5 * omega2[5];
      omega1[1] = 0.5 * omega2[6];
      omega1[2] = 0.5 * omega2[1];
      omega2[0] = 0.0;
      omega2[3] = -omega1[2];
      omega2[6] = omega1[1];
      omega2[1] = omega1[2];
      omega2[4] = 0.0;
      omega2[7] = -omega1[0];
      omega2[2] = -omega1[1];
      omega2[5] = omega1[0];
      omega2[8] = 0.0;
      a = ((sigma - 1.0) * S123[15] + 1.0) / a1;
      a1 = (((0.5 * a1 - sigma) + 1.0) * S123[15] - 1.0) / (a1 * sigma);
    } else {
      double th;
      th = std::acos(d);
      a = th / (2.0 * std::sqrt(1.0 - d * d));
      for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
        b_omega2_tmp = omega2_tmp << 2;
        omega2[3 * omega2_tmp] = S123[b_omega2_tmp] - S123[omega2_tmp];
        omega2[3 * omega2_tmp + 1] =
            S123[b_omega2_tmp + 1] - S123[omega2_tmp + 4];
        omega2[3 * omega2_tmp + 2] =
            S123[b_omega2_tmp + 2] - S123[omega2_tmp + 8];
      }
      double thSquare;
      omega1[0] = a * omega2[5];
      omega1[1] = a * omega2[6];
      omega1[2] = a * omega2[1];
      omega2[0] = 0.0;
      omega2[3] = -omega1[2];
      omega2[6] = omega1[1];
      omega2[1] = omega1[2];
      omega2[4] = 0.0;
      omega2[7] = -omega1[0];
      omega2[2] = -omega1[1];
      omega2[5] = omega1[0];
      omega2[8] = 0.0;
      thSquare = th * th;
      a1 = S123[15] * std::sin(th);
      b1 = S123[15] * std::cos(th);
      d = thSquare + sigma * sigma;
      a = (a1 * sigma + (1.0 - b1) * th) / (th * d);
      a1 = (c - ((b1 - 1.0) * sigma + a1 * th) / d) / thSquare;
    }
  }
  for (omega2_tmp = 0; omega2_tmp < 3; omega2_tmp++) {
    for (int k{0}; k < 3; k++) {
      rtemp = omega2_tmp + 3 * k;
      w[rtemp] = (a * omega2[rtemp] +
                  ((a1 * omega2[omega2_tmp] * omega2[3 * k] +
                    a1 * omega2[omega2_tmp + 3] * omega2[3 * k + 1]) +
                   a1 * omega2[omega2_tmp + 6] * omega2[3 * k + 2])) +
                 c * static_cast<double>(iv[rtemp]);
    }
  }
  r1 = 0;
  r2 = 1;
  r3 = 2;
  a1 = std::abs(w[0]);
  d = std::abs(w[1]);
  if (d > a1) {
    a1 = d;
    r1 = 1;
    r2 = 0;
  }
  if (std::abs(w[2]) > a1) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }
  w[r2] /= w[r1];
  w[r3] /= w[r1];
  w[r2 + 3] -= w[r2] * w[r1 + 3];
  w[r3 + 3] -= w[r3] * w[r1 + 3];
  w[r2 + 6] -= w[r2] * w[r1 + 6];
  w[r3 + 6] -= w[r3] * w[r1 + 6];
  if (std::abs(w[r3 + 3]) > std::abs(w[r2 + 3])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }
  w[r3 + 3] /= w[r2 + 3];
  w[r3 + 6] -= w[r3 + 3] * w[r2 + 6];
  a1 = S123[12];
  d = S123[13];
  b1 = S123[14];
  for (int k{0}; k < 3; k++) {
    b_omega2_tmp = k + 3 * r1;
    omega2[b_omega2_tmp] = static_cast<double>(iv[k]) / w[r1];
    rtemp = k + 3 * r2;
    omega2[rtemp] =
        static_cast<double>(iv[k + 3]) - omega2[b_omega2_tmp] * w[r1 + 3];
    omega2_tmp = k + 3 * r3;
    omega2[omega2_tmp] =
        static_cast<double>(iv[k + 6]) - omega2[b_omega2_tmp] * w[r1 + 6];
    omega2[rtemp] /= w[r2 + 3];
    omega2[omega2_tmp] -= omega2[rtemp] * w[r2 + 6];
    omega2[omega2_tmp] /= w[r3 + 6];
    omega2[rtemp] -= omega2[omega2_tmp] * w[r3 + 3];
    omega2[b_omega2_tmp] -= omega2[omega2_tmp] * w[r3];
    omega2[b_omega2_tmp] -= omega2[rtemp] * w[r2];
    e[k] = (omega2[k] * a1 + omega2[k + 3] * d) + omega2[k + 6] * b1;
    e[k + 3] = omega1[k];
  }
  e[6] = sigma;
}

void BlockMatrix::replaceBlock(double i, const double blockij[9])
{
  double colStart;
  double d;
  double rowStart;
  int b_i;
  int i1;
  int loop_ub;
  int unnamed_idx_0;
  rowStart = BlockSize[0] * (i - 1.0) + 1.0;
  colStart = BlockSize[1] * 0.0 + 1.0;
  d = (rowStart + BlockSize[0]) - 1.0;
  if (rowStart > d) {
    b_i = 0;
    i1 = 0;
  } else {
    b_i = static_cast<int>(rowStart) - 1;
    i1 = static_cast<int>(d);
  }
  d = (colStart + BlockSize[1]) - 1.0;
  unnamed_idx_0 = i1 - b_i;
  if (colStart > d) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(d);
  }
  for (i1 = 0; i1 < loop_ub; i1++) {
    for (int i2{0}; i2 < unnamed_idx_0; i2++) {
      Matrix[(b_i + i2) + Matrix.size(0) * i1] =
          blockij[i2 + unnamed_idx_0 * i1];
    }
  }
}

void BlockMatrix::replaceBlock()
{
  double colStart;
  double d;
  double rowStart;
  int b_loop_ub;
  int loop_ub;
  rowStart = BlockSize[0] * 0.0 + 1.0;
  colStart = BlockSize[1] * 0.0 + 1.0;
  d = (rowStart + BlockSize[0]) - 1.0;
  if (rowStart > d) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(d);
  }
  d = (colStart + BlockSize[1]) - 1.0;
  if (colStart > d) {
    b_loop_ub = 0;
  } else {
    b_loop_ub = static_cast<int>(d);
  }
  for (int i{0}; i < b_loop_ub; i++) {
    for (int i1{0}; i1 < loop_ub; i1++) {
      Matrix[i1 + Matrix.size(0) * i] = iv[i1 + loop_ub * i];
    }
  }
}

void Sim3Helpers::sim3ToSform(const double minVecSim3[7], double S[16])
{
  double R[9];
  double w[9];
  double wSquare[9];
  double a_tmp;
  double absxk;
  double c;
  double s;
  double scale;
  double t;
  double th;
  int S_tmp;
  w[0] = 0.0;
  w[3] = -minVecSim3[5];
  w[6] = minVecSim3[4];
  w[1] = minVecSim3[5];
  w[4] = 0.0;
  w[7] = -minVecSim3[3];
  w[2] = -minVecSim3[4];
  w[5] = minVecSim3[3];
  w[8] = 0.0;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(minVecSim3[3]);
  if (absxk > 3.3121686421112381E-170) {
    th = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    th = t * t;
  }
  absxk = std::abs(minVecSim3[4]);
  if (absxk > scale) {
    t = scale / absxk;
    th = th * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    th += t * t;
  }
  absxk = std::abs(minVecSim3[5]);
  if (absxk > scale) {
    t = scale / absxk;
    th = th * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    th += t * t;
  }
  th = scale * std::sqrt(th);
  s = std::exp(minVecSim3[6]);
  for (int i{0}; i < 3; i++) {
    for (S_tmp = 0; S_tmp < 3; S_tmp++) {
      wSquare[i + 3 * S_tmp] =
          (w[i] * w[3 * S_tmp] + w[i + 3] * w[3 * S_tmp + 1]) +
          w[i + 6] * w[3 * S_tmp + 2];
    }
  }
  if (std::abs(minVecSim3[6]) < 1.0E-5) {
    c = 1.0;
    if (th < 1.0E-5) {
      a_tmp = 0.5;
      t = 0.16666666666666666;
      std::memset(&R[0], 0, 9U * sizeof(double));
      R[0] = 1.0;
      R[4] = 1.0;
      R[8] = 1.0;
      for (int i{0}; i < 9; i++) {
        R[i] = (R[i] + w[i]) + wSquare[i] / 2.0;
      }
    } else {
      double thSquare;
      thSquare = th * th;
      a_tmp = (1.0 - std::cos(th)) / thSquare;
      scale = std::sin(th);
      t = (th - scale) / (thSquare * th);
      absxk = scale / th;
      std::memset(&R[0], 0, 9U * sizeof(double));
      R[0] = 1.0;
      R[4] = 1.0;
      R[8] = 1.0;
      for (int i{0}; i < 9; i++) {
        R[i] = (R[i] + absxk * w[i]) + a_tmp * wSquare[i];
      }
    }
  } else {
    c = (s - 1.0) / minVecSim3[6];
    if (th < 1.0E-5) {
      scale = minVecSim3[6] * minVecSim3[6];
      a_tmp = ((minVecSim3[6] - 1.0) * s + 1.0) / scale;
      t = (((0.5 * scale - minVecSim3[6]) + 1.0) * s - 1.0) /
          (scale * minVecSim3[6]);
      std::memset(&R[0], 0, 9U * sizeof(double));
      R[0] = 1.0;
      R[4] = 1.0;
      R[8] = 1.0;
      for (int i{0}; i < 9; i++) {
        R[i] = (R[i] + w[i]) + wSquare[i] / 2.0;
      }
    } else {
      double a1_tmp;
      double b1_tmp;
      double thSquare;
      a1_tmp = std::sin(th);
      scale = s * a1_tmp;
      b1_tmp = std::cos(th);
      absxk = s * b1_tmp;
      thSquare = th * th;
      t = thSquare + minVecSim3[6] * minVecSim3[6];
      a_tmp = (scale * minVecSim3[6] + (1.0 - absxk) * th) / (th * t);
      t = (c - ((absxk - 1.0) * minVecSim3[6] + scale * th) / t) / thSquare;
      absxk = a1_tmp / th;
      scale = (1.0 - b1_tmp) / thSquare;
      std::memset(&R[0], 0, 9U * sizeof(double));
      R[0] = 1.0;
      R[4] = 1.0;
      R[8] = 1.0;
      for (int i{0}; i < 9; i++) {
        R[i] = (R[i] + absxk * w[i]) + scale * wSquare[i];
      }
    }
  }
  for (int i{0}; i < 9; i++) {
    w[i] = (a_tmp * w[i] + t * wSquare[i]) + c * static_cast<double>(iv[i]);
  }
  scale = minVecSim3[0];
  absxk = minVecSim3[1];
  t = minVecSim3[2];
  for (int i{0}; i < 3; i++) {
    S_tmp = i << 2;
    S[S_tmp] = R[3 * i];
    S[S_tmp + 1] = R[3 * i + 1];
    S[S_tmp + 2] = R[3 * i + 2];
    S[i + 12] = (w[i] * scale + w[i + 3] * absxk) + w[i + 6] * t;
    S[S_tmp + 3] = 0.0;
  }
  S[15] = s;
}

double
TrustRegionIndefiniteDogLegSE2::solve(const ::coder::array<double, 2U> &seed,
                                      BlockMatrix &iobj_0, BlockMatrix **xSol,
                                      sparse &hess, double &solutionInfo_Error,
                                      double &solutionInfo_ExitFlag)
{
  BlockMatrix *b_xSol;
  sparse HessianTrial;
  sparse c;
  ::coder::array<double, 2U> x;
  ::coder::array<double, 2U> xn;
  ::coder::array<double, 1U> grad;
  ::coder::array<double, 1U> gradTrial;
  ::coder::array<double, 1U> stepDL_;
  ::coder::array<double, 1U> stepGN;
  ::coder::array<double, 1U> stepSD;
  ::coder::array<double, 1U> tmpd;
  ::coder::array<boolean_T, 1U> b_x;
  double delta;
  double solutionInfo_Iterations;
  int i;
  int nx;
  boolean_T localMin;
  NLPSolverExitFlags exitFlag;
  MaxNumIterationInternal = MaxNumIteration;
  MaxTimeInternal = MaxTime;
  SeedInternal.set_size(seed.size(0), 3);
  nx = seed.size(0) * 3;
  for (i = 0; i < nx; i++) {
    SeedInternal[i] = seed[i];
  }
  TimeObj.StartTime.tv_sec = tic(TimeObj.StartTime.tv_nsec);
  x.set_size(SeedInternal.size(0), 3);
  nx = SeedInternal.size(0) * 3;
  for (i = 0; i < nx; i++) {
    x[i] = SeedInternal[i];
  }
  (&(&iobj_0)[0])[0].Matrix.set_size(SeedInternal.size(0), 3);
  nx = SeedInternal.size(0) * 3;
  for (i = 0; i < nx; i++) {
    (&(&iobj_0)[0])[0].Matrix[i] = SeedInternal[i];
  }
  (&(&iobj_0)[0])[0].BlockSize[0] = 3.0;
  (&(&iobj_0)[0])[0].BlockSize[1] = 3.0;
  (&(&iobj_0)[0])[0].NumRowBlocks =
      static_cast<double>(SeedInternal.size(0)) / 3.0;
  (&(&iobj_0)[0])[0].NumColBlocks = 1.0;
  b_xSol = &(&(&iobj_0)[0])[0];
  solutionInfo_Iterations = 0.0;
  exitFlag = NLPSolverExitFlags::IterationLimitExceeded;
  solutionInfo_Error = nav::algs::internal::PoseGraphHelpers::poseGraphCost(
      SeedInternal, ExtraArgs.edgeNodePairs, ExtraArgs.edgeMeasurements,
      ExtraArgs.edgeInfoMats, ExtraArgs.tformSize, ExtraArgs.infoMatSize,
      ExtraArgs.poseDeltaLength, ExtraArgs.nodeMap, ExtraArgs.nodeDims,
      ExtraArgs.IsLandmarkNode, grad, hess);
  delta = InitialTrustRegionRadius;
  tmpd.set_size(grad.size(0));
  nx = grad.size(0);
  for (i = 0; i < nx; i++) {
    tmpd[i] = grad[i];
  }
  localMin = computeBasicSteps(tmpd, hess, stepSD, stepGN);
  if (localMin) {
    exitFlag = NLPSolverExitFlags::LocalMinimumFound;
  } else {
    double d;
    int b_i;
    boolean_T exitg1;
    boolean_T terminated;
    terminated = false;
    d = MaxNumIterationInternal;
    b_i = 0;
    exitg1 = false;
    while ((!exitg1) && (b_i <= static_cast<int>(d) - 1)) {
      double val;
      val = toc(TimeObj.StartTime.tv_sec, TimeObj.StartTime.tv_nsec);
      if (val > MaxTimeInternal) {
        exitFlag = NLPSolverExitFlags::TimeLimitExceeded;
        terminated = true;
        exitg1 = true;
      } else {
        double b_stepDL_;
        double bc;
        int currRowIdx;
        boolean_T exitg2;
        val = b_norm(stepSD);
        if (val >= delta) {
          val = delta / val;
          stepDL_.set_size(stepSD.size(0));
          nx = stepSD.size(0);
          for (i = 0; i < nx; i++) {
            stepDL_[i] = val * stepSD[i];
          }
        } else if (b_norm(stepGN) <= delta) {
          stepDL_.set_size(stepGN.size(0));
          nx = stepGN.size(0);
          for (i = 0; i < nx; i++) {
            stepDL_[i] = stepGN[i];
          }
        } else {
          if (stepGN.size(0) == stepSD.size(0)) {
            stepDL_.set_size(stepGN.size(0));
            nx = stepGN.size(0);
            for (i = 0; i < nx; i++) {
              stepDL_[i] = stepGN[i] - stepSD[i];
            }
          } else {
            minus(stepDL_, stepGN, stepSD);
          }
          bc = 0.0;
          nx = stepSD.size(0);
          for (i = 0; i < nx; i++) {
            bc += stepSD[i] * stepDL_[i];
          }
          b_stepDL_ = 0.0;
          nx = stepDL_.size(0);
          for (i = 0; i < nx; i++) {
            b_stepDL_ += stepDL_[i] * stepDL_[i];
          }
          val = (-bc +
                 std::sqrt(bc * bc + b_stepDL_ * (delta * delta - val * val))) /
                b_stepDL_;
          if (stepSD.size(0) == stepDL_.size(0)) {
            stepDL_.set_size(stepSD.size(0));
            nx = stepSD.size(0);
            for (i = 0; i < nx; i++) {
              stepDL_[i] = stepSD[i] + val * stepDL_[i];
            }
          } else {
            binary_expand_op(stepDL_, stepSD, val);
          }
        }
        nx = stepDL_.size(0);
        tmpd.set_size(stepDL_.size(0));
        for (currRowIdx = 0; currRowIdx < nx; currRowIdx++) {
          tmpd[currRowIdx] = std::abs(stepDL_[currRowIdx]);
        }
        b_x.set_size(tmpd.size(0));
        nx = tmpd.size(0);
        for (i = 0; i < nx; i++) {
          b_x[i] = (tmpd[i] < StepTolerance);
        }
        localMin = true;
        nx = 1;
        exitg2 = false;
        while ((!exitg2) && (nx <= b_x.size(0))) {
          if (!b_x[nx - 1]) {
            localMin = false;
            exitg2 = true;
          } else {
            nx++;
          }
        }
        if (localMin) {
          exitFlag = NLPSolverExitFlags::StepSizeBelowMinimum;
          terminated = true;
          exitg1 = true;
        } else {
          double costTrial;
          double d1;
          incrementX(x, stepDL_, xn);
          costTrial = nav::algs::internal::PoseGraphHelpers::poseGraphCost(
              xn, ExtraArgs.edgeNodePairs, ExtraArgs.edgeMeasurements,
              ExtraArgs.edgeInfoMats, ExtraArgs.tformSize,
              ExtraArgs.infoMatSize, ExtraArgs.poseDeltaLength,
              ExtraArgs.nodeMap, ExtraArgs.nodeDims, ExtraArgs.IsLandmarkNode,
              gradTrial, HessianTrial);
          d1 = solutionInfo_Error - costTrial;
          if (std::abs(d1) < FunctionTolerance) {
            exitFlag = NLPSolverExitFlags::ChangeInErrorBelowMinimum;
            terminated = true;
            exitg1 = true;
          } else {
            int b_c;
            int ridx;
            boolean_T guard1{false};
            ridx = hess.colidx[hess.colidx.size(0) - 1];
            if (ridx - 1 < 1) {
              nx = 0;
            } else {
              nx = ridx - 1;
            }
            tmpd.set_size(nx);
            for (i = 0; i < nx; i++) {
              tmpd[i] = 0.5 * hess.d[i];
            }
            if (ridx - 1 >= 1) {
              nx = ridx - 2;
            } else {
              nx = 0;
            }
            c.d.set_size(nx + 1);
            c.rowidx.set_size(nx + 1);
            for (i = 0; i <= nx; i++) {
              c.d[i] = 0.0;
              c.rowidx[i] = 0;
            }
            if (ridx - 1 < 1) {
              nx = 1;
            } else {
              nx = ridx;
            }
            for (i = 0; i <= nx - 2; i++) {
              c.rowidx[i] = hess.rowidx[i];
            }
            c.colidx.set_size(hess.colidx.size(0));
            nx = hess.colidx.size(0);
            for (i = 0; i < nx; i++) {
              c.colidx[i] = hess.colidx[i];
            }
            for (currRowIdx = 0; currRowIdx <= ridx - 2; currRowIdx++) {
              c.d[currRowIdx] = tmpd[currRowIdx];
            }
            nx = 1;
            i = hess.colidx.size(0);
            for (b_c = 0; b_c <= i - 2; b_c++) {
              ridx = c.colidx[b_c];
              c.colidx[b_c] = nx;
              while (ridx < c.colidx[b_c + 1]) {
                currRowIdx = c.rowidx[ridx - 1];
                val = c.d[ridx - 1];
                ridx++;
                if (val != 0.0) {
                  c.d[nx - 1] = val;
                  c.rowidx[nx - 1] = currRowIdx;
                  nx++;
                }
              }
            }
            c.colidx[c.colidx.size(0) - 1] = nx;
            tmpd.set_size(hess.m);
            nx = hess.m;
            for (i = 0; i < nx; i++) {
              tmpd[i] = 0.0;
            }
            if ((hess.n != 0) && (hess.m != 0) &&
                (c.colidx[c.colidx.size(0) - 1] - 1 != 0)) {
              i = hess.n;
              for (int acol{0}; acol < i; acol++) {
                bc = stepDL_[acol];
                b_c = c.colidx[acol];
                ridx = c.colidx[acol + 1];
                nx = ridx - c.colidx[acol];
                if (nx >= 4) {
                  int tmpd_tmp;
                  currRowIdx = (ridx - nx) + ((nx / 4) << 2);
                  for (int ap{b_c}; ap <= currRowIdx - 1; ap += 4) {
                    tmpd_tmp = c.rowidx[ap - 1] - 1;
                    tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap - 1] * bc;
                    tmpd[c.rowidx[ap] - 1] =
                        tmpd[c.rowidx[ap] - 1] + c.d[ap] * bc;
                    tmpd_tmp = c.rowidx[ap + 1] - 1;
                    tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap + 1] * bc;
                    tmpd_tmp = c.rowidx[ap + 2] - 1;
                    tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap + 2] * bc;
                  }
                  nx = ridx - 1;
                  for (int ap{currRowIdx}; ap <= nx; ap++) {
                    tmpd_tmp = c.rowidx[ap - 1] - 1;
                    tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap - 1] * bc;
                  }
                } else {
                  nx = ridx - 1;
                  for (int ap{b_c}; ap <= nx; ap++) {
                    int tmpd_tmp;
                    tmpd_tmp = c.rowidx[ap - 1] - 1;
                    tmpd[tmpd_tmp] = tmpd[tmpd_tmp] + c.d[ap - 1] * bc;
                  }
                }
              }
            }
            if (grad.size(0) == tmpd.size(0)) {
              b_stepDL_ = 0.0;
              nx = stepDL_.size(0);
              for (i = 0; i < nx; i++) {
                b_stepDL_ += -stepDL_[i] * (grad[i] + tmpd[i]);
              }
              val = d1 / b_stepDL_;
            } else {
              val = binary_expand_op(solutionInfo_Error, costTrial, stepDL_,
                                     grad, tmpd);
            }
            guard1 = false;
            if (val > 0.0) {
              xn.set_size(x.size(0), 3);
              nx = x.size(0) * x.size(1) - 1;
              for (i = 0; i <= nx; i++) {
                xn[i] = x[i];
              }
              incrementX(xn, stepDL_, x);
              solutionInfo_Iterations = static_cast<double>(b_i) + 1.0;
              solutionInfo_Error = costTrial;
              grad.set_size(gradTrial.size(0));
              nx = gradTrial.size(0);
              for (i = 0; i < nx; i++) {
                grad[i] = gradTrial[i];
              }
              hess = HessianTrial;
              localMin =
                  computeBasicSteps(gradTrial, HessianTrial, stepSD, stepGN);
              if (localMin) {
                exitFlag = NLPSolverExitFlags::LocalMinimumFound;
                terminated = true;
                exitg1 = true;
              } else {
                guard1 = true;
              }
            } else {
              guard1 = true;
            }
            if (guard1) {
              boolean_T b_guard1{false};
              localMin = false;
              b_guard1 = false;
              if (val > 0.75) {
                d1 = b_norm(stepDL_);
                if (d1 > 0.9 * delta) {
                  delta = 2.0 * d1;
                } else {
                  b_guard1 = true;
                }
              } else {
                b_guard1 = true;
              }
              if (b_guard1 && (val < 0.25)) {
                delta /= 4.0;
                if (delta < TrustRegionRadiusTolerance) {
                  localMin = true;
                }
              }
              if (localMin) {
                exitFlag = NLPSolverExitFlags::TrustRegionRadiusBelowMinimum;
                terminated = true;
                exitg1 = true;
              } else {
                b_i++;
              }
            }
          }
        }
      }
    }
    b_xSol = &(&(&iobj_0)[0])[1];
    (&(&iobj_0)[0])[1].Matrix.set_size(x.size(0), 3);
    nx = x.size(0) * 3;
    for (i = 0; i < nx; i++) {
      (&(&iobj_0)[0])[1].Matrix[i] = x[i];
    }
    (&(&iobj_0)[0])[1].BlockSize[0] = 3.0;
    (&(&iobj_0)[0])[1].BlockSize[1] = 3.0;
    (&(&iobj_0)[0])[1].NumRowBlocks = static_cast<double>(x.size(0)) / 3.0;
    (&(&iobj_0)[0])[1].NumColBlocks = 1.0;
    if (!terminated) {
      solutionInfo_Iterations = d;
    }
  }
  *xSol = b_xSol;
  solutionInfo_ExitFlag = static_cast<double>(exitFlag);
  return solutionInfo_Iterations;
}

void SEHelpers::veelogmSE3(const double T[16], double vec[6])
{
  creal_T u;
  creal_T v;
  double Vinv[9];
  double b_I[9];
  double wx[9];
  double wv[3];
  double a;
  double theta;
  double thetaSq;
  int wx_tmp;
  boolean_T guard1{false};
  boolean_T y;
  a = 0.5 * (((T[0] + T[5]) + T[10]) - 1.0);
  if (!(std::abs(a) > 1.0)) {
    u.re = std::acos(a);
  } else {
    v.re = a + 1.0;
    v.im = 0.0;
    ::poseGraphOptimize::coder::internal::scalar::b_sqrt(v);
    u.re = 1.0 - a;
    u.im = 0.0;
    ::poseGraphOptimize::coder::internal::scalar::b_sqrt(u);
    a = u.re;
    u.re = 2.0 * rt_atan2d_snf(a, v.re);
  }
  a = u.re / std::sin(u.re);
  for (int i{0}; i < 3; i++) {
    wx_tmp = i << 2;
    wx[3 * i] = T[wx_tmp] - T[i];
    wx[3 * i + 1] = T[wx_tmp + 1] - T[i + 4];
    wx[3 * i + 2] = T[wx_tmp + 2] - T[i + 8];
  }
  wv[0] = wx[5];
  wv[1] = wx[6];
  wv[2] = wx[1];
  guard1 = false;
  if ((!std::isinf(a)) && (!std::isnan(a))) {
    boolean_T exitg1;
    y = true;
    wx_tmp = 0;
    exitg1 = false;
    while ((!exitg1) && (wx_tmp < 3)) {
      if (!(wv[wx_tmp] == 0.0)) {
        y = false;
        exitg1 = true;
      } else {
        wx_tmp++;
      }
    }
    if (!y) {
      wv[0] = wx[5] * a / 2.0;
      wv[1] = wx[6] * a / 2.0;
      wv[2] = wx[1] * a / 2.0;
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }
  if (guard1) {
    std::memset(&b_I[0], 0, 9U * sizeof(double));
    b_I[0] = 1.0;
    b_I[4] = 1.0;
    b_I[8] = 1.0;
    for (int i{0}; i < 3; i++) {
      int I_tmp;
      wx_tmp = i << 2;
      b_I[3 * i] -= T[wx_tmp];
      I_tmp = 3 * i + 1;
      b_I[I_tmp] -= T[wx_tmp + 1];
      I_tmp = 3 * i + 2;
      b_I[I_tmp] -= T[wx_tmp + 2];
    }
    y = true;
    for (wx_tmp = 0; wx_tmp < 9; wx_tmp++) {
      if (y) {
        a = b_I[wx_tmp];
        if (std::isinf(a) || std::isnan(a)) {
          y = false;
        }
      } else {
        y = false;
      }
    }
    if (y) {
      ::poseGraphOptimize::coder::internal::svd(b_I, Vinv, wv, wx);
    } else {
      for (int i{0}; i < 9; i++) {
        wx[i] = rtNaN;
      }
    }
    a = 1.0 / std::sqrt((wx[6] * wx[6] + wx[7] * wx[7]) + wx[8] * wx[8]);
    wv[0] = wx[6] * a * u.re;
    wv[1] = wx[7] * a * u.re;
    wv[2] = wx[8] * a * u.re;
  }
  theta = std::sqrt((wv[0] * wv[0] + wv[1] * wv[1]) + wv[2] * wv[2]);
  thetaSq = theta * theta;
  a = (1.0 - std::cos(theta)) / thetaSq;
  if ((a < 2.2204460492503131E-16) || (std::isinf(a) || std::isnan(a))) {
    std::memset(&Vinv[0], 0, 9U * sizeof(double));
    Vinv[0] = 1.0;
    Vinv[4] = 1.0;
    Vinv[8] = 1.0;
  } else {
    wx[0] = 0.0;
    wx[3] = -wv[2];
    wx[6] = wv[1];
    wx[1] = wv[2];
    wx[4] = 0.0;
    wx[7] = -wv[0];
    wx[2] = -wv[1];
    wx[5] = wv[0];
    wx[8] = 0.0;
    a = 1.0 / thetaSq * (1.0 - std::sin(theta) / theta / (2.0 * a));
    std::memset(&b_I[0], 0, 9U * sizeof(double));
    for (wx_tmp = 0; wx_tmp < 3; wx_tmp++) {
      b_I[wx_tmp + 3 * wx_tmp] = 1.0;
      for (int i{0}; i < 3; i++) {
        Vinv[wx_tmp + 3 * i] =
            (wx[wx_tmp] * wx[3 * i] + wx[wx_tmp + 3] * wx[3 * i + 1]) +
            wx[wx_tmp + 6] * wx[3 * i + 2];
      }
    }
    for (int i{0}; i < 9; i++) {
      Vinv[i] = (b_I[i] - 0.5 * wx[i]) + a * Vinv[i];
    }
  }
  a = T[12];
  theta = T[13];
  thetaSq = T[14];
  for (int i{0}; i < 3; i++) {
    vec[i] = (Vinv[i] * a + Vinv[i + 3] * theta) + Vinv[i + 6] * thetaSq;
    vec[i + 3] = wv[i];
  }
}

} // namespace internal
} // namespace core
} // namespace robotics
} // namespace coder
static void binary_expand_op(double in1[3],
                             const ::coder::array<double, 2U> &in2,
                             const ::coder::array<double, 1U> &in3, int in4,
                             int in5)
{
  int stride_0_0;
  in1[0] = in2[in2.size(0) * 2];
  in1[1] = in2[in2.size(0) * 2 + 1];
  in1[2] = rt_atan2d_snf(in2[1], in2[0]);
  stride_0_0 = ((in5 - in4) + 1 != 1);
  in1[0] += in3[in4];
  in1[1] += in3[in4 + stride_0_0];
  in1[2] += in3[in4 + (stride_0_0 << 1)];
}

static void binary_expand_op(::coder::array<double, 2U> &in1,
                             const ::coder::array<double, 1U> &in2, int in3,
                             int in4)
{
  in1[in1.size(0) * 2] = in1[in1.size(0) * 2] + in2[in3];
  in1[in1.size(0) * 2 + 1] =
      in1[in1.size(0) * 2 + 1] + in2[in3 + ((in4 - in3) + 1 != 1)];
}

static void binary_expand_op(coder::nav::algs::internal::BlockInserter2 *in1,
                             int in2, int in4, int in5, const double in6_data[],
                             const int &in6_size)
{
  ::coder::array<double, 1U> b_in1;
  int stride_0_0;
  b_in1.set_size(in6_size);
  stride_0_0 = ((in5 - in4) + 1 != 1);
  for (int i{0}; i < in6_size; i++) {
    b_in1[i] = in1->Gradient[in4 + i * stride_0_0] + in6_data[i];
  }
  stride_0_0 = b_in1.size(0);
  for (int i{0}; i < stride_0_0; i++) {
    in1->Gradient[in2 + i] = b_in1[i];
  }
}

static double binary_expand_op(double in1, double in2,
                               const ::coder::array<double, 1U> &in3,
                               const ::coder::array<double, 1U> &in4,
                               const ::coder::array<double, 1U> &in5)
{
  ::coder::array<double, 1U> b_in4;
  double b_in3;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in5.size(0) == 1) {
    loop_ub = in4.size(0);
  } else {
    loop_ub = in5.size(0);
  }
  b_in4.set_size(loop_ub);
  stride_0_0 = (in4.size(0) != 1);
  stride_1_0 = (in5.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    b_in4[i] = in4[i * stride_0_0] + in5[i * stride_1_0];
  }
  b_in3 = 0.0;
  loop_ub = in3.size(0);
  for (int i{0}; i < loop_ub; i++) {
    b_in3 += -in3[i] * b_in4[i];
  }
  return (in1 - in2) / b_in3;
}

static void binary_expand_op(::coder::array<double, 1U> &in1,
                             const ::coder::array<double, 1U> &in2, double in3)
{
  ::coder::array<double, 1U> b_in2;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in1.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in1.size(0);
  }
  b_in2.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in1.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    b_in2[i] = in2[i * stride_0_0] + in3 * in1[i * stride_1_0];
  }
  in1.set_size(b_in2.size(0));
  loop_ub = b_in2.size(0);
  for (int i{0}; i < loop_ub; i++) {
    in1[i] = b_in2[i];
  }
}

namespace coder {
static double b_norm(const ::coder::array<double, 1U> &x)
{
  double y;
  if (x.size(0) == 0) {
    y = 0.0;
  } else {
    y = 0.0;
    if (x.size(0) == 1) {
      y = std::abs(x[0]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = x.size(0);
      for (int k{0}; k < kend; k++) {
        double absxk;
        absxk = std::abs(x[k]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }
      y = scale * std::sqrt(y);
    }
  }
  return y;
}

static void b_sparse(const ::coder::array<double, 1U> &varargin_1,
                     const ::coder::array<double, 1U> &varargin_2,
                     const ::coder::array<double, 1U> &varargin_3, sparse &y)
{
  anonymous_function b_this;
  ::coder::array<int, 1U> sortedIndices;
  ::coder::array<int, 1U> t;
  int i;
  int nc;
  int ns;
  int ny;
  nc = varargin_2.size(0);
  ns = varargin_1.size(0);
  b_this.workspace.b.set_size(varargin_1.size(0));
  for (int k{0}; k < ns; k++) {
    b_this.workspace.b[k] = static_cast<int>(varargin_1[k]);
  }
  ns = varargin_2.size(0);
  b_this.workspace.a.set_size(varargin_2.size(0));
  sortedIndices.set_size(varargin_2.size(0));
  for (int k{0}; k < ns; k++) {
    b_this.workspace.a[k] = static_cast<int>(varargin_2[k]);
    sortedIndices[k] = k + 1;
  }
  internal::introsort(sortedIndices, b_this.workspace.a.size(0), b_this);
  ny = b_this.workspace.a.size(0);
  t.set_size(b_this.workspace.a.size(0));
  ns = b_this.workspace.a.size(0);
  for (i = 0; i < ns; i++) {
    t[i] = b_this.workspace.a[i];
  }
  for (int k{0}; k < ny; k++) {
    b_this.workspace.a[k] = t[sortedIndices[k] - 1];
  }
  ny = b_this.workspace.b.size(0);
  t.set_size(b_this.workspace.b.size(0));
  ns = b_this.workspace.b.size(0);
  for (i = 0; i < ns; i++) {
    t[i] = b_this.workspace.b[i];
  }
  for (int k{0}; k < ny; k++) {
    b_this.workspace.b[k] = t[sortedIndices[k] - 1];
  }
  if ((b_this.workspace.b.size(0) == 0) || (b_this.workspace.a.size(0) == 0)) {
    ny = 0;
    y.n = 0;
  } else {
    ns = b_this.workspace.b.size(0);
    ny = b_this.workspace.b[0];
    for (int k{2}; k <= ns; k++) {
      i = b_this.workspace.b[k - 1];
      if (ny < i) {
        ny = i;
      }
    }
    y.n = b_this.workspace.a[b_this.workspace.a.size(0) - 1];
  }
  y.m = ny;
  ns = varargin_2.size(0);
  if (ns < 1) {
    ns = 1;
  }
  y.d.set_size(ns);
  y.colidx.set_size(y.n + 1);
  y.colidx[0] = 1;
  y.rowidx.set_size(ns);
  for (i = 0; i < ns; i++) {
    y.d[i] = 0.0;
    y.rowidx[i] = 0;
  }
  ns = 0;
  i = y.n;
  for (ny = 0; ny < i; ny++) {
    while ((ns + 1 <= nc) && (b_this.workspace.a[ns] == ny + 1)) {
      y.rowidx[ns] = b_this.workspace.b[ns];
      ns++;
    }
    y.colidx[ny + 1] = ns + 1;
  }
  for (int k{0}; k < nc; k++) {
    y.d[k] = varargin_3[sortedIndices[k] - 1];
  }
  y.fillIn();
}

namespace internal {
static void b_heapsort(::coder::array<int, 1U> &x, int xstart, int xend,
                       const anonymous_function &cmp)
{
  int idx;
  int n;
  n = (xend - xstart) - 1;
  for (idx = n + 2; idx >= 1; idx--) {
    heapify(x, idx, xstart, xend, cmp);
  }
  for (int k{0}; k <= n; k++) {
    int t;
    idx = (xend - k) - 1;
    t = x[idx];
    x[idx] = x[xstart - 1];
    x[xstart - 1] = t;
    heapify(x, 1, xstart, idx, cmp);
  }
}

namespace blas {
static void mtimes(const double A_data[], const int A_size[2],
                   const double B_data[], const int B_size[2], double C_data[],
                   int C_size[2])
{
  int inner;
  int mc;
  int nc;
  mc = A_size[0];
  inner = A_size[1];
  nc = B_size[1];
  C_size[0] = A_size[0];
  C_size[1] = B_size[1];
  for (int j{0}; j < nc; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * B_size[0];
    std::memset(&C_data[coffset], 0,
                static_cast<unsigned int>((mc + coffset) - coffset) *
                    sizeof(double));
    for (int k{0}; k < inner; k++) {
      double bkj;
      int aoffset;
      aoffset = k * A_size[0];
      bkj = B_data[boffset + k];
      for (int i{0}; i < mc; i++) {
        int b_i;
        b_i = coffset + i;
        C_data[b_i] += A_data[aoffset + i] * bkj;
      }
    }
  }
}

static void mtimes(const double A_data[], const int A_size[2],
                   const ::coder::array<double, 2U> &B, double C_data[],
                   int C_size[2])
{
  int inner;
  int mc;
  int nc;
  mc = A_size[1];
  inner = A_size[0];
  nc = B.size(1);
  C_size[0] = A_size[1];
  C_size[1] = B.size(1);
  for (int j{0}; j < nc; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * B.size(0);
    std::memset(&C_data[coffset], 0,
                static_cast<unsigned int>((mc + coffset) - coffset) *
                    sizeof(double));
    for (int k{0}; k < inner; k++) {
      double bkj;
      bkj = B[boffset + k];
      for (int i{0}; i < mc; i++) {
        int b_i;
        b_i = coffset + i;
        C_data[b_i] += A_data[i * A_size[0] + k] * bkj;
      }
    }
  }
}

static void xaxpy(double a, const double x[9], int ix0, double y[3])
{
  if (!(a == 0.0)) {
    for (int k{0}; k < 2; k++) {
      y[k + 1] += a * x[(ix0 + k) - 1];
    }
  }
}

static void xaxpy(double a, const double x[3], double y[9], int iy0)
{
  if (!(a == 0.0)) {
    for (int k{0}; k < 2; k++) {
      int i;
      i = (iy0 + k) - 1;
      y[i] += a * x[k + 1];
    }
  }
}

static void xaxpy(int n, double a, int ix0, double y[9], int iy0)
{
  if (!(a == 0.0)) {
    int i;
    i = n - 1;
    for (int k{0}; k <= i; k++) {
      int i1;
      i1 = (iy0 + k) - 1;
      y[i1] += a * y[(ix0 + k) - 1];
    }
  }
}

static double xdotc(int n, const double x[9], int ix0, const double y[9],
                    int iy0)
{
  double d;
  int i;
  d = 0.0;
  i = static_cast<unsigned char>(n);
  for (int k{0}; k < i; k++) {
    d += x[(ix0 + k) - 1] * y[(iy0 + k) - 1];
  }
  return d;
}

static double xnrm2(const double x[3])
{
  double scale;
  double y;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  for (int k{2}; k < 4; k++) {
    double absxk;
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      double t;
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      double t;
      t = absxk / scale;
      y += t * t;
    }
  }
  return scale * std::sqrt(y);
}

static double xnrm2(int n, const double x[9], int ix0)
{
  double scale;
  double y;
  int kend;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = (ix0 + n) - 1;
  for (int k{ix0}; k <= kend; k++) {
    double absxk;
    absxk = std::abs(x[k - 1]);
    if (absxk > scale) {
      double t;
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      double t;
      t = absxk / scale;
      y += t * t;
    }
  }
  return scale * std::sqrt(y);
}

static void xrot(double x[9], int ix0, int iy0, double c, double s)
{
  double temp;
  double temp_tmp;
  temp = x[iy0 - 1];
  temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = c * temp - s * temp_tmp;
  x[ix0 - 1] = c * temp_tmp + s * temp;
  temp = c * x[ix0] + s * x[iy0];
  x[iy0] = c * x[iy0] - s * x[ix0];
  x[ix0] = temp;
  temp = x[iy0 + 1];
  temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = c * temp - s * temp_tmp;
  x[ix0 + 1] = c * temp_tmp + s * temp;
}

static double xrotg(double &a, double &b, double &s)
{
  double absa;
  double absb;
  double c;
  double roe;
  double scale;
  roe = b;
  absa = std::abs(a);
  absb = std::abs(b);
  if (absa > absb) {
    roe = a;
  }
  scale = absa + absb;
  if (scale == 0.0) {
    s = 0.0;
    c = 1.0;
    a = 0.0;
    b = 0.0;
  } else {
    double ads;
    double bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= std::sqrt(ads * ads + bds * bds);
    if (roe < 0.0) {
      scale = -scale;
    }
    c = a / scale;
    s = b / scale;
    if (absa > absb) {
      b = s;
    } else if (c != 0.0) {
      b = 1.0 / c;
    } else {
      b = 1.0;
    }
    a = scale;
  }
  return c;
}

static void xswap(double x[9], int ix0, int iy0)
{
  double temp;
  temp = x[ix0 - 1];
  x[ix0 - 1] = x[iy0 - 1];
  x[iy0 - 1] = temp;
  temp = x[ix0];
  x[ix0] = x[iy0];
  x[iy0] = temp;
  temp = x[ix0 + 1];
  x[ix0 + 1] = x[iy0 + 1];
  x[iy0 + 1] = temp;
}

} // namespace blas
static void heapify(::coder::array<int, 1U> &x, int idx, int xstart, int xend,
                    const anonymous_function &cmp)
{
  int extremum;
  int extremumIdx;
  int i;
  int i1;
  int leftIdx;
  boolean_T changed;
  boolean_T exitg1;
  boolean_T varargout_1;
  changed = true;
  extremumIdx = (idx + xstart) - 2;
  leftIdx = ((idx << 1) + xstart) - 2;
  exitg1 = false;
  while ((!exitg1) && (leftIdx + 1 < xend)) {
    int cmpIdx;
    int i2;
    int xcmp;
    changed = false;
    extremum = x[extremumIdx];
    cmpIdx = leftIdx;
    xcmp = x[leftIdx];
    i = cmp.workspace.a[x[leftIdx] - 1];
    i1 = x[leftIdx + 1] - 1;
    i2 = cmp.workspace.a[i1];
    if (i < i2) {
      varargout_1 = true;
    } else if (i == i2) {
      varargout_1 = (cmp.workspace.b[x[leftIdx] - 1] < cmp.workspace.b[i1]);
    } else {
      varargout_1 = false;
    }
    if (varargout_1) {
      cmpIdx = leftIdx + 1;
      xcmp = x[leftIdx + 1];
    }
    i = cmp.workspace.a[x[extremumIdx] - 1];
    i1 = cmp.workspace.a[xcmp - 1];
    if (i < i1) {
      varargout_1 = true;
    } else if (i == i1) {
      varargout_1 =
          (cmp.workspace.b[x[extremumIdx] - 1] < cmp.workspace.b[xcmp - 1]);
    } else {
      varargout_1 = false;
    }
    if (varargout_1) {
      x[extremumIdx] = xcmp;
      x[cmpIdx] = extremum;
      extremumIdx = cmpIdx;
      leftIdx = ((((cmpIdx - xstart) + 2) << 1) + xstart) - 2;
      changed = true;
    } else {
      exitg1 = true;
    }
  }
  if (changed && (leftIdx + 1 <= xend)) {
    extremum = x[extremumIdx];
    i = cmp.workspace.a[x[extremumIdx] - 1];
    i1 = cmp.workspace.a[x[leftIdx] - 1];
    if (i < i1) {
      varargout_1 = true;
    } else if (i == i1) {
      varargout_1 = (cmp.workspace.b[x[extremumIdx] - 1] <
                     cmp.workspace.b[x[leftIdx] - 1]);
    } else {
      varargout_1 = false;
    }
    if (varargout_1) {
      x[extremumIdx] = x[leftIdx];
      x[leftIdx] = extremum;
    }
  }
}

static void insertionsort(::coder::array<int, 1U> &x, int xstart, int xend,
                          const anonymous_function &cmp)
{
  int i;
  i = xstart + 1;
  for (int k{i}; k <= xend; k++) {
    int idx;
    int xc;
    boolean_T exitg1;
    xc = x[k - 1] - 1;
    idx = k - 2;
    exitg1 = false;
    while ((!exitg1) && (idx + 1 >= xstart)) {
      int i1;
      boolean_T varargout_1;
      i1 = cmp.workspace.a[x[idx] - 1];
      if (cmp.workspace.a[xc] < i1) {
        varargout_1 = true;
      } else if (cmp.workspace.a[xc] == i1) {
        varargout_1 = (cmp.workspace.b[xc] < cmp.workspace.b[x[idx] - 1]);
      } else {
        varargout_1 = false;
      }
      if (varargout_1) {
        x[idx + 1] = x[idx];
        idx--;
      } else {
        exitg1 = true;
      }
    }
    x[idx + 1] = xc + 1;
  }
}

static void introsort(::coder::array<int, 1U> &x, int xend,
                      const anonymous_function &cmp)
{
  struct_T frame;
  if (xend > 1) {
    if (xend <= 32) {
      insertionsort(x, 1, xend, cmp);
    } else {
      stack st;
      int MAXDEPTH;
      int i;
      int pmax;
      int pmin;
      int pow2p;
      int t;
      boolean_T exitg1;
      pmax = 31;
      pmin = 0;
      exitg1 = false;
      while ((!exitg1) && (pmax - pmin > 1)) {
        t = (pmin + pmax) >> 1;
        pow2p = 1 << t;
        if (pow2p == xend) {
          pmax = t;
          exitg1 = true;
        } else if (pow2p > xend) {
          pmax = t;
        } else {
          pmin = t;
        }
      }
      MAXDEPTH = (pmax - 1) << 1;
      frame.xstart = 1;
      frame.xend = xend;
      frame.depth = 0;
      pmax = MAXDEPTH << 1;
      st.d.size[0] = pmax;
      for (i = 0; i < pmax; i++) {
        st.d.data[i] = frame;
      }
      st.d.data[0] = frame;
      st.n = 1;
      while (st.n > 0) {
        int frame_tmp_tmp;
        frame_tmp_tmp = st.n - 1;
        frame = st.d.data[st.n - 1];
        st.n--;
        i = frame.xend - frame.xstart;
        if (i + 1 <= 32) {
          insertionsort(x, frame.xstart, frame.xend, cmp);
        } else if (frame.depth == MAXDEPTH) {
          b_heapsort(x, frame.xstart, frame.xend, cmp);
        } else {
          int xmid;
          boolean_T varargout_1;
          xmid = (frame.xstart + i / 2) - 1;
          i = cmp.workspace.a[x[xmid] - 1];
          pmax = x[frame.xstart - 1];
          pmin = cmp.workspace.a[pmax - 1];
          if (i < pmin) {
            varargout_1 = true;
          } else if (i == pmin) {
            varargout_1 =
                (cmp.workspace.b[x[xmid] - 1] < cmp.workspace.b[pmax - 1]);
          } else {
            varargout_1 = false;
          }
          if (varargout_1) {
            x[frame.xstart - 1] = x[xmid];
            x[xmid] = pmax;
          }
          i = x[frame.xend - 1];
          pmax = cmp.workspace.a[i - 1];
          pmin = x[frame.xstart - 1];
          t = cmp.workspace.a[pmin - 1];
          if (pmax < t) {
            varargout_1 = true;
          } else if (pmax == t) {
            varargout_1 = (cmp.workspace.b[i - 1] < cmp.workspace.b[pmin - 1]);
          } else {
            varargout_1 = false;
          }
          if (varargout_1) {
            x[frame.xstart - 1] = i;
            x[frame.xend - 1] = pmin;
          }
          i = x[frame.xend - 1];
          pmax = cmp.workspace.a[i - 1];
          pmin = cmp.workspace.a[x[xmid] - 1];
          if (pmax < pmin) {
            varargout_1 = true;
          } else if (pmax == pmin) {
            varargout_1 =
                (cmp.workspace.b[i - 1] < cmp.workspace.b[x[xmid] - 1]);
          } else {
            varargout_1 = false;
          }
          if (varargout_1) {
            t = x[xmid];
            x[xmid] = i;
            x[frame.xend - 1] = t;
          }
          pow2p = x[xmid] - 1;
          x[xmid] = x[frame.xend - 2];
          x[frame.xend - 2] = pow2p + 1;
          pmax = frame.xstart - 1;
          pmin = frame.xend - 2;
          int exitg2;
          do {
            int exitg3;
            exitg2 = 0;
            pmax++;
            do {
              exitg3 = 0;
              i = cmp.workspace.a[x[pmax] - 1];
              if (i < cmp.workspace.a[pow2p]) {
                varargout_1 = true;
              } else if (i == cmp.workspace.a[pow2p]) {
                varargout_1 =
                    (cmp.workspace.b[x[pmax] - 1] < cmp.workspace.b[pow2p]);
              } else {
                varargout_1 = false;
              }
              if (varargout_1) {
                pmax++;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);
            pmin--;
            do {
              exitg3 = 0;
              i = cmp.workspace.a[x[pmin] - 1];
              if (cmp.workspace.a[pow2p] < i) {
                varargout_1 = true;
              } else if (cmp.workspace.a[pow2p] == i) {
                varargout_1 =
                    (cmp.workspace.b[pow2p] < cmp.workspace.b[x[pmin] - 1]);
              } else {
                varargout_1 = false;
              }
              if (varargout_1) {
                pmin--;
              } else {
                exitg3 = 1;
              }
            } while (exitg3 == 0);
            if (pmax + 1 >= pmin + 1) {
              exitg2 = 1;
            } else {
              t = x[pmax];
              x[pmax] = x[pmin];
              x[pmin] = t;
            }
          } while (exitg2 == 0);
          x[frame.xend - 2] = x[pmax];
          x[pmax] = pow2p + 1;
          if (pmax + 2 < frame.xend) {
            st.d.data[frame_tmp_tmp].xstart = pmax + 2;
            st.d.data[frame_tmp_tmp].xend = frame.xend;
            st.d.data[frame_tmp_tmp].depth = frame.depth + 1;
            st.n = frame_tmp_tmp + 1;
          }
          if (frame.xstart < pmax + 1) {
            st.d.data[st.n].xstart = frame.xstart;
            st.d.data[st.n].xend = pmax + 1;
            st.d.data[st.n].depth = frame.depth + 1;
            st.n++;
          }
        }
      }
    }
  }
}

namespace scalar {
static void b_sqrt(creal_T &x)
{
  double absxi;
  double absxr;
  double xi;
  double xr;
  xr = x.re;
  xi = x.im;
  if (xi == 0.0) {
    if (xr < 0.0) {
      absxr = 0.0;
      absxi = std::sqrt(-xr);
    } else {
      absxr = std::sqrt(xr);
      absxi = 0.0;
    }
  } else if (xr == 0.0) {
    if (xi < 0.0) {
      absxr = std::sqrt(-xi / 2.0);
      absxi = -absxr;
    } else {
      absxr = std::sqrt(xi / 2.0);
      absxi = absxr;
    }
  } else if (std::isnan(xr)) {
    absxr = rtNaN;
    absxi = rtNaN;
  } else if (std::isnan(xi)) {
    absxr = rtNaN;
    absxi = rtNaN;
  } else if (std::isinf(xi)) {
    absxr = std::abs(xi);
    absxi = xi;
  } else if (std::isinf(xr)) {
    if (xr < 0.0) {
      absxr = 0.0;
      absxi = xi * -xr;
    } else {
      absxr = xr;
      absxi = 0.0;
    }
  } else {
    absxr = std::abs(xr);
    absxi = std::abs(xi);
    if ((absxr > 4.4942328371557893E+307) ||
        (absxi > 4.4942328371557893E+307)) {
      absxr *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi * 0.5);
      if (absxi > absxr) {
        absxr = std::sqrt(absxi) * std::sqrt(absxr / absxi + 1.0);
      } else {
        absxr = std::sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxr = std::sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }
    if (xr > 0.0) {
      absxi = 0.5 * (xi / absxr);
    } else {
      if (xi < 0.0) {
        absxi = -absxr;
      } else {
        absxi = absxr;
      }
      absxr = 0.5 * (xi / absxi);
    }
  }
  x.re = absxr;
  x.im = absxi;
}

} // namespace scalar
static void svd(const double A[9], double U[9], double s[3], double V[9])
{
  double b_A[9];
  double b_s[3];
  double e[3];
  double work[3];
  double nrm;
  double rt;
  double snorm;
  double sqds;
  int ii;
  int kase;
  int m;
  int qjj;
  int qp1;
  int qq;
  int qq_tmp;
  b_s[0] = 0.0;
  e[0] = 0.0;
  work[0] = 0.0;
  b_s[1] = 0.0;
  e[1] = 0.0;
  work[1] = 0.0;
  b_s[2] = 0.0;
  e[2] = 0.0;
  work[2] = 0.0;
  for (qjj = 0; qjj < 9; qjj++) {
    b_A[qjj] = A[qjj];
    U[qjj] = 0.0;
    V[qjj] = 0.0;
  }
  for (int q{0}; q < 2; q++) {
    boolean_T apply_transform;
    qp1 = q + 2;
    qq_tmp = q + 3 * q;
    qq = qq_tmp + 1;
    apply_transform = false;
    nrm = blas::xnrm2(3 - q, b_A, qq_tmp + 1);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qq_tmp] < 0.0) {
        nrm = -nrm;
      }
      b_s[q] = nrm;
      if (std::abs(nrm) >= 1.0020841800044864E-292) {
        nrm = 1.0 / nrm;
        qjj = (qq_tmp - q) + 3;
        for (int k{qq}; k <= qjj; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        qjj = (qq_tmp - q) + 3;
        for (int k{qq}; k <= qjj; k++) {
          b_A[k - 1] /= b_s[q];
        }
      }
      b_A[qq_tmp]++;
      b_s[q] = -b_s[q];
    } else {
      b_s[q] = 0.0;
    }
    for (kase = qp1; kase < 4; kase++) {
      qjj = q + 3 * (kase - 1);
      if (apply_transform) {
        blas::xaxpy(
            3 - q,
            -(blas::xdotc(3 - q, b_A, qq_tmp + 1, b_A, qjj + 1) / b_A[qq_tmp]),
            qq_tmp + 1, b_A, qjj + 1);
      }
      e[kase - 1] = b_A[qjj];
    }
    for (ii = q + 1; ii < 4; ii++) {
      kase = (ii + 3 * q) - 1;
      U[kase] = b_A[kase];
    }
    if (q + 1 <= 1) {
      nrm = blas::xnrm2(e);
      if (nrm == 0.0) {
        e[0] = 0.0;
      } else {
        if (e[1] < 0.0) {
          e[0] = -nrm;
        } else {
          e[0] = nrm;
        }
        nrm = e[0];
        if (std::abs(e[0]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[0];
          for (int k{qp1}; k < 4; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          for (int k{qp1}; k < 4; k++) {
            e[k - 1] /= nrm;
          }
        }
        e[1]++;
        e[0] = -e[0];
        for (ii = qp1; ii < 4; ii++) {
          work[ii - 1] = 0.0;
        }
        for (kase = qp1; kase < 4; kase++) {
          blas::xaxpy(e[kase - 1], b_A, 3 * (kase - 1) + 2, work);
        }
        for (kase = qp1; kase < 4; kase++) {
          blas::xaxpy(-e[kase - 1] / e[1], work, b_A, 3 * (kase - 1) + 2);
        }
      }
      for (ii = qp1; ii < 4; ii++) {
        V[ii - 1] = e[ii - 1];
      }
    }
  }
  m = 1;
  b_s[2] = b_A[8];
  e[1] = b_A[7];
  e[2] = 0.0;
  U[6] = 0.0;
  U[7] = 0.0;
  U[8] = 1.0;
  for (int q{1}; q >= 0; q--) {
    qp1 = q + 2;
    qq = q + 3 * q;
    if (b_s[q] != 0.0) {
      for (kase = qp1; kase < 4; kase++) {
        qjj = (q + 3 * (kase - 1)) + 1;
        blas::xaxpy(3 - q, -(blas::xdotc(3 - q, U, qq + 1, U, qjj) / U[qq]),
                    qq + 1, U, qjj);
      }
      for (ii = q + 1; ii < 4; ii++) {
        kase = (ii + 3 * q) - 1;
        U[kase] = -U[kase];
      }
      U[qq]++;
      if (q - 1 >= 0) {
        U[3 * q] = 0.0;
      }
    } else {
      U[3 * q] = 0.0;
      U[3 * q + 1] = 0.0;
      U[3 * q + 2] = 0.0;
      U[qq] = 1.0;
    }
  }
  for (int q{2}; q >= 0; q--) {
    if ((q + 1 <= 1) && (e[0] != 0.0)) {
      blas::xaxpy(2, -(blas::xdotc(2, V, 2, V, 5) / V[1]), 2, V, 5);
      blas::xaxpy(2, -(blas::xdotc(2, V, 2, V, 8) / V[1]), 2, V, 8);
    }
    V[3 * q] = 0.0;
    V[3 * q + 1] = 0.0;
    V[3 * q + 2] = 0.0;
    V[q + 3 * q] = 1.0;
  }
  qq = 0;
  snorm = 0.0;
  for (int q{0}; q < 3; q++) {
    nrm = b_s[q];
    if (nrm != 0.0) {
      rt = std::abs(nrm);
      nrm /= rt;
      b_s[q] = rt;
      if (q + 1 < 3) {
        e[q] /= nrm;
      }
      kase = 3 * q;
      qjj = kase + 3;
      for (int k{kase + 1}; k <= qjj; k++) {
        U[k - 1] *= nrm;
      }
    }
    if (q + 1 < 3) {
      nrm = e[q];
      if (nrm != 0.0) {
        rt = std::abs(nrm);
        nrm = rt / nrm;
        e[q] = rt;
        b_s[q + 1] *= nrm;
        kase = 3 * (q + 1);
        qjj = kase + 3;
        for (int k{kase + 1}; k <= qjj; k++) {
          V[k - 1] *= nrm;
        }
      }
    }
    snorm = std::fmax(snorm, std::fmax(std::abs(b_s[q]), std::abs(e[q])));
  }
  while ((m + 2 > 0) && (qq < 75)) {
    boolean_T exitg1;
    qq_tmp = m + 1;
    ii = m + 1;
    exitg1 = false;
    while (!(exitg1 || (ii == 0))) {
      nrm = std::abs(e[ii - 1]);
      if ((nrm <= 2.2204460492503131E-16 *
                      (std::abs(b_s[ii - 1]) + std::abs(b_s[ii]))) ||
          (nrm <= 1.0020841800044864E-292) ||
          ((qq > 20) && (nrm <= 2.2204460492503131E-16 * snorm))) {
        e[ii - 1] = 0.0;
        exitg1 = true;
      } else {
        ii--;
      }
    }
    if (ii == m + 1) {
      kase = 4;
    } else {
      qjj = m + 2;
      kase = m + 2;
      exitg1 = false;
      while ((!exitg1) && (kase >= ii)) {
        qjj = kase;
        if (kase == ii) {
          exitg1 = true;
        } else {
          nrm = 0.0;
          if (kase < m + 2) {
            nrm = std::abs(e[kase - 1]);
          }
          if (kase > ii + 1) {
            nrm += std::abs(e[kase - 2]);
          }
          rt = std::abs(b_s[kase - 1]);
          if ((rt <= 2.2204460492503131E-16 * nrm) ||
              (rt <= 1.0020841800044864E-292)) {
            b_s[kase - 1] = 0.0;
            exitg1 = true;
          } else {
            kase--;
          }
        }
      }
      if (qjj == ii) {
        kase = 3;
      } else if (qjj == m + 2) {
        kase = 1;
      } else {
        kase = 2;
        ii = qjj;
      }
    }
    switch (kase) {
    case 1: {
      rt = e[m];
      e[m] = 0.0;
      for (int k{qq_tmp}; k >= ii + 1; k--) {
        double sm;
        sm = blas::xrotg(b_s[k - 1], rt, sqds);
        if (k > ii + 1) {
          rt = -sqds * e[0];
          e[0] *= sm;
        }
        blas::xrot(V, 3 * (k - 1) + 1, 3 * (m + 1) + 1, sm, sqds);
      }
    } break;
    case 2: {
      rt = e[ii - 1];
      e[ii - 1] = 0.0;
      for (int k{ii + 1}; k <= m + 2; k++) {
        double b;
        double sm;
        sm = blas::xrotg(b_s[k - 1], rt, sqds);
        b = e[k - 1];
        rt = -sqds * b;
        e[k - 1] = b * sm;
        blas::xrot(U, 3 * (k - 1) + 1, 3 * (ii - 1) + 1, sm, sqds);
      }
    } break;
    case 3: {
      double b;
      double scale;
      double sm;
      nrm = b_s[m + 1];
      scale = std::fmax(
          std::fmax(std::fmax(std::fmax(std::abs(nrm), std::abs(b_s[m])),
                              std::abs(e[m])),
                    std::abs(b_s[ii])),
          std::abs(e[ii]));
      sm = nrm / scale;
      nrm = b_s[m] / scale;
      rt = e[m] / scale;
      sqds = b_s[ii] / scale;
      b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0;
      nrm = sm * rt;
      nrm *= nrm;
      if ((b != 0.0) || (nrm != 0.0)) {
        rt = std::sqrt(b * b + nrm);
        if (b < 0.0) {
          rt = -rt;
        }
        rt = nrm / (b + rt);
      } else {
        rt = 0.0;
      }
      rt += (sqds + sm) * (sqds - sm);
      nrm = sqds * (e[ii] / scale);
      for (int k{ii + 1}; k <= qq_tmp; k++) {
        sm = blas::xrotg(rt, nrm, sqds);
        if (k > ii + 1) {
          e[0] = rt;
        }
        nrm = e[k - 1];
        b = b_s[k - 1];
        e[k - 1] = sm * nrm - sqds * b;
        rt = sqds * b_s[k];
        b_s[k] *= sm;
        qjj = 3 * (k - 1) + 1;
        kase = 3 * k + 1;
        blas::xrot(V, qjj, kase, sm, sqds);
        b_s[k - 1] = sm * b + sqds * nrm;
        sm = blas::xrotg(b_s[k - 1], rt, sqds);
        b = e[k - 1];
        rt = sm * b + sqds * b_s[k];
        b_s[k] = -sqds * b + sm * b_s[k];
        nrm = sqds * e[k];
        e[k] *= sm;
        blas::xrot(U, qjj, kase, sm, sqds);
      }
      e[m] = rt;
      qq++;
    } break;
    default:
      if (b_s[ii] < 0.0) {
        b_s[ii] = -b_s[ii];
        kase = 3 * ii;
        qjj = kase + 3;
        for (int k{kase + 1}; k <= qjj; k++) {
          V[k - 1] = -V[k - 1];
        }
      }
      qp1 = ii + 1;
      while ((ii + 1 < 3) && (b_s[ii] < b_s[qp1])) {
        rt = b_s[ii];
        b_s[ii] = b_s[qp1];
        b_s[qp1] = rt;
        qjj = 3 * ii + 1;
        kase = 3 * (ii + 1) + 1;
        blas::xswap(V, qjj, kase);
        blas::xswap(U, qjj, kase);
        ii = qp1;
        qp1++;
      }
      qq = 0;
      m--;
      break;
    }
  }
  s[0] = b_s[0];
  s[1] = b_s[1];
  s[2] = b_s[2];
}

} // namespace internal
static poseGraph *
optimizePoseGraph(poseGraph &b_poseGraph,
                  robotics::core::internal::BlockMatrix &iobj_0,
                  poseGraph &iobj_1)
{
  robotics::core::internal::BlockMatrix lobj_1[2];
  robotics::core::internal::BlockMatrix *obj;
  robotics::core::internal::BlockMatrix *posesUpdated;
  robotics::core::internal::TrustRegionIndefiniteDogLegSE2 solver;
  sparse hessian;
  ::coder::array<double, 2U> T1;
  ::coder::array<double, 2U> varargin_1;
  ::coder::array<double, 1U> c_poseGraph;
  ::coder::array<boolean_T, 1U> d_poseGraph;
  double R[9];
  double T1Offset[9];
  double c_T_tmp[9];
  double paramStruct_FirstNodePose[3];
  double R_idx_0;
  double R_idx_3;
  double T_tmp;
  double b_T_tmp;
  double paramStruct_FunctionTolerance;
  double paramStruct_GradientTolerance;
  double paramStruct_InitialTrustRegionRadius;
  double paramStruct_StepTolerance;
  double paramStruct_TrustRegionRadiusTolerance;
  double parsedResults_idx_0;
  double parsedResults_idx_1;
  int boffset;
  int coffset;
  int i;
  int loop_ub;
  boolean_T expl_temp;
  R_idx_0 =
      nav::algs::internal::PoseGraphOptimizer::parseOptimizePoseGraphInputs(
          R_idx_3, paramStruct_FunctionTolerance, expl_temp,
          paramStruct_GradientTolerance, paramStruct_StepTolerance,
          paramStruct_InitialTrustRegionRadius, paramStruct_FirstNodePose,
          paramStruct_TrustRegionRadiusTolerance, parsedResults_idx_0);
  parsedResults_idx_0 = b_poseGraph.MaxNumEdges;
  parsedResults_idx_1 = b_poseGraph.MaxNumNodes;
  iobj_1.MaxNumEdges = parsedResults_idx_0;
  iobj_1.MaxNumNodes = parsedResults_idx_1;
  iobj_1.NumEdges = 0.0;
  iobj_1.NumNodes = 1.0;
  iobj_1.NumLoopClosureEdges = 0.0;
  parsedResults_idx_0 = iobj_1.MaxNumNodes;
  i = static_cast<int>(parsedResults_idx_0 * 3.0);
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0].Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0].Matrix[i] = 0.0;
  }
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0].BlockSize[0] = 3.0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0].BlockSize[1] = 3.0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0].NumRowBlocks = parsedResults_idx_0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0].NumColBlocks = 1.0;
  iobj_1.NodeEstimates = &(&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[0];
  iobj_1.NodeEstimates->replaceBlock();
  loop_ub = static_cast<int>(iobj_1.MaxNumNodes);
  iobj_1.NodeMap.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.NodeMap[i] = 0.0;
  }
  loop_ub = static_cast<int>(iobj_1.MaxNumNodes);
  iobj_1.NodeDims.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.NodeDims[i] = 0.0;
  }
  loop_ub = static_cast<int>(iobj_1.MaxNumNodes);
  iobj_1.IsLandmarkNode.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.IsLandmarkNode[i] = false;
  }
  iobj_1.NodeMap[0] = 1.0;
  iobj_1.NodeDims[0] = 3.0;
  loop_ub = static_cast<int>(iobj_1.MaxNumEdges);
  iobj_1.EdgeNodePairs.set_size(loop_ub, 2);
  loop_ub <<= 1;
  for (i = 0; i < loop_ub; i++) {
    iobj_1.EdgeNodePairs[i] = 0.0;
  }
  loop_ub = static_cast<int>(iobj_1.MaxNumEdges);
  iobj_1.LoopClosureEdgeNodePairs.set_size(loop_ub, 2);
  loop_ub <<= 1;
  for (i = 0; i < loop_ub; i++) {
    iobj_1.LoopClosureEdgeNodePairs[i] = 0.0;
  }
  loop_ub = static_cast<int>(iobj_1.MaxNumEdges);
  iobj_1.LoopClosureEdgeIDsInternal.set_size(1, loop_ub);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.LoopClosureEdgeIDsInternal[i] = 0.0;
  }
  parsedResults_idx_0 = iobj_1.MaxNumEdges;
  i = static_cast<int>(parsedResults_idx_0 * 3.0);
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1].Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1].Matrix[i] = 0.0;
  }
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1].BlockSize[0] = 3.0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1].BlockSize[1] = 3.0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1].NumRowBlocks = parsedResults_idx_0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1].NumColBlocks = 1.0;
  iobj_1.EdgeMeasurements = &(&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[1];
  parsedResults_idx_0 = iobj_1.MaxNumEdges;
  i = static_cast<int>(parsedResults_idx_0 * 3.0);
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2].Matrix.set_size(i, 3);
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2].Matrix[i] = 0.0;
  }
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2].BlockSize[0] = 3.0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2].BlockSize[1] = 3.0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2].NumRowBlocks = parsedResults_idx_0;
  (&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2].NumColBlocks = 1.0;
  iobj_1.EdgeInfoMatrices = &(&(&(&(&(&(&iobj_0)[0])[0])[0])[0])[0])[2];
  loop_ub = b_poseGraph.EdgeNodePairs.size(0) << 1;
  iobj_1.EdgeNodePairs.set_size(b_poseGraph.EdgeNodePairs.size(0), 2);
  c_poseGraph.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    c_poseGraph[i] = b_poseGraph.EdgeNodePairs[i];
  }
  loop_ub = c_poseGraph.size(0);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.EdgeNodePairs[i] = c_poseGraph[i];
  }
  loop_ub = b_poseGraph.LoopClosureEdgeNodePairs.size(0) << 1;
  iobj_1.LoopClosureEdgeNodePairs.set_size(
      b_poseGraph.LoopClosureEdgeNodePairs.size(0), 2);
  c_poseGraph.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    c_poseGraph[i] = b_poseGraph.LoopClosureEdgeNodePairs[i];
  }
  loop_ub = c_poseGraph.size(0);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.LoopClosureEdgeNodePairs[i] = c_poseGraph[i];
  }
  loop_ub = b_poseGraph.LoopClosureEdgeIDsInternal.size(1);
  iobj_1.LoopClosureEdgeIDsInternal.set_size(
      1, b_poseGraph.LoopClosureEdgeIDsInternal.size(1));
  c_poseGraph.set_size(loop_ub);
  for (i = 0; i < loop_ub; i++) {
    c_poseGraph[i] = b_poseGraph.LoopClosureEdgeIDsInternal[i];
  }
  loop_ub = c_poseGraph.size(0);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.LoopClosureEdgeIDsInternal[i] = c_poseGraph[i];
  }
  obj = b_poseGraph.NodeEstimates;
  varargin_1.set_size(obj->Matrix.size(0), 3);
  loop_ub = obj->Matrix.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    varargin_1[i] = obj->Matrix[i];
  }
  parsedResults_idx_0 = obj->BlockSize[0];
  parsedResults_idx_1 = obj->BlockSize[1];
  (&(&(&iobj_0)[0])[0])[3].Matrix.set_size(varargin_1.size(0), 3);
  loop_ub = varargin_1.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    (&(&(&iobj_0)[0])[0])[3].Matrix[i] = varargin_1[i];
  }
  (&(&(&iobj_0)[0])[0])[3].BlockSize[0] = parsedResults_idx_0;
  (&(&(&iobj_0)[0])[0])[3].BlockSize[1] = parsedResults_idx_1;
  (&(&(&iobj_0)[0])[0])[3].NumRowBlocks =
      static_cast<double>(varargin_1.size(0)) / parsedResults_idx_0;
  (&(&(&iobj_0)[0])[0])[3].NumColBlocks = 3.0 / parsedResults_idx_1;
  iobj_1.NodeEstimates = &(&(&(&iobj_0)[0])[0])[3];
  obj = b_poseGraph.EdgeMeasurements;
  varargin_1.set_size(obj->Matrix.size(0), 3);
  loop_ub = obj->Matrix.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    varargin_1[i] = obj->Matrix[i];
  }
  parsedResults_idx_0 = obj->BlockSize[0];
  parsedResults_idx_1 = obj->BlockSize[1];
  (&(&(&iobj_0)[0])[0])[4].Matrix.set_size(varargin_1.size(0), 3);
  loop_ub = varargin_1.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    (&(&(&iobj_0)[0])[0])[4].Matrix[i] = varargin_1[i];
  }
  (&(&(&iobj_0)[0])[0])[4].BlockSize[0] = parsedResults_idx_0;
  (&(&(&iobj_0)[0])[0])[4].BlockSize[1] = parsedResults_idx_1;
  (&(&(&iobj_0)[0])[0])[4].NumRowBlocks =
      static_cast<double>(varargin_1.size(0)) / parsedResults_idx_0;
  (&(&(&iobj_0)[0])[0])[4].NumColBlocks = 3.0 / parsedResults_idx_1;
  iobj_1.EdgeMeasurements = &(&(&(&iobj_0)[0])[0])[4];
  obj = b_poseGraph.EdgeInfoMatrices;
  varargin_1.set_size(obj->Matrix.size(0), 3);
  loop_ub = obj->Matrix.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    varargin_1[i] = obj->Matrix[i];
  }
  parsedResults_idx_0 = obj->BlockSize[0];
  parsedResults_idx_1 = obj->BlockSize[1];
  (&(&(&iobj_0)[0])[0])[5].Matrix.set_size(varargin_1.size(0), 3);
  loop_ub = varargin_1.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    (&(&(&iobj_0)[0])[0])[5].Matrix[i] = varargin_1[i];
  }
  (&(&(&iobj_0)[0])[0])[5].BlockSize[0] = parsedResults_idx_0;
  (&(&(&iobj_0)[0])[0])[5].BlockSize[1] = parsedResults_idx_1;
  (&(&(&iobj_0)[0])[0])[5].NumRowBlocks =
      static_cast<double>(varargin_1.size(0)) / parsedResults_idx_0;
  (&(&(&iobj_0)[0])[0])[5].NumColBlocks = 3.0 / parsedResults_idx_1;
  iobj_1.EdgeInfoMatrices = &(&(&(&iobj_0)[0])[0])[5];
  iobj_1.NumNodes = b_poseGraph.NumNodes;
  c_poseGraph.set_size(b_poseGraph.NodeMap.size(0));
  loop_ub = b_poseGraph.NodeMap.size(0);
  for (i = 0; i < loop_ub; i++) {
    c_poseGraph[i] = b_poseGraph.NodeMap[i];
  }
  iobj_1.NodeMap.set_size(c_poseGraph.size(0));
  loop_ub = c_poseGraph.size(0);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.NodeMap[i] = c_poseGraph[i];
  }
  c_poseGraph.set_size(b_poseGraph.NodeDims.size(0));
  loop_ub = b_poseGraph.NodeDims.size(0);
  for (i = 0; i < loop_ub; i++) {
    c_poseGraph[i] = b_poseGraph.NodeDims[i];
  }
  iobj_1.NodeDims.set_size(c_poseGraph.size(0));
  loop_ub = c_poseGraph.size(0);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.NodeDims[i] = c_poseGraph[i];
  }
  d_poseGraph.set_size(b_poseGraph.IsLandmarkNode.size(0));
  loop_ub = b_poseGraph.IsLandmarkNode.size(0);
  for (i = 0; i < loop_ub; i++) {
    d_poseGraph[i] = b_poseGraph.IsLandmarkNode[i];
  }
  iobj_1.IsLandmarkNode.set_size(d_poseGraph.size(0));
  loop_ub = d_poseGraph.size(0);
  for (i = 0; i < loop_ub; i++) {
    iobj_1.IsLandmarkNode[i] = d_poseGraph[i];
  }
  iobj_1.NumEdges = b_poseGraph.NumEdges;
  iobj_1.NumLoopClosureEdges = b_poseGraph.NumLoopClosureEdges;
  T_tmp = std::sin(paramStruct_FirstNodePose[2]);
  b_T_tmp = std::cos(paramStruct_FirstNodePose[2]);
  solver.TimeObj.StartTime.tv_sec = 0.0;
  solver.TimeObj.StartTime.tv_nsec = 0.0;
  solver.MaxNumIteration = R_idx_0;
  solver.MaxTime = R_idx_3;
  solver.GradientTolerance = paramStruct_GradientTolerance;
  solver.StepTolerance = paramStruct_StepTolerance;
  solver.FunctionTolerance = paramStruct_FunctionTolerance;
  solver.InitialTrustRegionRadius = paramStruct_InitialTrustRegionRadius;
  solver.TrustRegionRadiusTolerance = paramStruct_TrustRegionRadiusTolerance;
  R_idx_0 = iobj_1.NumEdges;
  if (R_idx_0 < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(R_idx_0);
  }
  solver.ExtraArgs.edgeNodePairs.set_size(loop_ub, 2);
  for (i = 0; i < 2; i++) {
    for (boffset = 0; boffset < loop_ub; boffset++) {
      solver.ExtraArgs
          .edgeNodePairs[boffset + solver.ExtraArgs.edgeNodePairs.size(0) * i] =
          iobj_1.EdgeNodePairs[boffset + iobj_1.EdgeNodePairs.size(0) * i];
    }
  }
  R_idx_0 = iobj_1.NumNodes * 3.0;
  if (R_idx_0 < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(R_idx_0);
  }
  varargin_1.set_size(loop_ub, 3);
  for (i = 0; i < 3; i++) {
    for (boffset = 0; boffset < loop_ub; boffset++) {
      varargin_1[boffset + varargin_1.size(0) * i] =
          iobj_1.NodeEstimates
              ->Matrix[boffset + iobj_1.NodeEstimates->Matrix.size(0) * i];
    }
  }
  R_idx_0 = iobj_1.NumEdges * 3.0;
  if (R_idx_0 < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(R_idx_0);
  }
  solver.ExtraArgs.edgeInfoMats.set_size(loop_ub, 3);
  for (i = 0; i < 3; i++) {
    for (boffset = 0; boffset < loop_ub; boffset++) {
      solver.ExtraArgs
          .edgeInfoMats[boffset + solver.ExtraArgs.edgeInfoMats.size(0) * i] =
          iobj_1.EdgeInfoMatrices
              ->Matrix[boffset + iobj_1.EdgeInfoMatrices->Matrix.size(0) * i];
    }
  }
  R_idx_0 = iobj_1.NumEdges * 3.0;
  if (R_idx_0 < 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(R_idx_0);
  }
  solver.ExtraArgs.edgeMeasurements.set_size(loop_ub, 3);
  for (i = 0; i < 3; i++) {
    for (boffset = 0; boffset < loop_ub; boffset++) {
      solver.ExtraArgs
          .edgeMeasurements[boffset +
                            solver.ExtraArgs.edgeMeasurements.size(0) * i] =
          iobj_1.EdgeMeasurements
              ->Matrix[boffset + iobj_1.EdgeMeasurements->Matrix.size(0) * i];
    }
  }
  solver.ExtraArgs.tformSize[0] = 3.0;
  solver.ExtraArgs.infoMatSize[0] = 3.0;
  solver.ExtraArgs.tformSize[1] = 3.0;
  solver.ExtraArgs.infoMatSize[1] = 3.0;
  solver.ExtraArgs.poseDeltaLength = 3.0;
  solver.ExtraArgs.nodeMap.set_size(iobj_1.NodeMap.size(0));
  loop_ub = iobj_1.NodeMap.size(0);
  for (i = 0; i < loop_ub; i++) {
    solver.ExtraArgs.nodeMap[i] = iobj_1.NodeMap[i];
  }
  solver.ExtraArgs.nodeDims.set_size(iobj_1.NodeDims.size(0));
  loop_ub = iobj_1.NodeDims.size(0);
  for (i = 0; i < loop_ub; i++) {
    solver.ExtraArgs.nodeDims[i] = iobj_1.NodeDims[i];
  }
  solver.ExtraArgs.IsLandmarkNode.set_size(iobj_1.IsLandmarkNode.size(0));
  loop_ub = iobj_1.IsLandmarkNode.size(0);
  for (i = 0; i < loop_ub; i++) {
    solver.ExtraArgs.IsLandmarkNode[i] = iobj_1.IsLandmarkNode[i];
  }
  solver.solve(varargin_1, lobj_1[0], &posesUpdated, hessian,
               parsedResults_idx_1, R_idx_0);
  parsedResults_idx_0 = posesUpdated->BlockSize[0] * 0.0 + 1.0;
  parsedResults_idx_1 = posesUpdated->BlockSize[1] * 0.0 + 1.0;
  R_idx_0 = (parsedResults_idx_0 + posesUpdated->BlockSize[0]) - 1.0;
  if (parsedResults_idx_0 > R_idx_0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(R_idx_0);
  }
  R_idx_0 = (parsedResults_idx_1 + posesUpdated->BlockSize[1]) - 1.0;
  if (parsedResults_idx_1 > R_idx_0) {
    coffset = 0;
  } else {
    coffset = static_cast<int>(R_idx_0);
  }
  T1.set_size(loop_ub, coffset);
  for (i = 0; i < coffset; i++) {
    for (boffset = 0; boffset < loop_ub; boffset++) {
      T1[boffset + T1.size(0) * i] =
          posesUpdated->Matrix[boffset + posesUpdated->Matrix.size(0) * i];
    }
  }
  R_idx_0 = T1[0];
  parsedResults_idx_0 = T1[T1.size(0)];
  parsedResults_idx_1 = T1[1];
  R_idx_3 = T1[T1.size(0) + 1];
  c_T_tmp[0] = b_T_tmp;
  c_T_tmp[3] = -T_tmp;
  c_T_tmp[6] = paramStruct_FirstNodePose[0];
  c_T_tmp[1] = T_tmp;
  c_T_tmp[4] = b_T_tmp;
  c_T_tmp[7] = paramStruct_FirstNodePose[1];
  c_T_tmp[2] = 0.0;
  c_T_tmp[5] = 0.0;
  c_T_tmp[8] = 1.0;
  R[0] = R_idx_0;
  R[1] = parsedResults_idx_0;
  R[6] = -R_idx_0 * T1[T1.size(0) * 2] +
         -parsedResults_idx_1 * T1[T1.size(0) * 2 + 1];
  R[3] = parsedResults_idx_1;
  R[4] = R_idx_3;
  R[7] = -parsedResults_idx_0 * T1[T1.size(0) * 2] +
         -R_idx_3 * T1[T1.size(0) * 2 + 1];
  R[2] = 0.0;
  R[5] = 0.0;
  R[8] = 1.0;
  for (i = 0; i < 3; i++) {
    R_idx_0 = c_T_tmp[i];
    parsedResults_idx_0 = c_T_tmp[i + 3];
    parsedResults_idx_1 = c_T_tmp[i + 6];
    for (boffset = 0; boffset < 3; boffset++) {
      T1Offset[i + 3 * boffset] = (R_idx_0 * R[3 * boffset] +
                                   parsedResults_idx_0 * R[3 * boffset + 1]) +
                                  parsedResults_idx_1 * R[3 * boffset + 2];
    }
  }
  R_idx_0 = iobj_1.NumNodes;
  i = static_cast<int>(R_idx_0);
  for (int b_i{0}; b_i < i; b_i++) {
    int y_size_idx_1;
    posesUpdated->extractBlock(static_cast<double>(b_i) + 1.0, T1);
    loop_ub = T1.size(1);
    y_size_idx_1 = T1.size(1);
    for (int j{0}; j < loop_ub; j++) {
      coffset = j * 3;
      boffset = j * T1.size(0);
      for (int c_i{0}; c_i < 3; c_i++) {
        c_T_tmp[coffset + c_i] = (T1Offset[c_i] * T1[boffset] +
                                  T1Offset[c_i + 3] * T1[boffset + 1]) +
                                 T1Offset[c_i + 6] * T1[boffset + 2];
      }
    }
    parsedResults_idx_0 =
        posesUpdated->BlockSize[0] * ((static_cast<double>(b_i) + 1.0) - 1.0) +
        1.0;
    R_idx_0 = (parsedResults_idx_0 + posesUpdated->BlockSize[0]) - 1.0;
    if (parsedResults_idx_0 > R_idx_0) {
      boffset = 1;
    } else {
      boffset = static_cast<int>(parsedResults_idx_0);
    }
    for (loop_ub = 0; loop_ub < y_size_idx_1; loop_ub++) {
      posesUpdated
          ->Matrix[(boffset + posesUpdated->Matrix.size(0) * loop_ub) - 1] =
          c_T_tmp[3 * loop_ub];
      posesUpdated->Matrix[boffset + posesUpdated->Matrix.size(0) * loop_ub] =
          c_T_tmp[3 * loop_ub + 1];
      posesUpdated
          ->Matrix[(boffset + posesUpdated->Matrix.size(0) * loop_ub) + 1] =
          c_T_tmp[3 * loop_ub + 2];
    }
  }
  varargin_1.set_size(posesUpdated->Matrix.size(0), 3);
  loop_ub = posesUpdated->Matrix.size(0) * 3;
  for (i = 0; i < loop_ub; i++) {
    varargin_1[i] = posesUpdated->Matrix[i];
  }
  loop_ub = varargin_1.size(0);
  for (i = 0; i < 3; i++) {
    for (boffset = 0; boffset < loop_ub; boffset++) {
      iobj_1.NodeEstimates
          ->Matrix[boffset + iobj_1.NodeEstimates->Matrix.size(0) * i] =
          varargin_1[boffset + varargin_1.size(0) * i];
    }
  }
  return &iobj_1;
}

static double tic(double &tstart_tv_nsec)
{
  coderTimespec b_timespec;
  double tstart_tv_sec;
  if (!freq_not_empty) {
    freq_not_empty = true;
    coderInitTimeFunctions(&freq);
  }
  coderTimeClockGettimeMonotonic(&b_timespec, freq);
  tstart_tv_sec = b_timespec.tv_sec;
  tstart_tv_nsec = b_timespec.tv_nsec;
  return tstart_tv_sec;
}

static double toc(double tstart_tv_sec, double tstart_tv_nsec)
{
  coderTimespec b_timespec;
  if (!freq_not_empty) {
    freq_not_empty = true;
    coderInitTimeFunctions(&freq);
  }
  coderTimeClockGettimeMonotonic(&b_timespec, freq);
  return (b_timespec.tv_sec - tstart_tv_sec) +
         (b_timespec.tv_nsec - tstart_tv_nsec) / 1.0E+9;
}

} // namespace coder
static int div_s32(int numerator, int denominator)
{
  int quotient;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    unsigned int tempAbsQuotient;
    unsigned int u;
    if (numerator < 0) {
      tempAbsQuotient = ~static_cast<unsigned int>(numerator) + 1U;
    } else {
      tempAbsQuotient = static_cast<unsigned int>(numerator);
    }
    if (denominator < 0) {
      u = ~static_cast<unsigned int>(denominator) + 1U;
    } else {
      u = static_cast<unsigned int>(denominator);
    }
    tempAbsQuotient /= u;
    if ((numerator < 0) != (denominator < 0)) {
      quotient = -static_cast<int>(tempAbsQuotient);
    } else {
      quotient = static_cast<int>(tempAbsQuotient);
    }
  }
  return quotient;
}

static void minus(::coder::array<double, 1U> &in1,
                  const ::coder::array<double, 1U> &in2,
                  const ::coder::array<double, 1U> &in3)
{
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  if (in3.size(0) == 1) {
    loop_ub = in2.size(0);
  } else {
    loop_ub = in3.size(0);
  }
  in1.set_size(loop_ub);
  stride_0_0 = (in2.size(0) != 1);
  stride_1_0 = (in3.size(0) != 1);
  for (int i{0}; i < loop_ub; i++) {
    in1[i] = in2[i * stride_0_0] - in3[i * stride_1_0];
  }
}

static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else if (std::isinf(u0) && std::isinf(u1)) {
    int i;
    int i1;
    if (u0 > 0.0) {
      i = 1;
    } else {
      i = -1;
    }
    if (u1 > 0.0) {
      i1 = 1;
    } else {
      i1 = -1;
    }
    y = std::atan2(static_cast<double>(i), static_cast<double>(i1));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
  }
  return y;
}

static double rt_hypotd_snf(double u0, double u1)
{
  double a;
  double b;
  double y;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else if (std::isnan(b)) {
    y = rtNaN;
  } else {
    y = a * 1.4142135623730951;
  }
  return y;
}

void poseGraphOptimize(const ::coder::array<double, 2U> &absposes,
                       const ::coder::array<double, 2U> &loopNodePairs,
                       const ::coder::array<double, 2U> &loopPoses,
                       ::coder::array<double, 2U> &updatedPoses)
{
  coder::poseGraph b_pg;
  coder::poseGraph pg;
  coder::robotics::core::internal::BlockMatrix lobj_2[9];
  coder::robotics::core::internal::BlockMatrix *iobj_0;
  coder::robotics::core::internal::BlockMatrix *obj;
  ::coder::array<double, 2U> T;
  ::coder::array<double, 2U> b_y1;
  ::coder::array<double, 2U> nodeIds;
  double relPose[3];
  double colStart;
  double d;
  double rowStart;
  int dimSize;
  int iyStart;
  int loop_ub;
  int r;
  if (!isInitialized_poseGraphOptimize) {
    poseGraphOptimize_initialize();
  }
  iobj_0 = &lobj_2[0];
  pg.MaxNumEdges = 20000.0;
  pg.MaxNumNodes = 10000.0;
  pg.NumEdges = 0.0;
  pg.NumNodes = 1.0;
  pg.NumLoopClosureEdges = 0.0;
  rowStart = pg.MaxNumNodes;
  r = static_cast<int>(rowStart * 3.0);
  (&iobj_0[0])[0].Matrix.set_size(r, 3);
  iyStart = r * 3;
  for (r = 0; r < iyStart; r++) {
    (&iobj_0[0])[0].Matrix[r] = 0.0;
  }
  (&iobj_0[0])[0].BlockSize[0] = 3.0;
  (&iobj_0[0])[0].BlockSize[1] = 3.0;
  (&iobj_0[0])[0].NumRowBlocks = rowStart;
  (&iobj_0[0])[0].NumColBlocks = 1.0;
  pg.NodeEstimates = &(&iobj_0[0])[0];
  obj = pg.NodeEstimates;
  rowStart = obj->BlockSize[0] * 0.0 + 1.0;
  colStart = obj->BlockSize[1] * 0.0 + 1.0;
  d = (rowStart + obj->BlockSize[0]) - 1.0;
  if (rowStart > d) {
    iyStart = 0;
  } else {
    iyStart = static_cast<int>(d);
  }
  d = (colStart + obj->BlockSize[1]) - 1.0;
  if (colStart > d) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>(d);
  }
  for (r = 0; r < loop_ub; r++) {
    for (dimSize = 0; dimSize < iyStart; dimSize++) {
      obj->Matrix[dimSize + obj->Matrix.size(0) * r] =
          iv[dimSize + iyStart * r];
    }
  }
  iyStart = static_cast<int>(pg.MaxNumNodes);
  pg.NodeMap.set_size(iyStart);
  for (r = 0; r < iyStart; r++) {
    pg.NodeMap[r] = 0.0;
  }
  iyStart = static_cast<int>(pg.MaxNumNodes);
  pg.NodeDims.set_size(iyStart);
  for (r = 0; r < iyStart; r++) {
    pg.NodeDims[r] = 0.0;
  }
  iyStart = static_cast<int>(pg.MaxNumNodes);
  pg.IsLandmarkNode.set_size(iyStart);
  for (r = 0; r < iyStart; r++) {
    pg.IsLandmarkNode[r] = false;
  }
  pg.NodeMap[0] = 1.0;
  pg.NodeDims[0] = 3.0;
  iyStart = static_cast<int>(pg.MaxNumEdges);
  pg.EdgeNodePairs.set_size(iyStart, 2);
  iyStart <<= 1;
  for (r = 0; r < iyStart; r++) {
    pg.EdgeNodePairs[r] = 0.0;
  }
  iyStart = static_cast<int>(pg.MaxNumEdges);
  pg.LoopClosureEdgeNodePairs.set_size(iyStart, 2);
  iyStart <<= 1;
  for (r = 0; r < iyStart; r++) {
    pg.LoopClosureEdgeNodePairs[r] = 0.0;
  }
  iyStart = static_cast<int>(pg.MaxNumEdges);
  pg.LoopClosureEdgeIDsInternal.set_size(1, iyStart);
  for (r = 0; r < iyStart; r++) {
    pg.LoopClosureEdgeIDsInternal[r] = 0.0;
  }
  rowStart = pg.MaxNumEdges;
  r = static_cast<int>(rowStart * 3.0);
  (&iobj_0[0])[1].Matrix.set_size(r, 3);
  iyStart = r * 3;
  for (r = 0; r < iyStart; r++) {
    (&iobj_0[0])[1].Matrix[r] = 0.0;
  }
  (&iobj_0[0])[1].BlockSize[0] = 3.0;
  (&iobj_0[0])[1].BlockSize[1] = 3.0;
  (&iobj_0[0])[1].NumRowBlocks = rowStart;
  (&iobj_0[0])[1].NumColBlocks = 1.0;
  pg.EdgeMeasurements = &(&iobj_0[0])[1];
  rowStart = pg.MaxNumEdges;
  r = static_cast<int>(rowStart * 3.0);
  (&iobj_0[0])[2].Matrix.set_size(r, 3);
  iyStart = r * 3;
  for (r = 0; r < iyStart; r++) {
    (&iobj_0[0])[2].Matrix[r] = 0.0;
  }
  (&iobj_0[0])[2].BlockSize[0] = 3.0;
  (&iobj_0[0])[2].BlockSize[1] = 3.0;
  (&iobj_0[0])[2].NumRowBlocks = rowStart;
  (&iobj_0[0])[2].NumColBlocks = 1.0;
  pg.EdgeInfoMatrices = &(&iobj_0[0])[2];
  dimSize = absposes.size(0);
  if (absposes.size(0) == 0) {
    b_y1.set_size(0, 3);
  } else {
    iyStart = absposes.size(0) - 1;
    if (iyStart > 1) {
      iyStart = 1;
    }
    if (iyStart < 1) {
      b_y1.set_size(0, 3);
    } else {
      b_y1.set_size(absposes.size(0) - 1, 3);
      if (absposes.size(0) - 1 != 0) {
        iyStart = 0;
        for (r = 0; r < 3; r++) {
          loop_ub = r * dimSize;
          rowStart = absposes[loop_ub];
          for (int m{2}; m <= dimSize; m++) {
            colStart = absposes[(loop_ub + m) - 1];
            d = colStart;
            colStart -= rowStart;
            rowStart = d;
            b_y1[(iyStart + m) - 2] = colStart;
          }
          iyStart = (iyStart + dimSize) - 1;
        }
      }
    }
  }
  r = b_y1.size(0);
  for (loop_ub = 0; loop_ub < r; loop_ub++) {
    relPose[2] = b_y1[loop_ub + b_y1.size(0) * 2];
    rowStart = absposes[loop_ub + absposes.size(0) * 2];
    colStart = std::sin(rowStart);
    rowStart = std::cos(rowStart);
    relPose[0] =
        rowStart * b_y1[loop_ub] + colStart * b_y1[loop_ub + b_y1.size(0)];
    relPose[1] =
        -colStart * b_y1[loop_ub] + rowStart * b_y1[loop_ub + b_y1.size(0)];
    pg.addRelativePose(relPose, static_cast<double>(loop_ub) + 1.0,
                       (static_cast<double>(loop_ub) + 1.0) + 1.0);
  }
  r = loopNodePairs.size(0);
  for (loop_ub = 0; loop_ub < r; loop_ub++) {
    relPose[0] = loopPoses[loop_ub];
    relPose[1] = loopPoses[loop_ub + loopPoses.size(0)];
    relPose[2] = loopPoses[loop_ub + loopPoses.size(0) * 2];
    pg.addRelativePose(relPose, loopNodePairs[loop_ub],
                       loopNodePairs[loop_ub + loopNodePairs.size(0)]);
  }
  coder::optimizePoseGraph(pg, lobj_2[3], b_pg);
  rowStart = b_pg.NumNodes;
  if (std::isnan(rowStart)) {
    nodeIds.set_size(1, 1);
    nodeIds[0] = rtNaN;
  } else if (rowStart < 1.0) {
    nodeIds.set_size(1, 0);
  } else {
    nodeIds.set_size(1, static_cast<int>(rowStart - 1.0) + 1);
    iyStart = static_cast<int>(rowStart - 1.0);
    for (r = 0; r <= iyStart; r++) {
      nodeIds[r] = static_cast<double>(r) + 1.0;
    }
  }
  r = static_cast<int>(rowStart);
  updatedPoses.set_size(r, 3);
  for (loop_ub = 0; loop_ub < r; loop_ub++) {
    d = nodeIds[loop_ub];
    if (b_pg.IsLandmarkNode[static_cast<int>(d) - 1]) {
      b_pg.NodeEstimates->extractBlock(d, T);
      updatedPoses[loop_ub] = T[T.size(0) * 2];
      updatedPoses[loop_ub + updatedPoses.size(0)] = T[T.size(0) * 2 + 1];
      updatedPoses[loop_ub + updatedPoses.size(0) * 2] = rtNaN;
    } else {
      b_pg.NodeEstimates->extractBlock(d, T);
      updatedPoses[loop_ub] = T[T.size(0) * 2];
      updatedPoses[loop_ub + updatedPoses.size(0)] = T[T.size(0) * 2 + 1];
      updatedPoses[loop_ub + updatedPoses.size(0) * 2] =
          rt_atan2d_snf(T[1], T[0]);
    }
  }
}

void poseGraphOptimize_initialize()
{
  freq_not_empty = false;
  isInitialized_poseGraphOptimize = true;
}

void poseGraphOptimize_terminate()
{
  isInitialized_poseGraphOptimize = false;
}

} // namespace poseGraphOptimize
