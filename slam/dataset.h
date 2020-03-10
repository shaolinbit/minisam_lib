#ifndef DATASET_H_INCLUDED
#define DATASET_H_INCLUDED

#include "../mat/MatCal.h"
#include "../mat/Matrix.h"
#include "../nonlinear/NonlinearFactorGraph.h"


namespace minisam
{

typedef std::pair<NonlinearFactorGraph*, std::map<int,minimatrix*>*> GraphAndValues;


/// Indicates how noise parameters are stored in file
enum NoiseFormat {
  NoiseFormatG2O, ///< Information matrix I11, I12, I13, I22, I23, I33
  NoiseFormatTORO, ///< Information matrix, but inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr
  NoiseFormatGRAPH, ///< default: toro-style order, but covariance matrix !
  NoiseFormatCOV, ///< Covariance matrix C11, C12, C13, C22, C23, C33
  NoiseFormatAUTO  ///< Try to guess covariance matrix layout
};

/// Robust kernel type to wrap around quadratic noise model
enum KernelFunctionType {
  KernelFunctionTypeNONE, KernelFunctionTypeHUBER, KernelFunctionTypeTUKEY
};

/**
 * Load TORO 3D Graph*/

GraphAndValues load3D(const std::string& filename);

/**
 * Load TORO/G2O style graph files
 * @param filename
 * @param model optional noise model to use instead of one specified by file
 * @param maxID if non-zero cut out vertices >= maxID
 * @param addNoise add noise to the edges
 * @param smart try to reduce complexity of covariance to cheapest model
 * @param noiseFormat how noise parameters are stored
 * @param kernelFunctionType whether to wrap the noise model in a robust kernel
 * @return graph and initial values
 */
GraphAndValues load2D(const std::string& filename,
    GaussianNoiseModel* model =NULL, int maxID = 0, bool addNoise =
        false, bool smart = true, NoiseFormat noiseFormat = NoiseFormatAUTO, //
    KernelFunctionType kernelFunctionType = KernelFunctionTypeNONE);
/**
 * @brief This function parses a g2o file and stores the measurements into a
 * NonlinearFactorGraph and the initial guess in a Values structure
 * @param filename The name of the g2o file\
 * @param is3D indicates if the file describes a 2D or 3D problem
 * @param kernelFunctionType whether to wrap the noise model in a robust kernel
 * @return graph and initial values
 */
GraphAndValues readG2o(const std::string& g2oFile, const bool is3D = false,
    KernelFunctionType kernelFunctionType = KernelFunctionTypeNONE);
/**
 * @brief This function writes a g2o file from
 * NonlinearFactorGraph and a Values structure
 * @param filename The name of the g2o file to write
 * @param graph NonlinearFactor graph storing the measurements
 * @param estimate Values */

void writeG2o(const NonlinearFactorGraph& graph,
              const std::map<int,minimatrix*>& estimatep,const std::string& filename);
};
#endif // DATASET_H_INCLUDED
