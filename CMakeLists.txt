cmake_minimum_required (VERSION 2.6)
project (minisam)

SET(CMAKE_INSTALL_PREFIX /usr/include)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fpic -std=c11 ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fpic -std=c++11 ")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic -fvisibility=hidden -fvisibility-inlines-hidden -O3 -std=c11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -fvisibility=hidden -fvisibility-inlines-hidden -O3 -std=c++11")


############# # for pppbayestree
#[[
find_package(Boost REQUIRED COMPONENTS
filesystem  
)
if(NOT Boost_FOUND)
    message("Not found Boost")
endif()

include_directories(${Boost_INCLUDE_DIRS})
message("${Boost_INCLUDE_DIRS}")
message("${Boost_LIBRARIES}")
]]

#############  # for pppbayestree


set(HEAD_FILES_root
	gmfconfig.h
)

set(HEAD_FILES_3rdparty
	./3rdparty/ccolamd.h
	./3rdparty/SuiteSparse_config.h
)
set(HEAD_FILES_mat
	./mat/MatCal.h
	./mat/Matrix.h
	./mat/GaussianBlockMatrix.h
)
set(HEAD_FILES_geometry
	./geometry/Cal3_S2.h
	./geometry/Cal3_S2Stereo.h
	./geometry/Cal3Unified.h
	./geometry/Cal3Bundler.h
	./geometry/CalibratedCamera.h
	./geometry/EssentialMatrix.h
	./geometry/OrientedPlane3.h
	./geometry/PinholeCameraCal3S2.h
	./geometry/PinholePoseCal3S2.h
	./geometry/Pose2.h
	./geometry/Pose3.h
	./geometry/Quaternion.h
	./geometry/Rot2.h
	./geometry/Rot3.h
	./geometry/SimpleCamera.h
	./geometry/SO3.h
	./geometry/StereoCamera.h
	./geometry/StereoPoint2.h
	./geometry/Unit3.h
)
set(HEAD_FILES_inference
	./inference/BayesTree.h
        ./inference/BayesTreeCliqueBase.h
	./inference/ClusterTree.h
        ./inference/EliminateableFactorGraph.h
	./inference/EliminationTree.h
	./inference/Factor.h
	./inference/FactorGraph.h
	./inference/JunctionTree.h
	./inference/Ordering.h
	./inference/Symbol.h
	./inference/VariableIndex.h
	./inference/SymbolicConditional.h
)
set(HEAD_FILES_linear
	./linear/GaussianBayesNet.h
	./linear/GaussianConditional.h
	./linear/GaussianFactorGraph.h
	./linear/HessianFactor.h
	./linear/JacobianFactor.h
	./linear/NoiseModel.h
	./linear/RealGaussianFactor.h
	./linear/Scatter.h
)
set(HEAD_FILES_miniblas
        ./miniblas/memorypool.h
	./miniblas/miniblas.h
	./miniblas/minilinalg.h
	./miniblas/minimatrix_double.h
	./miniblas/minivector_double.h
)
set(HEAD_FILES_navigation
	./navigation/GPSFactor.h
	./navigation/ImuBias.h
	./navigation/ImuFactor.h
	./navigation/ManifoldPreintegration.h
	./navigation/NavState.h
	./navigation/PreintegratedRotation.h
	./navigation/PreintegrationParams.h
	./navigation/TangentPreintegration.h
)
set(HEAD_FILES_nonlinear
	./nonlinear/DoglegOptimizer.h
        ./nonlinear/DoglegOptimizerImpl.h
	./nonlinear/GaussNewtonOptimizer.h
	./nonlinear/ISAM2.h
	./nonlinear/ISAM2Clique.h
	./nonlinear/ISAM2Data.h
	./nonlinear/ISAM2-impl.h
	./nonlinear/ISAM2JunctionTree.h
	./nonlinear/ISAM2Params.h
        ./nonlinear/LevenbergMarquardtOptimizer.h
        ./nonlinear/LevenbergMarquardtParams.h
        ./nonlinear/LevenbergMarquardtState.h
        ./nonlinear/Marginals.h
	./nonlinear/NonlinearFactor.h
	./nonlinear/NonlinearFactorGraph.h
        ./nonlinear/NonlinearOptimizer.h
        ./nonlinear/NonlinearOptimizerParams.h
        ./nonlinear/NonlinearOptimizerState.h
)
set(HEAD_FILES_slam
	./slam/BetweenFactor.h
	./slam/PriorFactor.h
	./slam/ProjectionFactor.h
        ./slam/UnaryFactor.h
)

file(GLOB imukittiexamplegps_Dogleg "examples/imugpskitti/imukittiexamplegps_Dogleg_Cholesky.cpp")
file(GLOB imukittiexamplegps_gaussiannewton "examples/imugpskitti/imukittiexamplegps_gaussiannewton_CHOLESKY.cpp")
file(GLOB pppbayestree "examples/pppbayestree/pppbayestree.cpp")
file(GLOB visualisam2dogleg "examples/visualisam2/visualisam2dogleg.cpp")
file(GLOB visualisam2_gaussiannewton "examples/visualisam2/visualisam2_gaussiannewton.cpp")


aux_source_directory(examples/pppbayestree/configReader configReader)
aux_source_directory(examples/pppbayestree/gnssNavigation gnssNavigation)
aux_source_directory(examples/pppbayestree/gpstk gpstk)
aux_source_directory(examples/pppbayestree/robustModels robustModels)
aux_source_directory(examples/pppbayestree/slam slam)
aux_source_directory(examples/pppbayestree/utils utils)

#include_directories(/usr/local/include/eigen3/)
include_directories(.)
include_directories(${PROJECT_SOURCE_DIR}/examples)
link_directories(${PROJECT_SOURCE_DIR})

add_executable(imukittiexamplegps_Dogleg ${imukittiexamplegps_Dogleg})
target_link_libraries(imukittiexamplegps_Dogleg minisam)
add_executable(imukittiexamplegps_gaussiannewton ${imukittiexamplegps_gaussiannewton})
target_link_libraries(imukittiexamplegps_gaussiannewton minisam)

############################for pppbayestree
#[[
add_executable(pppbayestree  ${utils} ${configReader} ${gnssNavigation} ${gpstk} ${robustModels} ${slam} ${pppbayestree} )
target_link_libraries(pppbayestree minisam ${Boost_LIBRARIES})
]]
############################for pppbayestree


#add_executable(pppbayestree_dogleg ${pppbayestree_dogleg} ${utils} ${configReader} ${gnssNavigation} ${gpstk} ${robustModels} ${slam} ${pppbayestree} )
#target_link_libraries(pppbayestree_dogleg minisam)
add_executable(visualisam2dogleg ${visualisam2dogleg})
target_link_libraries(visualisam2dogleg minisam)
add_executable(visualisam2_gaussiannewton ${visualisam2_gaussiannewton})
target_link_libraries(visualisam2_gaussiannewton minisam)


install(FILES ${HEAD_FILES_root} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam)
install(FILES ${HEAD_FILES_3rdparty} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/3rdparty)
install(FILES ${HEAD_FILES_mat} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/mat)
install(FILES ${HEAD_FILES_geometry} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/geometry)
install(FILES ${HEAD_FILES_inference} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/inference)
install(FILES ${HEAD_FILES_linear} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/linear)
install(FILES ${HEAD_FILES_miniblas} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/miniblas)
install(FILES ${HEAD_FILES_navigation} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/navigation)
install(FILES ${HEAD_FILES_nonlinear} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/nonlinear)
install(FILES ${HEAD_FILES_slam} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/slam)
#install(FILES ${HEAD_FILES_symbolic} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/symbolic)
install(FILES libminisam.so DESTINATION /usr/lib)


