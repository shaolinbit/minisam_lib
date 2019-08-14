cmake_minimum_required (VERSION 2.6)
project (minisam)

SET(CMAKE_INSTALL_PREFIX /usr/include)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_CXX_COMPILER g++)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fpic -std=c11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fpic -std=c++11")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpic -fvisibility=hidden -fvisibility-inlines-hidden -O3 -std=c11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -fvisibility=hidden -fvisibility-inlines-hidden -O3 -std=c++11")

set(HEAD_FILES_root
	gmfconfig.h
)

set(HEAD_FILES_3rdparty
	./3rdparty/ccolamd.h
	./3rdparty/SuiteSparse_config.h
)
set(HEAD_FILES_base
	./base/MatCal.h
	./base/Matrix.h
	./base/SVBlockMatrix.h
)
set(HEAD_FILES_geometry
	./geometry/Cal3_S2.h
	./geometry/Cal3_S2Stereo.h
	./geometry/CalibratedCamera.h
	./geometry/PinholeCameraCal3S2.h
	./geometry/PinholePoseCal3S2.h
	./geometry/Pose2.h
	./geometry/Pose3.h
	./geometry/Quaternion.h
	./geometry/Rot2.h
	./geometry/Rot3.h
	./geometry/SimpleCamera.h
	./geometry/SO3.h
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
	./inference/VariableSlots.h
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
set(HEAD_FILES_navigation
	./navigation/ImuBias.h
	./navigation/ImuFactor.h
	./navigation/ManifoldPreintegration.h
	./navigation/NavState.h
	./navigation/PreintegratedRotation.h
	./navigation/PreintegrationBase.h
	./navigation/PreintegrationParams.h
	./navigation/TangentPreintegration.h
)
set(HEAD_FILES_nonlinear
	./nonlinear/DoglegOptimizer.h
        ./nonlinear/DoglegOptimizerImpl.h
	./nonlinear/ISAM2.h
	./nonlinear/ISAM2Data.h
	./nonlinear/ISAM2Clique.h
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
	./slam/BetweenFactorPose2.h
	./slam/BetweenFactorPose3.h
	./slam/PriorFactor.h
	./slam/PriorFactorPose2.h
	./slam/PriorFactorPose3.h
	./slam/ProjectionFactor.h
        ./slam/UnaryFactor.h
)
set(HEAD_FILES_symbolic
	./symbolic/SymbolicConditional.h
 )

install(FILES ${HEAD_FILES_root} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam)
install(FILES ${HEAD_FILES_3rdparty} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/3rdparty)
install(FILES ${HEAD_FILES_base} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/base)
install(FILES ${HEAD_FILES_geometry} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/geometry)
install(FILES ${HEAD_FILES_inference} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/inference)
install(FILES ${HEAD_FILES_linear} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/linear)
install(FILES ${HEAD_FILES_navigation} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/navigation)
install(FILES ${HEAD_FILES_nonlinear} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/nonlinear)
install(FILES ${HEAD_FILES_slam} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/slam)
install(FILES ${HEAD_FILES_symbolic} DESTINATION  ${CMAKE_INSTALL_PREFIX}/minisam/symbolic)
install(FILES libminisam.so DESTINATION /usr/lib)


