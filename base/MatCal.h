#ifndef MATCAL_H_INCLUDED
#define MATCAL_H_INCLUDED

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <list>
#include "../gmfconfig.h"
#include "../base/Matrix.h"



/*
This file was added by Shaolin for matrix function in isam.
*/




namespace minisam
{

double dot3d(const Eigen::Vector3d &p, const Eigen::Vector3d &q,
             Eigen::MatrixXd* H1,
             Eigen::MatrixXd* H2);
double dot3d(const Eigen::Vector3d &p, const Eigen::Vector3d &q);

double norm3d(const Eigen::Vector3d &p);
double norm3d(const Eigen::Vector3d &p, Eigen::MatrixXd* H);




Eigen::Vector3d cross3d(const Eigen::Vector3d &p,
                        const Eigen::Vector3d &q,
                        Eigen::Matrix3d* H1,
                        Eigen::Matrix3d* H2);
Eigen::Vector3d cross3d(const Eigen::Vector3d &p,
                        const Eigen::Vector3d &q,
                        Eigen::Matrix3d* H1);
Eigen::Vector3d cross3d(const Eigen::Vector3d &p,
                        const Eigen::Vector3d &q);
Eigen::Vector3d normalize3d(const Eigen::Vector3d &p,
                            Eigen::Matrix3d* H);
Eigen::Vector3d normalize3d(const Eigen::Vector3d &p);



Eigen::VectorXd VectorValuesgetvector(const std::map<int,Eigen::VectorXd>& VectorValues,const std::vector<int>& keys);


//template <class TPFactor>
//std::vector<int>& getvectorint_fromiterator(TPFactor& nrParents_);//no use in minisam


double VectorValuesDot(std::map<int,Eigen::VectorXd>& vectorvalues);

double VectorValuesDot(const std::map<int,Eigen::VectorXd>& vectorvalues1,const std::map<int,Eigen::VectorXd>& vectorvalues2);

void VectorValuesScaleMultiply(std::map<int,Eigen::VectorXd> *vectorvalues,double a);

void VectorValuesaxpy(std::map<int,Eigen::VectorXd> *y,
                      const std::map<int,Eigen::VectorXd>& x,double a);

double ListVectorDot(std::list<Eigen::VectorXd>& Errors);


std::map<int,Eigen::VectorXd> VectorValuesZero(const std::map<int,Eigen::VectorXd>& vvz);


std::map<int,Eigen::VectorXd> VectorValuesRetract(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);
std::map<int,Eigen::VectorXd> DVectorValuesRetract(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);

void DVectorValuesRetract_atmap(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2,std::map<int,Eigen::VectorXd>& vectorvalues3);

std::map<int,Eigen::VectorXd> VectorValuesSub(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);

std::map<int,Eigen::VectorXd> VectorValuesAdd(const std::map<int,Eigen::VectorXd>& vectorvalues1,
        const std::map<int,Eigen::VectorXd>& vectorvalues2);



void VectorValuesUpdate(std::map<int,Eigen::VectorXd>* vectorvalues1,const std::map<int,Eigen::VectorXd>& vectorvalues2);

double VectorValuesNorm(const std::map<int,Eigen::VectorXd>& vvz);


double VectorValuesSquareNorm(const std::map<int,Eigen::VectorXd>& vvz);



std::vector<int> getvectorfrom3keys(int key1,int key2,int key3);
std::vector<int> getvectorfrom5keys(int key1,int key2,int key3,int key4,int key5);
std::vector<int> getvectorfrom6keys(int key1,int key2,int key3,int key4,int key5,int key6);

Eigen::VectorXd ComposeTwoVec3(const Eigen::Vector3d& v1,const Eigen::Vector3d& v2);

Eigen::VectorXd Vec3ToVecX3(const Eigen::Vector3d& v);

Eigen::MatrixXd Mat3ToMatX3(const Eigen::Matrix3d& v);

};
#endif // MATRIX_H_INCLUDED
