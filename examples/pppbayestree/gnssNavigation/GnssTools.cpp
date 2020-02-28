/**
 * @file   GnssTools.cpp
 * @brief  Tools required to process GNSS data -- (i.e. ECEF to ENU transformation)
 * @author Ryan Watson & Jason Gross
 */

#include "GnssTools.h"

minivector obsMap(const minivector &p1, const minivector &p2, const int &Trop)
{
    /*
       inputs ::
       p1 --> ECEF xyz coordinates of satellite [meter]
       p2 --> ECEF xyz coordinates of receiver [meter]
       Trop --> Troposphere modeling switch
       outputs ::
       H --> measurement mapping
     */
    double r = sqrt( ((p1.data[0]-p2.data[0]))*(p1.data[0]-p2.data[0]) +
                     ((p1.data[1]-p2.data[1])*(p1.data[1]-p2.data[1])) + ((p1.data[2]-p2.data[2])*(p1.data[2]-p2.data[2])) );
    if (Trop == 1)
    {
        double el = calcEl(p1,p2);
        double mapT = tropMap(el);
        minivector H(5);
        // H << (p1.data[0]-p2.data[0])/r, (p1.data[1]- p2.data[1])/r,
        //(p1.data[2]-p2.data[2])/r, 1.0,mapT;
        H.data[0]=(p1.data[0]-p2.data[0])/r;
        H.data[1]=(p1.data[1]- p2.data[1])/r;
        H.data[2]=(p1.data[2]-p2.data[2])/r;
        H.data[3]=1.0;
        H.data[4]=mapT;
        return H;
    }
    else
    {
        minivector H(5);
        //H << (p1.data[0]-p2.data[0])/r, (p1.data[1]- p2.data[1])/r,
        //(p1.data[2]-p2.data[2])/r, 1.0, 0.0;
        H.data[0]=(p1.data[0]-p2.data[0])/r;
        H.data[1]=(p1.data[1]- p2.data[1])/r;
        H.data[2]=(p1.data[2]-p2.data[2])/r;
        H.data[3]=1.0;
        H.data[4]=0.0;
        return H;
    }
}

minivector  obsMapNED(const minivector &p1, const minivector &p2, const int &Trop)
{
    /*
       inputs ::
       p1 --> ECEF xyz coordinates of satellite [meter]
       p2 --> ECEF xyz coordinates of receiver [meter]
       Trop --> Troposphere modeling switch
       outputs ::
       H --> measurement mapping
     */
    double r = sqrt( ((p1.data[0]-p2.data[0]))*(p1.data[0]-p2.data[0]) + ((p1.data[1]-p2.data[1])*(p1.data[1]-p2.data[1])) + ((p1.data[2]-p2.data[2])*(p1.data[2]-p2.data[2])) );
    if (Trop == 1)
    {
        double el = calcElNed(p1);
        double mapT = tropMap(el);
        minivector H(5);
        //  H << (p1.data[0]-p2.data[0])/r, (p1.data[1]- p2.data[1])/r,
        // (p1.data[2]-p2.data[2])/r, 1.0,mapT;
        H.data[0]=(p1.data[0]-p2.data[0])/r;
        H.data[1]=(p1.data[1]-p2.data[1])/r;
        H.data[2]=(p1.data[2]-p2.data[2])/r;
        H.data[3]=1.0;
        H.data[4]=mapT;

        return H;
    }
    else
    {
        minivector H(5);
        //H << (p1.data[0]-p2.data[0])/r, (p1.data[1]- p2.data[1])/r,
        // (p1.data[2]-p2.data[2])/r, 1.0, 0.0;
        H.data[0]=(p1.data[0]-p2.data[0])/r;
        H.data[1]=(p1.data[1]-p2.data[1])/r;
        H.data[2]=(p1.data[2]-p2.data[2])/r;
        H.data[3]=1.0;
        H.data[4]=0.0;

        return H;
    }
}

mini_int_vector getPRN(const minimatrix &p)
{
    /*
       inputs ::
       p --> GNSS data packed in pre-specified format
       outputs ::
       prnVec --> Vector of prn numbers corresponding to visible satellites at given epoch
     */

    int nSats = static_cast<int>(p.data[0]);
    mini_int_vector prnVec(nSats);
    int count = 2;
    for (int i=0; i<nSats; i++)
    {
        prnVec.data[i]= static_cast<int>(p.data[count]);
        count = count+7;
    }
    return prnVec;
}

bool checkPRN(const mini_int_vector &p, const int &n)
{
    /*
       input ::
       p --> prn vector
       n --> prn value to check
       output ::
       true --> if n is present in current prn vec.
       false --> otherwise
     */

    for (int i=0; i<p.size; i++)
    {
        if (p.data[i] == n )
        {
            return true;
        }
    }
    return false;
}

minimatrix earthToNavTrans( const minivector &p1)
{
    /*
       inputs ::
        p1 --> ECEF xyz coordinates of the platform
       outputs ::
        R --> rotation matrix from Earth fram to Navigation frame
     */
    minivector llh=xyz2llh(p1);
    double sLat = sin(llh.data[0]);
    double sLon = sin(llh.data[1]);
    double cLat = cos(llh.data[0]);
    double cLon = cos(llh.data[1]);

    //Matrix R = (Matrix(3,3) << -1*sLat*cLon, -1*sLon, -1*cLat*cLon,
    //            -1*sLat*sLon, cLon, -1*cLat*sLon,
    //             cLat, 0.0, -1*sLat ).finished();
    minimatrix R(3,3);
    //  R<< -1*sLat*cLon, -1*sLon, -1*cLat*cLon,
    //             -1*sLat*sLon, cLon, -1*cLat*sLon,
    //             cLat, 0.0, -1*sLat;
    R.data[0]=-1*sLat*cLon;
    R.data[1]=-1*sLon;
    R.data[2]=-1*cLat*cLon;
    R.data[3]=-1*sLat*sLon;
    R.data[4]=cLon;
    R.data[5]=-1*cLat*sLon;
    R.data[6]=cLat;
    R.data[7]=0.0;
    R.data[8]= -1*sLat;

    minimatrix Rinv(3,3);
    minimatrix tempcov(R);
    minilinalg_luMatInverse(&Rinv,&tempcov);

    return Rinv;

    // return R.inverse();

}

double deltaObs(const minivector &p1, const minivector &p2, const double &meas)
{
    /*
       inputs ::
       p1 --> ECEF xyz coordinates of satellite [meter]
       p2 --> ECEF xyz coordinates of receiver [meter]
       meas --> observed range between satellite and receiver [meter]
       outputs ::
       deltaR --> difference observed pseudorange and modeled pseudorange [meter]
     */
    //double r = minisam::norm3d(p1-p2) + tropMap(calcEl(p1,p2))*tropDry(p2);
    minivector subp1p2(3);
    minivector_sub(&subp1p2,p1,p2);
    double r = minisam::norm3d(subp1p2);
    r += tropMap(calcEl(p1,p2))*tropDry(p2);
    return meas - r;
}

double deltaTrop(const minivector &p1, const minivector &p2)
{
    /*
       inputs ::
       p1 --> ECEF xyz coordinates of satellite [meter]
       p2 --> ECEF xyz coordinates of receiver [meter]
       meas --> observed range between satellite and receiver [meter]
     */
    return tropMap(calcEl(p1,p2))*tropDry(p2);
}

minivector xyz2llh(const minivector &p1)
{
    /*
       inputs ::
       p1 --> ECEF xyz receiver coordinates [meter]
       output ::
       posLLH --> latitude, longitude, height [rad,rad,meter]
     */
    double x2= pow(p1.data[0],2);
    double y2= pow(p1.data[1],2);
    double z2= pow(p1.data[2],2);

    double e= sqrt(1.0-((GNSS_semiMinor/GNSS_semiMajor)*(GNSS_semiMinor/GNSS_semiMajor)));
    double b2= pow(GNSS_semiMinor,2);
    double e2= pow(e,2);
    double ep = e*(GNSS_semiMajor/GNSS_semiMinor);
    double r = sqrt(x2+y2);
    double r2 = pow(r,2);
    double E2 = pow(GNSS_semiMajor,2)- pow(GNSS_semiMinor,2);
    double F = 54*b2*z2;
    double G = r2 + (1-e2)*z2-e2*E2;
    double c = (pow(e2,2)*F*r2)/(pow(G,3));
    double s = pow((1.0 + c + sqrt(c*c +2*c )),1.0/3.0);

    double P = F /(3.0 * (s+1/s+1)*(s+1/s+1)* G*G);
    double Q = sqrt(1.0 + 2*e2*e2*P);
    double ro = -(P*e2*r)/(1+Q) + sqrt((GNSS_semiMajor*GNSS_semiMajor/2)*(1+1/Q)-(P*(1-e2)*z2)/(Q*(1+Q))-P*r2/2);
    double tmp = pow((r-e2*ro),2);
    double U = sqrt( tmp +z2 );
    double V = sqrt( tmp + (1-e2)*z2 );
    double zo = (b2*p1.data[2])/(GNSS_semiMajor*V);

    double height = U*(1.0 - b2/(GNSS_semiMajor*V));
    double lat  = std::atan(( p1.data[2] + ep*ep*zo)/r);
    double temp= std::atan(p1.data[1]/p1.data[0]);
    double longitude;
    if(p1.data[0] >= 0.0)
    {
        longitude= temp;
    }
    else if(( p1.data[0] <0.0 )&&(p1.data[1]>=0.0))
    {
        longitude = temp + M_PI;
    }
    else
    {
        longitude = temp - M_PI;

    }
    minivector llh(lat,longitude,height);
    return llh;
}

minivector inertialToECEF( const minivector& inertialPosition, const double t, const double t0)
{
    /*
     * inputs ::
     *  inertialPos -- > ECI position vector
     *  t --> current time [sec]
     *  t0 --> time when coordinate frames aligned [sec]
     * ouputs ::
     *  ecefPos --> ECEF position vector [meters]
     */
    double cT = cos(GNSS_earthRot*(t-t0));
    double sT = sin(GNSS_earthRot*(t-t0));
    //Matrix Rie = (Matrix(3,3) << cT, sT, 0, -sT, cT, 0, 0, 0, 1 ).finished();
    minimatrix Rie(3,3);
    Rie.data[0]=cT;
    Rie.data[1]=sT;
    Rie.data[2]=0;
    Rie.data[3]=-sT;
    Rie.data[4]=cT;
    Rie.data[5]=0;
    Rie.data[6]=0;
    Rie.data[7]=0;
    Rie.data[8]=1;
    //Rie<<cT, sT, 0, -sT, cT, 0, 0, 0, 1 ;
    // minivector ecefPosition= Rie*inertialPosition;
    minivector ecefPosition(3);
    miniblas_dgemv(blasNoTrans,1.0,Rie,inertialPosition,0.0,&ecefPosition);

    return ecefPosition;
}

minivector enu2xyz(const minivector& p1, const minivector& p2)
{
    /*
       inputs ::
       p1 --> enu coordinates [meter]
       p2 --> ECEF xyz origin coordinates [meter]
       outputs ::
        posXYZ ---> ECEF XYZ position vector [meters]
     */
    minivector orgLLH = xyz2llh(p2);
    double sinPhi = sin(orgLLH.data[0]);
    double cosPhi = cos(orgLLH.data[0]);
    double sinLam = sin(orgLLH.data[1]);
    double cosLam = cos(orgLLH.data[1]);

    minimatrix R(3,3);
    // R << -1*sinLam, cosLam, 0,
    //              -1*sinPhi*cosLam, -1*sinPhi*sinLam, cosPhi,
    //              cosPhi*cosLam, cosPhi*sinLam, sinPhi;
    R.data[0]=-1*sinLam;
    R.data[1]=cosLam;
    R.data[2]=0;
    R.data[3]=-1*sinPhi*cosLam;
    R.data[4]=-1*sinPhi*sinLam;
    R.data[5]=cosPhi;
    R.data[6]=cosPhi*cosLam;
    R.data[7]=cosPhi*sinLam;
    R.data[8]=sinPhi;

    minivector deltaXYZ(3);
    //deltaXYZ = R.inverse()*p1;
    minimatrix Rinv(3,3);
    minimatrix Rtemp(R);
    minilinalg_luMatInverse(&Rinv,&Rtemp);

    miniblas_dgemv(blasNoTrans,1.0,Rinv,p1,0.0,&deltaXYZ);


    minivector_add(&deltaXYZ,p2);
    return deltaXYZ;

    // return p2 + deltaXYZ;
}

minivector ned2enu(const minivector& p1)
{
    //enuConv<< 0, 1, 0, 1, 0, 0, 0, 0, -1;
    // enuConv.data[0]=0;enuConv.data[1]=1;enuConv.data[2]=0;
    // enuConv.data[3]=1;enuConv.data[4]=0;enuConv.data[5]=0;
    //  enuConv.data[6]=0;enuConv.data[7]=0;enuConv.data[8]=-1;
    //return enuConv*p1;
    minivector result(p1.data[1],p1.data[0],-p1.data[2]);
    return result;
}

minivector xyz2enu(const minivector &p1, const minivector &p2)
{
    /*
       inputs ::
       p1 --> ECEF xyz coordinates [meter]
       p2 --> ECEF xyz origin coordinates [meter]
          outputs ::
              posENU --> ENU position coordinates [meters]
     */

    minivector posDiff(3);
    minivector_sub(&posDiff,p1,p2);

    minivector orgLLH = xyz2llh(p2);
    double sinPhi = sin(orgLLH.data[0]);
    double cosPhi = cos(orgLLH.data[0]);
    double sinLam = sin(orgLLH.data[1]);
    double cosLam = cos(orgLLH.data[1]);

    minimatrix R(3,3);
//    R<< (-1*sinLam), cosLam, 0, ((-1*sinPhi)*cosLam), ((-1*sinPhi)*sinLam), cosPhi,
  //  (cosPhi*cosLam), (cosPhi*sinLam), sinPhi;
    R.data[0]=-sinLam;
    R.data[1]=cosLam;
    R.data[2]=0;
    R.data[3]=((-1*sinPhi)*cosLam);
    R.data[4]=((-1*sinPhi)*sinLam);
    R.data[5]=cosPhi;
    R.data[6]=(cosPhi*cosLam);
    R.data[7]=(cosPhi*sinLam);
    R.data[8]=sinPhi;

    minivector pos(3);
    // pos = R*posDiff;
    miniblas_dgemv(blasNoTrans,1.0,R,posDiff,0.0,&pos);

    return pos;
}

double calcElNed(const minivector& p1)
{
    /*
     * inputs ::
     *  p1 --> ECEF xyz location [meters]
     * outpus ::
     *  EL --> elevation angle [rad]
     */
    minivector posENU = ned2enu(p1);
    double dnrm2=miniblas_dnrm2(&posENU);
    double El = std::atan2(posENU.data[2], dnrm2);
    //double El = std::atan2(posENU.data[2], gsl_blas_sqrt_dnrm2(posENU));
    return El;
}

double calcEl(const minivector& p1, const minivector& p2)
{
    /*
       inputs ::
       p1 --> ECEF xyz satellite coordinates [meter]
       p2 ---> ECEF xyz receiver coordinates [meter]
       output ::
       El --> elevation angle [rad]
     */
    minivector posENU = xyz2enu(p1,p2);
    double dnrm2=miniblas_dnrm2(&posENU);
    double El = std::atan2(posENU.data[2], dnrm2);
    //double El = std::atan2(posENU.data[2], gsl_blas_sqrt_dnrm2(posENU));
    return El;
}

double tropMap(const double& El)
{
    /*
       inputs ::
       El --> receiver to satellite elevation angle [rad]
       output ::
       m --> troposphere delay map
     */
    double m = 1.001/sqrt(0.002001+( pow(sin(El), 2)));
    return m;
}

double tropDry(const minivector& p1)
{
    /*
       inputs ::
       p1 --> ECEF xyz receiver coordinated [meter]
       output ::
       tropDelta --> delay associated with troposphere [meter]
     */
    minivector posLLH = xyz2llh(p1);
    double recHeight = posLLH.data[2]*0.001;         // receiver height [km]
    double pSea = GNSS_stdPressure* pow( ( (GNSS_stdTemp-(6.5*recHeight))/(GNSS_stdTemp+6.5)), 5.2459587 );
    double tropDelta =  ( (1e-6/5) * ( (77.624*(pSea/(GNSS_stdTemp+6.5*recHeight))) * (40136+148.72*((GNSS_stdTemp-(6.5*recHeight))- 288.16)) ) );
    return tropDelta;
}

double dopplerObs(double r1, double r2)
{
    /*
       inputs ::
       r1 --> carrier phase obs. timestep n
       r2 --> carrier phase obs. timestep n-1
       output ::
       tdcp --> time differenced carrier phase observable
     */
    double tdcp = r1 - r2;
    return tdcp;
}

double elDepWeight(const minivector& p1, const minivector& p2, double measWeight)
{
    /*
       inputs ::
       p1 --> ECEF xyz satellite coordinates [meter]
       p2 ---> ECEF xyz receiver coordinates [meter]
                          measNoise ---> initial noise applied to observable [meters]
       output ::
       r --> elevation angle dep. GNSS obs. weigth [rad]
     */
    double el = calcEl(p1,p2);
    double r = (1 / ( sin(el) )) * measWeight;
    return r;
}

minivector satVelocity(const minivector& p1, const minivector& p2, double ts)
{
    /*
       inputs ::
       p1 --> sat xyz at timestep n
       p2 --> sat xyz at timestep n-1
       ts --> timestep
       outputs ::
       satVel --> sat velocity
     */

    //minivector satVel = ( (p1 - p2)/ (2*ts) );
    minivector satVel(3);
    minivector_sub(&satVel,p1,p2);
    minivector_scale(&satVel,0.5/ts);

    return satVel;
}
