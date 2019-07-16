/**
 * @file   GnssTools.cpp
 * @brief  Tools required to process GNSS data -- (i.e. ECEF to ENU transformation)
 * @author Ryan Watson & Jason Gross
 */

#include "GnssTools.h"


Eigen::VectorXd obsMap(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const int &Trop) {
        /*
           inputs ::
           p1 --> ECEF xyz coordinates of satellite [meter]
           p2 --> ECEF xyz coordinates of receiver [meter]
           Trop --> Troposphere modeling switch
           outputs ::
           H --> measurement mapping
         */
        double r = sqrt( ((p1(0)-p2(0)))*(p1(0)-p2(0)) + ((p1(1)-p2(1))*(p1(1)-p2(1))) + ((p1(2)-p2(2))*(p1(2)-p2(2))) );
        if (Trop == 1) {
                double el = calcEl(p1,p2);
                double mapT = tropMap(el);
                Eigen::VectorXd H(5); H << (p1(0)-p2(0))/r, (p1(1)- p2(1))/r,
                (p1(2)-p2(2))/r, 1.0,mapT;
                return H;
        }
        else {
                Eigen::VectorXd H(5); H << (p1(0)-p2(0))/r, (p1(1)- p2(1))/r,
                (p1(2)-p2(2))/r, 1.0, 0.0;
                return H;
        }
}

Eigen::VectorXd  obsMapNED(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const int &Trop) {
        /*
           inputs ::
           p1 --> ECEF xyz coordinates of satellite [meter]
           p2 --> ECEF xyz coordinates of receiver [meter]
           Trop --> Troposphere modeling switch
           outputs ::
           H --> measurement mapping
         */
        double r = sqrt( ((p1(0)-p2(0)))*(p1(0)-p2(0)) + ((p1(1)-p2(1))*(p1(1)-p2(1))) + ((p1(2)-p2(2))*(p1(2)-p2(2))) );
        if (Trop == 1) {
                double el = calcElNed(p1);
                double mapT = tropMap(el);
                Eigen::VectorXd H(5); H << (p1(0)-p2(0))/r, (p1(1)- p2(1))/r,
                (p1(2)-p2(2))/r, 1.0,mapT;
                return H;
        }
        else {
                Eigen::VectorXd H(5); H << (p1(0)-p2(0))/r, (p1(1)- p2(1))/r,
                (p1(2)-p2(2))/r, 1.0, 0.0;
                return H;
        }
}

Eigen::VectorXi getPRN(const Eigen::MatrixXd &p) {
        /*
           inputs ::
           p --> GNSS data packed in pre-specified format
           outputs ::
           prnVec --> Vector of prn numbers corresponding to visible satellites at given epoch
         */

        int nSats = static_cast<int>(p(0));
        Eigen::VectorXi prnVec(nSats);
        int count = 2;
        for (int i=0; i<nSats; i++) {
                prnVec(i) = static_cast<int>(p(count));
                count = count+7;
        }
        return prnVec;
}

bool checkPRN(const Eigen::VectorXi &p, const int &n) {
        /*
           input ::
           p --> prn vector
           n --> prn value to check
           output ::
           true --> if n is present in current prn vec.
           false --> otherwise
         */

        for (int i=0; i<p.size(); i++) {
                if (p(i) == n ) {return true; }
        }
        return false;
}

Eigen::MatrixXd earthToNavTrans( const Eigen::Vector3d &p1){
        /*
           inputs ::
            p1 --> ECEF xyz coordinates of the platform
           outputs ::
            R --> rotation matrix from Earth fram to Navigation frame
         */
        Eigen::Vector3d llh=xyz2llh(p1);
        double sLat = sin(llh(0));
        double sLon = sin(llh(1));
        double cLat = cos(llh(0));
        double cLon = cos(llh(1));

        //Matrix R = (Matrix(3,3) << -1*sLat*cLon, -1*sLon, -1*cLat*cLon,
        //            -1*sLat*sLon, cLon, -1*cLat*sLon,
       //             cLat, 0.0, -1*sLat ).finished();
         Eigen::MatrixXd R(3,3);
         R<< -1*sLat*cLon, -1*sLon, -1*cLat*cLon,
                    -1*sLat*sLon, cLon, -1*cLat*sLon,
                    cLat, 0.0, -1*sLat;

        return R.inverse();

}

double deltaObs(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const double &meas){
        /*
           inputs ::
           p1 --> ECEF xyz coordinates of satellite [meter]
           p2 --> ECEF xyz coordinates of receiver [meter]
           meas --> observed range between satellite and receiver [meter]
           outputs ::
           deltaR --> difference observed pseudorange and modeled pseudorange [meter]
         */
        double r = minisam::norm3d(p1-p2) + tropMap(calcEl(p1,p2))*tropDry(p2);
        return meas - r;
}

double deltaTrop(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2){
        /*
           inputs ::
           p1 --> ECEF xyz coordinates of satellite [meter]
           p2 --> ECEF xyz coordinates of receiver [meter]
           meas --> observed range between satellite and receiver [meter]
         */
        return tropMap(calcEl(p1,p2))*tropDry(p2);
}

Eigen::Vector3d xyz2llh(const Eigen::Vector3d &p1){
        /*
           inputs ::
           p1 --> ECEF xyz receiver coordinates [meter]
           output ::
           posLLH --> latitude, longitude, height [rad,rad,meter]
         */
        double x2= pow(p1(0),2);
        double y2= pow(p1(1),2);
        double z2= pow(p1(2),2);

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
        double zo = (b2*p1(2))/(GNSS_semiMajor*V);

        double height = U*(1.0 - b2/(GNSS_semiMajor*V));
        double lat  = std::atan(( p1(2) + ep*ep*zo)/r);
        double temp= std::atan(p1(1)/p1(0));
        double longitude;
        if(p1(0) >= 0.0) {
                longitude= temp;
        }else if(( p1(0) <0.0 )&&(p1(1)>=0.0)) {
                longitude = temp + M_PI;
        } else{
                longitude = temp - M_PI;

        }
        Eigen::Vector3d llh(lat,longitude,height);
        return llh;
}

Eigen::Vector3d inertialToECEF( const Eigen::Vector3d& inertialPosition, const double t, const double t0){
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
        Eigen::MatrixXd Rie(3,3);
        Rie<<cT, sT, 0, -sT, cT, 0, 0, 0, 1 ;
        Eigen::Vector3d ecefPosition= Rie*inertialPosition;
        return ecefPosition;
}

Eigen::Vector3d enu2xyz(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
        /*
           inputs ::
           p1 --> enu coordinates [meter]
           p2 --> ECEF xyz origin coordinates [meter]
           outputs ::
            posXYZ ---> ECEF XYZ position vector [meters]
         */
        Eigen::Vector3d orgLLH = xyz2llh(p2);
        double sinPhi = sin(orgLLH(0));
        double cosPhi = cos(orgLLH(0));
        double sinLam = sin(orgLLH(1));
        double cosLam = cos(orgLLH(1));

         Eigen::MatrixXd R(3,3);
         R << -1*sinLam, cosLam, 0,
                     -1*sinPhi*cosLam, -1*sinPhi*sinLam, cosPhi,
                     cosPhi*cosLam, cosPhi*sinLam, sinPhi;

        Eigen::Vector3d deltaXYZ;
        deltaXYZ = R.inverse()*p1;
        return p2 + deltaXYZ;
}

Eigen::Vector3d ned2enu(const Eigen::Vector3d& p1) {

        Eigen::MatrixXd enuConv(3,3);
         enuConv<< 0, 1, 0, 1, 0, 0, 0, 0, -1;
        return enuConv*p1;
}

Eigen::Vector3d xyz2enu(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2){
        /*
           inputs ::
           p1 --> ECEF xyz coordinates [meter]
           p2 --> ECEF xyz origin coordinates [meter]
              outputs ::
                  posENU --> ENU position coordinates [meters]
         */

        Eigen::Vector3d posDiff = p1 - p2;
        Eigen::Vector3d orgLLH = xyz2llh(p2);
        double sinPhi = sin(orgLLH(0));
        double cosPhi = cos(orgLLH(0));
        double sinLam = sin(orgLLH(1));
        double cosLam = cos(orgLLH(1));

        Eigen::MatrixXd R(3,3);
        R<< (-1*sinLam), cosLam, 0, ((-1*sinPhi)*cosLam), ((-1*sinPhi)*sinLam), cosPhi, (cosPhi*cosLam), (cosPhi*sinLam), sinPhi;
        Eigen::Vector3d pos;
        pos = R*posDiff;
        Eigen::Vector3d posENU;
        return posENU = Eigen::Vector3d(pos(0), pos(1), pos(2));
}

double calcElNed(const Eigen::Vector3d& p1){
        /*
         * inputs ::
         *  p1 --> ECEF xyz location [meters]
         * outpus ::
         *  EL --> elevation angle [rad]
         */
        Eigen::Vector3d posENU = ned2enu(p1);
        double El = std::atan2(posENU(2), posENU.norm());
        return El;
}

double calcEl(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2){
        /*
           inputs ::
           p1 --> ECEF xyz satellite coordinates [meter]
           p2 ---> ECEF xyz receiver coordinates [meter]
           output ::
           El --> elevation angle [rad]
         */
        Eigen::Vector3d posENU = xyz2enu(p1,p2);
        double El = std::atan2(posENU(2), posENU.norm());
        return El;
}

double tropMap(const double& El){
        /*
           inputs ::
           El --> receiver to satellite elevation angle [rad]
           output ::
           m --> troposphere delay map
         */
        double m = 1.001/sqrt(0.002001+( pow(sin(El), 2)));
        return m;
}

double tropDry(const Eigen::Vector3d& p1){
        /*
           inputs ::
           p1 --> ECEF xyz receiver coordinated [meter]
           output ::
           tropDelta --> delay associated with troposphere [meter]
         */
        Eigen::Vector3d posLLH = xyz2llh(p1);
        double recHeight = posLLH(2)/1000;         // receiver height [km]
        double pSea = GNSS_stdPressure* pow( ( (GNSS_stdTemp-(6.5*recHeight))/(GNSS_stdTemp+6.5)), 5.2459587 );
        double tropDelta =  ( (1e-6/5) * ( (77.624*(pSea/(GNSS_stdTemp+6.5*recHeight))) * (40136+148.72*((GNSS_stdTemp-(6.5*recHeight))- 288.16)) ) );
        return tropDelta;
}

double dopplerObs(double r1, double r2){
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

double elDepWeight(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, double measWeight) {
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

Eigen::Vector3d satVelocity(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, double ts){
        /*
           inputs ::
           p1 --> sat xyz at timestep n
           p2 --> sat xyz at timestep n-1
           ts --> timestep
           outputs ::
           satVel --> sat velocity
         */

        Eigen::Vector3d satVel = ( (p1 - p2)/ (2*ts) );
        return satVel;
}


