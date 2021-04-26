import math

""" computations based on 
Pierre Billoir Nucl.Instrum.Meth. A225 (1984) 352-366
information matrix is the inverse of the covariance matrix
"""

class PFilter:
    def __init__(self, erry, errs, pbeta):
        self.erry = erry
        self.errs = errs
        self.pbeta = pbeta
        self.opt   = [0.,0.]
        self.info  = [0., 0., 0.]
        self.chi2  = 0.

    def reset(self):
        self.opt   = [0.,0.]
        self.info  = [0., 0., 0.]
        self.chi2  = 0.

    def initiate(self,y,a):
        self.reset()
        self.opt[0] = y
        self.opt[1] = a
        self.info[0] = 1./pow(2*self.erry,2)
        self.info[1] = 2. if a == 0. else 1./pow(0.5*a,2)
        self.info[2] = 0. # off-diagonal element of the 2x2 matrix

    def det(self,matrix):
        return matrix[0]*matrix[1] - matrix[2]*matrix[2]

    def invert(self, matrix):
        det = self.det(matrix)
        inverse = [0., 0., 0.]
        inverse[0] = matrix[1]/det
        inverse[1] = matrix[0]/det
        inverse[2] = -matrix[2]/det
        return inverse


    def getY(self):
        return self.opt[0]

    def getSlope(self):
        return self.opt[1]

    def getYerr(self):
        det = self.det(self.info)
        return math.sqrt(self.info[1]/det)

    def getSlopeErr(self):
        det = self.det(self.info)
        return 0. if det == 0. else math.sqrt(math.fabs(self.info[0]/det))
        
    def getCorr(self):
        erry = self.getYerr()
        errs = self.getSlopeErr()
        det = self.det(self.info)
        return (self.info[2]/det)/(erry*errs)
        
    def predict(self, step):
        return self.opt[0] + step * self.opt[1]

    def delta_y(self, ymeas, step):
        return math.fabs(self.predict(step)-ymeas)

    def computeChi2(self, ymeas, step):
        """ this actually is a fancy way to computer |ypred-ymeas| as errcov->0"""
        res = self.predict(step)-ymeas
        det = self.det(self.info)

        if(det != 0):
            errcov = self.info[1]/det
        else:
            print("det of I is 0")
            errcov=10.
        return res*res/(pow(self.erry,2)+errcov)

    def chi2_if_update(self, ymeas, step):
        """ what would be the chi2 if this point was added to the filter """
        ypred = self.predict(step)
        apred = self.opt[1]

        cov   = self.multScatt(step)
        info_inv = self.invert(self.info)
        info_inv = [i+a for i, a in zip(info_inv, cov)]
        info = self.invert(info_inv)

        info[2] = -info[0]*step + info[2]
        info[1] =  info[0]*pow(step,2) - 2.*info[2]*step + info[1]

        M = [1./pow(self.erry, 2), 0., 0.]
        i_m = [i+m for i,m in zip(info, M)]

        det = self.det(i_m)

        ymeas_err = ymeas/(self.erry*self.erry)

        yopt = (info[1]/det)*(ymeas_err + info[0]*ypred + info[2]*apred) - (info[2]/det)*(info[2]*ypred + info[1]*apred)
        aopt = -(info[2]/det)*(ymeas_err + info[0]*ypred + info[2]*apred) + (i_m[0])/det * (info[2]*ypred + info[1]*apred)
        
        ydelta = yopt - ypred
        adelta = aopt - apred

        """ chi2 is (opt-pred).T x I x (opt-pred) """
        chi2meas = i_m[0]*pow(ydelta,2) + info[1]*pow(adelta,2) + 2.*info[2]*ydelta*adelta

        return chi2meas




    def update(self, ymeas, step):
        """ update estimators currently for point n to point n+1 """

        ypred = self.predict(step)
        apred = self.opt[1]

        """1. Scattering : I*[n] = (I[n]^-1 + A[n])^-1 """
        """ <-> add the MS covariance mtx (A) to the information matrix """

        cov   = self.multScatt(step)
        info_inv = self.invert(self.info)
        info_inv = [i+a for i, a in zip(info_inv, cov)]
        self.info = self.invert(info_inv)
        
        
        """2. Propagate I*[n+1] = D[n].T^-1 x I*[n] x D[n]^-1 """
        """D[n] is the propagation matrix """
        self.info[2] = -self.info[0]*step + self.info[2]
        self.info[1] =  self.info[0]*pow(step,2) - 2.*self.info[2]*step + self.info[1]
        
        """3. Measurement I[n+1] = I*[n+1] + M[n] """
        """ M: measurement error M(0,0)=1./(erry*erry) """
        M = [1./pow(self.erry, 2), 0., 0.]
        i_m = [i+m for i,m in zip(self.info, M)]

        """4. Get new estimators (opt') given measurements """
        """ solve (I+M)*(opt'-pred)=M(meas-pred)"""
        det = self.det(i_m)

        ymeas_err = ymeas/(self.erry*self.erry)

        yopt = (self.info[1]/det)*(ymeas_err + self.info[0]*ypred + self.info[2]*apred) - (self.info[2]/det)*(self.info[2]*ypred + self.info[1]*apred)
        aopt = -(self.info[2]/det)*(ymeas_err + self.info[0]*ypred + self.info[2]*apred) + (i_m[0])/det * (self.info[2]*ypred + self.info[1]*apred)

        self.opt[0] = yopt
        self.opt[1] = aopt

        """add measurement error to the info matrix"""
        self.info[0] += 1./(self.erry*self.erry)
        
        ydelta = yopt - ypred
        adelta = aopt - apred

        """ chi2 is (opt-pred).T x I x (opt-pred) """
        chi2meas = self.info[0]*pow(ydelta,2) + self.info[1]*pow(adelta,2) + 2.*self.info[2]*ydelta*adelta
        self.chi2 += chi2meas
        return chi2meas

    def multScatt(self, step):
        """ from 
        The Kalman Filter Technique applied to Track Fitting in GLAST
        by Jose Hernando 
        """

        cov = [0., 0., 0.]
        if(self.pbeta == 0. or step == 0.):
            return cov

        X0 = 14. #LAr radiation length in cm
        theta_ms_square = pow(0.0136/self.pbeta, 2) * math.fabs(step)/X0

        ### this takes into account the projected incident angle
        incFac = pow( 1. + pow(self.opt[1], 2), 2.5)
        if(math.isnan(incFac)) : 
            incFac = 1.
        err = theta_ms_square * incFac


        cov[0] = err * step * step / 3.
        cov[1] = err
        cov[2] = err * math.fabs(step) / 2.
        
        #this is the naive case with no correlations
        """
        cov[0] = 0
        cov[1] = err
        cov[2] = 0.
        """
        return cov
        
    def getChi2(self):
        return self.chi2
        
