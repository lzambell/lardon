import math

class PFilter:
    def __init__(self, erry, errs, pbeta):
        self.erry = erry
        self.errs = errs
        self.pbeta = pbeta
        self.opt   = [0.,0.]
        self.info  = [0., 0., 0.]
        self.chi2  = 0.

    def reset(self):
        #self.erry = 0.32
        #self.errs = 0.05
        #self.pbeta = 1.
        self.opt   = [0.,0.]
        self.info  = [0., 0., 0.]
        self.chi2  = 0.

    def initiate(self,y,a):
        self.reset()
        self.opt[0] = y
        self.opt[1] = a
        self.info[0] = 1./pow(2*self.erry,2)
        self.info[1] = 2. if a == 0. else 1./pow(0.5*a,2)
        self.info[2] = 0.

    def det(self,matrix):
        return matrix[0]*matrix[1] - matrix[2]*matrix[2]
        
    def getY(self):
        return self.opt[0]

    def getSlope(self):
        return self.opt[1]

    def getYerr(self):
        det = self.det(self.info)
        return math.sqrt(self.info[1]/det)

    def getSlopeErr(self):
        det = self.det(self.info)
        return math.sqrt(self.info[0]/det)
        
    def getCorr(self):
        erry = self.getYerr()
        errs = self.getSlopeErr()
        det = self.det(self.info)
        return (self.info[2]/det)/(erry*errs)
        
    def predict(self, step):
        return self.opt[0] + step * self.opt[1]

    def computeChi2(self, ymeas, step):
        res = self.predict(step)-ymeas
        det = self.det(self.info)
        #self.info[0]*self.info[1]-self.info[2]*self.info[2]
        if(det != 0):
            errcov = self.info[1]/det
        else:
            errcov=10.
        return res*res/(pow(self.erry,2)+errcov)

    def update(self, ymeas, step):
        ypred = self.predict(step)
        apred = self.opt[1]
        cov   = self.multScatt(step)
        det_cov = self.det(cov)
        det_info = self.det(self.info)
        if(det_info==0):
            return 999.
        det_cross = 1./det_info + self.info[1] * cov[1]/det_info - 2.* self.info[2]*cov[2]/det_info
        det_prod = det_info * det_cross
        self.info[0] = self.info[0]/det_prod + cov[1]/det_cross
        self.info[1] = self.info[1]/det_prod + cov[0]/det_cross
        self.info[2] = self.info[2]/det_prod + cov[2]/det_cross
        
        """update the info matrix"""
        """info = propagation.T^-1 * info * propagation^-1"""
        self.info[2] = -self.info[0]*step + self.info[2]
        self.info[1] = self.info[0]*pow(step,2) - 2.*self.info[2]*step + self.info[1]
        
        """get opt', new best estimator given the measurement """
        """ solve (I+M)*(opt'-pred)=M(meas-pred)"""
        """ I: info, M: measurement error M(0,0)=1./(erry*erry)"""
        
        det = (self.info[0] + 1./(self.erry*self.erry))*self.info[1] - self.info[2]*self.info[2]
        ymeas_err = ymeas/(self.erry*self.erry)

        yopt = (self.info[1]/det)*(ymeas_err + self.info[0]*ypred + self.info[2]*apred) - (self.info[2]/det)*(self.info[2]*ypred + self.info[1]*apred)
        aopt = -(self.info[2]/det)*(ymeas_err + self.info[0]*ypred + self.info[2]*apred) + (self.info[0]+1./(self.erry*self.erry))/det * (self.info[2]*ypred + self.info[1]*apred)

        self.opt[0] = yopt
        self.opt[1] = aopt

        """add measurement error to the info matrix"""
        self.info[0] += 1./(self.erry*self.erry)
        
        ydelta = yopt - ypred
        adelta = aopt - apred
        chi2meas = self.info[0]*pow(ydelta,2) + self.info[1]*pow(adelta,2) + self.info[2]*ydelta*adelta
        self.chi2 += chi2meas
        return chi2meas

    def multScatt(self, step):
        cov = [0., 0., 0.]
        if(self.pbeta == 0. or step == 0.):
            return cov

        X0 = 14. #LAr radiation length in cm
        theta_ms_square = pow(0.0136/self.pbeta, 2) * math.fabs(step)/X0
        incFac = pow( 1. + pow(self.opt[1], 2), 2.5)
        if(math.isnan(incFac)) : 
            incFac = 1.
        err = theta_ms_square * incFac
        cov[0] = err * step * step
        cov[1] = err
        cov[2] = err * math.fabs(step)
        return cov
        
    def getChi2(self):
        return self.chi2
        
