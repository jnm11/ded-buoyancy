#!/usr/bin/python
""" 
PID Controller. see http://en.wikipedia.org/wiki/PID_controller
"""

import numpy as np 

class pid:
    """PID Controller
    """

    # fix initial conditions to exactly give the current velocity as output
    def __init__(self, Kp, Ki, Kd, Kt, XT, IT, DD, U, t, X, minU, maxU, absU, sclX):
        
        e = sclX*(X-XT)

        if Kp is None: Kp=0
        if Ki is None: Ki=0
        if Kd is None: Kd=0
        if Kt is None: Kt=0
        self.t      = np.float64(t)
        self.e      = np.float64(e)
        self.Kp     = np.float64(Kp)
        self.Ki     = np.float64(Ki)
        self.Kd     = np.float64(Kd)
        self.Kt     = np.float64(Kt)       # Low pass filter on output
        self.sclX   = np.float64(sclX)
        self.minU   = np.float64(minU)
        self.maxU   = np.float64(maxU)
        self.absU   = np.float64(absU)
        self.XT     = np.float64(XT)    # setpoint 
        self.ITerm  = np.float64(IT)
        self.PTerm  = np.float64(Kp*e)
        self.DTerm  = np.float64(DD)
        self.output = np.float64(U)
        self.print()

    def update_XT(self, XT):
        self.XT     = np.float64(XT)    # setpoint
                
    def print(self):
        print('    t      {} '.format(self.t))
        print('    e      {} '.format(self.e))
        print('    Kp     {} '.format(self.Kp))
        print('    Ki     {} '.format(self.Ki))
        print('    Kd     {} '.format(self.Kd))
        print('    Kt     {} '.format(self.Kt))
        print('    XT     {} '.format(self.XT))
        print('    sclX   {} '.format(self.sclX))
        print('    minU   {} '.format(self.minU))
        print('    maxU   {} '.format(self.maxU))
        print('    absU   {} '.format(self.absU))
        print('    ITerm  {} '.format(self.ITerm))
        print('    PTerm  {} '.format(self.PTerm))
        print('    DTerm  {} '.format(self.DTerm))
        print('    output {} '.format(self.output))

    def update(self, X, t):
        """Calculates PID value for given reference feedback
        
        .. math::
        u(t) = K_p [ e(t) + K_i \int_{0}^{t} e(t)dt + K_d {de}/{dt} ]
        
        .. figure:: images/pid_1.png
           :align:   center

           Test PID with Kp=1.2, Ki=1, Kd=0.001 (test_pid.py)

        """
        if not np.isfinite(X): return(self.output,self.ITerm,self.DTerm)
        
        dt = t - self.t
        
        if (dt > 0):
            e  = self.sclX*(X-self.XT)
            de = e - self.e
            
            DTerm=self.Kd*de/dt
            
            self.t = t
            self.e = e

            self.ITerm += self.Ki * e * dt
            self.DTerm = np.exp(-dt)*(self.DTerm-DTerm)+DTerm
            self.PTerm = self.Kp * e
            
            ITerm=min(self.maxU,max( self.minU,self.ITerm))
            DTerm=min(self.absU,max(-self.absU,self.DTerm))
            PTerm=min(self.absU,max(-self.absU,self.PTerm))
            if np.isfinite(ITerm): self.ITerm=ITerm
            if np.isfinite(DTerm): self.DTerm=DTerm
            if np.isfinite(PTerm): self.PTerm=PTerm
            f=min(self.maxU,max(self.minU,self.PTerm + self.ITerm + self.DTerm))
            self.output = (dt*f + self.Kt*self.output)/(dt+self.Kt)
        return(self.output,self.ITerm,self.DTerm)


  
 
