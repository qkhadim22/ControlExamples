#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Details:  It includes PD controller on spring-damper, single pendulum, double 
            # pendulum and lift boom.
#
# Author:   Qasim Khadim
# Contact: qasim.khadim@outlook.com,qkhadim22 (Github)
# Date:     2025-02-02
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

from Simulation import Models


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                         #MODEL SETTINGS#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

timeStep    = 1e-3                  #Simulation time step: Change it as desired.
T           = 10                    #Time period
ns          = int(T/timeStep)       #Number of steps

# Initial conditions for all models
theta1      = 14.2414062401257      #For single pendulum and liftboom
theta2      = -58.5351564492275     #For Two pendulums
dtheta1     = 0                     #For single pendulum and liftboom
dtheta2     = 0                     #For Two pendulums
s           = 800e-3                #For spring-damper
dots        = -0.2                  #For spring-damper

Vec0        = [s,dots,theta1,theta2,dtheta1,dtheta2 ]

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 'SpringDamper','Pendulum','TwoPendulum'[Kp=5000*100,Kd=0.5*5000*100], 'LiftBoom'[Kp=7e3,Kd=20]

model       = Models(nStepsTotal=ns, endTime=T, verboseMode=1,Kp=5e5,Kd=5e5*0.5,
                         Ki=0.02*5000*100, Case='SpringDamper')
inputVec    = model.CreateInputVector(Vec0)

data        = model.ComputeModel(inputVec, solutionViewer = True) #solutionViewer: for visualization

