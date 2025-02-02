#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Details:  It includes PD controller on spring-damper, single pendulum, double 
            # pendulum and lift boom.
#
# Author:   Qasim Khadim
# Contact: qasim.khadim@outlook.com,qkhadim22 (Github)
# Date:     2025-02-02
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import exudyn as exu
from exudyn.utilities import * #includes itemInterface and rigidBodyUtilities
import exudyn.graphics as graphics #only import if it does not conflict
from math import sin, cos, pi
import math as mt
from exudyn.utilities import SensorBody, SensorObject


from ControllerLibrary import Control


fileName1       = 'LowQualityPATU/Pillar.stl'
fileName2       = 'LowQualityPATU/LiftBoom.stl'

class Models():
  #initialize class 
  def __init__(self, nStepsTotal=100, endTime=0.5, Kp=10, Kd=1, Ki=0.20,
               Case = False,
               visualization = False,
               verboseMode = 0):

      self.nStepsTotal  = nStepsTotal
      self.endTime      = endTime
      self.Kp           = Kp
      self.Kd           = Kd
      self.Ki           = Ki
      self.Case         = Case
      
      self.timeVecOut   = np.arange(1,self.nStepsTotal+1)/self.nStepsTotal*self.endTime
      
  def CreateModel(self):
      self.SC  = exu.SystemContainer()
      self.mbs = self.SC.AddSystem()
      
  def GetOutputXAxisVector(self):
        return self.timeVecOut

  #create a randomized input vector
  #relCnt can be used to create different kinds of input vectors (sinoid, noise, ...)
  #isTest is True in case of test data creation
  def CreateInputVector(self, Vec0):

    if self.Case == 'TwoPendulum':
          vec   = np.zeros(4*self.nStepsTotal)
          vec[0]                     = Vec0[2]     
          vec[1*(self.nStepsTotal)]  = Vec0[3] 
          vec[2*(self.nStepsTotal)]  = Vec0[4] 
          vec[3*(self.nStepsTotal)]  = Vec0[5]         
    else:
          vec   = np.zeros(2*self.nStepsTotal)
          
          if self.Case == 'SpringDamper':
                 vec[0]                     = Vec0[0]     
                 vec[1*(self.nStepsTotal)]  = Vec0[1] 
          else:
                 vec[0]                     = Vec0[2]     
                 vec[1*(self.nStepsTotal)]  = Vec0[3] 
    return vec

  def SplitInputData(self, inputData):
        rv = {}
        if self.Case == 'TwoPendulum':
              rv['theta1']   = inputData[0*(self.nStepsTotal):1*(self.nStepsTotal)] 
              rv['theta2']   = inputData[1*(self.nStepsTotal):2*(self.nStepsTotal)]  
              rv['dtheta1']  = inputData[2*(self.nStepsTotal):3*(self.nStepsTotal)]  
              rv['dtheta2']  = inputData[3*(self.nStepsTotal):4*(self.nStepsTotal)]       
              
        else:
              vec   = np.zeros(2*self.nStepsTotal)
              
              if self.Case == 'SpringDamper':
                     rv['s1']         = inputData[0*(self.nStepsTotal):1*(self.nStepsTotal)]      
                     rv['dots1']      = inputData[1*(self.nStepsTotal):2*(self.nStepsTotal)] 
              else:
                     rv['theta1']     = inputData[0*(self.nStepsTotal):1*(self.nStepsTotal)]      
                     rv['dtheta1']    = inputData[1*(self.nStepsTotal):2*(self.nStepsTotal)]  
                     
        return rv

            
    #get number of simulation steps
  def GetNSimulationSteps(self):
        return self.nStepsTotal # x finer simulation than output
       
       
  #initialState contains position and velocity states as list of two np.arrays 
  def ComputeModel(self, inputData, verboseMode = 0, solutionViewer = False):
        self.CreateModel()
        # print('compute model')
        self.verboseMode = verboseMode
        output = []
        #set input data ...
        inputDict = self.SplitInputData(np.array(inputData))
        
        
        if self.Case == 'TwoPendulum':
              
              self.mbs.variables['theta1']  = inputData[self.nStepsTotal*0]
              self.mbs.variables['theta2']  = inputData[self.nStepsTotal*1]
              self.mbs.variables['dtheta1'] = inputData[self.nStepsTotal*2]
              self.mbs.variables['dtheta2'] = inputData[self.nStepsTotal*3]
              
              self.DoublePendulum(self.mbs.variables['theta1'],self.mbs.variables['theta2'],
                                  self.mbs.variables['dtheta1'], self.mbs.variables['dtheta2'])
        else:              
              if self.Case == 'SpringDamper':
                  self.mbs.variables['s1']      = inputData[self.nStepsTotal*0]
                  self.mbs.variables['dots1']   = inputData[self.nStepsTotal*1]
                  
                  self.SpringDamper(self.mbs.variables['s1'], self.mbs.variables['dots1'])
              else:
                     self.mbs.variables['theta1']      = inputData[self.nStepsTotal*0]
                     self.mbs.variables['dtheta1']     = inputData[self.nStepsTotal*1]
                     
                     if self.Case == 'Pendulum':
                             self.SinglePendulum(self.mbs.variables['theta1'], self.mbs.variables['dtheta1'])
                     else:
                             self.LiftBoom(self.mbs.variables['theta1'], self.mbs.variables['dtheta1'])
        
            
        return [inputDict] 
    



  def SpringDamper(self, x, dotx):

        self.x           = x
        self.dotx        = dotx
        self.dictSensors = {}
        
        oMass            = self.mbs.CreateMassPoint(referencePosition=[0,1.355,0],
                            initialDisplacement = [0,self.x,0],
                            initialVelocity= [0,self.dotx,0],
                            gravity = [0,-9.806,0],
                            physicsMass=1.6) 

        oGround = self.mbs.CreateGround()
        Controller          = Control(self.Kp,self.Ki,self.Kd)

        #create spring damper with reference length computed from reference positions (=L)
        #load user function, second joint:
        def ControlForce(mbs, t, itemIndex, u, v, k, d, f0):
           omega = 0.5
           phiDesired = pi*0.5*sin(omega*pi*t)
           dt = 1e-3
           P = Controller.P(phiDesired, u, dt)
           phiDesired_t = pi*0.5*omega*cos(omega*pi*t)

           D = Controller.D(phiDesired_t, v, dt)
           Force = P + D
           
           #torque = Control.PID(phiDesired, phi0, dt)
           
           return Force 
        
        Force         = self.mbs.CreateSpringDamper(bodyOrNodeList=[oMass, oGround],stiffness=0,damping=0,force=0,
                                                springForceUserFunction=ControlForce)
        
        # Force = self.mbs.AddLoad(LoadCoordinate(markerNumber = mRotNode0,
        #                            load = 1, loadUserFunction = ControlForce))

        #add load on body:
        #self.mbs.CreateForce(bodyNumber = oMass, loadVector = [0,-80,0])

        # add time integration scheme:
        self.mbs.Assemble()
        
        #add dependencies of load user functions on nodal coordinates:
        ltgN0 = self.mbs.systemData.GetNodeLTGODE2(oMass)
        
        #both loads depend on both bodies; this information needs to be added in order that
        #  the solver knows dependencies when computing Jacobians
        
        #self.mbs.systemData.AddODE2LoadDependencies(Force, list(ltgN0))
        
        self.simulationSettings = exu.SimulationSettings()   
        self.simulationSettings.timeIntegration.numberOfSteps  = self.GetNSimulationSteps()
        self.simulationSettings.timeIntegration.endTime        = self.endTime

        self.simulationSettings.timeIntegration.newton.numericalDifferentiation.doSystemWideDifferentiation = True
        self.simulationSettings.timeIntegration.computeLoadsJacobian = 2
        # simulationSettings.timeIntegration.newton.useModifiedNewton = True

        self.simulationSettings.solutionSettings.writeSolutionToFile = True
        self.simulationSettings.solutionSettings.solutionWritePeriod = 0.01

        self.simulationSettings.timeIntegration.verboseMode = 1
        self.simulationSettings.displayStatistics = True
        self.simulationSettings.displayComputationTime = True

        self.mbs.SolveDynamic(simulationSettings=self.simulationSettings)          
        self.mbs.SolutionViewer()
                                         
  def SinglePendulum(self, theta, dtheta):

         self.theta         = theta
         self.dtheta        = dtheta
         self.dictSensors   = {}
         
         oGround            = self.mbs.AddObject(ObjectGround())
         mGround            = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, localPosition=[-0.5, 0, 0]))
         
         # create nodes:
         n0                 = self.mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[0,0,0], initialVelocities=[0,0,0]))

         # create bodies:
         b0                 = self.mbs.AddObject(RigidBody2D(physicsMass=100, physicsInertia=1,nodeNumber=n0,
                                                              visualization=VRigidBody2D(graphicsData
                                                            =[graphics.Lines([[-0.5,0,0],[0.5,0,0]])])))

         # add markers and loads:
         m00                = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b0, localPosition=[-0.5, 0., 0.]))
         m01                = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b0, localPosition=[ 0.5, 0., 0.]))

         # add joints:
         jg0                = self.mbs.AddObject(RevoluteJoint2D(markerNumbers=[mGround,m00]))
          
         # add loads:
         mLoad0             = self.mbs.AddMarker(MarkerNodePosition(nodeNumber=n0))
         self.mbs.AddLoad(Force(markerNumber=mLoad0, loadVector=[0,-9.81*1,0]))
         
         #%%++++++++++++++++++++++++++++++++++++++++++
         #add controller using loads
         nGround             = self.mbs.AddNode(NodePointGround(visualization=VNodePointGround(show=False)))
         mCGround            = self.mbs.AddMarker(MarkerNodeCoordinate(nodeNumber=nGround, coordinate=0))
         mRotNode0           = self.mbs.AddMarker(MarkerNodeCoordinate(nodeNumber=n0, coordinate=2))
         Controller          = Control(self.Kp,self.Ki,self.Kd)


         #load user function, second joint:
         def UFtorque0(mbs, t, load):
            phi0 = mbs.GetMarkerOutput(mRotNode0, exu.OutputVariableType.Coordinates)[0]
            omega = 0.5
            phiDesired = pi*0.5*sin(omega*pi*t)
            dt = 1e-3
            P = Controller.P(phiDesired, phi0, dt)
            phiDesired_t = pi*0.5*omega*cos(omega*pi*t)
            phi0_t = mbs.GetMarkerOutput(mRotNode0, exu.OutputVariableType.Coordinates_t)[0]

            D = Controller.D(phiDesired_t, phi0_t, dt)
            torque = P + D
            
            #torque = Control.PID(phiDesired, phi0, dt)
            
            return torque 

         lTorque0 = self.mbs.AddLoad(LoadCoordinate(markerNumber = mRotNode0,
                                    load = 1, loadUserFunction = UFtorque0))

         # add time integration scheme:
         self.mbs.Assemble()
         
         #add dependencies of load user functions on nodal coordinates:
         ltgN0 = self.mbs.systemData.GetNodeLTGODE2(n0)
         #both loads depend on both bodies; this information needs to be added in order that
         #  the solver knows dependencies when computing Jacobians
         self.mbs.systemData.AddODE2LoadDependencies(lTorque0, list(ltgN0))
         
         self.simulationSettings = exu.SimulationSettings()   
         self.simulationSettings.timeIntegration.numberOfSteps  = self.GetNSimulationSteps()
         self.simulationSettings.timeIntegration.endTime        = self.endTime

         self.simulationSettings.timeIntegration.newton.numericalDifferentiation.doSystemWideDifferentiation = True
         self.simulationSettings.timeIntegration.computeLoadsJacobian = 2
         # simulationSettings.timeIntegration.newton.useModifiedNewton = True

         self.simulationSettings.solutionSettings.writeSolutionToFile = True
         self.simulationSettings.solutionSettings.solutionWritePeriod = 0.01

         self.simulationSettings.timeIntegration.verboseMode = 1
         self.simulationSettings.displayStatistics = True
         self.simulationSettings.displayComputationTime = True

         self.mbs.SolveDynamic(simulationSettings=self.simulationSettings)          
         self.mbs.SolutionViewer()
         
  def DoublePendulum(self, theta1,theta2, dtheta1,dtheta2):

       self.theta1         = theta1       
       self.theta2         = theta2       
       self.dtheta1        = dtheta1
       self.dtheta2        = dtheta2
       
       self.dictSensors    = {}
       
       oGround             = self.mbs.AddObject(ObjectGround())
       mGround             = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, localPosition=[-0.5, 0, 0]))
       
       # create nodes:
       n0                  = self.mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[0,0,0], initialVelocities=[0,0,0]))
       n1                  = self.mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[1,0,0], initialVelocities=[0,0,0]))

       # create bodies:
       b0                  = self.mbs.AddObject(RigidBody2D(physicsMass=100, physicsInertia=1,nodeNumber=n0,
                                                            visualization=VRigidBody2D(graphicsData=[graphics.Lines([[-0.5,0,0],[0.5,0,0]])])))

       b1                  = self.mbs.AddObject(RigidBody2D(physicsMass=1, physicsInertia=1,nodeNumber=n1,
                                                       visualization=VRigidBody2D(graphicsData=[graphics.Lines([[-0.5,0,0],[0.5,0,0]])])))
       # add markers and loads:
       m00                 = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b0, localPosition=[-0.5, 0., 0.]))
       m01                 = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b0, localPosition=[ 0.5, 0., 0.]))
       m10                 = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b1, localPosition=[-0.5, 0., 0.]))

       # add joints:
       jg0                 = self.mbs.AddObject(RevoluteJoint2D(markerNumbers=[mGround,m00]))
       j01                 = self.mbs.AddObject(RevoluteJoint2D(markerNumbers=[m01,m10]))
       # add loads:
       mLoad0              = self.mbs.AddMarker(MarkerNodePosition(nodeNumber=n0))
       mLoad1              = self.mbs.AddMarker(MarkerNodePosition(nodeNumber=n1))
       self.mbs.AddLoad(Force(markerNumber=mLoad0, loadVector=[0,-9.81*1,0]))
       self.mbs.AddLoad(Force(markerNumber=mLoad1, loadVector=[0,-9.81*1,0]))
       
       #%%++++++++++++++++++++++++++++++++++++++++++
       #add controller using loads
       nGround             = self.mbs.AddNode(NodePointGround(visualization=VNodePointGround(show=False)))
       mCGround            = self.mbs.AddMarker(MarkerNodeCoordinate(nodeNumber=nGround, coordinate=0))
       mRotNode0           = self.mbs.AddMarker(MarkerNodeCoordinate(nodeNumber=n0, coordinate=2))
       mRotNode1           = self.mbs.AddMarker(MarkerNodeCoordinate(nodeNumber=n1, coordinate=2))
       
       Controller          = Control(self.Kp,self.Ki,self.Kd)
       
       #load user function, first joint:
       def UFtorque1(mbs, t, load):
           omega = 0.5
           phi0 = mbs.GetMarkerOutput(mRotNode0, exu.OutputVariableType.Coordinates)[0]
           deltaPhi = mbs.GetMarkerOutput(mRotNode1, exu.OutputVariableType.Coordinates)[0] - phi0

           phi0_t = mbs.GetMarkerOutput(mRotNode0, exu.OutputVariableType.Coordinates_t)[0]
           deltaPhi_t = mbs.GetMarkerOutput(mRotNode1, exu.OutputVariableType.Coordinates_t)[0] - phi0_t

           phiDesired = -pi*0.5*sin(omega*pi*t)
           phiDesired_t = -pi*0.5*omega*cos(omega*pi*t)
           
           dt = 1e-3
           P = Controller.P(phiDesired, deltaPhi, dt)
           D = Controller.D(phiDesired_t, deltaPhi_t, dt)
           torque = P + D
           return torque

       #load user function, second joint:
       def UFtorque0(mbs, t, load):
          phi0 = mbs.GetMarkerOutput(mRotNode0, exu.OutputVariableType.Coordinates)[0]
          omega = 0.5
          phiDesired = pi*0.5*sin(omega*pi*t)
          dt = 1e-3
          P = Controller.P(phiDesired, phi0, dt)
          phiDesired_t = pi*0.5*omega*cos(omega*pi*t)
          phi0_t = mbs.GetMarkerOutput(mRotNode0, exu.OutputVariableType.Coordinates_t)[0]

          D = Controller.D(phiDesired_t, phi0_t, dt)
          torque = P + D
          
          #torque = Control.PID(phiDesired, phi0, dt)
          
          return torque - UFtorque1(mbs,t,0)

       lTorque0 = self.mbs.AddLoad(LoadCoordinate(markerNumber = mRotNode0,
                                  load = 1, loadUserFunction = UFtorque0))
       
       lTorque1 = self.mbs.AddLoad(LoadCoordinate(markerNumber = mRotNode1,
                            load = 1, loadUserFunction = UFtorque1))

       # add time integration scheme:
       self.mbs.Assemble()
       
       #add dependencies of load user functions on nodal coordinates:
       ltgN0 = self.mbs.systemData.GetNodeLTGODE2(n0)
       ltgN1 = self.mbs.systemData.GetNodeLTGODE2(n1)
       #both loads depend on both bodies; this information needs to be added in order that
       #  the solver knows dependencies when computing Jacobians
       self.mbs.systemData.AddODE2LoadDependencies(lTorque0, list(ltgN0))
       self.mbs.systemData.AddODE2LoadDependencies(lTorque1, list(ltgN0)+list(ltgN1))
       
       self.simulationSettings = exu.SimulationSettings()   
       self.simulationSettings.timeIntegration.numberOfSteps  = self.GetNSimulationSteps()
       self.simulationSettings.timeIntegration.endTime        = self.endTime

       self.simulationSettings.timeIntegration.newton.numericalDifferentiation.doSystemWideDifferentiation = True
       self.simulationSettings.timeIntegration.computeLoadsJacobian = 2
       # simulationSettings.timeIntegration.newton.useModifiedNewton = True

       self.simulationSettings.solutionSettings.writeSolutionToFile = True
       self.simulationSettings.solutionSettings.solutionWritePeriod = 0.01

       self.simulationSettings.timeIntegration.verboseMode = 1
       self.simulationSettings.displayStatistics = True
       self.simulationSettings.displayComputationTime = True

       self.mbs.SolveDynamic(simulationSettings=self.simulationSettings)          
       self.mbs.SolutionViewer()

  def LiftBoom(self, theta, dtheta):

        self.theta1         = theta
        self.dtheta1        = dtheta
        self.dictSensors    = {}
        
        oGround            = self.mbs.AddObject(ObjectGround())
        
        #0.17, 0.386113249-11.17e-3, 00
        
        #markerGround       = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, localPosition=[-1.229248, -0.055596, 0]))
        Marker1            = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, localPosition=[-1.029248, 0.04, 0]))
        Marker4            = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=oGround, localPosition=[-1.289248,1.060504, 0])) 
       
        #Rigid pillar        
        #pMid2           = np.array([1.229248, 0.055596, 0])  
        iCube2          = RigidBodyInertia(mass=143.66, com=np.array([0, 0, 0]) ,
                                        inertiaTensor=np.array([[1.295612, 2.776103,  -0.000003],
                                                                [ 2.776103,  110.443667, 0],
                                                                [ -0.000003,              0  ,  110.452812]]),
                                        inertiaTensorAtCOM=True)
        
        graphicsBody2   = GraphicsDataFromSTLfile(fileName2, color4blue,
                                        verbose=False, invertNormals=True,
                                        invertTriangles=True)
        graphicsBody2   = AddEdgesAndSmoothenNormals(graphicsBody2, edgeAngle=0.25*pi,
                                            addEdges=True, smoothNormals=True)
        graphicsCOM2    = GraphicsDataBasis(origin=iCube2.com, length=2*0.263342)
        [n2, b2]        = AddRigidBody(mainSys=self.mbs,
                            inertia=iCube2,  # includes COM
                            nodeType=exu.NodeType.RotationEulerParameters,
                            position=np.array([-0.09, 1.4261, 0]),  # pMid2
                            rotationMatrix=RotationMatrixZ(mt.radians(self.theta1)),
                            gravity=[0,-9.806,0],
                            graphicsDataList=[graphicsCOM2, graphicsBody2])
        
        Marker7         = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[-1.229248, -0.055596, 0]))                         
        Marker8         = self.mbs.AddMarker(MarkerBodyRigid(bodyNumber=b2, localPosition=[ -0.92675, 
                                                                                             -0.158096, 0]))  

        self.mbs.AddObject(GenericJoint(markerNumbers=[Marker4, Marker7],constrainedAxes=[1,1,1,1,1,0],
                                  visualization=VObjectJointGeneric(axesRadius=0.18*0.263342,axesLength=0.18)))           
        
        Controller      = Control(self.Kp,self.Ki,self.Kd)


        def ControlTorque(mbs, t, itemIndex, theta, dtheta, k, d, offset):
                omega = 0.005
                phiDesired = sin(omega*t*pi)
                dt = 1e-3
                P = Controller.P(phiDesired, theta, dt)
                phiDesired_t = pi*omega*cos(omega*t)

                D = Controller.D(phiDesired_t, dtheta, dt)
                Torque = P + D
                return Torque
            
        # def ControlForce(mbs, t, itemIndex, u, v, k, d, f0):
        #         omega = 0.5
        #         phiDesired = pi*0.5*sin(omega*pi*t)
        #         dt = 1e-3
        #         P = Controller.P(phiDesired, u, dt)
        #         phiDesired_t = pi*0.5*omega*cos(omega*pi*t)

        #         D = Controller.D(phiDesired_t, v, dt)
        #         Force = P + D
        
        #         #torque = Control.PID(phiDesired, phi0, dt)
        
        #         return Force 
     
        # Force         = self.mbs.CreateSpringDamper(bodyOrNodeList=[b2, oGround],stiffness=0,damping=0,force=0,
        #                                      springForceUserFunction=ControlForce)
               
        Torque         = self.mbs.AddObject(TorsionalSpringDamper(markerNumbers=[Marker1, Marker8],
                                                                 stiffness=200,damping=0,offset=0,
                                                                 springTorqueUserFunction=ControlTorque))

        # add time integration scheme:
        self.mbs.Assemble()
        
        #add dependencies of load user functions on nodal coordinates:
        #ltgN0 = self.mbs.systemData.GetNodeLTGODE2(n0)
        #both loads depend on both bodies; this information needs to be added in order that
        #  the solver knows dependencies when computing Jacobians
        #self.mbs.systemData.AddODE2LoadDependencies(lTorque0, list(ltgN0))
        
        self.simulationSettings = exu.SimulationSettings()   
        self.simulationSettings.timeIntegration.numberOfSteps  = self.GetNSimulationSteps()
        self.simulationSettings.timeIntegration.endTime        = self.endTime

        self.simulationSettings.timeIntegration.newton.numericalDifferentiation.doSystemWideDifferentiation = True
        #self.simulationSettings.timeIntegration.computeLoadsJacobian = 2
        # simulationSettings.timeIntegration.newton.useModifiedNewton = True
        self.simulationSettings.timeIntegration.newton.useModifiedNewton = True
        self.simulationSettings.timeIntegration.generalizedAlpha.spectralRadius  = 0.7
        self.simulationSettings.timeIntegration.relativeTolerance                = 1e-6
        self.simulationSettings.staticSolver.newton.relativeTolerance            = 1e-11 

        self.simulationSettings.solutionSettings.writeSolutionToFile = True
        self.simulationSettings.solutionSettings.solutionWritePeriod = 0.01

        self.simulationSettings.timeIntegration.verboseMode = 1
        self.simulationSettings.displayStatistics = True
        self.simulationSettings.displayComputationTime = True

        solverType=exu.DynamicSolverType.TrapezoidalIndex2 
        self.mbs.SolveDynamic(solverType=solverType, simulationSettings=self.simulationSettings)
        
        self.mbs.SolutionViewer()           
        
    #%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    # def Plotting(self, data): 
    #     Time = self.timeVecOut
    #     end = Time.shape[0]
       
    #     if self.Flexible:
    #         Angle1                   = np.loadtxt(f"solution/Simulation_Flexible_f_modes_{self.nModes}_angle1.txt", delimiter=',')
    #         Angle2                   = np.loadtxt(f"solution/Simulation_Flexible_f_modes_{self.nModes}_angle2.txt", delimiter=',')
    #         angularVelocity1         = np.loadtxt(f"solution/Simulation_Flexible_f_modes_{self.nModes}_angularVelocity1.txt", delimiter=',')
    #         angularVelocity2         = np.loadtxt(f"solution/Simulation_Flexible_f_modes_{self.nModes}_angularVelocity2.txt", delimiter=',')
   
    #         sDistance1               = np.loadtxt('solution/sDistance1.txt', delimiter=',')
    #         sDistance2               = np.loadtxt('solution/sDistance2.txt', delimiter=',')
    #         #sVelocity1               = np.loadtxt('solution/sVelocity1.txt', delimiter=',')
    #         #sVelocity2               = np.loadtxt('solution/sVelocity2.txt', delimiter=',')    
    #         sPressuresL             = np.loadtxt('solution/sPressuresL.txt', delimiter=',')
    #         sPressuresT             = np.loadtxt('solution/sPressuresT.txt', delimiter=',')
        
    #     else:
    #         Angle1                   = np.loadtxt('solution/Simulation_Rigid_angle1.txt', delimiter=',')
    #         Angle2                   = np.loadtxt('solution/Simulation_Rigid_angle2.txt', delimiter=',')
    #         angularVelocity1         = np.loadtxt('solution/Simulation_Rigid_angularVelocity1.txt', delimiter=',')
    #         angularVelocity2         = np.loadtxt('solution/Simulation_Rigid_angularVelocity2.txt', delimiter=',')
            
    #         sDistance1               = np.loadtxt('solution/sDistance1.txt', delimiter=',')
    #         sDistance2               = np.loadtxt('solution/sDistance2.txt', delimiter=',')
    #         #sVelocity1               = np.loadtxt('solution/sVelocity1.txt', delimiter=',')
    #         #sVelocity2               = np.loadtxt('solution/sVelocity2.txt', delimiter=',')    
    #         sPressuresL             = np.loadtxt('solution/sPressuresL.txt', delimiter=',')
    #         sPressuresT             = np.loadtxt('solution/sPressuresT.txt', delimiter=',')
        
        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,16], label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, np.rad2deg(Angle1[0:end, 3]), label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
        
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('Angle 1, degree')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(0, 35)
    #     plt.tight_layout()
    #     plt.savefig('solution/angle1.png')
    #     plt.show()


    #      # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,17], label='Experimental', marker='x', 
    #       linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, np.rad2deg(Angle2[0:end, 3])+7.96, label='Simulation', linestyle='--', marker='x', 
    #       linewidth=1, markersize=2, color='black')
         
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('Angle 2, degree')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(-60, 0)
    #     plt.tight_layout()
    #     plt.savefig('solution/angle2.png')
    #     plt.show()
        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,12], label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, np.rad2deg(angularVelocity1[0:end, 3]), label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('Angular velocity 1, deg/s')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(-10, 10)
    #     plt.tight_layout()
    #     plt.savefig('solution/angularVelocity1.png')
    #     plt.show()

        

        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,15], label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, np.rad2deg(angularVelocity2[0:end, 3])-np.rad2deg(angularVelocity1[0:end, 3]), label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('Angular velocity 2, deg/s')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(-25, 25)
    #     plt.tight_layout()
    #     plt.savefig('solution/angularVelocity2.png')
    #     plt.show()


    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,2]+1000*L_Cyl1, label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, (sDistance1[0:end, 1])*1000, label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('Actuator 1, mm')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(900, 1200)
    #     plt.tight_layout()
    #     plt.savefig('solution/sDistance1.png')
    #     plt.show()
        
        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,6]+1000*L_Cyl2+30, label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, (sDistance2[0:end, 1])*1000, label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('Actuator 2, mm')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(1250, 1600)
    #     plt.tight_layout()
    #     plt.savefig('solution/sDistance2.png')
    #     plt.show()


    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,3]*1e5, label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, sPressuresL[0:end,1], label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('p1,pa')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(100e5, 200e5)
    #     plt.tight_layout()
    #     plt.savefig('solution/p1.png')
    #     plt.show()
        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,4]*1e5, label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, sPressuresL[0:end,2], label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('p2, Pa')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(0e5, 100e5)
    #     plt.tight_layout()
    #     plt.savefig('solution/p2.png')
    #     plt.show()
        
        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,7]*1e5, label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, sPressuresT[0:end,1], label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('p3, Pa')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(0, 75e5)
    #     plt.tight_layout()
    #     plt.savefig('solution/p3.png')
    #     plt.show()
        
    #     # Lift actuator
    #     plt.figure(figsize=(10, 5))
    #     plt.plot(Time, ExpData[0:end,8]*1e5, label='Experimental', marker='x', 
    #      linewidth=1, markersize=2, color='green')
    #     plt.plot(Time, sPressuresT[0:end,2], label='Simulation', linestyle='--', marker='x', 
    #      linewidth=1, markersize=2, color='black')
    #     plt.xlabel('Time, s')  # Adjust label as appropriate
    #     plt.ylabel('p4, Pa')  # Adjust label as appropriate
    #     plt.legend(loc='upper left')  # Specify legend location
    #     plt.grid(True)  # Add grid
    #     # Set axis limits
    #     plt.xlim(0, self.endTime)
    #     plt.ylim(50e5, 150e5)
    #     plt.tight_layout()
    #     plt.savefig('solution/p4.png')
    #     plt.show()

        
    #     # data_dict = data[1]
        
    #     # for key, value in data_dict.items(): 
    #     #         title = f' {key}'  # Customize your title
    #     #         plotData = np.column_stack([Time, value])
    #     #         self.mbs.PlotSensor(plotData, title=title, newFigure = True, closeAll = False)

        
    #     return    
#%%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if __name__ == '__main__': #include this to enable parallel processing

    model = Models(nStepsTotal=250, endTime=1, 
                         verboseMode=1)
    
    inputData = model.CreateInputVector(0)
    output = model.ComputeModel(inputData, verboseMode=True, 
                                solutionViewer=False)
