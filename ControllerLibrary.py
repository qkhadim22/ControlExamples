#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Details:  It includes PD controller on spring-damper, single pendulum, double 
            # pendulum and lift boom.
#
# Author:   Qasim Khadim
# Contact: qasim.khadim@outlook.com,qkhadim22 (Github)
# Date:     2025-02-02
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
class Control:
    
    def __init__(self, Kp, Ki, Kd):
        self.Kp = Kp  # Proportional gain
        self.Ki = Ki  # Integral gain
        self.Kd = Kd  # Derivative gain
        self.prev_error = 0  # Previous error value
        self.integral = 0  # Integral of error
        
#+++++++++++++++++++++P-Controller++++++++++++++++++++++++++++++++++++        
    def P(self, setpoint, current_value, dt):
        # Compute the error
        error = setpoint - current_value
        
        # Proportional term
        P = self.Kp * error
        
        # Update previous error
        self.prev_error = error
        
        return P
    
    #+++++++++++++++++++++PI-Controller++++++++++++++++++++++++++++++++++++        
    def D(self, setpoint, current_value, dt):
           # Compute the error
           error = setpoint - current_value
           
           
           # Derivative term
           D = self.Kd * (error - self.prev_error)
           
           # Update previous error
           self.prev_error = error
          
           return D
        

#+++++++++++++++++++++PI-Controller++++++++++++++++++++++++++++++++++++        
    def PI(self, setpoint, current_value, dt):
        # Compute the error
        error = setpoint - current_value
       
        # Proportional term
        P = self.Kp * error
        
        # Integral term
        self.integral += error * dt
        
        I = self.Ki * self.integral
                
        # Update previous error
        self.prev_error = error
        
        return  P + I

#+++++++++++++++++++++PID-Controller++++++++++++++++++++++++++++++++++++        
    def PID(self, setpoint, current_value, 
            dt):
        # Compute the error
        error = setpoint - current_value
        
        # Proportional term
        P = self.Kp * error
        
        # Integral term
        self.integral += error * dt
        
        I = self.Ki * self.integral
        
        # Derivative term
        D = self.Kd * (error - self.prev_error) / dt
        
        # Update previous error
        self.prev_error = error
        
        return P + I + D

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++