#Casey Copeland
import os
import math
import tkinter as tk

class Vector4(object):
    def __init__(self,x,y,z,w=0.0):
        self.x=x
        self.y=y
        self.z=z
        self.w=w
    def __str__ (self):
        return "Vector4("+str(self.x)+","+str(self.y)+","+str(self.z)+","+str(self.w)+")"
    def __add__ (self, other):
        return Vector4(self.x+other.x, self.y+other.y, self.z+other.z, self.w+other.w)
    def __sub__ (self,other):
        return self+(-other)
    def __neg__ (self):
        return Vector4(-self.x, -self.y, -self.z)
    def __mul__ (self,other):
        hamilton = Vector4(self.w*other.x+self.x*other.w+self.y*other.z-self.z*other.y,
                self.w*other.y-self.x*other.z+self.y*other.w+self.z*other.x,
                self.w*other.z+self.x*other.y-self.y*other.x+self.z*other.w,
                self.w*other.w-self.x*other.x-self.y*other.y-self.z*other.z)
        return hamilton
    def distance (self):
        return math.pow(self.x*self.x+self.y*self.y+self.z*self.z+self.w*self.w,0.5)
    def normal (self):
        norm = self.distance()
        if norm == 0:
            return Vector4(0.0,0.0,0.0)
        unit = Vector4(self.x/norm,self.y/norm,self.z/norm,self.w/norm)
        #Set small values to zero so that the magnitude of the unit vector is closer to 1
        if unit.x*unit.x+unit.y*unit.y+unit.w*unit.w== 1.0:
            unit.z = 0.0
        if unit.x*unit.x+unit.z*unit.z+unit.w*unit.w == 1.0:
            unit.y = 0.0
        if unit.y*unit.y+unit.z*unit.z+unit.w*unit.w == 1.0:
            unit.x = 0.0
        if unit.x*unit.x+unit.y*unit.y+unit.z*unit.z == 1.0:
            unit.w = 0.0
        return unit
    def normalize (self):
        norm = self.normal()
        distance = self.distance()
        self.x = norm.x*distance
        self.y = norm.y*distance
        self.z = norm.z*distance
    def rotate (self,axis,theta):
        magnitude = self.distance()
        axis=axis.normal()
        rad = math.radians(theta)
        sineVal = math.sin(rad/2.0)
        quaternion = Vector4(axis.x*sineVal,axis.y*sineVal,
                axis.z*sineVal, math.cos(rad/2.0))
        conjugate = Vector4(-axis.x*sineVal,-axis.y*sineVal,
                -axis.z*sineVal, math.cos(rad/2.0))
        vec = (quaternion*self*conjugate).normal()
        vec.x*=magnitude
        vec.y*=magnitude
        vec.z*=magnitude
        return vec

def openMenuEntry (prevFrame,entryName,inputs,outputs,calcFunc):
    prevFrame.destroy()
    frame = tk.Tk()
    frame.title(entryName)

    rowCount = 0
    controls = []
    for inp in inputs:
        var = tk.StringVar()
        tk.Label(frame,text=inp).grid(row=rowCount,column=0)
        tk.Entry(textvariable=var).grid(row=rowCount,column=1)
        rowCount+=1
        controls.append(var)
    for output in outputs:
        var = tk.StringVar()
        tk.Label(frame,text=output).grid(row=rowCount,column=0)
        tk.Label(frame,textvariable=var).grid(row=rowCount,column=1)
        rowCount+=1
        controls.append(var)

    button = tk.Button(frame,text="Calculate",
            command= lambda : calcFunc(*controls))
    button.grid(row=rowCount,column=1)
    rowCount+=1
    button = tk.Button(frame,text="Return", command= lambda : returnMain(frame))
    button.grid(row=rowCount,column=1)
    frame.mainloop()

def createMenuEntry (frame,entryName,inputs,outputs,calcFunc):
    button = tk.Button(frame,text=entryName,
            command=lambda : openMenuEntry(frame,entryName,inputs,outputs,calcFunc))
    button.pack()

def returnMain (prevFrame):
    prevFrame.destroy()
    guiMain()

#treats circle as laying on the xy-plane
#Ex = (r*lambda)/(4*pi*epsilon)*integral(
#      (dtheta*(x2-rcos(theta)))/((x2-rcos(theta))^2+(y2-rsin(theta))^2+z^2)^(3/2) )
#Ey = (r*lambda)/(4*pi*epsilon)*integral(
#      (dtheta*(y2-rsin(theta)))/((x2-rcos(theta))^2+(y2-rsin(theta))^2+z^2)^(3/2) )
#Ez = (z*r*lambda)/(4*pi*epsilon)*integral(
#      dtheta/((x2-rcos(theta))^2+(y2-rsin(theta))^2+z^2)^(3/2) )
def circleIntegral (q,r, x, y, z, steps):
    k = 8.987551787*10**9
    density = q/(2.0*math.pi*r)
    constant = k*density*r
    field = Vector4(0.0,0.0,0.0)
    dt = (2.0*math.pi)/steps
    #for if dt can't split the length up into equal parts due to precision
    dtRemainder = 1.0
    theta = 0.0
    #uses midpoint summation for estimation
    while theta < 2.0*math.pi:
        if theta+dt > 2.0*math.pi:
            #divide by dt to account for outside multiplication by dt
            dtRemainder = (2.0*math.pi-theta)/dt
        theta+=dt/2.0
        denom = (math.pow((math.pow(x-r*math.cos(theta),2)+math.pow(y-r*math.sin(theta),2)+math.pow(z,2)),3.0/2.0))
        field.x+=dtRemainder*(x-r*math.cos(theta))/denom
        field.y+=dtRemainder*(y-r*math.sin(theta))/denom
        field.z+=dtRemainder/denom
        theta+=dt/2.0
    field.normalize()
    field.x*=constant*dt
    field.y*=constant*dt
    field.z*=z*constant*dt
    field.normalize()
    return field

#Disc laying on XY plane
#Uses double summation, so total number of parts to be summed split up
#between the two summations. This way both double integrals and single
#integrals take similar number of computations for same number of steps
def discIntegral (q, r, x, y, z,steps):
    k = 8.987551787*10**9
    density = q/(math.pi*r*r)
    constant = k*density
    field = Vector4(0.0,0.0,0.0)
    singleIntegrations = math.sqrt(steps)
    #Tries to perserve symmetry in calculations by minimizing dtRemainder
    dt = (2.0*math.pi)/(float(int(singleIntegrations)))
    dr = r/singleIntegrations
    #remainder for if delta can't split up into equal pieces
    dtRemainder = 1.0
    theta=0.0
    while theta < 2.0*math.pi:
        if theta+dt > 2.0*math.pi:
        #divide by dt to account for outside multiplication by dt
            dtRemainder = (2.0*math.pi-theta)/dt
        theta+=dt/2.0
        drRemainder = 1.0
        currentR = 0.0
        lineField = Vector4(0.0,0.0,0.0)
        #uses Milne's rule
        while currentR < r:
            if currentR+dr > r:
                #divide by dr to account for outside multiplication by dr
                drRemainder = (r-currentR)/dr
            currentR+=dr/4.0
            denom = (math.pow((math.pow(x-currentR*math.cos(theta),2)+math.pow(y-currentR*math.sin(theta),2)+math.pow(z,2)),3.0/2.0))
            lineField.x+=2*drRemainder*currentR*(x-currentR*math.cos(theta))/denom
            lineField.y+=2*drRemainder*currentR*(y-currentR*math.sin(theta))/denom
            lineField.z+=2*drRemainder*currentR/denom
            lineField.normalize()
            currentR+=dr/4.0
            denom = (math.pow((math.pow(x-currentR*math.cos(theta),2)+math.pow(y-currentR*math.sin(theta),2)+math.pow(z,2)),3.0/2.0))
            lineField.x-=drRemainder*currentR*(x-currentR*math.cos(theta))/denom
            lineField.y-=drRemainder*currentR*(y-currentR*math.sin(theta))/denom
            lineField.z-=drRemainder*currentR/denom
            lineField.normalize()
            currentR+=dr/4.0
            denom = (math.pow((math.pow(x-currentR*math.cos(theta),2)+math.pow(y-currentR*math.sin(theta),2)+math.pow(z,2)),3.0/2.0))
            lineField.x+=2*drRemainder*currentR*(x-currentR*math.cos(theta))/denom
            lineField.y+=2*drRemainder*currentR*(y-currentR*math.sin(theta))/denom
            lineField.z+=2*drRemainder*currentR/denom
            lineField.normalize()
            currentR+=dr/4.0
        lineField.normalize()
        lineField.x*=dtRemainder
        lineField.y*=dtRemainder
        lineField.z*=dtRemainder
        lineField.normalize()
        field+=lineField
        field.normalize()
        theta+=dt/2.0
    field.normalize()
    field.x*=constant*(dr/3.0)*dt
    field.y*=constant*(dr/3.0)*dt
    field.z*=z*constant*(dr/3.0)*dt
    field.normalize()
    return field

def guiCircleCalculate (q,r,x,y,z,steps,magnitude,fieldX,fieldY,fieldZ):
    efield = circleIntegral(float(q.get()),float(r.get()),
            float(x.get()),float(y.get()),float(z.get()),float(steps.get()))
    magnitude.set(str(efield.distance()))
    fieldX.set(str(efield.x))
    fieldY.set(str(efield.y))
    fieldZ.set(str(efield.z))

def guiDiscCalculate (q,r,x,y,z,steps,magnitude,fieldX,fieldY,fieldZ):
    efield = discIntegral(float(q.get()),float(r.get()),
            float(x.get()),float(y.get()),float(z.get()),float(steps.get()))
    magnitude.set(str(efield.distance()))
    fieldX.set(str(efield.x))
    fieldY.set(str(efield.y))
    fieldZ.set(str(efield.z))

def guiMain():
    menuCount=0
    frame = tk.Tk()
    frame.title("E field")
    createMenuEntry(frame,"Ring Charge",["Charge","Radius","Point x","Point y","Point z","Steps"],
            ["Magnitude","X Component","Y Component","Z Component"],guiCircleCalculate)
    createMenuEntry(frame,"Disc Charge",["Charge","Radius","Point x","Point y","Point z","Steps"],
            ["Magnitude","X Component","Y Component","Z Component"],guiDiscCalculate)
    frame.mainloop()

if __name__ == "__main__":
    guiMain()
