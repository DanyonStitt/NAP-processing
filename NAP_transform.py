# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 15:57:44 2022

@author: dst74

Transformation for the output files of the 
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sp
from Low_pass import butter_lowpass_filter
import pandas as pd

def transform(a, year):
    # Finding the expected impact velocity
    # Condition for heights of 7.5cm
    if "7.5" in a:
        height = 7.5
    # Condition for heights of 22.5cm
    elif "22.5" in a:
        height = 22.5
    else:
        b = a.find("cm")
        height = a[b-3:b-1]
        
    # height = float(height)/100
    # print(height)
    # expected_impact_velocity = np.sqrt(2*9.81*height) # Kinematic expected impact velocity
    
    ##########################################################################
    # Start the NAP transform
    file = np.loadtxt(a, delimiter=',')
    
    # Transforming voltages to accelerations (m/s^2)
    offset = file[:,:1000].mean(axis=0)
    # print(offset)
    # file[0:200,1:] = offset[1:]
    offset[0] = 0 # removing the offset for the time column
    lin_acc = file - offset

    # fig1, ax = plt.subplots(4)
    # ax[0].plot(lin_acc[:,0], lin_acc[:,4:7])
    # ax[0].set_title("COM")

    # ax[1].plot(lin_acc[:,0], lin_acc[:,1:4])
    # ax[1].set_title("x offset")

    # ax[2].plot(lin_acc[:,0], lin_acc[:,7:10])
    # ax[2].set_title("y offset")

    # ax[3].plot(lin_acc[:,0], lin_acc[:,10:])
    # ax[3].set_title("z offset")
    # plt.suptitle(a)

    # x-axis transform
    for i in range(1, 11, 3):
        lin_acc[:,i] = 9.81*lin_acc[:,i]/(3.3418*6.525e-3)
        
    # y-axis transform
    for i in range(2,12,3):
        lin_acc[:,i] = 9.81*lin_acc[:,i]/(3.3418*6.575e-3)

    # z-axis transform
    for i in range(3,13,3):
        lin_acc[:,i] = 9.81*lin_acc[:,i]/(3.3418*6.55e-3)
    
    ##########################################################################
    # Finding peak com acceleration and location of peak
    if year >= 2022:
        peak_acc = np.absolute(lin_acc[:,1:4]).max()
        pla_time = np.where(np.absolute(lin_acc[:,1:4]) == peak_acc)[0][0]
    else:
        peak_acc = np.absolute(lin_acc[:,4:7]).max()
        pla_time = np.where(np.absolute(lin_acc[:,4:7]) == peak_acc)[0][0]
    
    #  SECTION 1: FIXING THE PRE IMPACT DATA
    # Integrating every channel to find velocity and displacement
    lin_vel = (1/20000)*sp.cumtrapz(lin_acc, axis=0)
    lin_vel[:,0] = lin_acc[0:-1,0]
            
    # finding the dominant axis and it's maximum values
    peak_vel = np.absolute(lin_vel)[pla_time-2000:pla_time,:].max(0)
    dominant_vel = peak_vel[1:].max(0)
    dominant_axis = np.where(np.absolute(lin_vel) == dominant_vel)
    
    # Defining the x and y points for curve fitting
    time_1 = int(dominant_axis[0]) - 1000 # Start time
    impact_time = int(dominant_axis[0]) # End time
    axis = int(dominant_axis[1]) # Dominant axis
    l1 = [lin_vel[time_1,0], lin_vel[impact_time,0]] # x points
    l2 = [lin_vel[time_1, axis], lin_vel[impact_time, axis]] # y points
    
    # Fitting the 1st order polynomial to the data
    lin_fit_1 = np.polynomial.polynomial.polyfit(l1, l2, 1)
    x = np.linspace(l1[1]-0.5, l1[1], int((l1[1] - (l1[1]-0.5))*20000))
    y = lin_fit_1[1]*x + lin_fit_1[0]
    
    # Finding the zero crossing for the fitted curve
    if lin_vel[impact_time,axis] < 0:
        start = np.where(y < 0)[0][0]
    else:
        start = np.where(y > 0)[0][0]
    
    # Finding the index of the zero crossing in the accelerometer data set
    release_time = np.where(lin_acc[:,0] > x[start])[0][0]
    # Defing the end of the impact as 1 second after the start
    end_time = impact_time + 10000
    # Setting all post  impact acceleration values to zero
    lin_acc[end_time:-1,1:] = 0
    
    # Setting all acceleration values pre-release to zero
    lin_acc[0:release_time,1:] = 0
    lin_vel = (1/20000)*sp.cumtrapz(lin_acc, axis=0)
    lin_vel[:,0] = lin_acc[0:-1,0]

    
    # Checking for bugs
    fig1, ax = plt.subplots(4)
    ax[0].plot(lin_acc[0:-1,0], lin_vel[:,4:7])
    ax[0].set_title("COM")
    ax[0].set_xlabel("Time (s)")
    ax[0].set_ylabel("Velocity (m/s)")
    ax[0].minorticks_on()
    ax[0].grid(b=True, which='major')
    ax[0].grid(b=True, which="minor", linestyle='--')
    
    ax[1].plot(lin_acc[0:-1,0], lin_vel[:,1:4])
    ax[1].set_title("x offset")
    ax[1].set_xlabel("Time (s)")
    ax[1].set_ylabel("Velocity (m/s)")
    ax[1].minorticks_on()
    ax[1].grid(b=True, which='major')
    ax[1].grid(b=True, which="minor", linestyle='--')
    
    ax[2].plot(lin_acc[0:-1,0], lin_vel[:,7:10])
    ax[2].set_title("y offset")
    ax[2].set_xlabel("Time (s)")
    ax[2].set_ylabel("Velocity (m/s)")
    ax[2].minorticks_on()
    ax[2].grid(b=True, which='major')
    ax[2].grid(b=True, which="minor", linestyle='--')
    
    ax[3].plot(lin_acc[0:-1,0], lin_vel[:,10:])
    ax[3].set_title("z offset")
    ax[3].set_xlabel("Time (s)")
    ax[3].set_ylabel("Velocity (m/s)")
    ax[3].minorticks_on()
    ax[3].grid(b=True, which='major')
    ax[3].grid(b=True, which="minor", linestyle='--')
    plt.suptitle(a)
    
    
    ##########################################################################
    # Separating velocities into accelerameter arrays
    v_com = lin_vel[:,4:7]
    v_x = lin_vel[:,1:4]
    v_y = lin_vel[:,7:10]
    v_z = np.zeros(v_x.shape)
    v_z[:,0] = lin_vel[:,11] *-1
    v_z[:,1] = lin_vel[:,10]
    v_z[:,2] = lin_vel[:,12]
    
    # # Impact velocity check
    # impact_velocities = np.zeros(4)
    # impact_velocities[0] = (np.absolute(v_com).max()).max()
    # impact_velocities[1] = (np.absolute(v_x).max()).max()
    # impact_velocities[2] = (np.absolute(v_y).max()).max()
    # impact_velocities[3] = (np.absolute(v_z).max()).max()
    
    # impact_velocity_diff = np.absolute(impact_velocities.max() - impact_velocities.min())
    # impact_velocity_error = np.absolute(impact_velocities - expected_impact_velocity).max()
    
    # # If the velocities differ from each other too much
    # if impact_velocity_diff >= 0.2*impact_velocities.mean():
    #     # Find the xyz impact velocities for each accelerometer between t=0 and impact
    #     com_imp_vel = np.absolute(v_com)[0:impact_time,:].max(axis=0)
    #     x_imp_vel = np.absolute(v_x)[0:impact_time,:].max(axis=0)
    #     y_imp_vel = np.absolute(v_y)[0:impact_time,:].max(axis=0)
    #     z_imp_vel = np.absolute(v_z)[0:impact_time,:].max(axis=0)
        
    #     # Find the resultant impact velocities of each accelerometer
    #     res_imp_vel = np.zeros(4)
    #     res_imp_vel[0] = np.sqrt(com_imp_vel[0]**2 + com_imp_vel[1]**2 + com_imp_vel[2]**2)
    #     res_imp_vel[1] = np.sqrt(x_imp_vel[0]**2 + x_imp_vel[1]**2 + x_imp_vel[2]**2)
    #     res_imp_vel[2] = np.sqrt(y_imp_vel[0]**2 + y_imp_vel[1]**2 + y_imp_vel[2]**2)
    #     res_imp_vel[3] = np.sqrt(z_imp_vel[0]**2 + z_imp_vel[1]**2 + z_imp_vel[2]**2)
        
    #     # Find the difference between the mean and all accelerometers
    #     mean_res_imp_vel = res_imp_vel.mean()
    #     res_imp_vel_diff = np.absolute(res_imp_vel - mean_res_imp_vel)
        
    #     # If the resultant velocities still differ too much
    #     if res_imp_vel_diff.max() > 0.2*mean_res_imp_vel:
    #         print(a)
    #         print("Mean impact velocity: ", impact_velocities.mean())
    #         print("Measured velocities: ", impact_velocities)
    #         print("Mean resultant impact velocity: ", mean_res_imp_vel)
    #         print("Measured velocities: ", res_imp_vel)
    #         print(res_imp_vel_diff)
    #         # raise Exception("Bad impact. Impact velocity differs between accelerometers")
            
    # # If the velocities differ from the theoretical too much
    # elif impact_velocity_error > 0.2*expected_impact_velocity:
    #     # Find the resultant velocity from the three accelerometers
    #     com_imp_vel = np.absolute(v_com).max(axis=0)
    #     resultant_imp_vel = np.sqrt(com_imp_vel[0]**2 + com_imp_vel[1]**2 + com_imp_vel[2]**2)
    #     resultant_imp_vel_diff = np.absolute(resultant_imp_vel - expected_impact_velocity)
        
    #     if resultant_imp_vel_diff > 0.2*expected_impact_velocity:
    #             print(a)
    #             print("Expected velocity: ", expected_impact_velocity)
    #             print("Measured velocity: ", impact_velocities)
    #             raise Exception("Bad impact. Impact velocity differs from theoretical")
        
    ##########################################################################
    # Capturing the relevant part of the data
    lin_acc[end_time:-1,1:] = 0
    lin_acc[0:impact_time,1:] = 0
    lin_acc = lin_acc[impact_time-200:end_time,:]
    
    v_com = v_com[impact_time-200:end_time,:]
    
    # Separating each accelerometer into it's own array
    a_com = lin_acc[:,4:7]
    a_x = lin_acc[:,1:4]
    
    a_y = lin_acc[:,7:10]
    a_z = np.zeros(a_x.shape)
    a_z[:,0] = lin_acc[:,11] *-1
    a_z[:,1] = lin_acc[:,10]
    a_z[:,2] = lin_acc[:,12]
    
    # Accelerometer offset distances from the centre of mass accelerometer(m)
    rhoX = -50.8e-3
    rhoY = 48.3e-3
    rhoZ = 82.1e-3
    
    # Preallocating memory
    rot_acc = np.zeros(a_com.shape)
    
    # Acceleration calcs based on 'Measuring acceleration of a rigid body paper
    rot_acc[:,0] = ((a_y[:,2] - a_com[:,2])/(2*rhoY)) - ((a_z[:,1] - a_com[:,1])/(2*rhoZ))
    rot_acc[:,1] = ((a_z[:,0] - a_com[:,0])/(2*rhoZ)) - ((a_x[:,2] - a_com[:,2])/(2*rhoX))
    rot_acc[:,2] = ((a_x[:,1] - a_com[:,1])/(2*rhoX)) - ((a_y[:,0] - a_com[:,0])/(2*rhoY))
    
    time = lin_acc[:,0]
    
    # Integrating the rotational data
    rot_vel = (1/20000)*sp.cumtrapz(rot_acc, axis=0)
    
    # # Bug checking
    # fig2, ax4 = plt.subplots()
    # ax4.plot(lin_acc[0:-1,0], rot_vel, label="Adjusted")
    # ax4.set_title(a)
    # ax4.set_xlabel("Time (s)")
    # ax4.set_ylabel("Velocity (rad/s)")
    # ax4.minorticks_on()
    # ax4.grid(b=True, which='major')
    # ax4.grid(b=True, which="minor", linestyle='--')
    
    
    ##########################################################################
    # Filtering the data through a 8th order low pass filter
    fs = 20000
    cutoff = 300
    
    for i in range(0,3):
        rot_vel[:,i] = butter_lowpass_filter(rot_vel[:,i], cutoff, fs)
        rot_acc[:,i] = butter_lowpass_filter(rot_acc[:,i], cutoff, fs)
        a_com[:,i] = butter_lowpass_filter(a_com[:,i], cutoff, fs)
    
    # Converting linear acc to g's
    a_com = a_com/9.81
    
    # Fixing the lengths of the data sets
    a_com = a_com[0:-1,:]
    v_com = v_com[0:-1,:]
    rot_acc = rot_acc[0:-1,:]
    time = time[0:-1]
    time = np.linspace(-10, ((len(time)-200)/20), len(time))
    
    # Bug checking
    fig, ax = plt.subplots(3)
    ax[0].plot(time, a_com)
    ax[1].plot(time, rot_acc)
    ax[2].plot(time, rot_vel)
    ax[0].minorticks_on()
    ax[1].minorticks_on()
    ax[2].minorticks_on()
    ax[0].grid(b=True, which='major')
    ax[0].grid(b=True, which="minor", linestyle='--')
    ax[1].grid(b=True, which='major')
    ax[1].grid(b=True, which="minor", linestyle='--')
    ax[2].grid(b=True, which='major')
    ax[2].grid(b=True, which="minor", linestyle='--')
    plt.suptitle(a)
    plt.show()
    ##########################################################################
    # Exporting the data
    d = {"time" : time, "lin vel x" : v_com[:,0], "lin vel y" : v_com[:,1], "lin vel z" : v_com[:,2],
         "lin acc x" : a_com[:,0], "lin acc y" : a_com[:,1], "lin acc z" : a_com[:,2],
         "rot vel x" : rot_vel[:,0], "rot vel y" : rot_vel[:,1], "rot vel z" : rot_vel[:,2],
         "rot acc x" : rot_acc[:,0], "rot acc y" : rot_acc[:,1], "rot acc z" : rot_acc[:,2]}
    
    data_out = pd.DataFrame(data=d)
    
    return data_out

import os 
import glob

file_path = os.path.join("E:\\bug fix octobver 2\\", "*.txt")
filenames = glob.glob(file_path)

for file in filenames:
    print("\n", file, "\n")
    impact = transform(file, 2022)
    plt.show()