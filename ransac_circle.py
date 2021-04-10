camera_cx = 637.166
camera_cy = 369.465
camera_fx = 611.223;
camera_fy = 610.982;
import random
import numpy as np
def estimate_sphere(data):


    x1=data[0][0]
    y1=data[0][1]
    z1=data[0][2]
    x2 = data[1][ 0]
    y2 = data[1][ 1]
    z2 = data[1][ 2]
    x3 = data[2][ 0]
    y3 = data[2][ 1]
    z3 = data[2][ 2]
    x4 = data[3][ 0]
    y4 = data[3][1]
    z4 = data[3][ 2]
    a11 = 2 * (x2 - x1);
    a12 = 2 * (y2 - y1);
    a13 = 2 * (z2 - z1);
    a21 = 2 * (x3 - x2);
    a22 = 2 * (y3 - y2);
    a23 = 2 * (z3 - z2);
    a31 = 2 * (x4 - x3);
    a32 = 2 * (y4 - y3);
    a33 = 2 * (z4 - z3);
    b1 = x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + z2 * z2 - z1 * z1;
    b2 = x3 * x3 - x2 * x2 + y3 * y3 - y2 * y2 + z3 * z3 - z2 * z2;
    b3 = x4 * x4 - x3 * x3 + y4 * y4 - y3 * y3 + z4 * z4 - z3 * z3;
    d = a11 * a22 * a33 + a12 * a23 * a31 + a13 * a21 * a32 - a11 * a23 * a32 - a12 * a21 * a33 - a13 * a22 * a31;
    d1 = b1 * a22 * a33 + a12 * a23 * b3 + a13 * b2 * a32 - b1 * a23 * a32 - a12 * b2 * a33 - a13 * a22 * b3;
    d2 = a11 * b2 * a33 + b1 * a23 * a31 + a13 * a21 * b3 - a11 * a23 * b3 - b1 * a21 * a33 - a13 * b2 * a31;
    d3 = a11 * a22 * b3 + a12 * b2 * a31 + b1 * a21 * a32 - a11 * b2 * a32 - a12 * a21 * b3 - b1 * a22 * a31;
    x = d1 / d;
    y = d2 / d;
    z = d3 / d;

    radius=np.sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z))

    center=np.array([x,y,z])

    return center,radius

def verify_circle(data,center,radius,distance_threshold):
    x=data[:,0]
    y=data[:,1]
    z=data[:,2]
    dist=np.sqrt((center[0]-x)*(center[0]-x)+(center[1]-y)*(center[1]-y)+(center[2]-z)*(center[2]-z))
    error=np.abs(dist-radius)
    dist_right=np.sum(error<distance_threshold)

    cPerc=dist_right/len(data)
    return cPerc

def my_ransac(data,distance_threshold=10,P=0.95,
              sample_size=4,max_iterations=10000,radius_threshold=50):

    i=0
    K=100

    L_data=len(data)
    R_L=range(L_data)
    best_center = np.array([0, 0, 0])
    best_radius = 0
    if L_data<5:
        return best_center,best_radius

    best_perc=0
    # RANSAC-ING


    while i < K:

        s4=random.sample(R_L,sample_size)

        center,radius=estimate_sphere(data[s4])#求出球体参数，球心坐标和半径
        if radius>radius_threshold:
            continue
        cPerc=verify_circle(data,center,radius,distance_threshold)#内点占比
        if cPerc>best_perc :
            best_center=center
            best_radius=radius
            best_perc=cPerc
            wn=np.power(cPerc,4)
            K=(np.log(1-P)/np.log(1-wn))

        i+=1
        if i>max_iterations:
            break
    print("best_perc ,best_radius",best_perc,best_radius)
    return best_center,best_radius

