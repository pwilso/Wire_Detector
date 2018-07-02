# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 17:34:00 2017

@author: Paige
"""

from scipy import ndimage
import matplotlib.pyplot as plt
import numpy as np

### Functions #################################################################

def filter_pic(picture, sigma, threshold, scale, tilt):
    pic = ndimage.imread(picture, mode  = 'L', flatten = True)                  # import picture
    pic = ndimage.gaussian_filter(pic, sigma)                                   # apply Gaussian filter
                                 
    pic_grad = np.gradient(pic)                                                 # take gradient          
    Gx = pic_grad[0]
    Gy = pic_grad[1]
    mag = np.sqrt(Gx**2 + Gy**2)                                                # magnitude of gradient
    
    over = (mag > threshold).astype(int)
    thresh = mag*over 
    
    return thresh

def maximum(picture, row, column, L_or_R):                                      # find a local maximum to the pixel
    if L_or_R == "Right":
        while picture[row, column] <= picture[row, column + 1]:
            column += 1
    else:
        while picture[row, column] <= picture[row, column - 1]:
            column -= 1
    return column

def subpixel(picture, row, column, L_or_R, max_col):                            # adjust the maximum within the pixel
    g = []
    d = []
    while picture[row, column] != 0:
        g.append(picture[row, column])
        d.append(column - max_col)
        if L_or_R == "Right":
            column += 1
        else:
            column -= 1
    g = np.array(g)
    d = np.array(d)
    delta = sum(g*d)/sum(g)
    return delta 

def plotter(picture_ax, right_point, left_point, height_point):
    outline_pic = plt.figure(3)
    outline_pic_ax = outline_pic.add_subplot(111)
    outline_pic_ax.imshow(unfiltered)
    outline_pic_ax.axis("off")
    picture_ax.scatter(right_point, height_point, color = 'c', marker = '.')
    picture_ax.scatter(left_point, height_point, color = 'm', 
                           marker  = '.')
    plt.draw()       

def find_edge(m_value, n_value):
    m = m_value
    n = n_value
    value = pic_array[m, n]
    #print(value)
    
    while value == 0:                                                           # check right
        n += 1
        value = pic_array[m, n]
    n_nonzero = n
    n_max = maximum(pic_array, m, n, "Right")
    subpixel_n = subpixel(pic_array, m, n_nonzero, "Right", n_max) + n_max
    right = subpixel_n
    
    n = n_value
    value = pic_array[m, n]
    
    while value == 0:                                                           # check left
        n -= 1
        value = pic_array[m, n]
    n_nonzero = n
    n_max = maximum(pic_array, m, n, "Left")
    subpixel_n = subpixel(pic_array, m, n_nonzero, "Left", n_max) + n_max
    left = subpixel_n
    
    return (right, left, m)                                                     # returns the diameter and y coord in pixel values

def click(event):
    global diam, height, special, specialcount, diam_pix, height_pix
    n_orig = int(event.xdata)
    m_orig = int(event.ydata)
        
#    if len(diam_pix) < 2:
#        results = find_edge(m_orig, n_orig)
#        plotter(outline_pic_ax, results[0], results[1], results[2])
#        diam_pix.append(abs(results[0] - results[1]))
#        height_pix.append(results[2])
#    elif len(diam_pix) == 2:   
#        points = np.arange(height_pix[1], height_pix[0], -1*scale*step)
#        for p in points:
#            results = find_edge(p, n_orig)
#            plotter(outline_pic_ax, results[0], results[1], results[2])
#            diam_pix.insert(-1, abs(results[0] - results[1]))
#            height_pix.append(-1, p)
#    else:
#        if special == 0:
#            pass
#        else:
    if specialcount != 0:
        results = find_edge(m_orig, n_orig)
        plotter(outline_pic_ax, results[0], results[1], results[2])
        #plotter(thresh_pic_ax, results[0], results[1], results[2])
        diam_pix.append(abs(results[0] - results[1]))
        height_pix.append(results[2])
        specialcount -= 1
    else:
        results = find_edge(m_orig, n_orig)
        plotter(outline_pic_ax, results[0], results[1], results[2])
        diam_pix.append(abs(results[0] - results[1]))
        #print(diam_pix)
        height_pix.append(results[2])  
        #print(height_pix)
        diam = np.array(diam_pix)/scale
        height = np.array(height_pix)
        height = ((height - height[0])/scale)/np.sin(np.radians(tilt))
        #print("Height (nm)\tDiameter (nm)")
        for i in range(len(height)):
            print str(height[i]).ljust(15), '\t', diam[i]
        #print '\n'
        diam_pix = []
        height_pix = []
        diam = []
        height = []
        specialcount = special
    
### Parameters ################################################################

folder = "C:\Users\Paige\Desktop\Samples\Sample_1953\Oct30"
picture = "\\1953_1000D2_col2_x20_30deg.tif" 
img_name = folder + picture 

unfiltered = ndimage.imread(img_name)

tilt = 20.0                                                                     # image tilt in SEM
scale = 0.43                                                                 # SEM image scale

sigma = 1.5                                                                    # Gaussian smoothing
t = 4.0                                                                        # threshold 

step = 50                                                                       # distance between measurements in nm
special = 2
specialcount = special                                                          # number of special features to detect

diam_pix = []
height_pix = []

### Main ######################################################################

if __name__ == "__main__":
    press = 0
    number = 0
    
    outline_pic = plt.figure(3)
    outline_pic_ax = outline_pic.add_subplot(111)
    outline_pic_ax.imshow(unfiltered)
    
    pic_array = filter_pic(img_name, sigma, t, scale, tilt)
    
    thresh_pic = plt.figure(1)
    thresh_pic_ax = thresh_pic.add_subplot(111)
    thresh_pic_ax.imshow(pic_array, interpolation = "None", cmap = 'pink')
    thresh_pic_ax.axis("off")
    thresh_pic.subplots_adjust(left = 0, right = 1, top = 1, bottom = 0)
    plt.draw()
    
    cid = thresh_pic.canvas.mpl_connect('button_press_event', click)            # grab a point
    