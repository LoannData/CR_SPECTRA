#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 10:20:35 2019

@author: lbrahimi
"""

import numpy as np 
import matplotlib.pyplot as plt 


#function render(x0, y0, dir)
#{ //msg("render: "+x0+","+y0);
#  x_new=x0;
#  y_new=y0;
#  if ((x_old>0)&&(x_new>0))
#  { if (dir=="m") new Line(8*x_old, 8*y_old, 8*x_new, 8*y_new, "ff0000", 1);
#    else new Line(8*x_old, 8*y_old, 8*x_new, 8*y_new, "0000ff", 1);
#  }
#  //new Dot(8*x_new, 8*y_new, 6, 6, "#0000ff");
#  x_old=x_new;
#  y_old=y_new;
#}


def render(pts, x_old, y_old, x_new, y_new, di) : 
#    print pts, x_old, y_old, x_new, y_new
#    if (x_old > 0 and x_new > 0) : 
##        if (di == "m") : 
#        pts.append([x_new, y_new])
    same = False 
    for ii in range(len(pts)) : 
        if (x_new == pts[ii][0] and y_new == pts[ii][1]) : 
            same = True
    if (same == False) : 
#        if ((x_old >= ll and x_new >= ll) and 
#            (y_old >= tt and y_new >= tt) and 
#            (x_old <= ww-1 and x_new <= ww-1) and 
#            (y_old <= hh-1 and y_new <= hh-1)) :
        pts.append([x_new, y_new])
    return pts
#    if ((x_old > 0) and (x_new > 0)) : 
##        if (di == "m") : 
##        else : 
#        print pts
#        pts.append([x_new, y_new])
#        return pts
            

#, pts = pts

def go(x0, y0, dxl, dyl, dxr, dyr, di, pts = []) : 
#    print pts
#    pts = []
#  //x0, y0: start corner looking to the center of the rectangle
#  //dxl, dyl: vector from the start corner to the left corner of the rectangle
#  //dxr, dyr: vector from the start corner to the right corner of the rectangle
#  //dir: direction to go - "l"=left, "m"=middle, "r"=right
#  //msg("go: "+x0+", "+y0+", "+dxl+", "+dyl+", "+dxr+", "+dyr+", "+dir);
#    print abs((dxl+dyl)*(dxr+dyr))
    if (abs((dxl+dyl)*(dxr+dyr))<=6) : 
#        print "abs((dxl+dyl)*(dxr+dyr))<=6"
        if (abs(dxl + dyl) == 1) : 
#            print ("(abs(dxl + dyl) == 1) ")
            ddx = dxr/abs(dxr+dyr)
            ddy = dyr/abs(dxr+dyr)
            for ii in range(int(abs(dxr + dyr))) : 
                pts = render(pts, x0, y0, x0+ii*ddx+(dxl+ddx-1)/2, y0+ii*ddy+(dyl+ddy-1)/2, di)
            return pts
#            return pts
        if (abs(dxr + dyr) == 1) :
#            print " (abs(dxr + dyr) == 1) :"
            ddx = dxl/abs(dxl + dyl)
            ddy = dyl/abs(dxl + dyl)
            for ii in range(int(abs(dxl + dyl))) : 
                pts = render(pts, x0, y0, x0+ii*ddx+(dxr+ddx-1)/2, y0+ii*ddy+(dyr+ddy-1)/2, di)
            return pts
#            return pts
        if (di == "l") : 
#            print "(di == l) : "
            ddx = dxr/abs(dxr + dyr)
            ddy = dyr/abs(dxr + dyr)
            for ii in range(int(abs(dxr+dyr))) : 
                pts = render(pts, x0, y0, x0+ii*ddx+(dxl/2+ddx-1)/2, y0+ii*ddy+(dyl/2+ddy-1)/2, di)
            for ii in range(int(abs(dxr+dyr))-1, -1, -1) : 
                pts = render(pts, x0, y0, x0+ii*ddx+(dxl+dxl/2+ddx-1)/2, y0+ii*ddy+(dyl+dyl/2+ddy-1)/2, di) 
            return pts
#            return pts
        if (di == "r") : 
#            print "(di == r) : "
            ddx = dxl/abs(dxl + dyl)
            ddy = dyl/abs(dxl + dyl)
            for ii in range(0, int(abs(dxl + dyl))+1, 1) : 
#                pts = render(pts, x0, y0, x0+ii*ddx+(dxr/2+ddx-1)/2, y0+ii*ddy+(dyr/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+ii*ddx+(dxr/2+ddx-1)/2, y0+ii*ddy+(dyr/2+ddy-1)/2, di)
            for ii in range(int(abs(dxl + dyl)), -1, -1) : 
                pts = render(pts, x0, y0, x0+ii*ddx+(dxr+dxr/2+ddx-1)/2, y0+ii*ddy+(dyr+dyr/2+ddy-1)/2, di) 
            return pts
#            return pts
        if (di == "m") : 
#            print abs(dxr + dyr)
            if (abs(dxr + dyr) == 3) : 
                ddx=dxr/abs(dxr+dyr)
                ddy=dyr/abs(dxr+dyr)
                pts = render(pts, x0, y0, x0+(dxl/2+ddx-1)/2, y0+(dyl/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+(dxl+dxl/2+ddx-1)/2, y0+(dyl+dyl/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+ddx+(dxl+dxl/2+ddx-1)/2, y0+ddy+(dyl+dyl/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+ddx+(dxl/2+ddx-1)/2, y0+ddy+(dyl/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+2*ddx+(dxl/2+ddx-1)/2, y0+2*ddy+(dyl/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+2*ddx+(dxl+dxl/2+ddx-1)/2, y0+2*ddy+(dyl+dyl/2+ddy-1)/2, di)
                return pts
            if (abs(dxl + dyl) == 3) : 
                ddx=dxl/abs(dxl+dyl);
                ddy=dyl/abs(dxl+dyl);
                pts = render(pts, x0, y0, x0+(dxr/2+ddx-1)/2, y0+(dyr/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+(dxr+dxr/2+ddx-1)/2, y0+(dyr+dyr/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+ddx+(dxr+dxr/2+ddx-1)/2, y0+ddy+(dyr+dyr/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+ddx+(dxr/2+ddx-1)/2, y0+ddy+(dyr/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+2*ddx+(dxr/2+ddx-1)/2, y0+2*ddy+(dyr/2+ddy-1)/2, di)
                pts = render(pts, x0, y0, x0+2*ddx+(dxr+dxr/2+ddx-1)/2, y0+2*ddy+(dyr+dyr/2+ddy-1)/2, di)
                return pts
#            return pts
        
#        if (pts) : print pts
#        else : 
#            print "Problem"
        return pts
    
    # Divide into 2 parts if necessary
    if (2*(abs(dxl)+abs(dyl))>3*(abs(dxr)+abs(dyr))) : #left side much longer than right side
#        print "2*(abs(dxl)+abs(dyl))>3*(abs(dxr)+abs(dyr))) : #left side much longer than right side"
        dxl2 = np.round(dxl/2.)
        dyl2 = np.round(dyl/2.)
        if ((abs(dxr)+abs(dyr))%2==0) :  #right side is even
            if ((abs(dxl)+abs(dyl))%2==0) : #make 2 parts from even side
                if (di =="l") : 
                    if ((abs(dxl)+abs(dyl))%4==0) : #make 2 parts even-even from even side
                        pts = go(x0, y0, dxl2, dyl2, dxr, dyr, "l", pts = pts)
                        pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "l", pts = pts)
                        return pts
                    else  : #make 2 parts odd-odd from even side
                        pts = go(x0, y0, dxl2, dyl2, dxr, dyr, "m", pts = pts)
                        pts = go(x0+dxl2+dxr, y0+dyl2+dyr, -dxr, -dyr, dxl-dxl2, dyl-dyl2, "m", pts = pts)
                        return pts
            else : #make 2 parts from odd side
                if (di == "m") : 
                    if ((abs(dxl2)+abs(dyl2))%2==0) : 
                        pts = go(x0, y0, dxl2, dyl2, dxr, dyr, "l", pts = pts)
                        pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "m", pts = pts)
                        return pts
                    else : 
                        pts = go(x0, y0, dxl2, dyl2, dxr, dyr, "m", pts = pts)
                        pts = go(x0+dxl2+dxr, y0+dyl2+dyr, -dxr, -dyr, dxl-dxl2, dyl-dyl2, "r", pts = pts)
                        return pts
        else :  #right side is odd
            if (di=="l")  : 
                pts = go(x0, y0, dxl2, dyl2, dxr, dyr, "l", pts = pts)
                pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "l", pts = pts)
                return pts
            if (di=="m") : 
                pts = go(x0, y0, dxl2, dyl2, dxr, dyr, "l", pts = pts)
                pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr, dyr, "m", pts = pts)   
                return pts
    if (2*(abs(dxr)+abs(dyr))>3*(abs(dxl)+abs(dyl))) : #right side much longer than left side
#        print "2*(abs(dxr)+abs(dyr))>3*(abs(dxl)+abs(dyl))) : #right side much longer than left side"
        dxr2 = np.round(dxr/2.)
        dyr2 = np.round(dyr/2.)
        if ((abs(dxl)+abs(dyl))%2==0) :  #left side is even
#            print "left_side_even"
            if ((abs(dxr)+abs(dyr))%2==0) : #make 2 parts from even side
#                print "(abs(dxr)+abs(dyr))%2==0) : #make 2 parts from even side"
                if (di=="r") : 
                    if ((abs(dxr)+abs(dyr))%4==0) : #make 2 parts even-even from even side
#                        print "((abs(dxr)+abs(dyr))%4==0)"
                        pts = go(x0, y0, dxl, dyl, dxr2, dyr2, "r", pts = pts)
                        pts = go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "r", pts = pts)
                        return pts
                    else  :#make 2 parts odd-odd from even side
                        pts = go(x0, y0, dxl, dyl, dxr2, dyr2, "m", pts = pts)
                        pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxr-dxr2, dyr-dyr2, -dxl, -dyl, "m", pts = pts)  
                        return pts
            else : 
                if (di == "m") : 
                    if ((abs(dxr2)+abs(dyr2))%2==0) : 
                        pts = go(x0, y0, dxl, dyl, dxr2, dyr2, "r", pts = pts)
                        pts = go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "m", pts = pts)
                        return pts
          
                    else : 
                        pts = go(x0, y0, dxl, dyl, dxr2, dyr2, "m", pts = pts)
                        pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxr-dxr2, dyr-dyr2, -dxl, -dyl, "l", pts = pts)
                        return pts
        else : #left side is odd
#            print "left_side_odd"
            if (di == "r") : 
                pts = go(x0, y0, dxl, dyl, dxr2, dyr2, "r", pts = pts)
                pts = go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "r", pts = pts)
                return pts
            if (di == "m") : 
                pts = go(x0, y0, dxl, dyl, dxr2, dyr2, "r", pts = pts)
                pts = go(x0+dxr2, y0+dyr2, dxl, dyl, dxr-dxr2, dyr-dyr2, "m", pts = pts) 
                return pts
    # Divide into 2x2 parts
    if ((di=="l") or (di=="r")) : 
        dxl2 = np.round(dxl/2.)
        dyl2 = np.round(dyl/2.)
        dxr2 = np.round(dxr/2.)
        dyr2 = np.round(dyr/2.)
        if ((abs(dxl+dyl)%2==0) and (abs(dxr+dyr)%2==0)) : #even-even
            if (abs(dxl2+dyl2+dxr2+dyr2)%2==0) : #ee-ee or oo-oo
                if (di == "l") : 
#                    print pts
                    pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r", pts = pts)
#                    print pts
                    pts = go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "l", pts = pts)
#                    print pts
                    pts = go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "l", pts = pts)
                    pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r", pts = pts)
                    return pts
                else : 
                    pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l", pts = pts)
                    pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "r", pts = pts)
                    pts = go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "r", pts = pts)
                    pts = go(x0+dxr+dxl2, y0+dyr+dyl2, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l", pts = pts)   
                    return pts
            else :  #ee-oo or oo-ee
                if ((dxr2+dyr2)%2==0) : #ee-oo
                    if (di == "l") : 
                        pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r", pts = pts)
                        pts = go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m", pts = pts)
                        pts = go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m", pts = pts)
                        pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r", pts = pts)
                        return pts
                    else : #ee-oo for dir="r" not possible, so transforming into e-1,e+1-oo = oo-oo 
                        if (dxr2 != 0) : dxr2 += 1
                        else : dyr2 += 1
                        pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l", pts = pts)
                        pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m", pts = pts)
                        pts = go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m", pts = pts)
                        pts = go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l", pts = pts)
                        return pts
                else : #oo-ee
                    if (di == "r") : 
                        pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l", pts = pts)
                        pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m", pts = pts)
                        pts = go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m", pts = pts)
                        pts = go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l", pts = pts)
                        return pts
                    else : #oo-ee for dir="l" not possible, so transforming into oo-e-1,e+1 = oo-oo 
                        if (dxl2 != 0) : dxl2 += 1
                        else : dyl2 += 1
                        pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r", pts = pts)
                        pts = go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m", pts = pts)
                        pts = go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m", pts = pts)
                        pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r", pts = pts)
                        return pts
        else : #not even-even
            if ((abs(dxl+dyl)%2!=0) and (abs(dxr+dyr)%2!=0)) : #odd-odd
                if (dxl2%2!=0): dxl2=dxl-dxl2 #get it in a shape eo-eo 
                if (dyl2%2!=0): dyl2=dyl-dyl2
                if (dxr2%2!=0): dxr2=dxr-dxr2
                if (dyr2%2!=0): dyr2=dyr-dyr2
                if (di=="l") : 
                    pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r", pts = pts)
                    pts = go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m", pts = pts)
                    pts = go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m", pts = pts)
                    pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r", pts = pts)
                    return pts
                else : 
                    pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l", pts = pts)
                    pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m", pts = pts)
                    pts = go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m", pts = pts)
                    pts = go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l", pts = pts)
                    return pts
            else :  #even-odd or odd-even
                if (abs(dxl+dyl)%2==0) : #odd-even
                    if (di == "l") : 
                        if (dxr2%2!=0) :  dxr2=dxr-dxr2 #get it in a shape eo-xx 
                        if (dyr2%2!=0) :  dyr2=dyr-dyr2
                        if (abs(dxl+dyl)>2) : 
                            pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r", pts = pts)
                            pts = go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "l", pts = pts)
                            pts = go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "l", pts = pts)
                            pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r", pts = pts)
                            return pts
                        else : 
                            pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "r", pts = pts)
                            pts = go(x0+dxr2, y0+dyr2, dxl2, dyl2, dxr-dxr2, dyr-dyr2, "m", pts = pts)
                            pts = go(x0+dxr+dxl2, y0+dyr+dyl2, dxr2-dxr, dyr2-dyr, dxl-dxl2, dyl-dyl2, "m", pts = pts)
                            pts = go(x0+dxr2+dxl, y0+dyr2+dyl, dxl2-dxl, dyl2-dyl, -dxr2, -dyr2, "r", pts = pts)
                            return pts
                else : #even-odd
                    if (di == "r") : 
                        if (dxl2%2!=0): dxl2=dxl-dxl2  #get it in a shape xx-eo 
                        if (dyl2%2!=0): dyl2=dyl-dyl2
                        if (abs(dxr+dyr)>2) : 
                            pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l", pts = pts)
                            pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "r", pts = pts)
                            pts = go(x0+dxr2+dxl2, y0+dyr2+dyl2, dxl-dxl2, dyl-dyl2, dxr-dxr2, dyr-dyr2, "r", pts = pts)
                            pts = go(x0+dxr+dxl2, y0+dyr+dyl2, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l", pts = pts)
                            return pts
                        else : 
                            pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "l", pts = pts)
                            pts = go(x0+dxl2, y0+dyl2, dxl-dxl2, dyl-dyl2, dxr2, dyr2, "m", pts = pts)
                            pts = go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-dxr2, dyr-dyr2, dxl2-dxl, dyl2-dyl, "m", pts = pts)
                            pts = go(x0+dxl2+dxr, y0+dyl2+dyr, -dxl2, -dyl2, dxr2-dxr, dyr2-dyr, "l", pts = pts)
                            return pts
    else : #dir=="m" -> divide into 3x3 parts
        if ((abs(dxl+dyl)%2==0) and (abs(dxr+dyr)%2==0)) : 
            print "Error"
        if (abs(dxr+dyr)%2==0) :  #even-odd: oeo-ooo
            dxl2 = np.round(dxl/3.)
            dyl2 = np.round(dyl/3.)
            dxr2 = np.round(dxr/3.)
            dyr2 = np.round(dyr/3.)
            if ((dxl2+dyl2)%2==0) : #make it odd  
                dxl2=dxl-2*dxl2
                dyl2=dyl-2*dyl2
            if ((dxr2+dyr2)%2==0) : #make it odd (not necessary, however results are better for 12x30, 18x30 etc.) 
                if (abs(dxr2+dyr2)!=2) : 
                    if (dxr < 0) : dxr2 += 1
                    if (dxr > 0) : dxr2 -= 1
                    if (dyr < 0) : dyr2 +=1
                    if (dyr > 0) : dyr2 -=1
        else : 
#            print dxl
            dxl2 = np.round(dxl/3.)
            dyl2 = np.round(dyl/3.)
            dxr2 = np.round(dxr/3.)
            dyr2 = np.round(dyr/3.)
            if ((dxr2+dyr2)%2==0) : #make it odd  
                dxr2=dxr-2*dxr2
                dyr2=dyr-2*dyr2
            if ((dxl2+dyl2)%2==0) : #make it odd (not necessary, however results are better for 12x30, 18x30 etc.)
                if (abs(dxl2+dyl2)!=2) : 
                    if (dxl < 0) : dxl2 += 1
                    if (dxl > 0) : dxl2 -= 1
                    if (dyl < 0) : dyl2 += 1
                    if (dyl > 0) : dyl2 -= 1
        if (abs(dxl+dyl)<abs(dxr+dyr)) : 
            pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            pts = go(x0+dxl2+dxr2, y0+dyl2+dyr2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m", pts = pts)
            pts = go(x0+dxl-dxl2, y0+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            pts = go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m", pts = pts)
            pts = go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, 2*dxl2-dxl, 2*dyl2-dyl, 2*dxr2-dxr, 2*dyr2-dyr, "m", pts = pts)
            pts = go(x0+dxl2+dxr2, y0+dyl2+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m", pts = pts)
            pts = go(x0+dxr-dxr2, y0+dyr-dyr2, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            pts = go(x0+dxr+dxl2, y0+dyr+dyl2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m", pts = pts)
            pts = go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            return pts
        else : 
            pts = go(x0, y0, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            pts = go(x0+dxl2+dxr2, y0+dyl2+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m", pts = pts)
            pts = go(x0+dxr-dxr2, y0+dyr-dyr2, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            pts = go(x0+dxr+dxl2, y0+dyr+dyl2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m", pts = pts)
            pts = go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, 2*dxl2-dxl, 2*dyl2-dyl, 2*dxr2-dxr, 2*dyr2-dyr, "m", pts = pts)
            pts = go(x0+dxl2+dxr2, y0+dyl2+dyr2, -dxr2, -dyr2, dxl-2*dxl2, dyl-2*dyl2, "m", pts = pts)
            pts = go(x0+dxl-dxl2, y0+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            pts = go(x0+dxl+dxr2, y0+dyl+dyr2, dxr-2*dxr2, dyr-2*dyr2, -dxl2, -dyl2, "m", pts = pts)
            pts = go(x0+dxr-dxr2+dxl-dxl2, y0+dyr-dyr2+dyl-dyl2, dxl2, dyl2, dxr2, dyr2, "m", pts = pts)
            return pts
#    return pts
            
        

def spacefill(ll, tt, ww, hh) : 
    pts = []
    if (hh > ww) : # Go top -> Down
        if ((hh%2 == 1) and (ww%2 == 0)) : 
            pts = go(ll, tt, ww, 0, 0, hh, "m", pts = pts) # Go diagonal
#            print pts
        else : 
#            print "# Go top -> down"
            pts = go (ll, tt, ww, 0, 0, hh, "r", pts = pts) # Go top -> down
    else : # Go left -> right
        if ((ww % 2 == 1) and (hh % 2) == 0) : 
#            print "no"
            pts = go(ll, tt, ww, 0, 0, hh, "m", pts = pts) # Go diagonal
        else : 
#            print "hello"
            pts = go(ll, tt, ww, 0, 0, hh, "l", pts = pts) # Go left -> right
    return pts

def lineCorrector(pts, ll, tt, ww, hh) :
    ptnew = []
    for ii in range(len(pts)) : 
        x = pts[ii][0]
        y = pts[ii][1]
        if (x >= ll and y >= tt and x <= ww-ll-1 and y <= hh-tt-1) : 
            ptnew.append(pts[ii])
    return ptnew

def getRHilbert2D(ll, tt, ww, hh) : 
    pts = spacefill(ll, tt, ww, hh)
    pts = lineCorrector(pts, ll, tt, ww, hh)
    return pts

"""       
#h > w
ll = 0.
tt = 0.

ww = 512
hh = 64
pts = spacefill(ll, tt, ww, hh)
pts = lineCorrector(pts, ll, tt, ww, hh)





x = []
y = []
for ii in range(len(pts)) : 
    x.append(pts[ii][0])
    y.append(pts[ii][1])


fig_title = "./hilbert.png"
plt.figure(figsize=(int(ww/5.),int(hh/5.)))

plt.plot(x, y)

for ii in range(0, ww, 1) : 
    for jj in range(0, hh, 1) : 
        plt.plot(ll+ii, tt+jj, marker='o', c="black", ms=3.5)

plt.savefig(fig_title)
plt.show()
"""
