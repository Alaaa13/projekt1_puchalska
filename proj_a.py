# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:28:53 2022

@author: admin
"""

from math import sqrt, atan, sin, cos, degrees, radians, tan
import numpy as np
import statistics as st
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        #https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b) / self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
    
    def xyz2flh(self, X, Y, Z):
        '''
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (X, Y, Z)
        na współrzędne geodezyjne: długość szerokość i wysokośc elipsoidalna (fi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. fi można przeliczyć współrzędne z dokładnoscią ok 1 cm.  

        '''
        r = sqrt(X**2 + Y**2)
        fi_prev = atan(Z / (r * (1 - self.ecc2)))
        N = self.a / sqrt(1 - self.ecc2 * (sin(fi_prev)) ** 2)
        h = r / cos(fi_prev) - N
        fi_next = atan((Z / r) * (((1 - self.ecc2 * N / (N + h))**(-1))))
        epsilon = 0.0000001 / 206265
        while abs(fi_prev - fi_next) < epsilon:
            fi_prev = fi_next
            N = self.a / sqrt(1 - self.ecc2 * (sin(fi_prev)) ** 2)
            h = r / cos(fi_prev) - N
            fi_next = atan((Z / r) * (((1 - self.ecc2 * N / (N + h))**(-1))))
        fi = fi_prev
        lam = atan(Y / X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(fi)) ** 2)
        h = r / cos(fi) - N
        return degrees(fi), degrees(lam), h
    
    def flh2xyz(self, fi, lam, h):
        '''
        
        FI,LAMBDA H ->  X,Y,Z - algorytm transformacji współrzędnych geodezyjnych:
        długości, szerokości i wysokości elipsoidalnej (fi, lam, h) 
        na współrzedne ortokartezjańskie (X, Y, Z)
        '''
        
        N = self.a / sqrt(1 - self.ecc2 * (sin(fi)) ** 2)
        X = (N + h) * cos(fi) * cos(lam)
        Y = (N + h) * cos(fi)*sin(lam)
        Z = (N * (1 - self.ecc2) + h) * sin(fi)
        return X, Y, Z
    

    
    def s_A_z2neu(self,s, A, z):
        '''
         s, A, z ->  N,E,U - wyznaczenie wektora neu mając dane:4
         długosc, Azymut, oraz kąt zenitalny

        '''
     
        n = s*np.sin(z)*np.cos(A)
        e = s*np.sin(z)*np.sin(A)
        u = s*np.cos(z)
        return n, e, u
    
    def u2000(self, fi, lam, m_0, h):
        '''
        Układ 2000 - algorytm transformacji współrzędnych geodezyjnych:
        długości, szerokości i wysokości elipsoidalnej (fi, lam, h) 
        na współrzedne płaskie w układzie 2000 (X00, Y00, H)

        '''
        N = self.a /(sqrt(1-self.ecc2 * sin(fi)**2))
        t = tan(fi)
        n2 = self.ecc2 * cos(lam)**2

    
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24
        
        lam = radians(lam)
        lam_0 = radians(lam_0)
        l = lam - lam_0
    
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
    
    
        sigma = self.a* ((A_0*fi) - (A_2*sin(2*fi)) + (A_4*sin(4*fi)) - (A_6*sin(6*fi)))
    
        x = sigma + ((l**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))

        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (s*1000000) + 500000, 3)   
    
        return x00, y00, h 
    
    def u92(self, fi, lam, m_0):
        '''
        Układ 1992 - algorytm transformacji współrzędnych geodezyjnych:
        długości, szerokości i wysokości elipsoidalnej (fi, lam, h) 
        na współrzedne płaskie w układzie 1992 (X92, Y92, H)

        '''
      
        e_2 = self.ecc2/(1-self.ecc2)
        N = self.a/(sqrt(1-self.ecc2 * sin(fi)**2))
        t = tan(fi)
        n2 = e_2 * cos(lam)**2
        lam_0 = radians(19)
        l = lam - lam_0
        
        A_0 = 1 - (self.ecc2/4) - (3*(self.ecc2**2))/64 - (5*(self.ecc2**3))/256
        A_2 = 3/8 * (self.ecc2 + ((self.ecc2**2)/4) + ((15*self.ecc2**3)/128))
        A_4 = 15/256 * (self.ecc2**2 + (3*(self.ecc2**3))/4)
        A_6 = (35*(self.ecc2**3))/3072
        
        sigma = self.a* ((A_0*fi) - (A_2*sin(2*fi)) + (A_4*sin(4*fi)) - (A_6*sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*sin(fi)*cos(fi)) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
        
        return x92, y92 
