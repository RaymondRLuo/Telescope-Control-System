import sys
import numpy as np
import time

rad=np.pi/180.0
deg=180.0/np.pi

# Beijing
#lgt=116.47
#lat=39.92*rad

# Miyun
lgt=116.976
lat=40.558*rad

def start():
	print 1

def stop():
	print 0

def UTC2MJD(Y,M,D,h,m,s):
	ay=Y
	am=M
	if M<=2:
		ay=Y-1
		am=M+12
	jdn=np.floor(365.25*(ay+4716))+np.floor(30.6001*(am+1))+2.0-np.floor(ay/100)+np.floor(np.floor(ay/100)/4)+D-1524
	jd=jdn+(h-12)/24.0+m/1440.0+s/86400.0
	return jd-2400000.5

def UTC2LST(Y,M,D,h,m,s,lgt):
	mjd=UTC2MJD(Y,M,D,h,m,s)
	mjd_0=UTC2MJD(Y,M,D,0,0,0)
	t0=(mjd_0-51544.5)/36525.0
	t=(mjd-51544.5)/36525.0
	gmt=24110.54841+8640184.812866*t0+(0.093104-0.0000062*t)*t*t+1.0027379093*(h*3600+m*60+s)
	gmt=(gmt/3600)%24
	lst=gmt+lgt/15
	if lst>=24:
		lst=lst-24
	if lst<0:
		lst=lst+24
	return lst

def Equ2Hor(Y,M,D,h,m,s,ra,dec):
	#ut=h*3600+m*60+s
	lst=UTC2LST(Y,M,D,h,m,s,lgt)*15*rad
	tau=lst-ra
	#az=np.arctan(np.sin(tau)/(np.cos(tau)*np.sin(lat)-np.tan(dec)*np.cos(lat)))
	az=np.arctan2(np.sin(tau),np.cos(tau)*np.sin(lat)-np.tan(dec)*np.cos(lat))
	if az<0:
		az=az+2*np.pi
	el=np.arcsin(np.cos(lat)*np.cos(dec)*np.cos(tau)+np.sin(lat)*np.sin(dec))
	return az*deg,el*deg

def CalRStime(ra,dec):
	#if abs(np.tan(lat)*np.tan(dec))<=1:
	lst_r=ra-np.arccos(-np.tan(lat)*np.tan(dec))
	if lst_r<0:
		lst_r=lst_r+2*np.pi
	if lst_r>=2*np.pi:
		lst_r=lst_r-2*np.pi
	lst_s=ra+np.arccos(-np.tan(lat)*np.tan(dec))
	if lst_s<0:
		lst_s=lst_s+2*np.pi
	if lst_s>=2*np.pi:
		lst_s=lst_s-2*np.pi
	return lst_r*deg/15.0,lst_s*deg/15.0

def MJDdisplay(Y,M,D,h,m,s):
	mjd=UTC2MJD(Y,M,D,h,m,s)
	mjd='%.5f'%mjd
	return '  '+str(mjd)

def LSTdisplay(Y,M,D,h,m,s):
	#ut=h*3600+m*60+s
	lst=UTC2LST(Y,M,D,h,m,s,lgt)
	lst_h=int(np.floor(lst))
	lst_m=int(np.floor((lst-lst_h)*60))
	lst_s=int(np.floor(((lst-lst_h)*60-lst_m)*60))
	lst=str(lst_h)+':'+str(lst_m)+':'+str(lst_s)
	lst_struc=time.strptime(lst,'%H:%M:%S')
	lst=time.strftime('  %H:%M:%S',lst_struc)
	return lst

def Hordisplay(Y,M,D,h,m,s,ra,dec):
	#ra=ra*15*rad
	#dec=dec*rad
	az=Equ2Hor(Y,M,D,h,m,s,ra,dec)[0]
	el=Equ2Hor(Y,M,D,h,m,s,ra,dec)[1]
	return ' '+str('%.3f'%az), ' '+str('%.3f'%el)

def Posdisplay(Y,M,D,h,m,s,ra,dec):
	#ra=ra*15*rad
	#dec=dec*rad
	el=Equ2Hor(Y,M,D,h,m,s,ra,dec)[1]
	if el>=0:
		return 'Over the horizon.'
	else:
		return 'Under the horizon.'

def RStimedisplay(Y,M,D,h,m,s,ra,dec):
	#ra=ra*15*rad
	#dec=dec*rad
	if abs(np.tan(lat)*np.tan(dec))<=1:
		rt=CalRStime(ra,dec)[0]
		st=CalRStime(ra,dec)[1]
		rt_h=int(np.floor(rt))
		rt_m=int(np.floor((rt-rt_h)*60))
		rt_s=int(np.floor(((rt-rt_h)*60-rt_m)*60))
		rt=str(rt_h)+':'+str(rt_m)+':'+str(rt_s)
		st_h=int(np.floor(st))
		st_m=int(np.floor((st-st_h)*60))
		st_s=int(np.floor(((st-st_h)*60-st_m)*60))
		st=str(st_h)+':'+str(st_m)+':'+str(st_s)
		rt_struc=time.strptime(rt,'%H:%M:%S')
		rt=time.strftime('  %H:%M:%S',rt_struc)
		st_struc=time.strptime(st,'%H:%M:%S')
		st=time.strftime('  %H:%M:%S',st_struc)
		return rt,st
	else:
		if Equ2Hor(Y,M,D,h,m,s,ra,dec)[1]>=0:
			return 'Never set','Never set'
		if Equ2Hor(Y,M,D,h,m,s,ra,dec)[1]<0:
			return 'Never rise','Never rise'

def RStimerdisplay(Y,M,D,h,m,s,ra,dec):
	#ra=ra*15*rad
	#dec=dec*rad
	if abs(np.tan(lat)*np.tan(dec))<=1:
		#ut=h*3600+m*60+s
		lst=UTC2LST(Y,M,D,h,m,s,lgt)
		if CalRStime(ra,dec)[0]<lst:
			rt_rm=CalRStime(ra,dec)[0]+24.0-lst
		else:
			rt_rm=CalRStime(ra,dec)[0]-lst
		if CalRStime(ra,dec)[1]<lst:
			st_rm=CalRStime(ra,dec)[1]+24.0-lst
		else:
			st_rm=CalRStime(ra,dec)[1]-lst
		rt_rm_h=int(np.floor(rt_rm))
		rt_rm_m=int(np.floor((rt_rm-rt_rm_h)*60))
		rt_rm_s=int(np.floor(((rt_rm-rt_rm_h)*60-rt_rm_m)*60))
		rt_rm=str(rt_rm_h)+':'+str(rt_rm_m)+':'+str(rt_rm_s)
		st_rm_h=int(np.floor(st_rm))
		st_rm_m=int(np.floor((st_rm-st_rm_h)*60))
		st_rm_s=int(np.floor(((st_rm-st_rm_h)*60-st_rm_m)*60))
		st_rm=str(st_rm_h)+':'+str(st_rm_m)+':'+str(st_rm_s)
		rt_rm_struc=time.strptime(rt_rm,'%H:%M:%S')
		rt_rm=time.strftime('  %H:%M:%S',rt_rm_struc)
		st_rm_struc=time.strptime(st_rm,'%H:%M:%S')
		st_rm=time.strftime('  %H:%M:%S',st_rm_struc)
		return rt_rm,st_rm
	else:
		return 'None','None'

