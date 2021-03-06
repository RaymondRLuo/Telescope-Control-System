#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import datetime

rad=np.pi/180.0
deg=180.0/np.pi

lgt=116.976
lat=40.558*rad

cat=[]
psrname=[]
ra=[]
dec=[]
el=[]
lh=[]

utc=datetime.datetime.utcnow()
y=utc.year
m=utc.month
d=utc.day
h=np.arange(0,24,0.1)
m=0
s=0

f=open('miyun_list.txt')
fread=f.readlines()
for line in fread:
	dat=line.split()
	cat.append(dat)

for i in range(1,len(cat)):
	raj=float(cat[i][2])+float(cat[i][3])/60.0+float(cat[i][4])/3600.0
	if float(cat[i][5])<0:
		decj=float(cat[i][5])-float(cat[i][6])/60.0-float(cat[i][7])/3600.0
	else:
		decj=float(cat[i][5])+float(cat[i][6])/60.0+float(cat[i][7])/3600.0
	psrname.append(cat[i][1])
	ra.append(raj*15*rad)
	dec.append(decj*rad)

#print lh

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
	az=np.arctan2(np.sin(tau),np.cos(tau)*np.sin(lat)-np.tan(dec)*np.cos(lat))
	if az<0:
		az=az+2*np.pi
	el=np.arcsin(np.cos(lat)*np.cos(dec)*np.cos(tau)+np.sin(lat)*np.sin(dec))
	return az*deg,el*deg

for i in range(len(psrname)):
	ran=ra[i]
	decn=dec[i]
	elt=[]
	for j in range(len(h)):
		elt.append(Equ2Hor(y,m,d,h[j],m,s,ran,decn)[1])
	el.append(elt)

def bottomline(hr):
	elc=30+hr*0
	return elc

def UTC2LT(hr):
	localh=hr+8
	if localh>24:
		localh=localh-24
	return localh

for i in range(len(h)):
	lh.append(UTC2LT(h[i]))

lh=np.array(lh)

plt.title('Elevation Variations of Pulsars Today')
plt.xlabel('Universal Time(hour)')
plt.ylabel('Elevation(deg)')

#for i in range(0,15):
#	plt.plot(h,el[i],'r--',label=psrname[i])

plt.plot(lh,el[55],'r.',label=psrname[55])
plt.plot(lh,el[1],'b.',label=psrname[1])
plt.plot(lh,el[2],'g.',label=psrname[2])
plt.plot(lh,el[3],'m.',label=psrname[3])
plt.plot(lh,el[4],'y.',label=psrname[4])
plt.plot(lh,bottomline(lh),'k--',label='$30^o$ line')
plt.xlim([-1,25])
legend=plt.legend(loc='lower right',shadow=True,fontsize='x-small')
plt.show()