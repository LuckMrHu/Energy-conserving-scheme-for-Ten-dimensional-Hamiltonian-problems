program main
implicit none
character(len=20)::file1="error.dat"
character(len=20)::file2="xyz.dat"
character(len=20)::file3="r.dat"
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5
real*8 n,c1,c2,w,e,pe
real*8 h,endt
real*8::t=0d0,init=0d0,overt=0d0
real*8 inih,inil,temph,templ,r,fh,fl
real*8 inid,tempd,fli,el,num
integer i,j,k
integer::count=0d0
open(10,file=file1)
open(20,file=file2)
open(30,file=file3)
num=0d0       
el=0d0
pe=dacos(-1d0)
call initiala(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,e,inih,inil,h,endt,j,k,inid)
pause
call cpu_time(init)
1 do 2 i=1d0,j,k
  call step(t,h)
  !call rks(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
  !call iss(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
  !call ems(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
  call ecf(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
2 continue
  if(t.gt.1d2)then
  j=10d0
  else
  end if
  if(t.gt.1d3)then
  j=100d0
  else
  end if
  if(t.gt.1d4)then
  j=1000d0
  else
  end if
  if(t.gt.1d5)then
  j=10000d0
  else
  end if
  if(t.gt.1d6)then
  j=50000d0
  else
  end if
  if(t.gt.1d7)then
  j=100000d0
  else
  end if
  if(t.gt.1d8)then
  j=500000d0
  else
  end if
call hamiltonian(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,temph)
call angularmomentum(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,templ)
r=dsqrt(sx**2d0+sy**2d0+sz**2d0)
fh=dlog10(dabs(temph-inih)+1d-16)
fl=dlog10(dabs(templ-inil)+1d-16)
write(10,*)dlog10(t),fh,fl
write(20,*)sx,sy,sz
write(30,*)t,r
count=count+1d0
write(*,"(F7.4)")t/endt
if(t.lt.endt)goto 1
call cpu_time(overt)
write(*,*)dabs(overt-init)
write(*,*)count
write(*,*)num
pause
close(10)
close(20)
close(30)
stop
end

subroutine initiala(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,e,hami,mome,h,endt,j,k,inid)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,e,hami,mome,h,endt,inid
integer j,k
j=1d0
k=1d0
h=0.1d0
endt=1d7
c1=1d0
c2=1d0
w=1d0
n=1d0/(w+1d0/w+2d0)
inid=1d-8
call inivar(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,e)
call hamiltonian(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,hami)
call angularmomentum(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,mome)
write(*,*)hami
write(*,*)mome
return
end subroutine

subroutine inivar(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,e)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,e
real*8 pe
pe=dacos(-1d0)
e=0.0985d0
sx=40d0
sy=0d0
sz=0d0
cta=pe/4d0
ctb=pe/4d0
px=0d0
py=dsqrt((1d0-e)/sx)
!py=6.64731237968498d0/40d0
!py=12.4117721412935d0/150d0
pz=0d0
xia=0.1d0
xib=0.1d0
return
end subroutine

subroutine step(t,h)
implicit none
real*8 t,h
t=t+h
return
end subroutine

subroutine rks(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h
real*8 sx0,sy0,sz0,cta0,ctb0,px0,py0,pz0,xia0,xib0
real*8::a=1d0/2d0
real*8 k1sx,k1sy,k1sz,k1cta,k1ctb,k1px,k1py,k1pz,k1xia,k1xib
real*8 k2sx,k2sy,k2sz,k2cta,k2ctb,k2px,k2py,k2pz,k2xia,k2xib
real*8 dx,dy,dz,dcta,dctb,dpx,dpy,dpz,dxia,dxib
sx0=sx
sy0=sy
sz0=sz
cta0=cta
ctb0=ctb
px0=px
py0=py
pz0=pz
xia0=xia
xib0=xib
call xdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dx)
call ydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dy)
call zdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dz)
call ctadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dcta)
call ctbdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dctb)
call pxdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpx)
call pydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpy)
call pzdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpz)
call xiadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxia)
call xibdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxib)
k1sx=h*dx
k1sy=h*dy
k1sz=h*dz
k1cta=h*dcta
k1ctb=h*dctb
k1px=h*dpx
k1py=h*dpy
k1pz=h*dpz
k1xia=h*dxia
k1xib=h*dxib
sx=sx0+k1sx
sy=sy0+k1sy
sz=sz0+k1sz
cta=cta0+k1cta
ctb=ctb0+k1ctb
px=px0+k1px
py=py0+k1py
pz=pz0+k1pz
xia=xia0+k1xia
xib=xib0+k1xib
call xdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dx)
call ydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dy)
call zdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dz)
call ctadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dcta)
call ctbdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dctb)
call pxdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpx)
call pydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpy)
call pzdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpz)
call xiadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxia)
call xibdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxib)
k2sx=h*dx
k2sy=h*dy
k2sz=h*dz
k2cta=h*dcta
k2ctb=h*dctb
k2px=h*dpx
k2py=h*dpy
k2pz=h*dpz
k2xia=h*dxia
k2xib=h*dxib
sx=sx0+a*(k1sx+k2sx)
sy=sy0+a*(k1sy+k2sy)
sz=sz0+a*(k1sz+k2sz)
cta=cta0+a*(k1cta+k2cta)
ctb=ctb0+a*(k1ctb+k2ctb)
px=px0+a*(k1px+k2px)
py=py0+a*(k1py+k2py)
pz=pz0+a*(k1pz+k2pz)
xia=xia0+a*(k1xia+k2xia)
xib=xib0+a*(k1xib+k2xib)
return
end subroutine

subroutine iss(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h
real*8 sx0,sy0,sz0,cta0,ctb0,px0,py0,pz0,xia0,xib0
real*8::error=1d-15,d=0d0
real*8 sxm,sym,szm,ctam,ctbm,pxm,pym,pzm,xiam,xibm
real*8 z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
real*8 dx,dy,dz,dcta,dctb,dpx,dpy,dpz,dxia,dxib
integer::j=2000d0
integer i
i=0d0
sx0=sx
sy0=sy
sz0=sz
cta0=cta
ctb0=ctb
px0=px
py0=py
pz0=pz
xia0=xia
xib0=xib
sxm=sx
sym=sy
szm=sz
ctam=cta
ctbm=ctb
pxm=px
pym=py
pzm=pz
xiam=xia
xibm=xib
1 z1=sxm
  z2=sym
  z3=szm
  z4=ctam
  z5=ctbm
  z6=pxm
  z7=pym
  z8=pzm
  z9=xiam
  z10=xibm
  call xdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dx)
  sxm=sx0+h*dx
  sx=(sx0+sxm)/2d0
  call ydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dy)
  sym=sy0+h*dy
  sy=(sy0+sym)/2d0
  call zdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dz)
  szm=sz0+h*dz
  sz=(sz0+szm)/2d0
  call ctadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dcta)
  ctam=cta0+h*dcta
  cta=(cta0+ctam)/2d0
  call ctbdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dctb)
  ctbm=ctb0+h*dctb
  ctb=(ctb0+ctbm)/2d0
  call pxdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpx)
  pxm=px0+h*dpx
  px=(px0+pxm)/2d0
  call pydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpy)
  pym=py0+h*dpy
  py=(py0+pym)/2d0
  call pzdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpz)
  pzm=pz0+h*dpz
  pz=(pz0+pzm)/2d0
  call xiadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxia)
  xiam=xia0+h*dxia
  xia=(xia0+xiam)/2d0
  call xibdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxib)
  xibm=xib0+h*dxib
  xib=(xib0+xibm)/2d0
  i=i+1d0
  d=dsqrt((z1-sxm)**2d0+(z2-sym)**2d0+(z3-szm)**2d0+(z4-ctam)**2d0+(z5-ctbm)**2d0+(z6-pxm)**2d0+(z7-pym)**2d0+(z8-pzm)**2d0+(z9-xiam)**2d0+(z10-xibm)**2d0+1d-32)
  !write(*,*)d
  !pause
  if(i.gt.j)then
  goto 2
  else
  end if
  if(d.gt.error)goto 1
2 sx=sxm
  sy=sym
  sz=szm
  cta=ctam
  ctb=ctbm
  px=pxm
  py=pym
  pz=pzm
  xia=xiam
  xib=xibm
return
end subroutine

subroutine ems(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,h
real*8 sx0,sy0,sz0,cta0,ctb0,px0,py0,pz0,xia0,xib0
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5
real*8 q10,q20,q30,q40,q50,p10,p20,p30,p40,p50
real*8 dq1,dq2,dq3,dq4,dq5,dp1,dp2,dp3,dp4,dp5
real*8 dx,dy,dz,dcta,dctb,dpx,dpy,dpz,dxia,dxib
sx0=sx
sy0=sy
sz0=sz
cta0=cta
ctb0=ctb
px0=px
py0=py
pz0=pz
xia0=xia
xib0=xib
q1=sx
q2=sy
q3=sz
q4=cta
q5=ctb
p1=px
p2=py
p3=pz
p4=xia
p5=xib
q10=q1
q20=q2
q30=q3
q40=q4
q50=q5
p10=p1
p20=p2
p30=p3
p40=p4
p50=p5
call xdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dx)
call ydot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dy)
call zdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dz)
call ctadot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dcta)
call ctbdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dctb)
call pxdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp1)
call pydot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp2)
call pzdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp3)
call xiadot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp4)
call xibdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp5)
sx=sx0+(1d0/2d0)*h*dx
sy=sy0+(1d0/2d0)*h*dy
sz=sz0+(1d0/2d0)*h*dz
cta=cta0+(1d0/2d0)*h*dcta
ctb=ctb0+(1d0/2d0)*h*dctb
p1=p10+(1d0/2d0)*h*dp1
p2=p20+(1d0/2d0)*h*dp2
p3=p30+(1d0/2d0)*h*dp3
p4=p40+(1d0/2d0)*h*dp4
p5=p50+(1d0/2d0)*h*dp5
call xdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dq1)
call ydot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dq2)
call zdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dq3)
call ctadot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dq4)
call ctbdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dq5)
call pxdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dpx)
call pydot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dpy)
call pzdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dpz)
call xiadot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dxia)
call xibdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,n,c1,c2,w,dxib)
px=px0+h*dpx
py=py0+h*dpy
pz=pz0+h*dpz
xia=xia0+h*dxia
xib=xib0+h*dxib
q1=q10+h*dq1
q2=q20+h*dq2
q3=q30+h*dq3
q4=q40+h*dq4
q5=q50+h*dq5
call xdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dx)
call ydot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dy)
call zdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dz)
call ctadot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dcta)
call ctbdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dctb)
call pxdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp1)
call pydot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp2)
call pzdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp3)
call xiadot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp4)
call xibdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,n,c1,c2,w,dp5)
sx=sx+(1d0/2d0)*h*dx
sy=sy+(1d0/2d0)*h*dy
sz=sz+(1d0/2d0)*h*dz
cta=cta+(1d0/2d0)*h*dcta
ctb=ctb+(1d0/2d0)*h*dctb
p1=p1+(1d0/2d0)*h*dp1
p2=p2+(1d0/2d0)*h*dp2
p3=p3+(1d0/2d0)*h*dp3
p4=p4+(1d0/2d0)*h*dp4
p5=p5+(1d0/2d0)*h*dp5
sx=(sx+q1)/2d0
sy=(sy+q2)/2d0
sz=(sz+q3)/2d0
cta=(cta+q4)/2d0
ctb=(ctb+q5)/2d0
px=(px+p1)/2d0
py=(py+p2)/2d0
pz=(pz+p3)/2d0
xia=(xia+p4)/2d0
xib=(xib+p5)/2d0
return
end subroutine

subroutine hamiltonian(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,hami)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,hami
hami=((px**2d0+py**2d0+pz**2d0)/2d0-1d0/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))+((1d0/8d0)*(3d0*n-1d0)*((px**2d0+py**2d0+pz**2d0)**2d0)-(1d0/2d0)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*((3d0+n)*(px**2d0+py**2d0+pz**2d0)+n*(((sx/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*px+(sy/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*py+(sz/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*pz)**2d0))+(1d0/2d0)*(1d0/(sx**2d0+sy**2d0+sz**2d0)))+((1d0/16d0)*(1d0-5d0*n+5d0*n*n)*((px**2d0+py**2d0+pz**2d0)**3d0)+(1d0/8d0)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*&
((5d0-20d0*n-3d0*n*n)*((px**2d0+py**2d0+pz**2d0)**2d0)-2d0*n*n*(((sx/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*px+(sy/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*py+(sz/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*pz)**2d0)*(px**2d0+py**2d0+pz**2d0)-3d0*n*n*(((sx/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*px+(sy/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*py+(sz/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*pz)**4d0))+(1d0/2d0)*(1d0/(sx**2d0+sy**2d0+sz**2d0))*((5d0+8d0*n)*(px**2d0+py**2d0+pz**2d0)+3d0*n*(((sx/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*px+(sy/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*py+(sz/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*pz)**2d0))-(1d0/4d0)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))*(1d0+3d0*n))+(n/((sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))*((sy*pz-sz*py)*(((2d0+(3d0/2d0)*(1d0/w)))*((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dcos(cta))))+((2d0+(3d0/2d0)*w))*((((((c2/w)**2d0-xib**2d0)**&
(1d0/2d0))*dcos(ctb)))))+(sz*px-sx*pz)*(((2d0+(3d0/2d0)*(1d0/w)))*((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dsin(cta))))+((2d0+(3d0/2d0)*w))*((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dsin(ctb)))))+(sx*py-sy*px)*(((2d0+(3d0/2d0)*(1d0/w)))*((xia))+((2d0+(3d0/2d0)*w))*((xib))))+(n/2d0)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))*(3d0*(((((sx)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))))*((((1d0+1d0/w))*((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dcos(cta))))+((1d0+w))*((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dcos(ctb))))))+(((sy)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))))*((((1d0+1d0/w))*((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dsin(cta))))+((1d0+w))*((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dsin(ctb))))))+(((sz)*(1d0/((sx**2d0+sy**2d0+sz**2d0)**&
(1d0/2d0)))))*((((1d0+1d0/w))*((xia))+((1d0+w))*((xib))))))**2d0-((((((1d0+1d0/w))*((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dcos(cta))))+((1d0+w))*((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dcos(ctb))))))**2d0+((((1d0+1d0/w))*((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dsin(cta))))+((1d0+w))*((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dsin(ctb))))))**2d0+((((1d0+1d0/w))*((xia))+((1d0+w))*((xib))))**2d0)))
return
end subroutine

subroutine angularmomentum(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,lt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,lt
real*8 orbitl,totall
!orbitl=((sy*pz-sz*py)**2d0+(sz*px-sx*pz)**2d0+(sx*py-sy*px)**2d0)**(1d0/2d0)
totall=dsqrt((((sy*pz-sz*py))+((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dcos(cta))))+((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dcos(ctb)))))**2d0+(((sz*px-sx*pz))+((((((c1*w)**2d0-xia**2d0)**(1d0/2d0))*dsin(cta))))+((((((c2/w)**2d0-xib**2d0)**(1d0/2d0))*dsin(ctb)))))**2d0+(((sx*py-sy*px))+((xia))+((xib)))**2d0)
!lt=orbitl
lt=totall
return
end subroutine 

subroutine xdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxt
dxt=px-(4d0*n**2d0*px*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0+4d0*px*(3d0*n**2d0+20d0*n-5d0)*(px**2d0+py**2d0+pz**2d0)+(12d0*n**2d0*sx*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**3d0)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(4d0*n**2d0*sx*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))+(2d0*px*(8d0*n+5d0)+(6d0*n*sx*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0))-(2d0*px*(n+3d0)+(2d0*n*sx*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))+6d0*px*((5d0*n**2d0)/16d0-(5d0*n)/16d0+1d0/16d0)*(px**2d0+py**2d0+pz**2d0)**2d0+(n*(sz*(dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))-sy*(xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+4d0*px*((3d0*n)/8d0-1d0/8d0)*(px**2d0+py**2d0+pz**2d0)
return
end subroutine

subroutine ydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dyt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dyt
dyt=py-(4d0*n**2d0*py*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0+4d0*py*(3d0*n**2d0+20d0*n-5d0)*(px**2d0+py**2d0+pz**2d0)+(12d0*n**2d0*sy*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**3d0)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(4d0*n**2d0*sy*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))+(2d0*py*(8d0*n+5d0)+(6d0*n*sy*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0))-(2d0*py*(n+3d0)+(2d0*n*sy*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(n*(sz*(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*&
((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))-sx*(xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+6d0*py*((5d0*n**2d0)/16d0-(5d0*n)/16d0+1d0/16d0)*(px**2d0+py**2d0+pz**2d0)**2d0+4d0*py*((3d0*n)/8d0-1d0/8d0)*(px**2d0+py**2d0+pz**2d0)
return
end subroutine

subroutine zdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dzt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dzt
dzt=pz-(4d0*n**2d0*pz*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0+4d0*pz*(3d0*n**2d0+20d0*n-5d0)*(px**2d0+py**2d0+pz**2d0)+(12d0*n**2d0*sz*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**3d0)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(4d0*n**2d0*sz*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))+(2d0*pz*(8d0*n+5d0)+(6d0*n*sz*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0))-(2d0*pz*(n+3d0)+(2d0*n*sz*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))+6d0*pz*((5d0*n**2d0)/16d0-(5d0*n)/16d0+1d0/16d0)*&
(px**2d0+py**2d0+pz**2d0)**2d0+4d0*pz*((3d0*n)/8d0-1d0/8d0)*(px**2d0+py**2d0+pz**2d0)+(n*(sy*(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))-sx*(dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)
return
end subroutine

subroutine pxdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpxt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpxt
dpxt=sx/(sx**2d0+sy**2d0+sz**2d0)**2d0-sx/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(12d0*n**2d0*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**3d0*((px*sx**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-px/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))+4d0*n**2d0*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((px*sx**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-px/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sx*sy)/&
(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(3d0*sx*(3d0*n+1d0))/(4d0*(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0))+(n*(pz*(dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))-py*(xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(sx*((n+3d0)*(px**2d0+py**2d0+pz**2d0)+n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**&
(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))+(sx*((8d0*n+5d0)*(px**2d0+py**2d0+pz**2d0)+3d0*n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(sx**2d0+sy**2d0+sz**2d0)**2d0-(sx*(3d0*n**2d0*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**4d0+(3d0*n**2d0+20d0*n-5d0)*(px**2d0+py**2d0+pz**2d0)**2d0+2d0*n**2d0*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))-(3d0*n*sx*((dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))**2d0-3d0*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**&
(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0+(xia*(1d0/w+1d0)+xib*(w+1d0))**2d0+(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))**2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0))+(3d0*n*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**&
(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((sx**2d0*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w&
+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(sx*sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(3d0*n*sx*((xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))*(px*sy-py*sx)+(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))*&
(py*sz-pz*sy)-(dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))*(px*sz-pz*sx)))/(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0)+(3d0*n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((px*sx**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-px/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/&
(sx**2d0+sy**2d0+sz**2d0)-(n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((px*sx**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-px/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)
return
end subroutine

subroutine pydot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpyt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpyt
dpyt=sy/(sx**2d0+sy**2d0+sz**2d0)**2d0-sy/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(12d0*n**2d0*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**3d0*((py*sy**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-py/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))+4d0*n**2d0*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((py*sy**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-py/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(n*(pz*(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*&
((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))-px*(xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(3d0*sy*(3d0*n+1d0))/(4d0*(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0))-(sy*((n+3d0)*(px**2d0+py**2d0+pz**2d0)+n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))+(sy*((8d0*n+5d0)*(px**2d0+py**2d0+pz**2d0)+3d0*n*((px*sx)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(sx**2d0+sy**2d0+sz**2d0)**2d0-(sy*(3d0*n**2d0*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**4d0+(3d0*n**2d0+20d0*n-5d0)*(px**2d0+py**2d0+pz**2d0)**2d0+2d0*n**2d0*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))-(3d0*n*sy*((dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))**2d0-3d0*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*&
(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0+(xia*(1d0/w+1d0)+xib*(w+1d0))**2d0+(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))**2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0))-(3d0*n*sy*((xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))*(px*sy-py*sx)+(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))*(py*sz-pz*sy)-(dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*&
(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))*(px*sz-pz*sx)))/(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0)+(3d0*n*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((sy**2d0*(dsin(cta)*(c1**2d0*w**2d0-&
xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*sy*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(sy*sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(3d0*n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((py*sy**2d0)/&
(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-py/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)-(n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((py*sy**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-py/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sy)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(pz*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)
return
end subroutine

subroutine pzdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpzt)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dpzt
dpzt=sz/(sx**2d0+sy**2d0+sz**2d0)**2d0-sz/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(12d0*n**2d0*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**3d0*((pz*sz**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-pz/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(py*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))+4d0*n**2d0*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((pz*sz**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-pz/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sz)/&
(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(py*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(3d0*sz*(3d0*n+1d0))/(4d0*(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0))-(sz*((n+3d0)*(px**2d0+py**2d0+pz**2d0)+n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))+(sz*((8d0*n+5d0)*(px**2d0+py**2d0+pz**2d0)+3d0*n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(sx**2d0+sy**2d0+sz**2d0)**2d0+(n*(py*(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))-px*(dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(sz*(3d0*n**2d0*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**&
(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**4d0+(3d0*n**2d0+20d0*n-5d0)*(px**2d0+py**2d0+pz**2d0)**2d0+2d0*n**2d0*(px**2d0+py**2d0+pz**2d0)*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0))/(8d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))-(3d0*n*sz*((dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))**2d0-3d0*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**&
(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))**2d0+(xia*(1d0/w+1d0)+xib*(w+1d0))**2d0+(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))**2d0))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0))-(3d0*n*sz*((xib*((3d0*w)/2d0+2d0)+xia*(3d0/(2d0*w)+2d0))*(px*sy-py*sx)+(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))*(py*sz-pz*sy)-(dsin(ctb)*&
(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0))*(px*sz-pz*sx)))/(sx**2d0+sy**2d0+sz**2d0)**(5d0/2d0)+(3d0*n*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((sz**2d0*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(xia*(1d0/w+1d0)+xib*(w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*sz*(dcos(cta)*(c1**2d0*w**2d0&
-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(sy*sz*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(3d0*n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((pz*sz**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-pz/(sx**2d0+sy**2d0+sz**2d0)**&
(1d0/2d0)+(px*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(py*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)-(n*((px*sx)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(py*sy)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(pz*sz)/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((pz*sz**2d0)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-pz/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(px*sx*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)+(py*sy*sz)/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)
return
end subroutine

subroutine ctadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dcta)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dcta
dcta=-(n*((3d0/(2d0*w)+2d0)*(px*sy-py*sx)-(xia*dcos(cta)*(3d0/(2d0*w)+2d0)*(py*sz-pz*sy))/(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)+(xia*dsin(cta)*(3d0/(2d0*w)+2d0)*(px*sz-pz*sx))/(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(n*(2d0*(xia*(1d0/w+1d0)+xib*(w+1d0))*(1d0/w+1d0)+6d0*((sx*xia*dcos(cta)*(1d0/w+1d0))/((c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(sz*(1d0/w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*xia*dsin(cta)*(1d0/w+1d0))/((c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(2d0*xia*dcos(cta)*(1d0/w+1d0)*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**&
(1d0/2d0)*(w+1d0)))/(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)-(2d0*xia*dsin(cta)*(1d0/w+1d0)*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))
return
end subroutine

subroutine ctbdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dctb)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dctb
dctb=-(n*(((3d0*w)/2d0+2d0)*(px*sy-py*sx)-(xib*dcos(ctb)*((3d0*w)/2d0+2d0)*(py*sz-pz*sy))/(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)+(xib*dsin(ctb)*((3d0*w)/2d0+2d0)*(px*sz-pz*sx))/(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)-(n*(2d0*(xia*(1d0/w+1d0)+xib*(w+1d0))*(w+1d0)+6d0*((sx*xib*dcos(ctb)*(w+1d0))/((c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(sz*(w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*xib*dsin(ctb)*(w+1d0))/((c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*&
(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)))*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-(2d0*xib*dsin(ctb)*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w&
+1d0))*(w+1d0))/(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)-(2d0*xib*dcos(ctb)*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))*(w+1d0))/(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))
return
end subroutine

subroutine xiadot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxia)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxia
dxia=(n*(6d0*((sx*dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)-(sy*dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**&
(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))-2d0*dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))+2d0*dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*&
(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))-(n*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0)*(px*sz-pz*sx)+dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(3d0/(2d0*w)+2d0)*(py*sz-pz*sy)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)
return
end subroutine

subroutine xibdot(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxib)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,n,c1,c2,w,dxib
dxib=-(n*(6d0*((sy*dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)-(sx*dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))*((sz*(xia*(1d0/w+1d0)+xib*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sx*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*&
(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0)+(sy*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0)))/(sx**2d0+sy**2d0+sz**2d0)**(1d0/2d0))+2d0*dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(dcos(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))*(w+1d0)-2d0*dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(dsin(cta)*(c1**2d0*w**2d0-xia**2d0)**(1d0/2d0)*(1d0/w+1d0)+dsin(ctb)*&
(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*(w+1d0))*(w+1d0)))/(2d0*(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0))-(n*(dcos(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)*(px*sz-pz*sx)+dsin(ctb)*(c2**2d0/w**2d0-xib**2d0)**(1d0/2d0)*((3d0*w)/2d0+2d0)*(py*sz-pz*sy)))/(sx**2d0+sy**2d0+sz**2d0)**(3d0/2d0)
return
end subroutine

subroutine ecf(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,n,c1,c2,w,h)
implicit real*8 (a-z)
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,n,c1,c2,w,h
real*8 q10,q20,q30,q40,q50,p10,p20,p30,p40,p50
real*8 q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m
real*8 z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
real*8 h0000010000,h0000000000,h0000110000,h0001110000,h0011110000,h0111110000,h1111110000,h1111110001,h1111110011,h1111110111,h1111111111
real*8 h0000001000,h0000011000,h0000111000,h0001111000,h0011111000,h0111111000,h1111111000,h1111111001,h1111111011
real*8 h0000000100,h0000001100,h0000011100,h0000111100,h0001111100,h0011111100,h0111111100,h1111111100,h1111111101
real*8 h0000000010,h0000000110,h0000001110,h0000011110,h0000111110,h0001111110,h0011111110,h0111111110,h1111111110
real*8 h0000000001,h0000000011,h0000000111,h0000001111,h0000011111,h0000111111,h0001111111,h0011111111,h0111111111
real*8 h1000000000,h1000000001,h1000000011,h1000000111,h1000001111,h1000011111,h1000111111,h1001111111,h1011111111
real*8 h0100000000,h1100000000,h1100000001,h1100000011,h1100000111,h1100001111,h1100011111,h1100111111,h1101111111
real*8 h0010000000,h0110000000,h1110000000,h1110000001,h1110000011,h1110000111,h1110001111,h1110011111,h1110111111
real*8 h0001000000,h0011000000,h0111000000,h1111000000,h1111000001,h1111000011,h1111000111,h1111001111,h1111011111
real*8 h0000100000,h0001100000,h0011100000,h0111100000,h1111100000,h1111100001,h1111100011,h1111100111,h1111101111
real*8::error=1d-12
real*8::d=0d0
real*8 dist
integer::j=15d0
integer i
i=1d0
q10=q1
q20=q2
q30=q3
q40=q4
q50=q5
p10=p1
p20=p2
p30=p3
p40=p4
p50=p5
q1m=q1
q2m=q2
q3m=q3
q4m=q4
q5m=q5
p1m=p1
p2m=p2
p3m=p3
p4m=p4
p5m=p5
call rks(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h)
1 z1=q1m
  z2=q2m
  z3=q3m
  z4=q4m
  z5=q5m
  z6=p1m
  z7=p2m
  z8=p3m
  z9=p4m
  z10=p5m
call hamiltonian(q10,q20,q30,q40,q50,p1m,p20,p30,p40,p50,n,c1,c2,w,h0000010000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,n,c1,c2,w,h1000010000)
call hamiltonian(q10,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h0000110001)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1100011000)
call hamiltonian(q10,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h0001110011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1110011100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h0011110111)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1111011110)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q1m,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h1000000000)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,n,c1,c2,w,h0000100001)
call hamiltonian(q1m,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h1100001000)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0001100011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h1110001100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0011100111)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h1111001110)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111101111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111101111)
if(dabs(p1m-p10).ne.0d0)then
q1m=h*(1d0/(p1m-p10))*(1d0/10d0)*((h0000010000-h0000000000+h1000010000-h1000000000+h0000110001-h0000100001+h1100011000-h1100001000+h0001110011-h0001100011+h1110011100-h1110001100+h0011110111-h0011100111+h1111011110-h1111001110+h0111111111-h0111101111+h1111111111-h1111101111))+q10
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h0000001000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h0100001000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1000011000)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h0110001100)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1000111001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h0111001110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1001111011)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111101111)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1011111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0100000000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,n,c1,c2,w,h1000010000)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0110000100)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h1000110001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0111000110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h1001110011)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0111100111)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h1011110111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h1111110111)
if(dabs(p2m-p20).ne.0d0)then
q2m=h*(1d0/(p2m-p20))*(1d0/10d0)*((h0000001000-h0000000000+h0100001000-h0100000000+h1000011000-h1000010000+h0110001100-h0110000100+h1000111001-h1000110001+h0111001110-h0111000110+h1001111011-h1001110011+h0111101111-h0111100111+h1011111111-h1011110111+h1111111111-h1111110111))+q20
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0000000100)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0010000100)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h0100001100)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0011000110)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1100011100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0011100111)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1100111101)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h1011110111)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1101111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0010000000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h0100001000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0011000010)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1100011000)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0011100011)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1100111001)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h1011110011)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1101111011)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1111111011)
if(dabs(p3m-p30).ne.0d0)then
q3m=h*(1d0/(p3m-p30))*(1d0/10d0)*((h0000000100-h0000000000+h0010000100-h0010000000+h0100001100-h0100001000+h0011000110-h0011000010+h1100011100-h1100011000+h0011100111-h0011100011+h1100111101-h1100111001+h1011110111-h1011110011+h1101111111-h1101111011+h1111111111-h1111111011))+q30
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0000000010)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0001000010)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0010000110)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0001100011)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h0110001110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h1001110011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1110011110)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1101111011)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1110111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0001000000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0010000100)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p40,p5m,n,c1,c2,w,h0001100001)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h0110001100)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h1001110001)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1110011100)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1101111001)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1110111101)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1111111101)
if(dabs(p4m-p40).ne.0d0)then
q4m=h*(1d0/(p4m-p40))*(1d0/10d0)*((h0000000010-h0000000000+h0001000010-h0001000000+h0010000110-h0010000100+h0001100011-h0001100001+h0110001110-h0110001100+h1001110011-h1001110001+h1110011110-h1110011100+h1101111011-h1101111001+h1110111111-h1110111101+h1111111111-h1111111101))+q40
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p5m,n,c1,c2,w,h0000000001)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,n,c1,c2,w,h0000100001)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0001000011)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h1000110001)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0011000111)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1100111001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111001111)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1110111101)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111011111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p50,n,c1,c2,w,h0000100000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0001000010)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p50,n,c1,c2,w,h1000110000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0011000110)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1100111000)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h0111001110)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1110111100)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1111011110)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1111111110)
if(dabs(p5m-p50).ne.0d0)then
q5m=h*(1d0/(p5m-p50))*(1d0/10d0)*((h0000000001-h0000000000+h0000100001-h0000100000+h0001000011-h0001000010+h1000110001-h1000110000+h0011000111-h0011000110+h1100111001-h1100111000+h0111001111-h0111001110+h1110111101-h1110111100+h1111011111-h1111011110+h1111111111-h1111111110))+q50
else
end if
call hamiltonian(q1m,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h1000000000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,n,c1,c2,w,h1000010000)
call hamiltonian(q1m,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h1100001000)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h1000110001)
call hamiltonian(q1m,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h1110001100)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h1001110011)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h1111001110)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h1011110111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111101111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p1m,p20,p30,p40,p50,n,c1,c2,w,h0000010000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h0100001000)
call hamiltonian(q10,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h0000110001)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h0110001100)
call hamiltonian(q10,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h0001110011)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h0111001110)
call hamiltonian(q10,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h0011110111)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111101111)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111111111)
if(dabs(q1m-q10).ne.0d0)then
p1m=(-h)*(1d0/(q1m-q10))*(1d0/10d0)*((h1000000000-h0000000000+h1000010000-h0000010000+h1100001000-h0100001000+h1000110001-h0000110001+h1110001100-h0110001100+h1001110011-h0001110011+h1111001110-h0111001110+h1011110111-h0011110111+h1111101111-h0111101111+h1111111111-h0111111111))+p10
else
end if
call hamiltonian(q10,q2m,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0100000000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h0100001000)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0110000100)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1100011000)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0111000110)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1100111001)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0111100111)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1101111011)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h1111110111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p2m,p30,p40,p50,n,c1,c2,w,h0000001000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0010000100)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1000011000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0011000110)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1000111001)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0011100111)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1001111011)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,n,c1,c2,w,h1011110111)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1011111111)
if(dabs(q2m-q20).ne.0d0)then
p2m=(-h)*(1d0/(q2m-q20))*(1d0/10d0)*((h0100000000-h0000000000+h0100001000-h0000001000+h0110000100-h0010000100+h1100011000-h1000011000+h0111000110-h0011000110+h1100111001-h1000111001+h0111100111-h0011100111+h1101111011-h1001111011+h1111110111-h1011110111+h1111111111-h1011111111))+p20
else
end if
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0010000000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0010000100)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0011000010)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h0110001100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0011100011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1110011100)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h1011110011)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1110111101)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1111111011)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p3m,p40,p50,n,c1,c2,w,h0000000100)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0001000010)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p3m,p40,p50,n,c1,c2,w,h0100001100)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0001100011)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1100011100)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,n,c1,c2,w,h1001110011)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1100111101)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,n,c1,c2,w,h1101111011)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1101111111)
if(dabs(q3m-q30).ne.0d0)then
p3m=(-h)*(1d0/(q3m-q30))*(1d0/10d0)*((h0010000000-h0000000000+h0010000100-h0000000100+h0011000010-h0001000010+h0110001100-h0100001100+h0011100011-h0001100011+h1110011100-h1100011100+h1011110011-h1001110011+h1110111101-h1100111101+h1111111011-h1101111011+h1111111111-h1101111111))+p30
else
end if
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0001000000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0001000010)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p40,p5m,n,c1,c2,w,h0001100001)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0011000110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h1001110001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h0111001110)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1101111001)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1111011110)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1111111101)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p4m,p50,n,c1,c2,w,h0000000010)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,n,c1,c2,w,h0000100001)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p4m,p50,n,c1,c2,w,h0010000110)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,n,c1,c2,w,h1000110001)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p4m,p50,n,c1,c2,w,h0110001110)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,n,c1,c2,w,h1100111001)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1110011110)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,n,c1,c2,w,h1110111101)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1110111111)
if(dabs(q4m-q40).ne.0d0)then
p4m=(-h)*(1d0/(q4m-q40))*(1d0/10d0)*((h0001000000-h0000000000+h0001000010-h0000000010+h0001100001-h0000100001+h0011000110-h0010000110+h1001110001-h1000110001+h0111001110-h0110001110+h1101111001-h1100111001+h1111011110-h1110011110+h1111111101-h1110111101+h1111111111-h1110111111))+p40
else
end if
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p50,n,c1,c2,w,h0000100000)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,n,c1,c2,w,h0000100001)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p50,n,c1,c2,w,h1000110000)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0001100011)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1100111000)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0011100111)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1110111100)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111101111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1111111110)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,n,c1,c2,w,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p5m,n,c1,c2,w,h0000000001)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,n,c1,c2,w,h1000010000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p5m,n,c1,c2,w,h0001000011)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,n,c1,c2,w,h1100011000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p5m,n,c1,c2,w,h0011000111)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,n,c1,c2,w,h1110011100)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p5m,n,c1,c2,w,h0111001111)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,n,c1,c2,w,h1111011110)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p5m,n,c1,c2,w,h1111011111)
if(dabs(q5m-q50).ne.0d0)then
p5m=(-h)*(1d0/(q5m-q50))*(1d0/10d0)*((h0000100000-h0000000000+h0000100001-h0000000001+h1000110000-h1000010000+h0001100011-h0001000011+h1100111000-h1100011000+h0011100111-h0011000111+h1110111100-h1110011100+h0111101111-h0111001111+h1111111110-h1111011110+h1111111111-h1111011111))+p50
else
end if  
  i=i+1d0
  d=dsqrt((z1-q1m)**2d0+(z2-q2m)**2d0+(z3-q3m)**2d0+(z4-q4m)**2d0+(z5-q5m)**2d0+(z6-p1m)**2d0+(z7-p2m)**2d0+(z8-p3m)**2d0+(z9-p4m)**2d0+(z10-p5m)**2d0+1d-32)
  !write(*,*)d
  !pause
  if(i.gt.j)then
  goto 2
  else
  end if
  if(d.gt.error)goto 1
2 q1=q1m
  q2=q2m
  q3=q3m
  q4=q4m
  q5=q5m
  p1=p1m
  p2=p2m
  p3=p3m
  p4=p4m
  p5=p5m
return
end subroutine