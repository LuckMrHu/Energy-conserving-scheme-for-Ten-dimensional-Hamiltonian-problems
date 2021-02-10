program main
implicit none
character(len=20)::file1="data1.dat"
character(len=20)::file2="data2.dat"
real*8 beta,epsi
common beta,epsi
real*8 ht,t1,tout1,t2,tout2,t3,tout3
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5
real*8 x1,x2,x3,x4,x5,v1,v2,v3,v4,v5
real*8 h,endt,reorbit
real*8::t=0d0,init=0d0,overt=0d0
real*8 inih,inis,temph,temps,fh,fs,r1,r2,r
real*8 inid,tempd,fli,el,num
integer i,j,k
integer::count=0d0
reorbit=0d0
open(10,file=file1)
open(20,file=file2)
call cpu_time(init)
call constant(beta,epsi,ht,endt,inid)
call initiala(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,inih,inis,fli,el,num,tout1,tout2,t1,t2)
call initialb(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,x1,x2,x3,x4,x5,v1,v2,v3,v4,v5,inid)
write(*,*)inih
write(*,*)inis
pause     
do 3 j=1d0,51d0,1d0     
do 2 i=1d0,101d0,1d0      
1       tout1=tout1+ht
       !call rks(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,ht)
	   !call iss(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,ht)
	   !call ems(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,ht)
	   call ecf(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,ht)
	   call ecf(x1,x2,x3,x4,x5,v1,v2,v3,v4,v5,beta,epsi,ht) 
	   tempd=dsqrt((q1-x1)**2d0+(q2-x2)**2d0+(q3-x3)**2d0+(p1-v1)**2d0+(p2-v2)**2d0+(p3-v3)**2d0+(q4-x4)**2d0+(q5-x5)**2d0+(p4-v4)**2d0+(p5-v5)**2d0)
	   if(tempd.ge.1d0)then
       x1=q1+inid*(x1-q1)/tempd
       x2=q2+inid*(x2-q2)/tempd
       x3=q3+inid*(x3-q3)/tempd
       x4=q4+inid*(x4-q4)/tempd
       x5=q5+inid*(x5-q5)/tempd
       v1=p1+inid*(v1-p1)/tempd
       v2=p2+inid*(v2-p2)/tempd
       v3=p3+inid*(v3-p3)/tempd
       v4=p4+inid*(v4-p4)/tempd
       v5=p5+inid*(v5-p5)/tempd
	   num=num+1d0
       tempd=inid
	   else
	   end if
	   fli=num*dlog10(1d0/inid)+dlog10(tempd/inid)
if(tout1.lt.3800d0)goto 1
call hamiltonian(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,temph)
call norm(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,temps)
fh=dlog10(dabs(temph-inih)+1d-16)
fs=dlog10(dabs(temps-inis)+1d-16)
write(10,*)epsi,beta,fli
write(20,*)fh,fs
count=count+1d0
write(*,*)count
beta=beta+0.1499d0
!epsi=epsi+0.2d0
call initiala(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,inih,inis,fli,el,num,tout1,tout2,t1,t2)
call initialb(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,x1,x2,x3,x4,x5,v1,v2,v3,v4,v5,inid)
2 continue
beta=0.01d0
!epsi=-10d0
!beta=beta+0.1798d0
!epsi=epsi+0.2d0
epsi=epsi+0.4d0
call initiala(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,inih,inis,fli,el,num,tout1,tout2,t1,t2)
call initialb(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,x1,x2,x3,x4,x5,v1,v2,v3,v4,v5,inid)
3 continue
call cpu_time(overt)
write(*,*)dabs(overt-init)
write(*,*)count
pause
close(10)
close(20)
stop
end

subroutine constant(beta,epsi,h,endt,inid)
implicit none
real*8 beta,epsi,h,endt,inid
beta=0.01d0
epsi=-10d0
h=0.01d0
endt=3800d0
inid=1d-8
return
end subroutine

subroutine initiala(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,inih,inis,fli,el,num,tout1,tout2,t1,t2)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,inih,inis,fli,el,num,tout1,tout2,t1,t2
el=0d0
fli=0d0
num=0d0
tout1=0d0
tout2=0d0
t1=0d0
t2=0d0
q1=0.5d0
q2=0.4d0
q3=0.3d0
q4=0.2d0
q5=0.1d0
p1=0.05d0
p2=0.04d0
p3=0.03d0
p4=0.02d0
p5=0.01d0
call hamiltonian(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,inih)
call norm(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,inis)
return
end subroutine

subroutine initialb(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,x1,x2,x3,x4,x5,v1,v2,v3,v4,v5,inid)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,x1,x2,x3,x4,x5,v1,v2,v3,v4,v5,inid
x1=q1+inid
x2=q2
x3=q3
x4=q4
x5=q5
v1=p1
v2=p2
v3=p3
v4=p4
v5=p5
return
end subroutine

subroutine timestep(t,h)
implicit none
real*8 t,h
t=t+h
return
end subroutine

subroutine norm(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,norms)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,norms
norms=(1d0/2d0)*(q1**2d0+q2**2d0+q3**2d0+q4**2d0+q5**2d0+p1**2d0+p2**2d0+p3**2d0+p4**2d0+p5**2d0)
return
end subroutine

subroutine hamiltonian(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,hami)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,hami
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
hami=(epsi1/2d0)*(q1**2d0+p1**2d0)+(beta/8d0)*((q1**2d0+p1**2d0)**2d0)-(p1*p2+q1*q2)&
&+(epsi2/2d0)*(q2**2d0+p2**2d0)+(beta/8d0)*((q2**2d0+p2**2d0)**2d0)-(p2*p3+q2*q3)&
&+(epsi3/2d0)*(q3**2d0+p3**2d0)+(beta/8d0)*((q3**2d0+p3**2d0)**2d0)-(p3*p4+q3*q4)&
&+(epsi4/2d0)*(q4**2d0+p4**2d0)+(beta/8d0)*((q4**2d0+p4**2d0)**2d0)-(p4*p5+q4*q5)&
&+(epsi5/2d0)*(q5**2d0+p5**2d0)+(beta/8d0)*((q5**2d0+p5**2d0)**2d0)
return
end subroutine

subroutine qadot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqa)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqa
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dqa=epsi1*p1-p2+(beta*p1*(p1**2d0+q1**2d0))/2d0
return
end subroutine

subroutine qbdot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqb)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqb
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dqb=epsi2*p2-p3-p1+(beta*p2*(p2**2d0+q2**2d0))/2d0
return
end subroutine

subroutine qcdot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqc)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqc
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dqc=epsi3*p3-p4-p2+(beta*p3*(p3**2d0+q3**2d0))/2d0
return
end subroutine

subroutine qddot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqd)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqd
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dqd=epsi4*p4-p5-p3+(beta*p4*(p4**2d0+q4**2d0))/2d0
return
end subroutine

subroutine qedot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqe)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dqe
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dqe=epsi5*p5-p4+(beta*p5*(p5**2d0+q5**2d0))/2d0
return
end subroutine

subroutine padot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpa)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpa
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dpa=q2-epsi1*q1-(beta*q1*(p1**2d0+q1**2d0))/2d0
return
end subroutine

subroutine pbdot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpb)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpb
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dpb=q1+q3-epsi2*q2-(beta*q2*(p2**2d0+q2**2d0))/2d0
return
end subroutine

subroutine pcdot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpc)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpc
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dpc=q2+q4-epsi3*q3-(beta*q3*(p3**2d0+q3**2d0))/2d0
return
end subroutine

subroutine pddot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpd)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpd
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dpd=q3+q5-epsi4*q4-(beta*q4*(p4**2d0+q4**2d0))/2d0
return
end subroutine

subroutine pedot(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpe)
implicit none
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,epsi,beta,dpe
real*8 epsi1,epsi2,epsi3,epsi4,epsi5
epsi1=epsi
epsi2=epsi
epsi3=epsi
epsi4=epsi
epsi5=epsi
dpe=q4-epsi5*q5-(beta*q5*(p5**2d0+q5**2d0))/2d0
return
end subroutine

subroutine iss(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,h)    
implicit none
real*8,dimension(5)::p,q,e
real*8 b,h,z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi
real*8,dimension(5)::pn,qn,pm,qm
real*8::r=0d0,error=1d-15
p(1)=p1
p(2)=p2
p(3)=p3
p(4)=p4
p(5)=p5
q(1)=q1
q(2)=q2
q(3)=q3
q(4)=q4
q(5)=q5
b=beta
e(1)=epsi
e(2)=epsi
e(3)=epsi
e(4)=epsi
e(5)=epsi
pn=p
pm=p
qn=q
qm=q
100 z1=qm(1)
    z2=qm(2)
    z3=qm(3)
    z4=qm(4)
    z5=qm(5)
    z6=pm(1)
    z7=pm(2)
    z8=pm(3)
    z9=pm(4)
    z10=pm(5)
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    qm(1)=qn(1)+h*(p(1)*(e(1)+b*((q(1)**2d0+p(1)**2d0)/2d0))-p(2))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    qm(2)=qn(2)+h*(p(2)*(e(2)+b*((q(2)**2d0+p(2)**2d0)/2d0))-p(3)-p(1))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    qm(3)=qn(3)+h*(p(3)*(e(3)+b*((q(3)**2d0+p(3)**2d0)/2d0))-p(4)-p(2))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    qm(4)=qn(4)+h*(p(4)*(e(4)+b*((q(4)**2d0+p(4)**2d0)/2d0))-p(5)-p(3))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    qm(5)=qn(5)+h*(p(5)*(e(5)+b*((q(5)**2d0+p(5)**2d0)/2d0))-p(4))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    pm(1)=pn(1)+h*(-q(1)*(e(1)+b*((q(1)**2d0+p(1)**2d0)/2d0))+q(2))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    pm(2)=pn(2)+h*(-q(2)*(e(2)+b*((q(2)**2d0+p(2)**2d0)/2d0))+q(3)+q(1))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    pm(3)=pn(3)+h*(-q(3)*(e(3)+b*((q(3)**2d0+p(3)**2d0)/2d0))+q(4)+q(2))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    pm(4)=pn(4)+h*(-q(4)*(e(4)+b*((q(4)**2d0+p(4)**2d0)/2d0))+q(5)+q(3))
    q(1)=(qm(1)+qn(1))/2d0
    q(2)=(qm(2)+qn(2))/2d0
    q(3)=(qm(3)+qn(3))/2d0
    q(4)=(qm(4)+qn(4))/2d0
    q(5)=(qm(5)+qn(5))/2d0
    p(1)=(pm(1)+pn(1))/2d0
    p(2)=(pm(2)+pn(2))/2d0
    p(3)=(pm(3)+pn(3))/2d0
    p(4)=(pm(4)+pn(4))/2d0
    p(5)=(pm(5)+pn(5))/2d0
    pm(5)=pn(5)+h*(-q(5)*(e(5)+b*((q(5)**2d0+p(5)**2d0)/2d0))+q(4))
    r=dsqrt((qm(1)-z1)**2d0+(qm(2)-z2)**2d0+(qm(3)-z3)**2d0+(qm(4)-z4)**2d0+(qm(5)-z5)**2d0+(pm(1)-z6)**2d0+(pm(2)-z7)**2d0+(pm(3)-z8)**2d0+(pm(4)-z9)**2d0+(pm(5)-z10)**2d0)
    q=qm
    p=pm
    if(r.gt.error)goto 100
	q1=q(1)
	q2=q(2)
	q3=q(3)
	q4=q(4)
	q5=q(5)
	p1=p(1)
	p2=p(2)
	p3=p(3)
	p4=p(4)
	p5=p(5)
return
end subroutine

subroutine rks(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,h)
implicit none
real*8,dimension(5)::p,q,e
real*8,dimension(2,10)::dk                        
real*8,dimension(5)::inip
real*8,dimension(5)::iniq
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi
real*8 b,h
p(1)=p1
p(2)=p2
p(3)=p3
p(4)=p4
p(5)=p5
q(1)=q1
q(2)=q2
q(3)=q3
q(4)=q4
q(5)=q5
b=beta
e(1)=epsi
e(2)=epsi
e(3)=epsi
e(4)=epsi
e(5)=epsi
iniq=q
inip=p
dk=0d0
dk(1,1)=h*(p(1)*(e(1)+b*(1d0/2d0)*(q(1)**2d0+p(1)**2d0))-p(2))
dk(1,2)=h*(p(2)*(e(2)+b*(1d0/2d0)*(q(2)**2d0+p(2)**2d0))-p(3)-p(1))
dk(1,3)=h*(p(3)*(e(3)+b*(1d0/2d0)*(q(3)**2d0+p(3)**2d0))-p(4)-p(2))
dk(1,4)=h*(p(4)*(e(4)+b*(1d0/2d0)*(q(4)**2d0+p(4)**2d0))-p(5)-p(3))
dk(1,5)=h*(p(5)*(e(5)+b*(1d0/2d0)*(q(5)**2d0+p(5)**2d0))-p(4))
dk(1,6)=h*(((-q(1))*(e(1)+b*(1d0/2d0)*(q(1)**2d0+p(1)**2d0)))+q(2))
dk(1,7)=h*(((-q(2))*(e(2)+b*(1d0/2d0)*(q(2)**2d0+p(2)**2d0)))+q(3)+q(1))
dk(1,8)=h*(((-q(3))*(e(3)+b*(1d0/2d0)*(q(3)**2d0+p(3)**2d0)))+q(4)+q(2))
dk(1,9)=h*(((-q(4))*(e(4)+b*(1d0/2d0)*(q(4)**2d0+p(4)**2d0)))+q(5)+q(3))
dk(1,10)=h*(((-q(5))*(e(5)+b*(1d0/2d0)*(q(5)**2d0+p(5)**2d0)))+q(4))
q(1)=iniq(1)+(1d0/2d0)*dk(1,1)
q(2)=iniq(2)+(1d0/2d0)*dk(1,2)
q(3)=iniq(3)+(1d0/2d0)*dk(1,3)
q(4)=iniq(4)+(1d0/2d0)*dk(1,4)
q(5)=iniq(5)+(1d0/2d0)*dk(1,5)
p(1)=inip(1)+(1d0/2d0)*dk(1,6)
p(2)=inip(2)+(1d0/2d0)*dk(1,7)
p(3)=inip(3)+(1d0/2d0)*dk(1,8)
p(4)=inip(4)+(1d0/2d0)*dk(1,9)
p(5)=inip(5)+(1d0/2d0)*dk(1,10)
dk(2,1)=h*(p(1)*(e(1)+b*(1d0/2d0)*(q(1)**2d0+p(1)**2d0))-p(2))
dk(2,2)=h*(p(2)*(e(2)+b*(1d0/2d0)*(q(2)**2d0+p(2)**2d0))-p(3)-p(1))
dk(2,3)=h*(p(3)*(e(3)+b*(1d0/2d0)*(q(3)**2d0+p(3)**2d0))-p(4)-p(2))
dk(2,4)=h*(p(4)*(e(4)+b*(1d0/2d0)*(q(4)**2d0+p(4)**2d0))-p(5)-p(3))
dk(2,5)=h*(p(5)*(e(5)+b*(1d0/2d0)*(q(5)**2d0+p(5)**2d0))-p(4))
dk(2,6)=h*(((-q(1))*(e(1)+b*(1d0/2d0)*(q(1)**2d0+p(1)**2d0)))+q(2))
dk(2,7)=h*(((-q(2))*(e(2)+b*(1d0/2d0)*(q(2)**2d0+p(2)**2d0)))+q(3)+q(1))
dk(2,8)=h*(((-q(3))*(e(3)+b*(1d0/2d0)*(q(3)**2d0+p(3)**2d0)))+q(4)+q(2))
dk(2,9)=h*(((-q(4))*(e(4)+b*(1d0/2d0)*(q(4)**2d0+p(4)**2d0)))+q(5)+q(3))
dk(2,10)=h*(((-q(5))*(e(5)+b*(1d0/2d0)*(q(5)**2d0+p(5)**2d0)))+q(4))
q(1)=iniq(1)+dk(2,1)
q(2)=iniq(2)+dk(2,2)
q(3)=iniq(3)+dk(2,3)
q(4)=iniq(4)+dk(2,4)
q(5)=iniq(5)+dk(2,5)
p(1)=inip(1)+dk(2,6)
p(2)=inip(2)+dk(2,7)
p(3)=inip(3)+dk(2,8)
p(4)=inip(4)+dk(2,9)
p(5)=inip(5)+dk(2,10)
q1=q(1)
q2=q(2)
q3=q(3)
q4=q(4)
q5=q(5)
p1=p(1)
p2=p(2)
p3=p(3)
p4=p(4)
p5=p(5)
return
end subroutine

subroutine ecf(q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,h)
implicit real*8 (a-z)
real*8 q1,q2,q3,q4,q5,p1,p2,p3,p4,p5,beta,epsi,h
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
real*8::error=1d-15
real*8::d=0d0
real*8 dist
integer::j=30d0
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
call rks(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,beta,epsi,h)
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
call hamiltonian(q10,q20,q30,q40,q50,p1m,p20,p30,p40,p50,epsi,beta,h0000010000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,epsi,beta,h1000010000)
call hamiltonian(q10,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h0000110001)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,epsi,beta,h1100011000)
call hamiltonian(q10,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h0001110011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,epsi,beta,h1110011100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h0011110111)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,epsi,beta,h1111011110)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h0111111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q1m,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h1000000000)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,epsi,beta,h0000100001)
call hamiltonian(q1m,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h1100001000)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,epsi,beta,h0001100011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h1110001100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,epsi,beta,h0011100111)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h1111001110)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,epsi,beta,h0111101111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,epsi,beta,h1111101111)
if(dabs(p1m-p10).ne.0d0)then
q1m=h*(1d0/(p1m-p10))*(1d0/10d0)*((h0000010000-h0000000000+h1000010000-h1000000000+h0000110001-h0000100001+h1100011000-h1100001000+h0001110011-h0001100011+h1110011100-h1110001100+h0011110111-h0011100111+h1111011110-h1111001110+h0111111111-h0111101111+h1111111111-h1111101111))+q10
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h0000001000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h0100001000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p2m,p30,p40,p50,epsi,beta,h1000011000)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h0110001100)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1000111001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h0111001110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1001111011)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,epsi,beta,h0111101111)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1011111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0100000000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,epsi,beta,h1000010000)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0110000100)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h1000110001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0111000110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h1001110011)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,epsi,beta,h0111100111)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h1011110111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h1111110111)
if(dabs(p2m-p20).ne.0d0)then
q2m=h*(1d0/(p2m-p20))*(1d0/10d0)*((h0000001000-h0000000000+h0100001000-h0100000000+h1000011000-h1000010000+h0110001100-h0110000100+h1000111001-h1000110001+h0111001110-h0111000110+h1001111011-h1001110011+h0111101111-h0111100111+h1011111111-h1011110111+h1111111111-h1111110111))+q20
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0000000100)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0010000100)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h0100001100)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0011000110)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p3m,p40,p50,epsi,beta,h1100011100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,epsi,beta,h0011100111)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1100111101)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h1011110111)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1101111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0010000000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h0100001000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p30,p4m,p50,epsi,beta,h0011000010)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,epsi,beta,h1100011000)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p30,p4m,p5m,epsi,beta,h0011100011)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1100111001)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h1011110011)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1101111011)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1111111011)
if(dabs(p3m-p30).ne.0d0)then
q3m=h*(1d0/(p3m-p30))*(1d0/10d0)*((h0000000100-h0000000000+h0010000100-h0010000000+h0100001100-h0100001000+h0011000110-h0011000010+h1100011100-h1100011000+h0011100111-h0011100011+h1100111101-h1100111001+h1011110111-h1011110011+h1101111111-h1101111011+h1111111111-h1111111011))+q30
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p4m,p50,epsi,beta,h0000000010)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,epsi,beta,h0001000010)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0010000110)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,epsi,beta,h0001100011)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h0110001110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h1001110011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p4m,p50,epsi,beta,h1110011110)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1101111011)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1110111111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p40,p50,epsi,beta,h0001000000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0010000100)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p40,p5m,epsi,beta,h0001100001)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h0110001100)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h1001110001)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,epsi,beta,h1110011100)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1101111001)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1110111101)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1111111101)
if(dabs(p4m-p40).ne.0d0)then
q4m=h*(1d0/(p4m-p40))*(1d0/10d0)*((h0000000010-h0000000000+h0001000010-h0001000000+h0010000110-h0010000100+h0001100011-h0001100001+h0110001110-h0110001100+h1001110011-h1001110001+h1110011110-h1110011100+h1101111011-h1101111001+h1110111111-h1110111101+h1111111111-h1111111101))+q40
else
end if
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p5m,epsi,beta,h0000000001)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,epsi,beta,h0000100001)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p5m,epsi,beta,h0001000011)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h1000110001)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p5m,epsi,beta,h0011000111)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1100111001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p5m,epsi,beta,h0111001111)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1110111101)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111011111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p50,epsi,beta,h0000100000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,epsi,beta,h0001000010)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p50,epsi,beta,h1000110000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0011000110)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p50,epsi,beta,h1100111000)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h0111001110)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p50,epsi,beta,h1110111100)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,epsi,beta,h1111011110)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p50,epsi,beta,h1111111110)
if(dabs(p5m-p50).ne.0d0)then
q5m=h*(1d0/(p5m-p50))*(1d0/10d0)*((h0000000001-h0000000000+h0000100001-h0000100000+h0001000011-h0001000010+h1000110001-h1000110000+h0011000111-h0011000110+h1100111001-h1100111000+h0111001111-h0111001110+h1110111101-h1110111100+h1111011111-h1111011110+h1111111111-h1111111110))+q50
else
end if
call hamiltonian(q1m,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h1000000000)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,epsi,beta,h1000010000)
call hamiltonian(q1m,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h1100001000)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h1000110001)
call hamiltonian(q1m,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h1110001100)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h1001110011)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h1111001110)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h1011110111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,epsi,beta,h1111101111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p1m,p20,p30,p40,p50,epsi,beta,h0000010000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h0100001000)
call hamiltonian(q10,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h0000110001)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h0110001100)
call hamiltonian(q10,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h0001110011)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h0111001110)
call hamiltonian(q10,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h0011110111)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,epsi,beta,h0111101111)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h0111111111)
if(dabs(q1m-q10).ne.0d0)then
p1m=(-h)*(1d0/(q1m-q10))*(1d0/10d0)*((h1000000000-h0000000000+h1000010000-h0000010000+h1100001000-h0100001000+h1000110001-h0000110001+h1110001100-h0110001100+h1001110011-h0001110011+h1111001110-h0111001110+h1011110111-h0011110111+h1111101111-h0111101111+h1111111111-h0111111111))+p10
else
end if
call hamiltonian(q10,q2m,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0100000000)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h0100001000)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0110000100)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,epsi,beta,h1100011000)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0111000110)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1100111001)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,epsi,beta,h0111100111)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1101111011)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h1111110111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p2m,p30,p40,p50,epsi,beta,h0000001000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0010000100)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p2m,p30,p40,p50,epsi,beta,h1000011000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0011000110)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1000111001)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,epsi,beta,h0011100111)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1001111011)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p3m,p4m,p5m,epsi,beta,h1011110111)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1011111111)
if(dabs(q2m-q20).ne.0d0)then
p2m=(-h)*(1d0/(q2m-q20))*(1d0/10d0)*((h0100000000-h0000000000+h0100001000-h0000001000+h0110000100-h0010000100+h1100011000-h1000011000+h0111000110-h0011000110+h1100111001-h1000111001+h0111100111-h0011100111+h1101111011-h1001111011+h1111110111-h1011110111+h1111111111-h1011111111))+p20
else
end if
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0010000000)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0010000100)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p30,p4m,p50,epsi,beta,h0011000010)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h0110001100)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p30,p4m,p5m,epsi,beta,h0011100011)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,epsi,beta,h1110011100)
call hamiltonian(q1m,q20,q3m,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h1011110011)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1110111101)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1111111011)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p3m,p40,p50,epsi,beta,h0000000100)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,epsi,beta,h0001000010)
call hamiltonian(q10,q2m,q30,q40,q50,p10,p2m,p3m,p40,p50,epsi,beta,h0100001100)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,epsi,beta,h0001100011)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p3m,p40,p50,epsi,beta,h1100011100)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p4m,p5m,epsi,beta,h1001110011)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1100111101)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p4m,p5m,epsi,beta,h1101111011)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1101111111)
if(dabs(q3m-q30).ne.0d0)then
p3m=(-h)*(1d0/(q3m-q30))*(1d0/10d0)*((h0010000000-h0000000000+h0010000100-h0000000100+h0011000010-h0001000010+h0110001100-h0100001100+h0011100011-h0001100011+h1110011100-h1100011100+h1011110011-h1001110011+h1110111101-h1100111101+h1111111011-h1101111011+h1111111111-h1101111111))+p30
else
end if
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p40,p50,epsi,beta,h0001000000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p50,epsi,beta,h0001000010)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p40,p5m,epsi,beta,h0001100001)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0011000110)
call hamiltonian(q1m,q20,q30,q4m,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h1001110001)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h0111001110)
call hamiltonian(q1m,q2m,q30,q4m,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1101111001)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,epsi,beta,h1111011110)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1111111101)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p4m,p50,epsi,beta,h0000000010)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,epsi,beta,h0000100001)
call hamiltonian(q10,q20,q3m,q40,q50,p10,p20,p3m,p4m,p50,epsi,beta,h0010000110)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p5m,epsi,beta,h1000110001)
call hamiltonian(q10,q2m,q3m,q40,q50,p10,p2m,p3m,p4m,p50,epsi,beta,h0110001110)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p5m,epsi,beta,h1100111001)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p4m,p50,epsi,beta,h1110011110)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p5m,epsi,beta,h1110111101)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1110111111)
if(dabs(q4m-q40).ne.0d0)then
p4m=(-h)*(1d0/(q4m-q40))*(1d0/10d0)*((h0001000000-h0000000000+h0001000010-h0000000010+h0001100001-h0000100001+h0011000110-h0010000110+h1001110001-h1000110001+h0111001110-h0110001110+h1101111001-h1100111001+h1111011110-h1110011110+h1111111101-h1110111101+h1111111111-h1110111111))+p40
else
end if
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p50,epsi,beta,h0000100000)
call hamiltonian(q10,q20,q30,q40,q5m,p10,p20,p30,p40,p5m,epsi,beta,h0000100001)
call hamiltonian(q1m,q20,q30,q40,q5m,p1m,p20,p30,p40,p50,epsi,beta,h1000110000)
call hamiltonian(q10,q20,q30,q4m,q5m,p10,p20,p30,p4m,p5m,epsi,beta,h0001100011)
call hamiltonian(q1m,q2m,q30,q40,q5m,p1m,p2m,p30,p40,p50,epsi,beta,h1100111000)
call hamiltonian(q10,q20,q3m,q4m,q5m,p10,p20,p3m,p4m,p5m,epsi,beta,h0011100111)
call hamiltonian(q1m,q2m,q3m,q40,q5m,p1m,p2m,p3m,p40,p50,epsi,beta,h1110111100)
call hamiltonian(q10,q2m,q3m,q4m,q5m,p10,p2m,p3m,p4m,p5m,epsi,beta,h0111101111)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p50,epsi,beta,h1111111110)
call hamiltonian(q1m,q2m,q3m,q4m,q5m,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111111111)
!call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p50,epsi,beta,h0000000000)
call hamiltonian(q10,q20,q30,q40,q50,p10,p20,p30,p40,p5m,epsi,beta,h0000000001)
call hamiltonian(q1m,q20,q30,q40,q50,p1m,p20,p30,p40,p50,epsi,beta,h1000010000)
call hamiltonian(q10,q20,q30,q4m,q50,p10,p20,p30,p4m,p5m,epsi,beta,h0001000011)
call hamiltonian(q1m,q2m,q30,q40,q50,p1m,p2m,p30,p40,p50,epsi,beta,h1100011000)
call hamiltonian(q10,q20,q3m,q4m,q50,p10,p20,p3m,p4m,p5m,epsi,beta,h0011000111)
call hamiltonian(q1m,q2m,q3m,q40,q50,p1m,p2m,p3m,p40,p50,epsi,beta,h1110011100)
call hamiltonian(q10,q2m,q3m,q4m,q50,p10,p2m,p3m,p4m,p5m,epsi,beta,h0111001111)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p50,epsi,beta,h1111011110)
call hamiltonian(q1m,q2m,q3m,q4m,q50,p1m,p2m,p3m,p4m,p5m,epsi,beta,h1111011111)
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

subroutine ems(sx,sy,sz,cta,ctb,px,py,pz,xia,xib,beta,epsi,h)
implicit none
real*8 sx,sy,sz,cta,ctb,px,py,pz,xia,xib,beta,epsi,h
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
call qadot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dx)
call qbdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dy)
call qcdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dz)
call qddot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dcta)
call qedot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dctb)
call padot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp1)
call pbdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp2)
call pcdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp3)
call pddot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp4)
call pedot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp5)
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
call qadot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dq1)
call qbdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dq2)
call qcdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dq3)
call qddot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dq4)
call qedot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dq5)
call padot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dpx)
call pbdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dpy)
call pcdot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dpz)
call pddot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dxia)
call pedot(sx,sy,sz,cta,ctb,p1,p2,p3,p4,p5,epsi,beta,dxib)
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
call qadot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dx)
call qbdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dy)
call qcdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dz)
call qddot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dcta)
call qedot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dctb)
call padot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp1)
call pbdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp2)
call pcdot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp3)
call pddot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp4)
call pedot(q1,q2,q3,q4,q5,px,py,pz,xia,xib,epsi,beta,dp5)
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