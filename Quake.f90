!---------------------------------------------
!Variables for use by all program units
!---------------------------------------------
Module constants
Implicit None
Integer,parameter :: N = 500                          !Number of blocks
Double precision,parameter :: dt = 1.0E-1            !step size
Integer,parameter :: steps = 15000                    !Number of time steps
Double precision, parameter :: tmax = dt*DBLE(steps) !Total time
Double precision, parameter :: vo = 1.0E-2           !Plate speed
Double precision, parameter :: F0 = 50.0E0           !Maximum friction force
Double precision, parameter :: vf = 5.0E-4           !Velocity dependence of friction
Double precision,parameter :: dto2 = dt/2.0E0
Double precision,parameter :: two = 2.0E0
Double precision,parameter :: six = 6.0E0
Double precision,parameter :: half = 0.5E0
End module constants
!---------------------------------------------
!Main Program
!---------------------------------------------
Program Quake
Use constants
Implicit None
Double precision, Dimension(:), Allocatable :: dx
Double precision, Dimension(:), Allocatable :: xp
Double precision, Dimension(:), Allocatable :: m
Double precision, Dimension(:), Allocatable :: kp
Double precision, Dimension(:), Allocatable :: kc
Double precision t
Integer I,J,seed
External rpm; Double precision rpm
Allocate(dx(N));Allocate(xp(N));Allocate(m(N));Allocate(kp(N));Allocate(kc(N))

!Initializations
t = 0.0E0
Call system_clock(seed)
Call srand(seed)
Do I = 1,N
  dx(I) = 0.0E0
  dx(I) = 5.0E0*rpm()
  xp(I) = 0.0E0
  xp(I) = 2.0E0*rpm()
  m(I) = 1.0E0 + 0.2*rpm()
  kp(I) = 40.0E0 !+ 10.0E0*rpm()
  kc(I) = 2.5E2 !+ 50.0E0*rpm()
End Do

!Main Loop
Do While(t.lt.tmax)
  t = t + dt
  Call RK(t,dx,xp,m,kp,kc)
End Do

Deallocate(dx);Deallocate(xp);Deallocate(m);Deallocate(kp);Deallocate(kc)
End program Quake
!---------------------------------------------
!Subroutines and Functions
!---------------------------------------------
Subroutine RK(t,dx,xp,m,kp,kc)
Use constants
Implicit None
External F1;Double precision F1
External F2;Double precision F2
External signf; Double precision signf
Double precision t,dx(N),xp(N),xnext,xprev,dxs,xps,m(N),kp(N),kc(N)
Double precision k1,k2,k3,k4,l1,l2,l3,l4
Integer I
Do I = 1,N
  If(I.eq.1) then
    xprev = 0.0E0
  Else
    xprev = dx(I-1)
  End If

  If (I.eq.N) then
    xnext = 0.0E0
  Else
    xnext = dx(I+1)
  End If

  k1 = dt*F1(xp(I))
  l1 = dt*F2(t,dx(I),xp(I),xprev,xnext,m(I),kp,kc,I)

  k2 = dt*F1(xp(I) + l1/two)
  l2 = dt*F2(t+dto2,dx(I)+k1/two,xp(I)+l1/two,xprev,xnext,m(I),kp,kc,I)

  k3 = dt*F1(xp(I) + l2/two)
  l3 = dt*F2(t+dto2,dx(I)+k2/two,xp(I)+l2/two,xprev,xnext,m(I),kp,kc,I)

  k4 = dt*F1(xp(I) + l3)
  l4 = dt*F2(t+dt,dx(I)+k3,xp(I)+l3,xprev,xnext,m(I),kp,kc,I)

  dxs = (k1 + two*k2 + two*k3  +k4)/six
  dx(I) = dx(I) + dxs
  xps = (l1 + two*l2 + two*l3 + l4)/six

  If(signf(xps)*signf(xp(I)).eq.1.0E0) then
    xp(I) = xps + xp(I)
  End If

  If(signf(xps)*signf(xp(I)).eq.-1.0E0) then
    xp(I) = 0.0E0
  End If

  If(dxs/dt.ge.1.0E-4)Write(*,*) t,dxs/dt
End Do
Return
End subroutine RK
!---------------------------------------------
Function F1(xp)
Use constants
Implicit None
Double precision F1,xp
F1 = xp
Return
End function F1
!---------------------------------------------
!Spring forces
Function F2(t,x,xp,xprev,xnext,m,kp,kc,I)
Use constants
Implicit None
Double precision F2,t,x,xp,xprev,xnext
Double precision m,kp(N),kc(N)
Double precision springs,f,limit
External friciton; Double precision friction
External signf; Double precision signf
Integer I
limit = 5.0E-4

!Free left end
If(xprev.eq.0.0E0) then
  springs = kc(I+1)*xnext - kc(I)*x + kp(I)*(vo*t-x)
End If

!Free right end
If(xnext.eq.0.0E0) then
  springs = kc(I-1)*xprev - kc(I)*x + kp(I)*(vo*t-x)
End If

!Middle positions
If(xprev.ne.0.0E0 .and. xnext.ne.0.0E0) then
  springs = kc(I+1)*xnext+kc(I-1)*xprev-kc(I)*2.0E0*x + kp(I)*(vo*t - x)
End If

If(abs(xp).lt.limit .and. springs.lt.F0) then
  F2 = 0.0E0
End If

If(abs(xp).lt.limit .and. springs.gt.F0) then
  F2 = springs + (-1.0E0*signf(springs)*F0)
End If

If(abs(xp).gt.limit) then
  f = friction(xp)
  F2 = springs + f
End If

F2 = m*F2

Return
End function F2
!---------------------------------------------
!Friction force
Function friction(xp)
Use constants
Implicit None
Double precision friction,xp,limit
External signf;Double precision signf
!friction = (-1.0E0 * F0 * signf(xp)) / (1.0E0 + abs(xp/vf))
limit = 5.0E-4
If(xp.lt.limit) friction = 0.0E0
If(xp.ge.limit) friction = -1.0E0*(signf(xp)*F0)/2.0E0
Return
End function friction
!---------------------------------------------
!Sign function
Function signf(x)
Implicit None
Double precision x,signf
If(x.eq.0.0E0) signf = 0.0E0
signf = x/abs(x)
Return
End function signf
!---------------------------------------------
!Random number between -1 and 1
Function rpm()
Use constants
Implicit None
Double precision rpm
rpm = rand()
If(rpm.lt.0.5E0) rpm = -1.0E0*rand()
If(rpm.gt.0.5E0) rpm = rand()
Return
End function rpm
!---------------------------------------------

                                                                 
