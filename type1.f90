program langevin

! This code simulates type-I active ring at Wi=370

     implicit none

     integer, parameter :: n= 200
     real(8), parameter :: box=200.0,pi=4.0*atan(1.0),sigma=2.0d0,rc=2**(1/6)*sigma,ft=1.0d0
     real(8), parameter    :: eps = 1.0d0, temp = 0.10d0, kbt = 0.10d0, mass = 1.0d0,beeta=1.0d0,shearr=0.050d0
     real(8)               :: ranfx(n),ranfy(n),ranfz(n),fsx(n), fsy(n), fsz(n),fax(n),fay(n),faz(n)
     real(8)               :: x(n), y(n), z(n),strain,v_le(n),x_le(n),xl(n),yl(n),zl(n)
     real(8)               :: vx(n),vy(n),vz(n),shear
     real(8)               :: v2, scale1, vxav, vyav, vzav
     real(8)               :: fx(n), fy(n), fz(n), oldtfx(n), oldtfy(n), oldtfz(n), tfx(n), tfy(n), tfz(n)
     real(8)               ::  dt, ke, pe,t, pe_sp,pa, vcm
     integer(8)            :: i, itr,it,k,sim_time,JSEED,ISEED,re_it
     integer(4)            :: timeArray(3),ierror
     character(len=50)     :: pos,vel,lpos
     
     
     pa=0.0
     dt = 5*0.0001d0
     sim_time = 50000
     itr = int(sim_time/dt)
     

     open(3,file='energy.dat')

     call itime(timeArray)
     IF( timeArray(3) .lt. 10 ) then
     if( timeArray(3) .eq. 0 ) then 
       JSEED = 100*timeArray(1) + timeArray(2) + 10000 
     else
       JSEED = 100*timeArray(1) + timeArray(2) + timeArray(3)*10000
      endif
     else 
       JSEED = 100*timeArray(1) + timeArray(2) + timeArray(3)*100 
     endif
     ISEED = - JSEED
     
!------------------  INITIALISATION  -------------------- 
     open(500,file='ini.dat')  ! initialize with ini.dat data file
!     open(600,file='re_vel.dat') 
!     open(700,file='restart_time.dat')
!      call init_pos(x, y, z)
     do i = 1,n
       read(500,*) x(i), y(i), z(i)
!       read(600,*) vx(i), vy(i), vz(i)
     end do
!     read(700,*) re_it
     close(500)!;close(600);close(700)

     call init_vel(vx, vy, vz)
     call spring(x,y,z,fsx,fsy,fsz,pe_sp)
     call lj_force(fx, fy, fz, x, y, z, pe)
     call random_force(ranfx,ranfy,ranfz)
     call active(x,y,z,fax,fay,faz)

     strain = 0.0; x_le = 0.0
     xl = x; yl = y; zl = z 
!------------------ MAIN TIME LOOP-----------------------
     do it = 1, itr

     if (it .gt. 1000000)then  ! shear is applied after small time period 
         shear = shearr
!         strain = shear*dt*it
!         strain = strain - anint(strain/box)*box
         v_le = shear*y
         x_le = x_le + v_le*(dt)
     else
             shear = 0.0d0
     end if

     tfx = fx + ranfx + fsx + fax 
     tfy = fy + ranfy + fsy + fay
     tfz = fz + ranfz + fsz + faz

! Friction incorporated here 
     x = x + vx * dt + 0.5d0 * (tfx - beeta*vx)/mass * (dt**2) + shear*(y)*dt
     y = y + vy * dt + 0.5d0 * (tfy - beeta*vy)/mass * (dt**2) 
     z = z + vz * dt + 0.5d0 * (tfz - beeta*vz)/mass * (dt**2)

       
     oldtfx = tfx
     oldtfy = tfy
     oldtfz = tfz

     call spring(x,y,z,fsx,fsy,fsz,pe_sp)
     call lj_force(fx, fy, fz, x, y, z, pe)
     call random_force(ranfx,ranfy,ranfz)
     call active(x,y,z,fax,fay,faz)

     tfx = fx + ranfx + fsx + fax 
     tfy = fy + ranfy + fsy + fay
     tfz = fz + ranfz + fsz + faz
 
     vx = vx + dt * 0.5d0 * ( oldtfx + tfx )/mass - dt*vx*beeta/mass + dt*shear*vy
     vy = vy + dt * 0.5d0 * ( oldtfy + tfy )/mass - dt*vy*beeta/mass     
     vz = vz + dt * 0.5d0 * ( oldtfz + tfz )/mass - dt*vz*beeta/mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! modified PBC with Lees-Edwards boundary conditions
     do i = 1, n
     
     if(x(i) .gt. box)then
        xl(i)=x(i)-anint(x(i)/box)*box
     end if

     if(y(i) .lt. 0.0)then
        xl(i) = x(i)-x_le(i)
        yl(i) = y(i) + box
     end if
     if(y(i) .gt. box)then
         xl(i) = x(i)+x_le(i)
         yl(i) = y(i) - box
     end if
     if(z(i) .gt. box)then
            zl(i) = z(i)-box
     end if
     if(z(i) .lt. 0.0)then
           zl(i) = z(i)+box
     end if
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

     if(mod(it,10000) == 0) then
       ke = 0.5d0 * sum(vx ** 2 + vy ** 2 + vz **  2 )  !KINETIC ENERGY
       ke = ke / dfloat(n)
       t = (sum(vx ** 2  + vy ** 2 + vz **  2))/(3*(float(n))) !TEMPERATURE

       write(3,*) it*dt, t, ke, pe+pe_sp, ke+pe+pe_sp ! time, temp, kinetic, pot, tot en
     end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!save for restart !!!!!!!!!!!!!!!!!!!!!!!!
      if (mod(it,50000) ==0) then
        
      open(500, file = 're_pos.dat')
        DO i = 1,n
         WRITE(500,*) x(i), y(i), z(i)
        END DO
        WRITE(500,*) x(1), y(1), z(1)
      close(500)
      open(600, file = 're_vel.dat')
        do i = 1,n
          write(600,*) vx(i), vy(i), vz(i)
        end do
      close(600)
      open(700, file='restart_time.dat')
       write(700,*) it
      close(700)

     end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   This section writes position data for every 20000^th iteration in separate file
     if(mod(it,20000)==0)then
        write(pos, '(I9)') it

        write(lpos, '(I9)') it
        open(100,file='pos_'//trim(adjustl(pos))//'.dat') ! global coordinates

        open(300,file='lp_'//trim(adjustl(lpos))//'.dat') ! local coordinates

        DO i = 1,n
         WRITE(100,*) x(i), y(i), z(i)
         WRITE(300,*) xl(i), yl(i), zl(i)

        END DO
        WRITE(100,*) x(1), y(1), z(1)
        WRITE(300,*) xl(1), yl(1), zl(1)
       close(100);close(300)

     end if


   
     end do

     contains

!--------------------------------------------------------------------------------
    

     subroutine init_pos(x, y, z)

     implicit none
      real(8), parameter :: r = 250.0d0/pi , theta = 1.80d0, m = 4.0*atan(1.0)/180.0

     real(8), intent(inout) :: x(n), y(n), z(n)
     integer :: i
     x = 100.0; y = 100.0; z = 100.0
     do i = 2,n
      x(i) = x(i) + r*cos(theta*m*(i-1))
      y(i) = y(i) + r*sin(theta*m*(i-1))
     end do
     x(1) = x(1) + r

     end subroutine
      
     subroutine init_vel(vx, vy, vz)

     implicit none
     real(8), intent(inout) :: vx(n), vy(n), vz(n)
!     real(8), intent(out) :: t
     real(8) :: v2, scale1, vxav, vyav, vzav,ran2
 
     do i = 1, n
          vx(i) = (ran2(ISEED) - 0.50d0)
          vy(i) = (ran2(ISEED) - 0.50d0)
          vz(i) = (ran2(ISEED) - 0.50d0)
     end do
     v2 = sum(vx ** 2 + vy ** 2 + vz ** 2)
     v2 = v2 / dfloat(n)
     scale1  = dsqrt(3 * kbt / v2)
 
     vxav = sum(vx) / dfloat(n)
     vyav = sum(vy) / dfloat(n)
     vzav = sum(vz) / dfloat(n)
 
     vx = (vxav-vx) * scale1
     vy = (vyav-vy) * scale1
     vz = (vzav-vz) * scale1

     return
     end subroutine

     SUBROUTINE PBC(x,y,z)
     implicit none
     REAL(kind=8), INTENT(INOUT)::x(n),y(n),z(n)
     do i =1,n
     if(x(i)>box) then 
      x(i)=x(i)-box
     endif
     if(y(i)>box) then 
      y(i)=y(i)-box
     endif
     if(z(i)>box) then 
      z(i)=z(i)-box
     endif
     if(x(i)<0.0) then 
      x(i)=x(i)+box
     endif
     if(y(i)<0.0) then 
      y(i)=y(i)+box
     endif
     if(z(i)<0.0) then 
      z(i)=z(i)+box
     endif
     end do
    END SUBROUTINE PBC 


     subroutine random_force(ranfx,ranfy,ranfz)
      
     implicit none
     integer(kind = 8)             :: i
     
     real(kind= 8) , parameter     :: e1=0.36787944d0
     real(kind= 8)                 :: ran2,p,u1,u2
     real(kind = 8)  ,intent(out)  :: ranfx(n), ranfy(n), ranfz(n)
   
      p = sqrt(2.0d0*beeta*kbt/dt)
     do i= 1,n
  1   u1 = ran2(ISEED)
     u2 = ran2(JSEED)
     u2 = (2.*u2-1.)*sqrt(2.*e1)
     if(-4.*u1*u1*log(u1).lt.u2*u2 ) go to 1
     ranfx(i) = u2/u1
     end do
     do i= 1,n
  2   u1 = ran2(ISEED)
     u2 = ran2(JSEED)
     u2 = (2.*u2-1.)*sqrt(2.*e1)
     if(-4.*u1*u1*log(u1).lt.u2*u2 ) go to 2
     ranfy(i) = u2/u1
     end do
     do i= 1,n
  3   u1 = ran2(ISEED)
     u2 = ran2(JSEED)
     u2 = (2.*u2-1.)*sqrt(2.*e1)
     if(-4.*u1*u1*log(u1).lt.u2*u2 ) go to 3
     ranfz(i) = u2/u1
     end do

     ranfx = ranfx*p
     ranfy = ranfy*p
     ranfz = ranfz*p

     return
     end subroutine
    

     subroutine lj_force(fx, fy, fz, x, y, z, pe)
       implicit none
       real(8), intent(inout) :: fx(n), fy(n), fz(n), x(n), y(n), z(n)
       real(8), intent(inout) :: pe
       real(8) :: dx, dy, dz, dr, fr,s6, s12
       integer(8) :: i, j

 
       fx = 0.0d0
       fy = 0.0d0
       fz = 0.0d0
       pe = 0.0d0
 
       do i =  1, n - 1
          do j = i + 1, n
             dx = x(i) - x(j)
             dy = y(i) - y(j)
             dz = z(i) - z(j)
 
             dx = dx - box * anint(dy / box) 
             dy = dy - box * anint(dy / box)
             dz = dz - box * anint(dz / box)
 
             dr = dsqrt(dx ** 2 + dy ** 2 + dz ** 2)
             s6 = (sigma/dr)**6
             s12 = s6**2
             if(dr .le. rc) then
                   fr = (24.0d0/dr) * (2.0d0*s12 - s6)  !sigma= eps =1
                   pe = pe + 4.0d0 *eps* (s12 - s6) + eps
 
                   fx(i) = fx(i) + fr * dx / dr
                   fy(i) = fy(i) + fr * dy / dr
                   fz(i) = fz(i) + fr * dz / dr
 
                   fx(j) = fx(j) - fr * dx / dr
                   fy(j) = fy(j) - fr * dy / dr
                   fz(j) = fz(j) - fr * dz / dr
 
             end if
      
          end do
       end do
    
 
       pe = pe / dfloat(n) !POTENTIAL ENERGY
       return
       end subroutine

       subroutine spring(x,y,z,fsx,fsy,fsz,pe_sp)
       implicit none

       real(8), parameter     :: ks = 2000.0, l0 = 2.50d0
       real(8), intent(inout) :: x(n),y(n),z(n),fsx(n),fsy(n),fsz(n)
       real(8), intent(inout) :: pe_sp
       real(8)                :: dx,dy,dz,dr,dx_end,dy_end,dz_end,dr_end,ff

       integer(8)             :: i,j

       fsx = 0.0d0
       fsy = 0.0d0
       fsz = 0.0d0
       pe_sp = 0.0d0

       dx_end = x(n) - x(1)
       dy_end = y(n) - y(1)
       dz_end = z(n) - z(1)
       dx_end = dx_end-(real(box)*nint(dx_end/real(box)))
       dy_end = dy_end-(real(box)*nint(dy_end/real(box)))
       dz_end = dz_end-(real(box)*nint(dz_end/real(box)))
       dr_end = sqrt((dx_end**2)+(dy_end**2)+(dz_end**2))
       ff = ks*(dr_end-l0)

       fsx(n)= fsx(n)- ff*(dx_end/dr_end)
       fsy(n)= fsy(n)- ff*(dy_end/dr_end)
       fsz(n)= fsz(n)- ff*(dz_end/dr_end)


       fsx(1)= fsx(1)+ ff*(dx_end/dr_end)
       fsy(1)= fsy(1)+ ff*(dy_end/dr_end)
       fsz(1)= fsz(1)+ ff*(dz_end/dr_end)
       pe_sp = pe_sp + 0.5*ks*(dr_end-l0)**2


       do i = 1, n-1
          j = i+1
          dx = x(i)-x(j)
          dy = y(i)-y(j)
          dz = z(i)-z(j)
          dx = dx-(real(box)*nint(dx/real(box)))
          dy = dy-(real(box)*nint(dy/real(box)))
          dz = dz-(real(box)*nint(dz/real(box)))
          dr = sqrt((dx**2)+(dy**2)+(dz**2))
          ff = ks*(dr-l0)
          
         fsx(i)= fsx(i)- ff*(dx/dr)
         fsy(i)= fsy(i)- ff*(dy/dr)
         fsz(i)= fsz(i)- ff*(dz/dr)
         fsx(j)= fsx(j)+ ff*(dx/dr)
         fsy(j)= fsy(j)+ ff*(dy/dr)
         fsz(j)= fsz(j)+ ff*(dz/dr)
         pe_sp = pe_sp + 0.5*ks*(dr-l0)**2
         pe_sp = pe_sp/dfloat(n)
      end do   

      return
      end subroutine

     subroutine active(x,y,z,fax,fay,faz)
     implicit none

     integer(8)            :: i,j
     real(8)  ,intent(inout) :: fax(n), fay(n), faz(n), x(n), y(n), z(n)
     real(8)               :: dx,dy,dz,dx1,dy1,dz1,dxn,dyn,dzn,dr,dr1,drn

     fax = 0.0d0; fay = 0.0d0; faz = 0.0d0

     do i = 1,n-1
 
        dx = x(i+1) - x(i-1)
        dy = y(i+1) - y(i-1)
        dz = z(i+1) - z(i-1)
        dx = dx - box * anint(dx / box)
        dy = dy - box * anint(dy / box)
        dz = dz - box * anint(dz / box)
        dr = sqrt(dx**2 + dy**2 + dz**2)
        fax(i) = fax(i) + ft*(dx/dr)
        fay(i) = fay(i) + ft*(dy/dr)
        faz(i) = faz(i) + ft*(dz/dr)
       
      end do

      dxn = x(1)-x(n-1); dyn = y(1)-y(n-1); dzn = z(1)-z(n-1)
      dxn = dxn - box * anint(dxn / box)
      dyn = dyn - box * anint(dyn / box)
      dzn = dzn - box * anint(dzn / box)
      drn = sqrt(dxn**2 + dyn**2 + dzn**2)
      fax(n) = fax(n) + ft*(dxn/drn)
      fay(n) = fay(n) + ft*(dyn/drn)
      faz(n) = faz(n) + ft*(dzn/drn)
  
      return
      end subroutine

      
 end program

 FUNCTION ran2(idum)
    INTEGER(8) idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL(8) ran2,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,&
    IQ2=52774,IR1=12211, IR2=3791,NTAB=32, &
    NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
         do j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
          enddo
       iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
    END function ran2
    
 
     

    
