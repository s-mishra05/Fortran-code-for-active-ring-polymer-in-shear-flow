program vmd

!This program makes vmd file from position data files in the COM frame 
! Small section (N=30) is named "O" while others are named "N" for visual aid

implicit none
real(8),dimension(200)::x,y,z
integer::i,j,k,n,step
character(len=50)::pos
real(8)::xcm

n=200

open(121,file="200.xyz") ! output vmd file name
do step =20000,200000000,20000

IF(step .lt. 10000) THEN

    WRITE(pos,'(a4,i4,a4)')"pos_",step,".dat"
           
   ELSE IF(step .lt. 100000) THEN
     WRITE(pos,'(a4,i5,a4)')"pos_",step,".dat"
           
   ELSE IF(step .lt. 1000000) THEN
     WRITE(pos,'(a4,i6,a4)')"pos_",step,".dat"
           
   ELSE IF(step .lt. 10000000) THEN
     WRITE(pos,'(a4,i7,a4)')"pos_",step,".dat"

   ELSE IF(step .lt. 100000000) THEN
     WRITE(pos,'(a4,i8,a4)')"pos_",step,".dat"

   ELSE IF(step .lt. 1000000000) THEN
     WRITE(pos,'(a4,i9,a4)')"pos_",step,".dat"
           

   Endif

 OPEN (130, FILE = pos)
        DO I=1,200
          READ(130,*)X(I),Y(I),Z(I)

        END DO

xcm = sum(x)/real(N) ! center of mass frame in x-direction (flow)
do i=1,N
x(i)=x(i)-xcm
end do



write(121,*)200

write(121,*)step
 DO i = 1,N
      if(i .le. 30)then
      write(121,*)"O",x(i),y(i),z(i)
      else
      write(121,*)"N",x(i),y(i),z(i)
      end if

ENDDO
 

 
 enddo

end program



