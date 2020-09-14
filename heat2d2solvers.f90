!  heat2d2solvers.f90 
!
!  FUNCTIONS:
!  heat2d2solvers - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: heat2d2solvers

!   by mehrdad baba hosein pour
!
!  PURPOSE:  Entry point for the console application.
! beyond sky!
!****************************************************************************


    program heat_equation

    implicit none
 integer :: nx,nt,ny,i,j,t,s,l,h        !nx= nodes in x , ny=nodes in y
 real,parameter:: alpha=0.5
 real :: dx,dy,p,st,dt,ttime   !ttime= total time 
 real,allocatable,dimension(:,:)::u_old,u_new,x,y,upre,unow,ufut    !upre=u previous !ufut=u in future
 
 print*,"inter your length = "
 read*,l
 print*,"inter your hieght = "
 read*,h
 print*,"inter your nodes in x = "
 read*,nx
 print*,"inter your nodes in y = "
 read*,ny
 nt=1500       !time stepping 
 dx=l/real(nx-1)  !space stepping in x
 dy=h/real(ny-1)   !space stepping in y
 print*,"inter your total time = "
 read*,ttime
 dt=ttime/real(nt-1)
 p=alpha*dt*((1/dx**2)+(1/dy**2))
 st=dt*nt
 
Write(*,*) 'please input your solver : 1 for explcit solver (FTCS) or 2 for alternative directional implcit (ADI)'
Read(*,*) s
 
 allocate(u_old(nx,ny),u_new(nx,ny),upre(nx,ny),unow(nx,ny),ufut(nx,ny),X(nx,ny),Y(nx,ny))
 
 !grid generation 

do i=1,nx
    do j=1,ny
    X(i,j)=(i-1)*dx
    Y(i,j)=(j-1)*dy
    end do
end do

 
 
! initial condition 
 u_old(:,:)=50
 
!boudary condition
u_old(2:,ny)=(0.33333)*(4*u_old(ny-1,2:ny-1)-u_old(ny-2,2:ny-1)) !up
do j=1,ny
    u_old(1,j)=28*0.66*dx*(0-100)+(0.33333)*(4*u_old(2,j)-u_old(3,j))    !left
end do 
u_old(nx,2:ny-1)=50 !right                       
u_old(2:,1)=100 !down
 

u_new=u_old

!show initial and boudary condition  
   Open(1,FIle='initial_condition.plt')      
   Write(1,*) 'Variables=X,Y,unew'
   write(1,*)'zone I=',nx,'J=',ny
do  J=1,ny
   Do I=1,nx
        Write(1,*) X(i,j),Y(i,j),u_old(i,j)
    END do 
end do 


if (s==1) then 
   
      write(*,*) "------------------------------------- FTCS--------------------------------------"
      print*,'solving time is = ',st
    write(*,*)'this solver is based on stablity number, generally stablity number must be lower than 0.5'
    write(*,*)'CAUTION'
    write(*,*)'if your DELTA X= DELTA Y the stability number must be lower than 0.25'
   write(*,*)'your stablility number is =  '
   write(*,*) p 
 
     print*,"please wait....."    
     !FTCS_solver
     
  do t=1,nt 
         do j=2,ny-1
             do i=2,nx-1
         
                 u_new(i,j)=(alpha*dt/dx**2)*(u_old(i+1,j)-2*u_old(i,j)+u_old(i-1,j))+(Alpha*dt/dy**2)*(u_old(i,j+1)-2*u_old(i,j)+u_old(i,j-1))+u_old(i,j)
                 
               
             end do
         end do 
  
        !show  result
   Open(2,FIle='FTCS_result.plt')      
   Write(2,*) 'Variables="X","Y","unew"'
   write(2,*)'zone I=',nx,'J=',ny
do  J=1,ny
   Do I=1,nx
        Write(2,*) X(i,j),Y(i,j),u_new(i,j)
    END do 
end do  

  u_old=u_new
  end do 
end if


     if (s==2) then 
         
        
         
   write(*,*)"------------------------------------- ADI--------------------------------------"
       write(*,*)'this solver is stabled under any circumstances'
       
        
 u_old(:,:)=50
 

u_old(2:,ny)=(0.33333)*(4*u_old(ny-1,2:ny-1)-u_old(ny-2,2:ny-1)) !up
do j=1,ny
    u_old(1,j)=28*0.66*dx*(0-100)+(0.33333)*(4*u_old(2,j)-u_old(3,j))   !left
end do 
u_old(nx,2:ny-1)=50 !right    !                   
u_old(2:,1)=100 !down


 unow(:,:)=50
 
unow(2:,ny)=(0.33333)*(4*u_old(ny-1,2:ny-1)-u_old(ny-2,2:ny-1)) !up
do j=1,ny
    u_old(1,j)=28*0.66*dx*(0-100)+(0.33333)*(4*u_old(2,j)-u_old(3,j))    !left 
end do 
unow(nx,2:ny-1)=50 !right    !                   
unow(2:,1)=100 !down

ufut=u_old

  print*,"please wait....."      
 do t=1,nt          
    ! X sweep 
    do i=2,nx-1
       do j=2,ny-1
    unow(i,j)=(alpha*dt/2)*((((unow(I+1,j)-2*unow(i,j)+unow(i-1,j))/dx**2))+(((u_old(i,j+1)-2*u_old(i,j)+u_old(i,j-1))/dy**2)))+u_old(i,j)
 
       end do 
    end do 
    
    !y sweep 
    do i=2,nx-1
        do j=2,ny-1
            ufut(i,j)=(alpha*dt/2)*((((unow(i+1,j)-2*unow(I,J)+unow(i-1,j))/dx**2))+(((ufut(i,j+1)-2*ufut(i,j)+ufut(i,j-1))/dy**2)))+unow(i,j)
         
        end do 
    end do 
    
            
  ! show result  
   Open(3,FIle='ADI_result.plt')
Write(3,*) 'Variables=X,Y,ufut'
write(3,*)'zone I=',nx,'J=',ny
do  J=1,ny
   Do I=1,nx
        Write(3,*) X(i,j),Y(i,j),ufut(i,j)
    END do 
end do  
    u_old=unow
    unow=ufut

 end do 
     end if 
   
     print*,"done!"
pause

end program heat_equation





