        program fghn

                implicit none 
                integer :: i,j,np,nx,k,NE,nvib,ntime
                real*8 ::dx,dp,xmin,xmax,pi,mass
                real*8 ::D,B,Xe,Kb,Tp1,Tp2,Tp3,sigma
                real*8 :: Tp4,omega,time
                real*8 :: Z1,Z2,Z3,Z4
                real*8 :: Xv1,Xv2,Xv3,Xv4 
                real*8, allocatable :: x(:),E(:),V(:),a(:),gauss(:)
                real*8, allocatable :: XvTp2(:),XvTp3(:),XvTp4(:)
                real*8, allocatable :: H(:,:),T(:,:),S(:,:),C(:,:)
                real*8, allocatable :: rho(:,:)


                integer :: LWORK,INFO
                real*8, allocatable :: WORK(:)

                nx=701
                NE=67
                np=(nx-1)/2
                xmin=0.0
                xmax=6.0
                dx=(xmax-xmin)/(nx-1)
                pi=2.0*asin(1.0)
                dp=(2.0*pi)/((nx-1)*dx)
                mass=12861.93
                D=0.3640
                B=1.4285
                Xe=2.075
                Kb=3.16685d-6
                Tp1=4
                Tp2=300
                Tp3=1000
                Tp4=4000
                sigma=0.25d0

                allocate(x(nx),V(nx))

                do i=1,nx
                x(i)=xmin+(i-1)*dx
                enddo

                allocate(H(nx,nx),S(nx,nx),T(nx,nx))

                do i=1,nx 
                do j=1,nx
                T(i,j)=0.0
                S(i,j)=0.0
                do k=1,np
                T(i,j)= T(i,j)+(((k*dp)**2)*cos(k*dp*(x(j)-x(i))))
                enddo
                T(i,j)=T(i,j)/((nx-1)*mass)
                H(i,j)=T(i,j)
                if(i.eq.j)then

                        V(i)=D*(1-exp(-B*(x(i)-Xe)))**2
                        H(i,j)=H(i,j)+V(i)
                        S(i,j)=1.0
                endif
                enddo
                enddo

                allocate(C(nx,nx),E(nx))

                LWORK=3*nx
                ALLOCATE(WORK(LWORK))
                CALL DSYGV (1,'V','U',nx,H,nx,S,nx,E,WORK,LWORK,INFO)
                DEALLOCATE(WORK)





                !do i=1,NE
                !write(6,*)E(5),V(5)
                !enddo
                !write(6,*)'INFO',INFO

                C=(H/dsqrt(dx))

                nvib=0

                do i=1,nx
                if(E(i).lt.D)then
                !write(6,*)i,E(i)
                nvib=nvib+1
                endif
                enddo

                write(6,*)'Número de estados vibracionales:',nvib
                
                allocate(a(nvib),gauss(nx))


                a=0.d0
                  do i=1,nvib
                    do j=1,nx
                     gauss(j)=dexp(-(x(j)-(Xe+1.d0))**2/(4.d0*sigma**2))
                     gauss(j)=gauss(j)/dsqrt(sigma*dsqrt(2.d0*pi))
                      a(i)=a(i)+gauss(j)*C(j,i)*dx
                    enddo
                    write(6,*)a(i)
                  enddo
               ! write(6,*)sum(a(:)*a(:))

                ntime=1

                allocate(rho(nx,ntime))

                time=0.d0

                do i=1,nx
                  do j=1,nvib
                    do k=1,nvib

                     omega = E(j) - E(i)
                     rho(i,1)=rho(i,1)+a(j)*a(k)*C(i,j)*C(i,k)
     $                        *dcos(omega*time)
                     enddo
                  enddo   
                write(77,*)x(i),rho(i,1)
               enddo
                

                stop'jaja'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





                open(1,file='wfn.dat')

                do i=1,nx
                 write(1,*)x(i),V(i),
     $           (C(i,1)**2+E(1)),
     $           (C(i,2)**2+E(2)),
     $           (C(i,3)**2+E(3)),
     $           (C(i,18)**2+E(18))
                enddo
                
                close(1)


                Z1=0.0
                do i=1,NE
                Z1=Z1+(dexp(-E(i)/(Kb*Tp1)))
                enddo

                Z2=0.0
                do i=1,NE
                Z2=Z2+(dexp(-E(i)/(Kb*Tp2)))
                enddo
                
                Z3=0.0
                do i=1,NE
                Z3=Z3+(dexp(-E(i)/(Kb*Tp3)))
                enddo

                Z4=0.0
                do i=1,NE
                Z4=Z4+(dexp(-E(i)/(Kb*Tp4)))
                enddo

                allocate(XvTp2(NE),XvTp3(NE),XvTp4(NE))

                
                write(6,*)'Z1',Z1
                write(6,*)'Z2',Z2
                write(6,*)'Z3',Z3
                write(6,*)'Z4',Z4

                Xv1=((dexp(-E(1)/(Kb*Tp1)))/Z1)
                Xv2=((dexp(-E(1)/(Kb*Tp2)))/Z2)
                Xv3=((dexp(-E(1)/(Kb*Tp3)))/Z3)
                Xv4=((dexp(-E(1)/(Kb*Tp4)))/Z4)

                write(6,*)'Xv1',Xv1
                write(6,*)'Xv2',Xv2
                write(6,*)'Xv3',Xv3
                write(6,*)'Xv4',Xv4
               

                do i=1,NE
                XvTp2(i)=((dexp(-E(i)/(Kb*Tp2)))/Z2)
                XvTp3(i)=((dexp(-E(i)/(Kb*Tp3)))/Z3)
                XvTp4(i)=((dexp(-E(i)/(Kb*Tp4)))/Z4)
                enddo


!                write(6,*)'Area'
!                Area=0.0
!                do i=1,nx
!                Area=Area+dx*C(i,1)*C(i,1)
!                enddo
!                write(6,*)'Area',Area,'dx',dx,dsqrt(dx)

!                C=H/dsqrt(dx) 

!                Area=0.0
!                do i=1,nx
!                Area=Area+C(i,1)*dexp(-(x(i)-4.0)**2)
!                enddo
!                write(6,*)'Area',Area

                open(2,file='Xv.dat')
                do i=1,NE
                write(2,*)i,XvTp2(i),XvTp3(i),XvTp4(i)
                enddo
                close(2)


        end program
