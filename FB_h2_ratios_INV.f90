        program BF
        use prob; use kinds; use denseop; use matrix_oper
!$ use omp_lib

        implicit none
        integer , allocatable :: ilev(:),nlev1(:),ped1(:,:),cum_nlev1(:)
        REAL (r8), allocatable :: ZAZt(:,:),V_inv(:,:),ZZt(:,:,:)
        REAL (r8), allocatable :: X(:,:),y(:),e(:),Z(:,:,:),Za(:,:)
        REAL (r8), allocatable :: dat(:,:)
        REAL (r8), allocatable :: A_mat(:,:)   
        REAL (r8), allocatable :: aux(:,:),XX(:,:),XY(:),b(:),integral(:),log_pcond_h2_0_muestras(:),ratios_fact(:),ratios_fact_a(:)
        REAL (r8), allocatable :: log_pcond_ratios_0(:),pmar_ratios_0(:),log_pcond_ratios_0_muestras(:,:)
        
        REAL (r8)            :: ys,vp,h2,h2a,unif,se,vp_orig,gl,vare,c,log_det,diff
        REAL (r8)            :: den1,nmet,imet,den2,acon,fact,den0,log_pcond_h2_0,pmar_h2_0,sut,suma,sum_ratios,zero1=0.0,t1,t2
        integer              :: nanim,ndat,neff1,nround,neq,vne_1,nintegral
        integer              :: k,id,id_s,id_d,i,j,io,iad,ii,kk,jj,ijk,ij
        integer              :: saved_rounds,nburnin
        integer              :: num_threads_original,nfijos,max_nlev,nproces, n_inversions
        character *50        :: dat_file,ped1_file
              
        
                                
        write (*,*) "FICHERO DE DATOS"
        read (*,*) dat_file
        print *,dat_file
        
        write (*,*) "FICHERO DE PEDIGREE"
        read (*,*) ped1_file
        print *,ped1_file
        
        write (*,*) "# registros en datos"
        read (*,*) ndat
        print *,ndat
        
        write (*,*) "# individuos en ped1igree"
        read (*,*) nanim
        print *,nanim
        
        write (*,*) "# de efectos fijos y aleatorios "
        read (*,*) neff1
        allocate(ilev(neff1), nlev1(0:(neff1)),cum_nlev1(0:(neff1)),ratios_fact(neff1),ratios_fact_a(neff1)  )
        nlev1=0
        print *,neff1
        
        write (*,*) "niveles de los efectos (excluyendo al genético) "
        read (*,*) nlev1(1:(neff1))
        print *,nlev1(1:(neff1))
        cum_nlev1=0 
        do i=1,neff1
         cum_nlev1(i)= cum_nlev1(i-1)+nlev1(i)
        enddo
        write (*,*) "niveles acumulados de los efectos (excluyendo al genético) "
        print *,cum_nlev1(1:(neff1))
        
        ratios_fact=0.0
        ratios_fact_a=0.0
        write (*,*) "Valor inicial para los ratios de las varianzas sobre la varianza fenotipica"
        read (*,*)ratios_fact(1:neff1)
        print *,ratios_fact

        write (*,*) "Valor inicial para varianza Fenotipica"
        read (*,*)vp
        print *,vp
        
        write (*,*) "Valor inicial para h2 "
        read (*,*)h2
        print *,h2
        
        write (*,*) "# de iteraciones "
        read (*,*)nround        
        print *,nround
        
        write (*,*) "# de burnin "
        read (*,*)nburnin               
        print *,nburnin
         
        write (*,*)"# de procesadores"
        read(*,*)nproces
        print *,nproces
                
                
                
        vp_orig=vp
        gl=-2

        
               
        nfijos=0
        neq=0
        max_nlev=0
        do i=1,neff1
         if(ratios_fact(i).eq.0) then
          nfijos=nfijos+1
          neq=neq+nlev1(i)
         endif 
         if(nlev1(i).gt.max_nlev)max_nlev=nlev1(i)
        enddo       
        print *,nfijos,neq,max_nlev
               
        nintegral=100
        saved_rounds=0

        allocate(ped1(nanim,3))
        allocate(ZAZt(ndat,ndat),V_inv(ndat,ndat),ZZt(neff1,ndat,ndat))
        allocate(X(neq,ndat),e(ndat),y(ndat),dat(ndat,nfijos))
        allocate(A_mat(nanim,nanim),Za(nanim,ndat) )
        allocate(aux(neq,ndat),XX(neq,neq),XY(neq),b(neq),Z(neff1,max_nlev,ndat) )
        allocate(integral(0:nintegral),log_pcond_h2_0_muestras(nround-nburnin))
        allocate(log_pcond_ratios_0(neff1),pmar_ratios_0(neff1),log_pcond_ratios_0_muestras(neff1,nround-nburnin) )
        
        log_pcond_h2_0_muestras=0
        log_pcond_ratios_0_muestras=0
        
        num_threads_pregs=nproces
                 
!$omp parallel
!$omp master
!$ num_threads_original = omp_get_num_threads()
!$omp end master
!$omp end parallel
!$
!$ if(num_threads_pregs > 0) then
!$    call omp_set_num_threads(num_threads_pregs)
!$ end if
   print *, num_threads_original, num_threads_pregs
   print*,''
   print '(a)'      ,' *--------------------------------------------------------------*'
   print '(3a)'     ,' *                  Dist Version ',version,'                    *'
   print '(a)'      ,' *--------------------------------------------------------------*'
!$ if (num_threads_pregs>0) then
!$   print '(a,i2,a)'      ,' *             Optimized OpenMP Version 1 - ',num_threads_pregs,' threads            *'
!$ else
!$   print '(a,i2,a)'      ,' *             Optimized OpenMP Version 2 - ',num_threads_original,' threads            *'
!$ endif

               
        
        open(68,file='salida')
!        open(876,file='Muestras_fijos') 
!        open(877,file='Muestras_VC') 

        OPEN(11,FILE=ped1_file)
        io=0
        DO 
          READ(11,*,iostat=io)ID,ID_s,ID_d
          if(io.ne.0)exit
          ped1(id,1)=id
          ped1(id,2)=id_s
          ped1(id,3)=id_d
        ENDDO
        CLOSE(11)
        print *,'Creando A...'
        A_mat=-99.

        DO i=1,nanim
         DO j=i,nanim
          A_mat(i,j)=parentesco(i,j)
          if(i.ne.j)A_mat(j,i)=A_mat(i,j)
         ENDDO
        ENDDO
       
        
        open(11,file=dat_file)
        i=0
        x=0
        z=0
        do 
         read(11,*,iostat=io)ys,ilev,iad
         if(io.ne.0)exit
         i=i+1
         y(i)=ys
         dat(i,:)=ilev
         do j=1,nfijos
          k=cum_nlev1(j-1)+ilev(j)
          x(k,i)=1
         enddo
         do j=nfijos+1,neff1
           k=ilev(j)
           Z(j,k,i)=1
         enddo          
         Za(iad,i)=1      
        enddo
        close(11)
        
!        open(878,file='X_matrix')
!        do i=1,neq
!         write(878,'(1000(f10.7,x))') (x(i,j),j=1,ndat)
!        enddo 
!        close(878)
        
        ZAZt=0
          ZAZt=mymatmul(transpose( Za),mymatmul(A_mat,Za) )
!        open (123,file='ZAZt')
!        do i=1,ndat
!          write(123,'(10000(f10.7,x))')(ZAZt(i,j),j=1,ndat)
!        enddo
!        close(123)

        ZZt=0
        do k=1,neff1
         if(ratios_fact(k).ne.0)then
          ZZt(k,:,:)=mymatmul(transpose( Z(k, 1:nlev1(k) ,:)),Z(k,1:nlev1(k),:) )
         endif
        enddo
                        
!        open (333,file='ZZt')
!        do k=1,neff1
!        if(ratios_fact(k).ne.0)then
!         do i=1,ndat
!           write(333,'(10000(f10.7,x))')(ZZt(k,i,j),j=1,ndat)
!         enddo
!         write(333,*)'#####'
!        endif 
!        enddo        
!        close(333)
      
!        b=(/10,11,13,0,2,-3,0,-2/)
         b=0.
        imet=1
        fact=0.01 !parametro para el random walk M-H
  
!        open(455,file='V_inv')
!        open(444,file='X.V_inv.y')
!        open(466,file='X.V_inv.X')
!        open(456,file='V')
          
   call cpu_time(t1)
   n_inversions=0          
   
   do ijk=1, nround

!       call cpu_time(t1)
      print *,'iteración ---->',ijk

!##################################

         call     crea_inv_V(h2=h2,ratios=ratios_fact,log_det=log_det)               
                       
!c       construccion de la matriz de coeficientes.
!c       XtV-1        
         aux=0.0
         aux=mymatmul(x, V_inv*(1./vp) )
!c       XtV-1X
         xx=0.0
         xx=mymatmul(aux,transpose(x))
!        construccion de rhs
!c       xtv-1y
         XY=0.0
         XY=mymatmul(aux,y)
              
!c       muestro de gibbs

        do i=1,neq
         b(i)=xy(i)
         do j=1,neq  
          if (i.ne.j) b(i)=b(i)- (xx(i,j)*b(j))
         enddo
         b(i)=b(i)/xx(i,i)
         vare=1.0/xx(i,i)
         b(i)=gen_normal(b(i),vare)
        enddo      
        b(4)=b(4)-b(4)
        b(5)=b(5)-b(4)
        b(6)=b(6)-b(4)
        b(7)=b(7)-b(7)
        b(8)=b(8)-b(7)

!        write(876,'(5000(f8.3,x) )' )b

!##############
!c       muestreo  vp ,POSTERIOR 
!##############
        e=0       
        do j=1,ndat
         e(j)=y(j)
         do k=1,nfijos
          i=cum_nlev1(k-1)+dat(j,k)
          e(j)=e(j)-b(i)
         enddo
        enddo       
        se=0
        do i=1,ndat
         do j=1,ndat
          se=se+e(i)*V_inv(i,j)*e(j)
         enddo
        enddo
       
        vne_1=ndat-2
        se=1.0/(se)
        vp=gen_invwishart(se,vne_1)
        print '(a50,f10.4 )'," vp:",vp
!##############
!c        muestreo  h2 y ratios , metropolis-hastings
!##############

        den1= eval_conditional(h2=h2 ,ratios=ratios_fact)  
        nmet=0
        j=0

       do while(nmet.lt.imet)       
         sum_ratios=2.
         do while (sum_ratios > 1.0)
          h2a=gen_uniform()        
          ratios_fact_a=0.0
          sum_ratios=h2a
          do i=1,neff1
           if(ratios_fact(i).ne.0.)then
             ratios_fact_a(i)=gen_uniform()
             sum_ratios=sum_ratios+ratios_fact_a(i)
           endif  
          enddo
         enddo
        den2=eval_conditional(h2=h2a,ratios=ratios_fact_a)
        if (den2.gt.den1) then
          h2=h2a
          ratios_fact=ratios_fact_a
          den1=den2
          nmet=nmet+1
        else
         unif=gen_uniform()
         acon=exp(den2-den1)
         if (unif.lt.acon)then
          h2=h2a
          ratios_fact=ratios_fact_a
          den1=den2
          nmet=nmet+1
         endif
        endif
        j=j+1
       enddo
  
       print '(a50,f10.4,x,i10 )'," h2/evaluaciones:",h2,j 
        do i=1,neff1
           if(ratios_fact(i).ne.0.)then
            print '(a50,i10,x,f10.4,x,i10 )'," ratio/valor/evaluaciones:",i,ratios_fact(i),j      
           endif  
        enddo
!        write(877,'(100(f10.4,x))' )vp,h2,ratios_fact

!############## 


       if(ijk>nburnin)then        

        saved_rounds=saved_rounds+1 
  
!c       calculo de la densidad de h2 en 0
        sum_ratios=0
        do k=1,neff1
            if(ratios_fact(k).ne.0)sum_ratios=sum_ratios+ratios_fact(k)
        enddo
        sum_ratios=1-sum_ratios
        
        h2a=0.
        integral=0
        den1=eval_conditional(h2=h2a,ratios=ratios_fact)  
        integral(0)=den1
        den0=den1
        
        diff=sum_ratios/nintegral
        
        do ij=1,nintegral
          h2a=h2a+diff
          if(h2a .gt. sum_ratios)exit
          den1=eval_conditional(h2=h2a,ratios=ratios_fact)            
          integral(ij)=den1
        enddo
        c=maxval(integral)
        suma=0
        do ij=0,nintegral
         suma=suma+exp(integral(ij)-c)
        enddo
        suma=log(suma/(nintegral+1))+c
        log_pcond_h2_0=den0-suma
      

        log_pcond_h2_0_muestras(saved_rounds)=log_pcond_h2_0
         c=maxval(log_pcond_h2_0_muestras(1:saved_rounds))
         sut=0
         pmar_h2_0=0
         do ij=1,saved_rounds
           sut=sut+exp(log_pcond_h2_0_muestras(ij)-c)
         enddo
         sut=log(sut/saved_rounds)+c
         pmar_h2_0=exp(sut)
                  
         print '(a80,x,g20.7,x,g20.7)',"p(h2=0|y) / log[ p(h2=0|y,...)] para la iteracion",pmar_h2_0,log_pcond_h2_0
        
!c       calculo de la densidad de ratios en 0        
        do k=1,neff1
         if(ratios_fact(k).ne.0)then
         
           sum_ratios=h2
           
           do i=1,neff1
            if(ratios_fact(i).ne.0 .and. i.ne.k)sum_ratios=sum_ratios+ratios_fact(i)
           enddo
           sum_ratios=1-sum_ratios  
           ratios_fact_a=ratios_fact
           ratios_fact_a(k)=0.
             
            integral=0
            den1=eval_conditional(h2=h2,ratios=ratios_fact_a)
            integral(0)=den1
            den0=den1

            diff=sum_ratios/nintegral

            do ij=1,nintegral
             ratios_fact_a(k)=ratios_fact_a(k)+diff
             if(ratios_fact_a(k) .gt. sum_ratios)exit
             den1=eval_conditional(h2=h2,ratios=ratios_fact_a)
             integral(ij)=den1
            enddo
         
            c=maxval(integral)
            suma=0
            do ij=0,nintegral
             suma=suma+exp(integral(ij)-c)
            enddo
            suma=log(suma/(nintegral+1))+c
            log_pcond_ratios_0(k)=den0-suma


            log_pcond_ratios_0_muestras(k,saved_rounds)=log_pcond_ratios_0(k)
            c=maxval(log_pcond_ratios_0_muestras(k,1:saved_rounds))
            sut=0
            pmar_ratios_0(k)=0
            do ij=1,saved_rounds
             sut=sut+exp(log_pcond_ratios_0_muestras(k,ij)-c)
            enddo
            sut=log(sut/saved_rounds)+c
            pmar_ratios_0(k)=exp(sut)
            
            print '(a80,x,i,x,g20.7,x,g20.7)',"ratio /p(ratio=0|y) / log[ p(ratio=0|y,...)] para la iteracion",k,pmar_ratios_0(k),log_pcond_ratios_0(k)          
         endif
        enddo
                 
        write(68,'(i10,40(g20.7,x))')ijk,vp,h2,ratios_fact,log_pcond_h2_0,log_pcond_ratios_0, pmar_h2_0, pmar_ratios_0      
               
       endif
              
        enddo

call cpu_time(t2)
print *,'tiempo TOTAL ',(t2-t1), ' segundos'
print *,'# de inversiones',  n_inversions

        
!        close(455)
!        close(444)
!        close(456)
!        close(457)
!        close(466)
        
        contains

!----------------------------------------------------------------
   subroutine crea_inv_V(h2,ratios,log_det)
   integer :: i,j
   real *8 :: sum_ratios,log_det,h2,ratios(:),t1,t2
   
   V_inv=ZAZt*h2
   sum_ratios=h2
   do i=1,neff1
    if(ratios(i).ne.0)then
     V_inv=V_inv+ZZt(i,:,:)*ratios(i)
     sum_ratios=sum_ratios+ratios(i)
    endif 
   enddo   
   
   do i=1,ndat
     V_inv(i,i)= V_inv(i,i)+(1.0-sum_ratios)
   enddo
           
!   call cpu_time(t1)        
    call BlasInverseSym_logdet(V_inv,ndat,log_det)
!   call cpu_time(t2)
   n_inversions=n_inversions+1
!   print *,'tiempo INVERSION ',n_inversions,(t2-t1), ' segundos'   
   end subroutine

!----------------------------------------------------------------
  function eval_conditional(h2,ratios) result(densidad)
  real *8 :: h2,ratios(:)
  real *8 :: densidad,log_det,densidad1,a,b,log_prior,b1,b2
  integer :: i,j
 
        densidad1=0
        densidad =0         
        call crea_inv_V(h2,ratios,log_det)
        
        do i=1,ndat
         do j=1,ndat
           densidad1=densidad1 +e(i) * V_inv(i,j) * e(j) 
         enddo  
        enddo
        
        
       ! a=1.0
       ! b=2.0
       ! b1=(a-1.0)*log(h2)
       ! b2=(b-1.0)*log(1.0-h2)
       ! log_prior=b1+b2 ! lo explicito para a=1 y b=2
       log_prior=log(1.0-h2)
        do i=1,neff1
         if(ratios(i).ne.0)then
          !b1=(a-1.0)*log(ratios(i))
          !b2=(b-1.0)*log(1.0-ratios(i))
          !log_prior=log_prior+( b1 + b2 ) 
          log_prior=log_prior+log(1.0-ratios(i))
        endif
       enddo        
        
      densidad=-0.5*log_det-(0.5/vp)*densidad1+log_prior
  end function
!----------------------------------------------------------------
   real function inbree(an)
   implicit none
   integer::an,idam,isir
   
   isir=ped1(an,2)
   idam=ped1(an,3)
   if(idam==0.or.isir==0) then
      inbree=0
   else
      inbree=.5*parentesco(idam,isir) 
   endif
   end function
!----------------------------------------------------------------
   recursive real function parentesco(an1,an2) result(r)
   implicit none
   integer::an1,an2
   if(an1.ne.0 .and. an2.ne.0)then
     if(A_mat(an1,an2) .ne. -99)then
      r=A_mat(an1,an2)
      return
     endif
   endif
   if(an1==0.or.an2==0) then
      r=0
   elseif(an1==an2) then
      r=1.0+inbree(an1)
   else
      if (an1 < an2)then 
           r=.5*(parentesco(an1,ped1(an2,2))+parentesco(an1,ped1(an2,3)))
      else
           r=.5*(parentesco(an2,ped1(an1,2))+parentesco(an2,ped1(an1,3)))
      endif 
   endif
   if(an1.ne.0 .and. an2.ne.0)then
    A_mat(an1,an2)=r
    A_mat(an2,an1)=r
   endif
   end function


       END PROGRAM    


