  
PROGRAM Distancias
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  !Criação da subrotina do cálculo de distância por Mahalanobis                 !
  !Orientador: Cosme Ferreira da Ponte Neto                                     !
  !Aluno: Victor Ribeiro Carreira                                               !
  !Este programa visa comparar as distâncias Euclidiana e  Mahalananobis        !
  !Cálculo de distâncias                                                        !
  !Subrotina Mahalanobis                                                        !
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
 IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! DECLARAÇÃO DAS VARIÁVEIS GLOBAIS !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
INTEGER, PARAMETER::SP = SELECTED_REAL_KIND(p=4, r=4)
INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(p=8, r=8)
    
REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:):: tr1, tr2, s, menor, g, rg, trT, & 
cov, md, mdT, alfa, d2, g1, g2

INTEGER(KIND=DP):: i, j, k, ie, np, ij, nt1, nt2, np2

REAL(KIND=DP)::a1, ctr1x, ctr2x, ctr1y, ctr2y, d_p1_c1, d_p1_c2, d_p2_c1, d_p2_c2, &
d_p3_c1, d_p3_c2, d_p4_c1, d_p4_c2, dist

CHARACTER(LEN=11):: l(4), lito(1000)
CHARACTER(LEN=80):: cab 


ALLOCATE(tr1(10,2), tr2(5,2), s(100,100), g(1000,20), rg(1000,2), trT(8,1000),&
cov(2,2), md(2,1), mdT(1,2), alfa(1,2), d2(1,1), g1(10,2), g2(1,2))

!implicit real*8(a-h,o-z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!ALOCANDO OS ARQUIVOS DE ENTRADA E DE SAÍDA !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


OPEN(1,file='grupo1')						
OPEN(2,file='entrada2.txt')
OPEN(3,file='saida2.txt')
OPEN(4,file='grupo2')	
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!LENDO OS ARQUIVOS DE ENTRADA!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



READ(2,15) cab    ! lê cabeçalho e armazena em cab
READ(2,15) cab    ! lê linha em branco abaixo do cabeçalho e armazena em cab

 ie=1
DO WHILE (.TRUE.)
  READ(2,*,END=8) rg(ie,1),g(ie,1),g(ie,2)
  ie=ie+1
END DO
8 CONTINUE
CLOSE(2)

np=ie-1  ! numero de pontos da entrada
WRITE(6,*) "n de dados de entrada",np
WRITE(6,*)'coordenada dos pontos'

DO i=1,np
  WRITE(6,*) g(i,1),'    ',g(i,2)
END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!LENDO OS ARQUIVOS DOS GRUPOS!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ij=1
DO WHILE (.TRUE.)
  READ(1,*,end=9) tr1(ij,1),tr1(ij,2),a1
  ij=ij+1
END DO
9 CONTINUE
CLOSE(1)

nt1=ij-1 !número de dados do grupo1
WRITE(6,*) "n de dados do grupo 1",nt1

ij=1
DO WHILE (.TRUE.)
  READ(4,*,END=10) tr2(ij,1),tr2(ij,2),a1
  ij=ij+1
END DO
10 CONTINUE
CLOSE(4)
 
 nt2=ij-1 !número de dados do grupo1

WRITE(6,*) "n de dados do grupo 2",nt2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! DISTÂNCIA DE EUCLIDIANA !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Calculo do centroide 1
  ctr1x=0d0
  ctr1y=0d0

 DO i=1,nt1
   ctr1x=ctr1x+tr1(i,1)
   ctr1y=ctr1y+tr1(i,2)
 END DO 

 ctr1x=ctr1x/dfloat(nt1)
 ctr1y=ctr1y/dfloat(nt1)

!Calculo do centroide 2
 ctr2x=0d0
 ctr2y=0d0

DO i=1,nt2
  ctr2x=ctr2x+tr2(i,1)
  ctr2y=ctr2y+tr2(i,2)
END DO 

 ctr2x=ctr2x/dfloat(nt2)
 ctr2y=ctr2y/dfloat(nt2)

WRITE(6,*)'centroide 1',ctr1x,'   ',ctr1y
WRITE(6,*)'centroide 2',ctr2x,'   ',ctr2y


! distancia do ponto 1 ao centroide da classe 1
	d_p1_c1=dsqrt((g(1,1)-ctr1x)**2+(g(1,2)-ctr1y)**2)
! distancia do ponto 1 ao centroide da classe 2
	d_p1_c2=dsqrt((g(1,1)-ctr2x)**2+(g(1,2)-ctr2y)**2)

! distancia do ponto 2 ao centroide da classe 1
	d_p2_c1=dsqrt((g(2,1)-ctr1x)**2+(g(2,2)-ctr1y)**2)
! distancia do ponto 2 ao centroide da classe 2
	d_p2_c2=dsqrt((g(2,1)-ctr2x)**2+(g(2,2)-ctr2y)**2)

! distancia do ponto 3 ao centroide da classe 1
	d_p3_c1=dsqrt((g(3,1)-ctr1x)**2+(g(3,2)-ctr1y)**2)
! distancia do ponto 3 ao centroide da classe 2
	d_p3_c2=dsqrt((g(3,1)-ctr2x)**2+(g(3,2)-ctr2y)**2)

! distancia do ponto 4 ao centroide da classe 1
	d_p4_c1=dsqrt((g(4,1)-ctr1x)**2+(g(4,2)-ctr1y)**2)
! distancia do ponto 4 ao centroide da classe 2
	d_p4_c2=dsqrt((g(4,1)-ctr2x)**2+(g(4,2)-ctr2y)**2)


WRITE(6,*) '======================='
WRITE(6,*) '=======EUCLIDES========'
WRITE(6,*) 'd_p1_c1=',d_p1_c1
WRITE(6,*) 'd_p2_c1=',d_p2_c1
WRITE(6,*) 'd_p3_c1=',d_p3_c1
WRITE(6,*) 'd_p4_c1=',d_p4_c1
WRITE(6,*) '==============='
WRITE(6,*) 'd_p1_c2=',d_p1_c2
WRITE(6,*) 'd_p2_c2=',d_p2_c2
WRITE(6,*) 'd_p3_c2=',d_p3_c2
WRITE(6,*) 'd_p4_c2=',d_p4_c2
WRITE(6,*) '==============='

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! DISTÂNCIA DE MAHALANOBIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PRINT*, 'antes',tr1(1,1),tr1(3,1)

!distancias entre os dois grupos
CALL maha(tr1,nt1,tr2,nt2,2,dist)
WRITE(6,*)'distancia entre os grupos=',dist
PRINT*, 'depois1',tr1(1,1),tr1(3,1)



!distancias entre os dois grupos
PAUSE
CALL maha(tr1,nt1,tr2,nt2,2,dist)
WRITE(6,*)'distancia entre os grupos=',dist

PRINT*, 'depois1',tr1(1,1),tr1(3,1)

WRITE(6,*) '======MAHA========='

!	classe 1
DO i=1,nt1
  g1(i,1)=tr1(i,1)
  g1(i,2)=tr1(i,2)
END DO

  g2(1,1)=g(1,1)
  g2(1,2)=g(1,2)

CALL maha(tr1,nt1,g2,1,2,dist)
WRITE(6,*)'d_p1_c1 maha=',dist

 g2(1,1)=g(2,1) 
 g2(1,2)=g(2,2)


CALL maha(tr1,nt1,g2,1,2,dist)
WRITE(6,*)'d_p2_c1 maha=',dist

 g2(1,1)=g(3,1)
 g2(1,2)=g(3,2)

CALL maha(tr1,nt1,g2,1,2,dist)
WRITE(6,*)'d_p3_c1 maha=',dist
 
  g2(1,1)=g(4,1)
  g2(1,2)=g(4,2)

CALL maha(tr1,nt1,g2,1,2,dist)
WRITE(6,*)'d_p4_c1 maha=',dist

!classe 2

WRITE(6,*) '==============='

DO i=1,nt2
  g1(i,1)=tr2(i,1)
  g1(i,2)=tr2(i,2)
END DO 

  g2(1,1)=g(1,1)
  g2(1,2)=g(1,2)

CALL maha(tr2,nt2,g2,1,2,dist)
WRITE(6,*)'d_p1_c2 maha=',dist

  g2(1,1)=g(2,1)
  g2(1,2)=g(2,2)

CALL maha(tr2,nt2,g2,1,2,dist)
WRITE(6,*)'d_p2_c2 maha=',dist

 g2(1,1)=g(3,1)
 g2(1,2)=g(3,2)

CALL maha(tr2,nt2,g2,1,2,dist)
WRITE(6,*)'d_p3_c2 maha=',dist

 g2(1,1)=g(4,1)
 g2(1,2)=g(4,2)

CALL maha(tr2,nt2,g2,1,2,dist)
WRITE(6,*)'d_p4_c2 maha=',dist

!distancias entre os dois grupos

CALL maha(tr1,nt1,tr2,nt2,2,dist)
WRITE(6,*)'distancia entre os grupos=',dist

11	FORMAT(4(ES12.4E3,2x))
12	FORMAT(I3,2x,3(f6.2,2x))
13	FORMAT(4(ES12.4E3,2x))
14	FORMAT(4(ES9.2E2,2x))
15	FORMAT(A71)
16	FORMAT(A11,8(ES9.2E3))
17	FORMAT(A30,2x,ES12.4E3)
18	FORMAT(1(f6.2,2x),2x,A11,2x,ES12.4E3)
19	FORMAT(A8,A2,f4.1,A2,f4.1)




PRINT*,' ************ FIM *************'
PRINT*,''

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SUBROTINAS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS

SUBROUTINE media_desvio(x,n,media,desvio_padrao)
  IMPLICIT NONE
  INTEGER(KIND=DP), INTENT(IN)::n
  REAL(KIND=DP), INTENT(IN), ALLOCATABLE, DIMENSION(:):: x
  REAL(KIND=DP), INTENT(OUT):: media, desvio_padrao
  INTEGER(KIND=DP):: i,j
  REAL(KIND=DP):: soma, soma2

n=30! voce escolhe

 ALLOCATE(x(n))
 !implicit real*8(a-h,o-z)
 !real*8 x(n),media,desvio_padrao
 soma = 0d0
 soma2 = 0d0

 DO j=1,n	
  soma = soma + x(j)
  soma2 = soma2 + x(j)**2
 END DO 

 media = soma/n
 desvio_padrao=dsqrt((soma2-soma**2/n)/(n-1))

END SUBROUTINE media_desvio


!#####################################################
	
SUBROUTINE soma(n,x,w)
 IMPLICIT NONE
 INTEGER(KIND=DP):: i,e
 REAL(KIND=DP), ALLOCATABLE, INTENT(IN)::n=30! pode mudar
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:), INTENT(IN)::x	
 REAL(KIND=DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT)::w   

ALLOCATE(x(e))
    !integer:: e,i
	!real*8 :: x(e),w
 
 w=0
 
 DO i=1,e
  w= w + x(i)
 END DO 

END SUBROUTINE soma

!#####################################################

SUBROUTINE maha(g11,np1,g22,np2,ndim,dist)      
	
!! subrotina que calcula a distância de mahalanobis entre
!! dois agrupamentos de elementos com dimensão ndim 	

 IMPLICIT NONE
 REAL(KIND=DP),DIMENSION(:,:), ALLOCATABLE:: g11, g22
 REAL(KIND=DP):: dist
 INTEGER(KIND=DP):: np1, np2, ndim
 REAL(KIND=DP),DIMENSION(:,:), ALLOCATABLE:: g1, g2, g1T, g2T, cov1, &
 cov2, covag, md, mdT, alfa, d2
 REAL(KIND=DP),DIMENSION(:), ALLOCATABLE::	soma, xm1, m2
 INTEGER(KIND=DP):: i, j, k  

 ALLOCATE(g1(np1,ndim),g2(np2,ndim),&
 g1T(ndim,np1),g2T(ndim,np2),cov1(ndim,ndim),cov2(ndim,ndim),&
 covag(ndim,ndim),soma(ndim),xm1(ndim),m2(ndim),g22(np2,ndim),&
 md(ndim,1),mdT(1,ndim),alfa(1,ndim),d2(1,1),g11(np1,ndim)) 

  
!implicit real*8(a-h,o-z)

!	real*8,intent(in)::g1(np1,ndim),g2(np2,ndim)
!	real*8,intent(out)::dist
!	integer,intent(in)::np1,np2,ndim 


!REAL*8 g1(np1,ndim),g2(np2,ndim),g1T(ndim,np1),g2T(ndim,np2),&
!cov1(ndim,ndim),cov2(ndim,ndim),covag(ndim,ndim),soma(ndim),xm1(ndim),&
!m2(ndim),g22(np2,ndim),md(ndim,1),mdT(1,ndim),alfa(1,ndim),d2(1,1),g11(np1,ndim)


	g1=g11
	g2=g22

!	grupo 1	

DO j=1,ndim
  soma(j)=0d0
   DO i=1,np1
    soma(j)=soma(j)+g1(i,j)
   END DO 
END DO

DO i=1,ndim
  xm1(i)=soma(i)/dfloat(np1)
END DO	

!	grupo 2	

DO j=1,ndim
  soma(j)=0d0
  DO i=1,np2
   soma(j)=soma(j)+g2(i,j)
  END DO
END DO

DO i=1,ndim
  xm2(i)=soma(i)/DFLOAT(np2)
  !xm2(i)=soma(i)/REAL(np2)
END DO	

!vetor das diferenças - será escrito sobre a matrizes g1 e g2

DO j=1,ndim
 DO i=1,np1
  g1(i,j)=g1(i,j)-xm1(j)
  END DO 
END DO 

DO j=1,ndim
 DO i=1,np2
  g2(i,j)=g2(i,j)-xm2(j)
 END DO 
END DO 	

!     --------GRUPO 1 ---------------------
!	criando a matriz transposta g1T
!	-------------- -------------------
DO i=1,np1    !107 ! número de equações 
 DO j=1,ndim   !2
  g1T(j,i)=g1(i,j)
 END DO 
END DO 
!----------------------------------------------------
!	 - multiplicação de matrizes
!	   multiplicação de g1T por g1 

DO k=1,ndim
 DO j=1,ndim
  cov1(j,k)=0.d0
  DO i=1,np1	
   cov1(j,k)=cov1(j,k)+g1T(j,i)*g1(i,k)
  END DO
 END DO
END DO 

DO i=1,ndim
 DO j=1,ndim
  cov1(i,j)=cov1(i,j)/dfloat(np1)
 END DO
END DO

!	write(6,*) '======covari�ncia 1 ======'
!	write(6,*) cov1(1,1),cov1(1,2)
!	write(6,*) cov1(2,1),cov1(2,2)

!     --------GRUPO 2 ---------------------
!	criando a matriz transposta g2T

DO i=1,np2
 DO j=1,ndim
  g2T(j,i)=g2(i,j)
 END DO 
END DO 

!---------------------------------------------------
!	 - multiplicação de matrizes
!	   multiplicação de g2T por g2 

DO k=1,ndim
 DO j=1,ndim
  cov2(j,k)=0.d0
   DO i=1,np2
    cov2(j,k)=cov2(j,k)+g2T(j,i)*g2(i,k)
   END DO
  END DO 
END DO 

DO i=1,ndim
  DO j=1,ndim
   cov2(i,j)=cov2(i,j)/dfloat(np2)
  END DO 
END DO 

!	write(6,*) '======covariância 2 ======'
!	write(6,*) cov2(1,1),cov2(1,2)
!	write(6,*) cov2(2,1),cov2(2,2)


!	-------- covariância agrupada------

DO i=1,ndim
  DO j=1,ndim
   covag(i,j)=dfloat(np1)*cov1(i,j)/(dfloat(np1+np2))+dfloat(np2)*cov2(i,j)/(dfloat(np1+np2))
  END DO
END DO 

!	write(6,*) '======covariância agrupada ======'
!	write(6,*) covag(1,1),covag(1,2)
!	write(6,*) covag(2,1),covag(2,2)	

!	inversão da matriz covag - usando subrotina

 CALL INVERT(covag,ndim)

!	write(6,*) '====== inv covariância agrupada ======'
!	write(6,*) covag(1,1),covag(1,2)
!	write(6,*) covag(2,1),covag(2,2)

!	diferenicas médias

 DO i=1,ndim
  md(i,1)=xm1(i)-xm2(i)
 END DO 

!	write(6,*) '====== diferencias medias ======'
!	write(6,*) md(1,1)
!	write(6,*) md(2,1)

!	criando a matriz transposta mdT
!	---------------------------
 DO i=1,ndim
  DO j=1,1
   mdT(j,i)=md(i,j)
  END DO 
END DO 

!----------------------------------------------------
!	multiplicaçãoo de mdT por cov^-1 
!	 - multiplicação de matrizes

 DO k=1,ndim	
  DO j=1,1	
   alfa(j,k)=0.d0
    DO i=1,ndim
	 alfa(j,k)=alfa(j,k)+mdT(j,i)*covag(i,k)
    END DO 
  END DO 
END DO 

!----------------------------------------------------
!	multiplica��o de alfa por md 
!	 - multiplica��o de matrizes

DO k=1,1
 DO j=1,1
  d2(j,k)=0.d0
   DO i=1,ndim  !2	!
    d2(j,k)=d2(j,k)+alfa(j,i)*md(i,k)
   END DO 
  END DO 
END DO 

	dist=dsqrt(d2(1,1))

      return
END SUBROUTINE maha

!#####################################################


SUBROUTINE INVERT(A,i)      
 real*8 A(i,i),B(i)
 integer i,im,j,k,l
   IM=I-1
   DO 5 K=1,I
    DO 2 J=1,IM
    2 B(J)=A(1,J+1)/A(1,1)
    B(I)=1.d0/A(1,1)
     DO 4 L=1,IM
      DO 3 J=1,IM
      3 A(L,J)=A(L+1,J+1)-A(L+1,1)*B(J)
      4 A(L,I)=-A(L+1,1)*B(I)
      DO 5 J=1,I
       5 A(I,J)=B(J)
RETURN

END SUBROUTINE INVERT

END PROGRAM Distancias





