PROGRAM teste_mahalanobis

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  !Criação da subrotina do cálculo de distância euclideana                     !
  !Orientador: Cosme Ferreira da Ponte Neto                                    !
  !Aluno: Victor Ribeiro Carreira                                              !
  !Este programa validar a subrotina mahalanobeana                             !
  !Categoria: classificador                                                    !
  !Subrotina mahalanobeana                                                     !
  !Para usar compilação com flags utilize:                                     !
  !gfortran -fbounds-check -fbacktrace -Wall -Wextra -pedantic                 !
  !"pasta/subpasta/nomedopragrama.f95" -o nomedoexecutável                     !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

IMPLICIT NONE

INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

INTEGER(KIND=SP):: i,j,k,x,y,z
REAL(KIND=DP)::maha
REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::lito1, lito2, ponto

21 FORMAT(15(F4.2,2x))

ALLOCATE(lito1(5,15),lito2(1,15))

lito1=0d0
lito2=0d0

DO i=1,5
 DO j=1,15
   IF(i>j)THEN
     lito1(i,j)=3.42
   ELSE IF(i<j)THEN
     lito1(i,j)=1.72
   ELSE
     lito1(i,j)=8.98
   ENDIF
   !lito1(i,j)=2.0
 END DO
END DO

DO i=1,1
 DO j=1,15
   IF(i>j)THEN
     lito2(i,j)=7.72
   ELSE IF(i<j)THEN
     lito2(i,j)=6.98
   ELSE
     lito2(i,j)=0.25
   ENDIF
   !lito2(i,j)=3.0
 END DO
END DO

PRINT*,'---------------------'
PRINT*,'Matriz lito1'
PRINT*,'---------------------'
WRITE(*,21) lito1
PRINT*,'---------------------'
PRINT*,'Matriz lito2'
PRINT*,'---------------------'
WRITE(*,21) lito2
PRINT*,'---------------------'


!Criando um Esferóide 

!Interface 2
DO j=1,Nx
  IF (j.ge.245 .and. j.le.345) THEN
    zc=150+NINT(SQRT(2500.0-(j-295.0)**2))
  ELSE
    zc=150
  ENDIF
  DO i=zc,Nz
    vel(i,j)=2100.0
  ENDDO
ENDDO


!*******************************************************************************************!


CONTAINS


SUBROUTINE mahalanobeana(g11,np1,g22,np2,ndim,dist)

!  	subrotina que calcula a distância de mahalanobis entre
!  	dois agrupamentos de elementos com dimensão ndim

   IMPLICIT NONE
    INTEGER, PARAMETER::SP = SELECTED_INT_KIND(r=8)
    INTEGER, PARAMETER::DP = SELECTED_REAL_KIND(12,100)

    INTEGER(KIND=DP), INTENT(IN):: np1
    INTEGER(KIND=SP), INTENT(IN):: np2, ndim
    REAL(KIND=DP),INTENT(OUT):: dist

    INTEGER(KIND=SP):: i,j,k
    REAL(KIND=DP),ALLOCATABLE,DIMENSION(:)::soma, xm1, xm2
    REAL(KIND=DP),ALLOCATABLE, DIMENSION(:,:):: g1, g2, g1T, g2T, cov1, cov2, &
    covag, g11, g22, md, mdT, alfa, d2

    ALLOCATE(soma(ndim),xm1(ndim),xm2(ndim))

    ALLOCATE(g1(np1,ndim),g2(np2,ndim),g1T(ndim,np1),g2T(ndim,np2),&
    cov1(ndim,ndim),cov2(ndim,ndim),covag(ndim,ndim),md(ndim,1),&
    mdT(1,ndim),alfa(1,ndim),d2(1,1))

    g1=g11
    g2=g22

!  	grupo 1

  DO j=1,ndim
    soma(j)=0d0
    DO i=1,np1
      soma(j)=soma(j)+g1(i,j)
    END DO
  END DO

  DO i=1,ndim
    xm1(i)=soma(i)/dfloat(np1)
  END DO

!  	grupo 2

  DO j=1,ndim
    soma(j)=0d0
    DO i=1,np2
      soma(j)=soma(j)+g2(i,j)
    END DO
  END DO

  DO i=1,ndim
    xm2(i)=soma(i)/dfloat(np2)
  END DO

!  	vetor das diferenças - será escrito sobre a matrizes g1 e g2


  DO j=1,ndim
    DO i=1,np1
      g1(i,j)=g1(i,j)-xm1(j)
    END DO
  END DO

  DO  j=1,ndim
    DO i=1,np2
      g2(i,j)=g2(i,j)-xm2(j)
    END DO
  END DO

!      --------GRUPO 1 ---------------------
!  	criando a matriz transposta g1T
!  	-------------- -------------------
  DO i=1,np1    !107 ! número de equações
    DO j=1,ndim   !2
      g1T(j,i)=g1(i,j)
    END DO
  END DO
!  ----------------------------------------------------
!  	 - multiplicação de matrizes
!  	   multiplicação de g1T por g1

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

!  	write(6,*) '======covariância 1 ======'
!  	write(6,*) cov1(1,1),cov1(1,2)
!  	write(6,*) cov1(2,1),cov1(2,2)

!      --------GRUPO 2 ---------------------
!  	criando a matriz transposta g2T

  DO i=1,np2
    DO j=1,ndim
      g2T(j,i)=g2(i,j)
    END DO
  END DO

!  ---------------------------------------------------
!  	 - multiplicação de matrizes
!  	   multiplicação de g2T por g2

   DO k=1,ndim
     DO j=1,ndim
       cov2(j,k)=0.d0
       DO i=1,np2
         cov2(j,k)=cov2(j,k)+g2T(j,i)*g2(i,k)
       END DO
     END DO
   END DO

   DO  i=1,ndim
     DO j=1,ndim
       cov2(i,j)=cov2(i,j)/dfloat(np2)
     END DO
   END DO

   ! WRITE(6,*) '======covariância 2 ======'
   ! WRITE(6,*) cov2(1,1),cov2(1,2)
   ! WRITE(6,*) cov2(2,1),cov2(2,2)


!  	-------- covariância agrupada------

   DO i=1,ndim
     DO j=1,ndim
       covag(i,j)=dfloat(np1)*cov1(i,j)/(dfloat(np1+np2))+ &
       dfloat(np2)*cov2(i,j)/(dfloat(np1+np2))
     END DO
   END DO

    !WRITE(6,*) '======covariância agrupada ======'
    !WRITE(6,*) covag(1,1),covag(1,2)
    !WRITE(6,*) covag(2,1),covag(2,2)

!  	inversao da matriz covag - usando subrotina

   CALL INVERT(covag,ndim)

    !WRITE(6,*) '====== inv covariância agrupada ======'
    !WRITE(6,*) covag(1,1),covag(1,2)
    !WRITE(6,*) covag(2,1),covag(2,2)

!  	diferenicas médias

   DO i=1,ndim
     md(i,1)=xm1(i)-xm2(i)
   END DO

!  	write(6,*) '====== diferencias medias ======'
!  	write(6,*) md(1,1)
!  	write(6,*) md(2,1)

!  	criando a matriz transposta mdT
!  	---------------------------

  DO i=1,ndim
    DO j=1,1
      mdT(j,i)=md(i,j)
    END DO
  END DO

!  ----------------------------------------------------
!  	multiplicação de mdT por cov^-1
!  	 - multiplicação de matrizes


  DO k=1,ndim
    DO j=1,1
      alfa(j,k)=0.d0
      DO i=1,ndim
        alfa(j,k)=alfa(j,k)+mdT(j,i)*covag(i,k)
      END DO
    END DO
  END DO

!  ----------------------------------------------------
!  	multiplicação de alfa por md
!  	 - multiplicação de matrizes

  DO k=1,1
    DO j=1,1
      d2(j,k)=0.d0
      DO i=1,ndim  !2	!
        d2(j,k)=d2(j,k)+alfa(j,i)*md(i,k)
      END DO
    END DO
  END DO

  dist=dsqrt(d2(1,1))


END SUBROUTINE mahalanobeana


!--------------------------------------------------------------------------


   SUBROUTINE INVERT(A,i)
      integer i,im,j,k,l
      real*8 A(i,i),B(i)

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

   END SUBROUTINE INVERT


END PROGRAM teste_mahalanobis 