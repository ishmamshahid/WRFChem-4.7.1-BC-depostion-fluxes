MODULE module_wetdep_ls

USE module_state_description, only: p_qv,p_qc,p_so2,p_sulf,                &
                                    p_bc1,p_bc2,p_oc1,p_oc2,p_seas_1,      &
                                    p_seas_2,p_seas_3,p_seas_4,p_dms

CONTAINS

SUBROUTINE wetdep_ls(dt,var,rain,moist,rho,num_moist,                     &
                     num_chem,numgas,dz8w,vvel,chem_opt,                   &
                     bc1wetdep_ls,bc2wetdep_ls,                            &
                     ids,ide,jds,jde,kds,kde,                              &
                     ims,ime,jms,jme,kms,kme,                              &
                     its,ite,jts,jte,kts,kte)

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: num_chem,numgas,num_moist,                       &
                          chem_opt,                                         &
                          ids,ide,jds,jde,kds,kde,                         &
                          ims,ime,jms,jme,kms,kme,                         &
                          its,ite,jts,jte,kts,kte

   REAL, INTENT(IN) :: dt

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),                    &
         INTENT(IN) :: moist

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme),                              &
         INTENT(IN) :: rho,dz8w,vvel

   REAL, DIMENSION(ims:ime,kms:kme,jms:jme,1:num_chem),                   &
         INTENT(INOUT) :: var

   REAL, DIMENSION(ims:ime,jms:jme),                                      &
         INTENT(IN) :: rain

   REAL, DIMENSION(ims:ime,jms:jme),                                      &
         INTENT(INOUT) :: bc1wetdep_ls,bc2wetdep_ls

!  REAL, DIMENSION(ims:ime,jms:jme),                                      &
!        INTENT(INOUT) :: var_rmv

   REAL, DIMENSION(its:ite,jts:jte) :: var_sum,var_rmv

   REAL, DIMENSION(its:ite,kts:kte,jts:jte) :: var_rmvl

   REAL, DIMENSION(its:ite,jts:jte) :: frc,var_sum_clw,rain_clw

   REAL :: dvar,factor,clsum,alpha,rho_water

   INTEGER :: nv,i,j,k,km,kb,kbeg

   rho_water = 1000.

!  write(0,*) 'in wetdepls, numgas,num_chem = ',numgas,num_chem,chem_opt
!  nv=p_bc1

   DO nv=1,num_chem

      IF (nv.le.numgas .and. nv.ne.p_sulf) CYCLE

      alpha = .5    ! scavenging factor

      IF (chem_opt >= 300) THEN

         ! Hydrophobic BC1 and OC1 are not wet scavenged
         IF (nv.eq.p_bc1 .or. nv.eq.p_oc1 .or. nv.eq.p_dms) CYCLE

         ! Hydrophilic BC2 and OC2 scavenging efficiency
         IF (nv.eq.p_bc2 .or. nv.eq.p_oc2) alpha = 0.8

      ENDIF

      IF (nv.eq.p_sulf   .or.                                               &
          nv.eq.p_seas_1 .or.                                               &
          nv.eq.p_seas_2 .or.                                               &
          nv.eq.p_seas_3 .or.                                               &
          nv.eq.p_seas_4) alpha = 1.

      DO i=its,ite
         DO j=jts,jte

            var_sum_clw(i,j) = 0.
            var_sum(i,j)     = 0.
            var_rmvl(i,:,j)  = 0.
            frc(i,j)         = 0.
            rain_clw(i,j)    = 0.

            IF (rain(i,j).gt.1.e-10) THEN

!              write(0,*)i,j,rain(i,j)

!
! Convert rain back to rate
!
               rain_clw(i,j) = rain(i,j)/dt

!
! Total cloud water
!
               DO k=1,kte-1

                  dvar = MAX(0.,moist(i,k,j,p_qc) *                         &
                                rho(i,k,j)          *                         &
                                vvel(i,k,j)         *                         &
                                dz8w(i,k,j))

                  var_sum_clw(i,j) = var_sum_clw(i,j) + dvar

                  var_sum(i,j) = var_sum(i,j) +                             &
                                 var(i,k,j,nv)*rho(i,k,j)

               ENDDO

               IF (var_sum(i,j).gt.1.e-10 .and.                            &
                   var_sum_clw(i,j).gt.1.e-10) THEN

!
! Assuming frc is constant, it is the conversion factor,
! similar to the convective parameterization.
!
                  frc(i,j) = rain_clw(i,j)/var_sum_clw(i,j)

!                 write(0,*)'frc ',frc(i,j),var_sum_clw(i,j),var_sum(i,j)

                  frc(i,j) = MAX(1.e-6,MIN(frc(i,j),.005))

               ENDIF

            ENDIF

         ENDDO
      ENDDO

!
! Wet scavenging
!
      DO i=its,ite
         DO j=jts,jte

            IF (rain(i,j).gt.1.e-10      .and.                             &
                var_sum(i,j).gt.1.e-10   .and.                             &
                var_sum_clw(i,j).gt.1.e-10) THEN

               DO k=kts,kte-2

                  IF (var(i,k,j,nv).gt.1.e-16 .and.                        &
                      moist(i,k,j,p_qc).gt.0.) THEN

                     factor = MAX(0.,frc(i,j)   *                           &
                                      rho(i,k,j) *                           &
                                      dz8w(i,k,j)*                           &
                                      vvel(i,k,j))

!                    write(0,*)'var before ',k,km,var(i,k,j,nv),factor

                     dvar = alpha*factor/(1.+factor)*var(i,k,j,nv)

!
! Accumulated grid-scale BC wet deposition
!
! dvar units:
!     ug BC / kg dry air
!
! dvar * rho * dz8w:
!     ug BC / m2
!
! Multiplication by 1.e-9 converts:
!     ug/m2 -> kg/m2
!
! BC1 is skipped above for chem_opt >= 300, so BC1 wet deposition
! remains zero. BC2 is hydrophilic and is wet scavenged.
!
                     IF (chem_opt >= 300 .and. nv.eq.p_bc2) THEN

                        bc2wetdep_ls(i,j) = bc2wetdep_ls(i,j) +             &
                                           dvar          *                   &
                                           rho(i,k,j)    *                   &
                                           dz8w(i,k,j)   *                   &
                                           1.e-9

                     ENDIF

                     var(i,k,j,nv) = MAX(1.e-16,                            &
                                         var(i,k,j,nv)-dvar)

!                    write(0,*)'var after ',km,var(i,k,j,nv),dvar

                  ENDIF

               ENDDO

!              var_rmv(i,j)=var_rmv(i,j)+var_rmvl(i,j)

            ENDIF

         ENDDO
      ENDDO

   ENDDO

END SUBROUTINE wetdep_ls

END MODULE module_wetdep_ls

