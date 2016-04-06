Subroutine calc_velocity_serial_2d
    use vpm_vars         
    use pmeshpar
    use parvar
    use pmgrid
    use openmpth
    Implicit None
    double precision ::  dphidx, dphidy,dphidz,dpsizdx, dpsizdy
    double precision ::  wdudx, wdvdy, wdwdz, velxp, velyp, velzp, velxm, velym, velzm
    double precision ::  upi,umi,upj,umj,upk,umk
    double precision ::  vpi,vmi,vpj,vmj,vpk,vmk
    double precision ::  wpi,wmi,wpj,wmj,wpk,wmk
    integer          :: i, j, k

    DXpm2 = 2 * DXpm
    DYpm2 = 2 * DYpm
    DZpm2 = 2 * DZpm
    k = 1
    !$omp do
    do j = NYs_bl(1) + 1, NYf_bl(1) - 1
       do i = NXs_bl(1) + 1, NXf_bl(1) - 1



            dpsizdx  = (SOL_pm(1,i + 1, j, k)   - SOL_pm(1,i-1,j,k)) / (2d0*DXpm)
            dpsizdy  = (SOL_pm(1,i , j+1, k)    - SOL_pm(1,i,j-1,k)) / (2d0*DYpm)


            velvrx_pm(i, j, k)  =  + dpsizdy 
            velvry_pm(i, j, k)  =  - dpsizdx 


        enddo
    enddo

End Subroutine calc_velocity_serial_2d
