PROGRAM average
!1 -46275.219200 -46275.2        24370.950472 24371
!2 -46273.688606 -46274.5        24349.836356 24360.4
!3 -46270.864119 -46273.3        24330.359232 24350.4
!4 -46267.069370 -46271.7        24312.507618 24340.9
 IMPLICIT NONE
 REAL(8)               ::   datos(1:4) = 0.0
 REAL(8)               ::   average_v  = 0.0
 REAL(8)               ::   sd_v       = 0.0
 REAL(8)               ::   partition  = 0.0
 REAL(8)               ::   ener0      = 99999999.9999999
 REAL(8)               ::   beta = 300.0
 REAL(8),parameter     ::   kb = 0.000086173324 ! eV K^-1
 INTEGER               ::   i,j,k,ierr
 CHARACTER (LEN=80)    ::   line
 beta  = 1.0/(kb*beta)
 k = 0
 do0: DO
      READ (5,'(A)',iostat=ierr) line
      IF (ierr /= 0) EXIT do0
      k = k + 1
      READ(line,*)i,( datos(j) ,j=1,4)
      IF(datos(2) <= ener0) ener0 = datos(2)
 END DO do0
 REWIND(5)
 do1: DO
      READ (5,'(A)',iostat=ierr) line
      IF (ierr /= 0) EXIT do1
      READ(line,*)i,( datos(j) ,j=1,4)
      partition = partition + exp(-beta*(datos(2)-ener0))
      average_v = average_v + datos(4)*exp(-beta*(datos(2)-ener0))
      sd_v = sd_v + datos(4)*datos(4)*exp(-beta*(datos(2)-ener0))
 END DO do1
 average_v = average_v / partition
 sd_v = sqrt((sd_v / partition) - average_v*average_v)
 WRITE(6,*)average_v,sd_v,(average_v/8.0)**(1.0/3.0),(sd_v/8.0)**(1/3.0)
END PROGRAM average
