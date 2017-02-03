SUBROUTINE PrtConn(IOut,N,LAnMMP,IAnMMP)

 INTEGER :: IOut, I, J, N, LAnMMP
 INTEGER :: IAnMMP(N,LAnMMP)

 2000 FORMAT('   ---- Gaussian-style connectivity ----------------------')
 2100 FORMAT(2x,9(I5,1x))

  WRITE(IOut,2000)
  DO I=1,N
    WRITE(IOut,2100) (IAnMMP(I,J),J=1,LAnMMP)
  ENDDO

END SUBROUTINE
