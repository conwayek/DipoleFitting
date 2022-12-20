subroutine poten(xp,muX, r1, r2, cos_theta)

implicit none
integer , parameter :: dp=kind(1.d0)
real(dp), parameter :: PI=3.14159265358979323846264338328_dp
double precision, intent(in) :: xp(1:200)
integer :: n_r_X, n_theta_X
integer :: n_r_Y, n_theta_Y
real(dp), intent(in) :: r1, r2, cos_theta
real(dp), intent(out) :: muX

integer, parameter :: n_max_coeffs=200       ! Statically allocate arrays for the fit coeffs.
integer, parameter :: max_exponent=9
real(dp) :: d_X(n_max_coeffs)                ! fit coefficients X
real(dp) :: d_Y(n_max_coeffs)                ! fit coefficients Y
integer  :: index_X(3,n_max_coeffs)          ! index table containing the exponents of the X fit
integer  :: index_Y(3,n_max_coeffs)          ! index table containing the exponents of the Y fit
real(dp) :: powers_x1(0:max_exponent), powers_x2(0:max_exponent), powers_x3(0:max_exponent)
integer :: ii, ncoeffs_X, ncoeffs_Y
real(dp) :: x1, x2, x3, theta
real(dp) :: damp1, damp2, value1, value2, cos_theta_half, sin_theta_half
real(dp) :: damp, muOH
!**********************************************************************************************


   n_r_X = 8; n_theta_X = 8; ncoeffs_X = 200

!              X1_POWER            X2_POWER                 X3_POWER           PARAMETER
   index_X(1,001)= 0; index_X(2,001)= 0; index_X(3,001)= 1; d_X(001)=+xp(  1) 
   index_X(1,002)= 0; index_X(2,002)= 0; index_X(3,002)= 2; d_X(002)=+xp(  2) 
   index_X(1,003)= 0; index_X(2,003)= 0; index_X(3,003)= 3; d_X(003)=+xp(  3) 
   index_X(1,004)= 0; index_X(2,004)= 0; index_X(3,004)= 4; d_X(004)=+xp(  4) 
   index_X(1,005)= 0; index_X(2,005)= 0; index_X(3,005)= 5; d_X(005)=+xp(  5) 
   index_X(1,006)= 0; index_X(2,006)= 0; index_X(3,006)= 6; d_X(006)=+xp(  6) 
   index_X(1,007)= 0; index_X(2,007)= 0; index_X(3,007)= 7; d_X(007)=+xp(  7) 
   index_X(1,008)= 0; index_X(2,008)= 0; index_X(3,008)= 8; d_X(008)=+xp(  8) 
   index_X(1,009)= 0; index_X(2,009)= 2; index_X(3,009)= 1; d_X(009)=+xp(  9) 
   index_X(1,010)= 0; index_X(2,010)= 2; index_X(3,010)= 2; d_X(010)=+xp( 10) 
   index_X(1,011)= 0; index_X(2,011)= 2; index_X(3,011)= 3; d_X(011)=+xp( 11) 
   index_X(1,012)= 0; index_X(2,012)= 2; index_X(3,012)= 4; d_X(012)=+xp( 12) 
   index_X(1,013)= 0; index_X(2,013)= 2; index_X(3,013)= 5; d_X(013)=+xp( 13) 
   index_X(1,014)= 0; index_X(2,014)= 2; index_X(3,014)= 6; d_X(014)=+xp( 14) 
   index_X(1,015)= 0; index_X(2,015)= 2; index_X(3,015)= 7; d_X(015)=+xp( 15) 
   index_X(1,016)= 0; index_X(2,016)= 2; index_X(3,016)= 8; d_X(016)=+xp( 16) 
   index_X(1,017)= 0; index_X(2,017)= 4; index_X(3,017)= 1; d_X(017)=+xp( 17) 
   index_X(1,018)= 0; index_X(2,018)= 4; index_X(3,018)= 2; d_X(018)=+xp( 18) 
   index_X(1,019)= 0; index_X(2,019)= 4; index_X(3,019)= 3; d_X(019)=+xp( 19) 
   index_X(1,020)= 0; index_X(2,020)= 4; index_X(3,020)= 4; d_X(020)=+xp( 20) 
   index_X(1,021)= 0; index_X(2,021)= 4; index_X(3,021)= 5; d_X(021)=+xp( 21) 
   index_X(1,022)= 0; index_X(2,022)= 4; index_X(3,022)= 6; d_X(022)=+xp( 22) 
   index_X(1,023)= 0; index_X(2,023)= 4; index_X(3,023)= 7; d_X(023)=+xp( 23) 
   index_X(1,024)= 0; index_X(2,024)= 4; index_X(3,024)= 8; d_X(024)=+xp( 24) 
   index_X(1,025)= 0; index_X(2,025)= 6; index_X(3,025)= 1; d_X(025)=+xp( 25) 
   index_X(1,026)= 0; index_X(2,026)= 6; index_X(3,026)= 2; d_X(026)=+xp( 26) 
   index_X(1,027)= 0; index_X(2,027)= 6; index_X(3,027)= 3; d_X(027)=+xp( 27) 
   index_X(1,028)= 0; index_X(2,028)= 6; index_X(3,028)= 4; d_X(028)=+xp( 28) 
   index_X(1,029)= 0; index_X(2,029)= 6; index_X(3,029)= 5; d_X(029)=+xp( 29) 
   index_X(1,030)= 0; index_X(2,030)= 6; index_X(3,030)= 6; d_X(030)=+xp( 30) 
   index_X(1,031)= 0; index_X(2,031)= 6; index_X(3,031)= 7; d_X(031)=+xp( 31) 
   index_X(1,032)= 0; index_X(2,032)= 6; index_X(3,032)= 8; d_X(032)=+xp( 32) 
   index_X(1,033)= 0; index_X(2,033)= 8; index_X(3,033)= 1; d_X(033)=+xp( 33) 
   index_X(1,034)= 0; index_X(2,034)= 8; index_X(3,034)= 2; d_X(034)=+xp( 34) 
   index_X(1,035)= 0; index_X(2,035)= 8; index_X(3,035)= 3; d_X(035)=+xp( 35) 
   index_X(1,036)= 0; index_X(2,036)= 8; index_X(3,036)= 4; d_X(036)=+xp( 36) 
   index_X(1,037)= 0; index_X(2,037)= 8; index_X(3,037)= 5; d_X(037)=+xp( 37) 
   index_X(1,038)= 0; index_X(2,038)= 8; index_X(3,038)= 6; d_X(038)=+xp( 38) 
   index_X(1,039)= 0; index_X(2,039)= 8; index_X(3,039)= 7; d_X(039)=+xp( 39) 
   index_X(1,040)= 0; index_X(2,040)= 8; index_X(3,040)= 8; d_X(040)=+xp( 40) 
   index_X(1,041)= 1; index_X(2,041)= 0; index_X(3,041)= 1; d_X(041)=+xp( 41) 
   index_X(1,042)= 1; index_X(2,042)= 0; index_X(3,042)= 2; d_X(042)=+xp( 42) 
   index_X(1,043)= 1; index_X(2,043)= 0; index_X(3,043)= 3; d_X(043)=+xp( 43) 
   index_X(1,044)= 1; index_X(2,044)= 0; index_X(3,044)= 4; d_X(044)=+xp( 44) 
   index_X(1,045)= 1; index_X(2,045)= 0; index_X(3,045)= 5; d_X(045)=+xp( 45) 
   index_X(1,046)= 1; index_X(2,046)= 0; index_X(3,046)= 6; d_X(046)=+xp( 46) 
   index_X(1,047)= 1; index_X(2,047)= 0; index_X(3,047)= 7; d_X(047)=+xp( 47) 
   index_X(1,048)= 1; index_X(2,048)= 0; index_X(3,048)= 8; d_X(048)=+xp( 48) 
   index_X(1,049)= 1; index_X(2,049)= 2; index_X(3,049)= 1; d_X(049)=+xp( 49) 
   index_X(1,050)= 1; index_X(2,050)= 2; index_X(3,050)= 2; d_X(050)=+xp( 50) 
   index_X(1,051)= 1; index_X(2,051)= 2; index_X(3,051)= 3; d_X(051)=+xp( 51) 
   index_X(1,052)= 1; index_X(2,052)= 2; index_X(3,052)= 4; d_X(052)=+xp( 52) 
   index_X(1,053)= 1; index_X(2,053)= 2; index_X(3,053)= 5; d_X(053)=+xp( 53) 
   index_X(1,054)= 1; index_X(2,054)= 2; index_X(3,054)= 6; d_X(054)=+xp( 54) 
   index_X(1,055)= 1; index_X(2,055)= 2; index_X(3,055)= 7; d_X(055)=+xp( 55) 
   index_X(1,056)= 1; index_X(2,056)= 2; index_X(3,056)= 8; d_X(056)=+xp( 56) 
   index_X(1,057)= 1; index_X(2,057)= 4; index_X(3,057)= 1; d_X(057)=+xp( 57) 
   index_X(1,058)= 1; index_X(2,058)= 4; index_X(3,058)= 2; d_X(058)=+xp( 58) 
   index_X(1,059)= 1; index_X(2,059)= 4; index_X(3,059)= 3; d_X(059)=+xp( 59) 
   index_X(1,060)= 1; index_X(2,060)= 4; index_X(3,060)= 4; d_X(060)=+xp( 60) 
   index_X(1,061)= 1; index_X(2,061)= 4; index_X(3,061)= 5; d_X(061)=+xp( 61) 
   index_X(1,062)= 1; index_X(2,062)= 4; index_X(3,062)= 6; d_X(062)=+xp( 62) 
   index_X(1,063)= 1; index_X(2,063)= 4; index_X(3,063)= 7; d_X(063)=+xp( 63) 
   index_X(1,064)= 1; index_X(2,064)= 4; index_X(3,064)= 8; d_X(064)=+xp( 64) 
   index_X(1,065)= 1; index_X(2,065)= 6; index_X(3,065)= 1; d_X(065)=+xp( 65) 
   index_X(1,066)= 1; index_X(2,066)= 6; index_X(3,066)= 2; d_X(066)=+xp( 66) 
   index_X(1,067)= 1; index_X(2,067)= 6; index_X(3,067)= 3; d_X(067)=+xp( 67) 
   index_X(1,068)= 1; index_X(2,068)= 6; index_X(3,068)= 4; d_X(068)=+xp( 68) 
   index_X(1,069)= 1; index_X(2,069)= 6; index_X(3,069)= 5; d_X(069)=+xp( 69) 
   index_X(1,070)= 1; index_X(2,070)= 6; index_X(3,070)= 6; d_X(070)=+xp( 70) 
   index_X(1,071)= 1; index_X(2,071)= 6; index_X(3,071)= 7; d_X(071)=+xp( 71) 
   index_X(1,072)= 1; index_X(2,072)= 6; index_X(3,072)= 8; d_X(072)=+xp( 72) 
   index_X(1,073)= 2; index_X(2,073)= 0; index_X(3,073)= 1; d_X(073)=+xp( 73) 
   index_X(1,074)= 2; index_X(2,074)= 0; index_X(3,074)= 2; d_X(074)=+xp( 74) 
   index_X(1,075)= 2; index_X(2,075)= 0; index_X(3,075)= 3; d_X(075)=+xp( 75) 
   index_X(1,076)= 2; index_X(2,076)= 0; index_X(3,076)= 4; d_X(076)=+xp( 76) 
   index_X(1,077)= 2; index_X(2,077)= 0; index_X(3,077)= 5; d_X(077)=+xp( 77) 
   index_X(1,078)= 2; index_X(2,078)= 0; index_X(3,078)= 6; d_X(078)=+xp( 78) 
   index_X(1,079)= 2; index_X(2,079)= 0; index_X(3,079)= 7; d_X(079)=+xp( 79) 
   index_X(1,080)= 2; index_X(2,080)= 0; index_X(3,080)= 8; d_X(080)=+xp( 80) 
   index_X(1,081)= 2; index_X(2,081)= 2; index_X(3,081)= 1; d_X(081)=+xp( 81) 
   index_X(1,082)= 2; index_X(2,082)= 2; index_X(3,082)= 2; d_X(082)=+xp( 82) 
   index_X(1,083)= 2; index_X(2,083)= 2; index_X(3,083)= 3; d_X(083)=+xp( 83) 
   index_X(1,084)= 2; index_X(2,084)= 2; index_X(3,084)= 4; d_X(084)=+xp( 84) 
   index_X(1,085)= 2; index_X(2,085)= 2; index_X(3,085)= 5; d_X(085)=+xp( 85) 
   index_X(1,086)= 2; index_X(2,086)= 2; index_X(3,086)= 6; d_X(086)=+xp( 86) 
   index_X(1,087)= 2; index_X(2,087)= 2; index_X(3,087)= 7; d_X(087)=+xp( 87) 
   index_X(1,088)= 2; index_X(2,088)= 2; index_X(3,088)= 8; d_X(088)=+xp( 88) 
   index_X(1,089)= 2; index_X(2,089)= 4; index_X(3,089)= 1; d_X(089)=+xp( 89) 
   index_X(1,090)= 2; index_X(2,090)= 4; index_X(3,090)= 2; d_X(090)=+xp( 90) 
   index_X(1,091)= 2; index_X(2,091)= 4; index_X(3,091)= 3; d_X(091)=+xp( 91) 
   index_X(1,092)= 2; index_X(2,092)= 4; index_X(3,092)= 4; d_X(092)=+xp( 92) 
   index_X(1,093)= 2; index_X(2,093)= 4; index_X(3,093)= 5; d_X(093)=+xp( 93) 
   index_X(1,094)= 2; index_X(2,094)= 4; index_X(3,094)= 6; d_X(094)=+xp( 94) 
   index_X(1,095)= 2; index_X(2,095)= 4; index_X(3,095)= 7; d_X(095)=+xp( 95) 
   index_X(1,096)= 2; index_X(2,096)= 4; index_X(3,096)= 8; d_X(096)=+xp( 96) 
   index_X(1,097)= 2; index_X(2,097)= 6; index_X(3,097)= 1; d_X(097)=+xp( 97) 
   index_X(1,098)= 2; index_X(2,098)= 6; index_X(3,098)= 2; d_X(098)=+xp( 98) 
   index_X(1,099)= 2; index_X(2,099)= 6; index_X(3,099)= 3; d_X(099)=+xp( 99) 
   index_X(1,100)= 2; index_X(2,100)= 6; index_X(3,100)= 4; d_X(100)=+xp(100) 
   index_X(1,101)= 2; index_X(2,101)= 6; index_X(3,101)= 5; d_X(101)=+xp(101) 
   index_X(1,102)= 2; index_X(2,102)= 6; index_X(3,102)= 6; d_X(102)=+xp(102) 
   index_X(1,103)= 2; index_X(2,103)= 6; index_X(3,103)= 7; d_X(103)=+xp(103) 
   index_X(1,104)= 2; index_X(2,104)= 6; index_X(3,104)= 8; d_X(104)=+xp(104) 
   index_X(1,105)= 3; index_X(2,105)= 0; index_X(3,105)= 1; d_X(105)=+xp(105) 
   index_X(1,106)= 3; index_X(2,106)= 0; index_X(3,106)= 2; d_X(106)=+xp(106) 
   index_X(1,107)= 3; index_X(2,107)= 0; index_X(3,107)= 3; d_X(107)=+xp(107) 
   index_X(1,108)= 3; index_X(2,108)= 0; index_X(3,108)= 4; d_X(108)=+xp(108) 
   index_X(1,109)= 3; index_X(2,109)= 0; index_X(3,109)= 5; d_X(109)=+xp(109) 
   index_X(1,110)= 3; index_X(2,110)= 0; index_X(3,110)= 6; d_X(110)=+xp(110) 
   index_X(1,111)= 3; index_X(2,111)= 0; index_X(3,111)= 7; d_X(111)=+xp(111) 
   index_X(1,112)= 3; index_X(2,112)= 0; index_X(3,112)= 8; d_X(112)=+xp(112) 
   index_X(1,113)= 3; index_X(2,113)= 2; index_X(3,113)= 1; d_X(113)=+xp(113) 
   index_X(1,114)= 3; index_X(2,114)= 2; index_X(3,114)= 2; d_X(114)=+xp(114) 
   index_X(1,115)= 3; index_X(2,115)= 2; index_X(3,115)= 3; d_X(115)=+xp(115) 
   index_X(1,116)= 3; index_X(2,116)= 2; index_X(3,116)= 4; d_X(116)=+xp(116) 
   index_X(1,117)= 3; index_X(2,117)= 2; index_X(3,117)= 5; d_X(117)=+xp(117) 
   index_X(1,118)= 3; index_X(2,118)= 2; index_X(3,118)= 6; d_X(118)=+xp(118) 
   index_X(1,119)= 3; index_X(2,119)= 2; index_X(3,119)= 7; d_X(119)=+xp(119) 
   index_X(1,120)= 3; index_X(2,120)= 2; index_X(3,120)= 8; d_X(120)=+xp(120) 
   index_X(1,121)= 3; index_X(2,121)= 4; index_X(3,121)= 1; d_X(121)=+xp(121) 
   index_X(1,122)= 3; index_X(2,122)= 4; index_X(3,122)= 2; d_X(122)=+xp(122) 
   index_X(1,123)= 3; index_X(2,123)= 4; index_X(3,123)= 3; d_X(123)=+xp(123) 
   index_X(1,124)= 3; index_X(2,124)= 4; index_X(3,124)= 4; d_X(124)=+xp(124) 
   index_X(1,125)= 3; index_X(2,125)= 4; index_X(3,125)= 5; d_X(125)=+xp(125) 
   index_X(1,126)= 3; index_X(2,126)= 4; index_X(3,126)= 6; d_X(126)=+xp(126) 
   index_X(1,127)= 3; index_X(2,127)= 4; index_X(3,127)= 7; d_X(127)=+xp(127) 
   index_X(1,128)= 3; index_X(2,128)= 4; index_X(3,128)= 8; d_X(128)=+xp(128) 
   index_X(1,129)= 4; index_X(2,129)= 0; index_X(3,129)= 1; d_X(129)=+xp(129) 
   index_X(1,130)= 4; index_X(2,130)= 0; index_X(3,130)= 2; d_X(130)=+xp(130) 
   index_X(1,131)= 4; index_X(2,131)= 0; index_X(3,131)= 3; d_X(131)=+xp(131) 
   index_X(1,132)= 4; index_X(2,132)= 0; index_X(3,132)= 4; d_X(132)=+xp(132) 
   index_X(1,133)= 4; index_X(2,133)= 0; index_X(3,133)= 5; d_X(133)=+xp(133) 
   index_X(1,134)= 4; index_X(2,134)= 0; index_X(3,134)= 6; d_X(134)=+xp(134) 
   index_X(1,135)= 4; index_X(2,135)= 0; index_X(3,135)= 7; d_X(135)=+xp(135) 
   index_X(1,136)= 4; index_X(2,136)= 0; index_X(3,136)= 8; d_X(136)=+xp(136) 
   index_X(1,137)= 4; index_X(2,137)= 2; index_X(3,137)= 1; d_X(137)=+xp(137) 
   index_X(1,138)= 4; index_X(2,138)= 2; index_X(3,138)= 2; d_X(138)=+xp(138) 
   index_X(1,139)= 4; index_X(2,139)= 2; index_X(3,139)= 3; d_X(139)=+xp(139) 
   index_X(1,140)= 4; index_X(2,140)= 2; index_X(3,140)= 4; d_X(140)=+xp(140) 
   index_X(1,141)= 4; index_X(2,141)= 2; index_X(3,141)= 5; d_X(141)=+xp(141) 
   index_X(1,142)= 4; index_X(2,142)= 2; index_X(3,142)= 6; d_X(142)=+xp(142) 
   index_X(1,143)= 4; index_X(2,143)= 2; index_X(3,143)= 7; d_X(143)=+xp(143) 
   index_X(1,144)= 4; index_X(2,144)= 2; index_X(3,144)= 8; d_X(144)=+xp(144) 
   index_X(1,145)= 4; index_X(2,145)= 4; index_X(3,145)= 1; d_X(145)=+xp(145) 
   index_X(1,146)= 4; index_X(2,146)= 4; index_X(3,146)= 2; d_X(146)=+xp(146) 
   index_X(1,147)= 4; index_X(2,147)= 4; index_X(3,147)= 3; d_X(147)=+xp(147) 
   index_X(1,148)= 4; index_X(2,148)= 4; index_X(3,148)= 4; d_X(148)=+xp(148) 
   index_X(1,149)= 4; index_X(2,149)= 4; index_X(3,149)= 5; d_X(149)=+xp(149) 
   index_X(1,150)= 4; index_X(2,150)= 4; index_X(3,150)= 6; d_X(150)=+xp(150) 
   index_X(1,151)= 4; index_X(2,151)= 4; index_X(3,151)= 7; d_X(151)=+xp(151) 
   index_X(1,152)= 4; index_X(2,152)= 4; index_X(3,152)= 8; d_X(152)=+xp(152) 
   index_X(1,153)= 5; index_X(2,153)= 0; index_X(3,153)= 1; d_X(153)=+xp(153) 
   index_X(1,154)= 5; index_X(2,154)= 0; index_X(3,154)= 2; d_X(154)=+xp(154) 
   index_X(1,155)= 5; index_X(2,155)= 0; index_X(3,155)= 3; d_X(155)=+xp(155) 
   index_X(1,156)= 5; index_X(2,156)= 0; index_X(3,156)= 4; d_X(156)=+xp(156) 
   index_X(1,157)= 5; index_X(2,157)= 0; index_X(3,157)= 5; d_X(157)=+xp(157) 
   index_X(1,158)= 5; index_X(2,158)= 0; index_X(3,158)= 6; d_X(158)=+xp(158) 
   index_X(1,159)= 5; index_X(2,159)= 0; index_X(3,159)= 7; d_X(159)=+xp(159) 
   index_X(1,160)= 5; index_X(2,160)= 0; index_X(3,160)= 8; d_X(160)=+xp(160) 
   index_X(1,161)= 5; index_X(2,161)= 2; index_X(3,161)= 1; d_X(161)=+xp(161) 
   index_X(1,162)= 5; index_X(2,162)= 2; index_X(3,162)= 2; d_X(162)=+xp(162) 
   index_X(1,163)= 5; index_X(2,163)= 2; index_X(3,163)= 3; d_X(163)=+xp(163) 
   index_X(1,164)= 5; index_X(2,164)= 2; index_X(3,164)= 4; d_X(164)=+xp(164) 
   index_X(1,165)= 5; index_X(2,165)= 2; index_X(3,165)= 5; d_X(165)=+xp(165) 
   index_X(1,166)= 5; index_X(2,166)= 2; index_X(3,166)= 6; d_X(166)=+xp(166) 
   index_X(1,167)= 5; index_X(2,167)= 2; index_X(3,167)= 7; d_X(167)=+xp(167) 
   index_X(1,168)= 5; index_X(2,168)= 2; index_X(3,168)= 8; d_X(168)=+xp(168) 
   index_X(1,169)= 6; index_X(2,169)= 0; index_X(3,169)= 1; d_X(169)=+xp(169) 
   index_X(1,170)= 6; index_X(2,170)= 0; index_X(3,170)= 2; d_X(170)=+xp(170) 
   index_X(1,171)= 6; index_X(2,171)= 0; index_X(3,171)= 3; d_X(171)=+xp(171) 
   index_X(1,172)= 6; index_X(2,172)= 0; index_X(3,172)= 4; d_X(172)=+xp(172) 
   index_X(1,173)= 6; index_X(2,173)= 0; index_X(3,173)= 5; d_X(173)=+xp(173) 
   index_X(1,174)= 6; index_X(2,174)= 0; index_X(3,174)= 6; d_X(174)=+xp(174) 
   index_X(1,175)= 6; index_X(2,175)= 0; index_X(3,175)= 7; d_X(175)=+xp(175) 
   index_X(1,176)= 6; index_X(2,176)= 0; index_X(3,176)= 8; d_X(176)=+xp(176) 
   index_X(1,177)= 6; index_X(2,177)= 2; index_X(3,177)= 1; d_X(177)=+xp(177) 
   index_X(1,178)= 6; index_X(2,178)= 2; index_X(3,178)= 2; d_X(178)=+xp(178) 
   index_X(1,179)= 6; index_X(2,179)= 2; index_X(3,179)= 3; d_X(179)=+xp(179) 
   index_X(1,180)= 6; index_X(2,180)= 2; index_X(3,180)= 4; d_X(180)=+xp(180) 
   index_X(1,181)= 6; index_X(2,181)= 2; index_X(3,181)= 5; d_X(181)=+xp(181) 
   index_X(1,182)= 6; index_X(2,182)= 2; index_X(3,182)= 6; d_X(182)=+xp(182) 
   index_X(1,183)= 6; index_X(2,183)= 2; index_X(3,183)= 7; d_X(183)=+xp(183) 
   index_X(1,184)= 6; index_X(2,184)= 2; index_X(3,184)= 8; d_X(184)=+xp(184) 
   index_X(1,185)= 7; index_X(2,185)= 0; index_X(3,185)= 1; d_X(185)=+xp(185) 
   index_X(1,186)= 7; index_X(2,186)= 0; index_X(3,186)= 2; d_X(186)=+xp(186) 
   index_X(1,187)= 7; index_X(2,187)= 0; index_X(3,187)= 3; d_X(187)=+xp(187) 
   index_X(1,188)= 7; index_X(2,188)= 0; index_X(3,188)= 4; d_X(188)=+xp(188) 
   index_X(1,189)= 7; index_X(2,189)= 0; index_X(3,189)= 5; d_X(189)=+xp(189) 
   index_X(1,190)= 7; index_X(2,190)= 0; index_X(3,190)= 6; d_X(190)=+xp(190) 
   index_X(1,191)= 7; index_X(2,191)= 0; index_X(3,191)= 7; d_X(191)=+xp(191) 
   index_X(1,192)= 7; index_X(2,192)= 0; index_X(3,192)= 8; d_X(192)=+xp(192) 
   index_X(1,193)= 8; index_X(2,193)= 0; index_X(3,193)= 1; d_X(193)=+xp(193) 
   index_X(1,194)= 8; index_X(2,194)= 0; index_X(3,194)= 2; d_X(194)=+xp(194) 
   index_X(1,195)= 8; index_X(2,195)= 0; index_X(3,195)= 3; d_X(195)=+xp(195) 
   index_X(1,196)= 8; index_X(2,196)= 0; index_X(3,196)= 4; d_X(196)=+xp(196) 
   index_X(1,197)= 8; index_X(2,197)= 0; index_X(3,197)= 5; d_X(197)=+xp(197) 
   index_X(1,198)= 8; index_X(2,198)= 0; index_X(3,198)= 6; d_X(198)=+xp(198) 
   index_X(1,199)= 8; index_X(2,199)= 0; index_X(3,199)= 7; d_X(199)=+xp(199) 
   index_X(1,200)= 8; index_X(2,200)= 0; index_X(3,200)= 8; d_X(200)=+xp(200) 

!***************                                                     
! FITTING VARIABLES                                                  
!theta = dacos(cos_theta)    
!theta=cos_theta

!***********
! Check that input makes sense
!if( abs(cos_theta) > 1._dp) write(*,*) 'WARNING: |cos(theta)| > 1 in DIPS'
if( r1 < 0._dp)             write(*,*) 'WARNING: r1 < 0 in DIPS'
if( r2 < 0._dp)             write(*,*) 'WARNING: r2 < 0 in DIPS'


!***************
! FITTING VARIABLES
theta = cos_theta!acos(cos_theta)
!***************
x1 = r1 + r2
x2 = r2 - r1
x3 = PI - theta
!***************


!***************
! PRE-COMPUTE POWERS OF THE VARIABLES
powers_x1(0) = 1._dp
powers_x2(0) = 1._dp
powers_x3(0) = 1._dp

do ii=1, max(n_r_X, n_r_Y)
 powers_x1(ii) =  x1*powers_x1(ii-1)
 powers_x2(ii) =  x2*powers_x2(ii-1)
enddo

do ii=1, max(n_theta_X, n_theta_Y)
 powers_x3(ii) =  x3*powers_x3(ii-1)
enddo

!***************
muX = 0._dp
do ii = 1, ncoeffs_X
  muX = muX + d_X(ii) * powers_x1(index_X(1,ii)) * powers_x2(index_X(2,ii)) * powers_x3(index_X(3,ii))
enddo
!***************

!add damping/asymptotic terms
damp1 = damp(r1)
damp2 = damp(r2)
value1 = muOH(r1)*( 1._dp - damp2 )
value2 = muOH(r2)*( 1._dp - damp1 )

cos_theta_half = sqrt( (1._dp + cos_theta)/2._dp  )
sin_theta_half = sqrt( (1._dp - cos_theta)/2._dp  )

muX = muX*damp1*damp2 + cos_theta_half*( value1 + value2  )


!************************************

end subroutine poten


!==================================
function damp(r)
! Damping function.
! It is exactly =0 for r< rMIN;
! It is exactly =1 for r> rMAX;
IMPLICIT NONE
integer, parameter :: fp =kind(1.d0)
real(fp), parameter :: PI=3.14159265358979323846264338328_fp
real(fp), intent(in) :: r
real(fp) :: damp
real(fp), parameter :: rMIN=3.5_fp, rMAX=5.5_fp
real(fp) :: t

t = ( ( (r-rMIN)/(rMAX-rMIN) ) -0.5_fp )

if( t <= -0.5_fp) then
   damp = 1._fp
   return
elseif(t >= 0.5_fp) then
   damp = 0._fp
   return
endif

damp = 1._fp / ( EXP(4._fp*tan(PI*t)) + 1._fp )

end function damp


!==================================
function muOH(r)
IMPLICIT NONE
! Lorenzo Lodi 7 April 2011
! It's a fit in even-tempered Gaussian functions of 59 AQCC[7,5]/XP aug-cc-pV5Z dipoles for the ground X state of OH
! between 1 and 4 bohrs. Dipoles should be accurate to about 2.e-3. Relativistic & core-corrections are NOT included.
! NB: There is a NASTY cancellation error with this fit which eats up about 5 digits of accuracy.
! double precision reals are compulsory. RMS of the fit is 5e-5 a.u. The largest residuals are about 0.7e-3 a.u.
! Note: Kahan summation algorithm (compensated summation) does not help here. The problem is ill-conditioned
!       as sum_i |x_i| / |sum_i x_i| = ~ 3e+7
! It might be solved by using as basis orthogonalised Gaussians, but I didn't deem it necessary at this time.
!
! Fit has reasonable behaviour for r-> +oo and also not too bad for r->0.
integer, parameter :: fp =kind(1.d0)
real(fp), intent(in) :: r
real(fp) :: muOH
integer, parameter :: n_coeffs = 8
real(fp), parameter :: a = 0.934_fp
real(fp) :: coeffs(n_coeffs)
real(fp) :: powers_a(n_coeffs)
real(fp) :: r2
integer :: i
logical :: zFirstRun  = .true.
save zFirstRun, coeffs, powers_a


 powers_a(1) = a
 do i=2, n_coeffs
   powers_a(i) = powers_a(i-1)*a
 enddo

 coeffs( 1) =   -51632.29191732509_fp
 coeffs( 2) =   336439.76799518_fp
 coeffs( 3) =  -948384.0949923932_fp
 coeffs( 4) =  1499393.4799029354_fp
 coeffs( 5) = -1436572.608839631_fp
 coeffs( 6) =   834824.6196200893_fp
 coeffs( 7) =  -272815.3271033254_fp
 coeffs( 8) =    38746.11983311285_fp


r2 = r**2
muOH = 0._fp
do i=1, n_coeffs
  muOH = muOH + coeffs(i)*exp( -powers_a(i)*r2 )
enddo

muOH = muOH*r2

end function muOH


