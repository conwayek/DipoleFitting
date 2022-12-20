subroutine poten(xp,muY, r1, r2, cos_theta)

implicit none
integer , parameter :: dp=kind(1.d0)
real(dp), parameter :: PI=3.14159265358979323846264338328_dp
double precision, intent(in) :: xp(1:225)
integer :: n_r_X, n_theta_X
integer :: n_r_Y, n_theta_Y
real(dp), intent(in) :: r1, r2, cos_theta
real(dp), intent(out) :: muY

integer, parameter :: n_max_coeffs=225       ! Statically allocate arrays for the fit coeffs.
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

   n_r_Y = 9; n_theta_Y = 8; ncoeffs_Y = 225
   index_Y(1,001)= 0; index_Y(2,001)= 1; index_Y(3,001)= 0; d_Y(001)=+xp(  1) 
   index_Y(1,002)= 0; index_Y(2,002)= 1; index_Y(3,002)= 1; d_Y(002)=+xp(  2) 
   index_Y(1,003)= 0; index_Y(2,003)= 1; index_Y(3,003)= 2; d_Y(003)=+xp(  3) 
   index_Y(1,004)= 0; index_Y(2,004)= 1; index_Y(3,004)= 3; d_Y(004)=+xp(  4) 
   index_Y(1,005)= 0; index_Y(2,005)= 1; index_Y(3,005)= 4; d_Y(005)=+xp(  5) 
   index_Y(1,006)= 0; index_Y(2,006)= 1; index_Y(3,006)= 5; d_Y(006)=+xp(  6) 
   index_Y(1,007)= 0; index_Y(2,007)= 1; index_Y(3,007)= 6; d_Y(007)=+xp(  7) 
   index_Y(1,008)= 0; index_Y(2,008)= 1; index_Y(3,008)= 7; d_Y(008)=+xp(  8) 
   index_Y(1,009)= 0; index_Y(2,009)= 1; index_Y(3,009)= 8; d_Y(009)=+xp(  9) 
   index_Y(1,010)= 0; index_Y(2,010)= 3; index_Y(3,010)= 0; d_Y(010)=+xp( 10) 
   index_Y(1,011)= 0; index_Y(2,011)= 3; index_Y(3,011)= 1; d_Y(011)=+xp( 11) 
   index_Y(1,012)= 0; index_Y(2,012)= 3; index_Y(3,012)= 2; d_Y(012)=+xp( 12) 
   index_Y(1,013)= 0; index_Y(2,013)= 3; index_Y(3,013)= 3; d_Y(013)=+xp( 13) 
   index_Y(1,014)= 0; index_Y(2,014)= 3; index_Y(3,014)= 4; d_Y(014)=+xp( 14) 
   index_Y(1,015)= 0; index_Y(2,015)= 3; index_Y(3,015)= 5; d_Y(015)=+xp( 15) 
   index_Y(1,016)= 0; index_Y(2,016)= 3; index_Y(3,016)= 6; d_Y(016)=+xp( 16) 
   index_Y(1,017)= 0; index_Y(2,017)= 3; index_Y(3,017)= 7; d_Y(017)=+xp( 17) 
   index_Y(1,018)= 0; index_Y(2,018)= 3; index_Y(3,018)= 8; d_Y(018)=+xp( 18) 
   index_Y(1,019)= 0; index_Y(2,019)= 5; index_Y(3,019)= 0; d_Y(019)=+xp( 19) 
   index_Y(1,020)= 0; index_Y(2,020)= 5; index_Y(3,020)= 1; d_Y(020)=+xp( 20) 
   index_Y(1,021)= 0; index_Y(2,021)= 5; index_Y(3,021)= 2; d_Y(021)=+xp( 21) 
   index_Y(1,022)= 0; index_Y(2,022)= 5; index_Y(3,022)= 3; d_Y(022)=+xp( 22) 
   index_Y(1,023)= 0; index_Y(2,023)= 5; index_Y(3,023)= 4; d_Y(023)=+xp( 23) 
   index_Y(1,024)= 0; index_Y(2,024)= 5; index_Y(3,024)= 5; d_Y(024)=+xp( 24) 
   index_Y(1,025)= 0; index_Y(2,025)= 5; index_Y(3,025)= 6; d_Y(025)=+xp( 25) 
   index_Y(1,026)= 0; index_Y(2,026)= 5; index_Y(3,026)= 7; d_Y(026)=+xp( 26) 
   index_Y(1,027)= 0; index_Y(2,027)= 5; index_Y(3,027)= 8; d_Y(027)=+xp( 27) 
   index_Y(1,028)= 0; index_Y(2,028)= 7; index_Y(3,028)= 0; d_Y(028)=+xp( 28) 
   index_Y(1,029)= 0; index_Y(2,029)= 7; index_Y(3,029)= 1; d_Y(029)=+xp( 29) 
   index_Y(1,030)= 0; index_Y(2,030)= 7; index_Y(3,030)= 2; d_Y(030)=+xp( 30) 
   index_Y(1,031)= 0; index_Y(2,031)= 7; index_Y(3,031)= 3; d_Y(031)=+xp( 31) 
   index_Y(1,032)= 0; index_Y(2,032)= 7; index_Y(3,032)= 4; d_Y(032)=+xp( 32) 
   index_Y(1,033)= 0; index_Y(2,033)= 7; index_Y(3,033)= 5; d_Y(033)=+xp( 33) 
   index_Y(1,034)= 0; index_Y(2,034)= 7; index_Y(3,034)= 6; d_Y(034)=+xp( 34) 
   index_Y(1,035)= 0; index_Y(2,035)= 7; index_Y(3,035)= 7; d_Y(035)=+xp( 35) 
   index_Y(1,036)= 0; index_Y(2,036)= 7; index_Y(3,036)= 8; d_Y(036)=+xp( 36) 
   index_Y(1,037)= 0; index_Y(2,037)= 9; index_Y(3,037)= 0; d_Y(037)=+xp( 37) 
   index_Y(1,038)= 0; index_Y(2,038)= 9; index_Y(3,038)= 1; d_Y(038)=+xp( 38) 
   index_Y(1,039)= 0; index_Y(2,039)= 9; index_Y(3,039)= 2; d_Y(039)=+xp( 39) 
   index_Y(1,040)= 0; index_Y(2,040)= 9; index_Y(3,040)= 3; d_Y(040)=+xp( 40) 
   index_Y(1,041)= 0; index_Y(2,041)= 9; index_Y(3,041)= 4; d_Y(041)=+xp( 41) 
   index_Y(1,042)= 0; index_Y(2,042)= 9; index_Y(3,042)= 5; d_Y(042)=+xp( 42) 
   index_Y(1,043)= 0; index_Y(2,043)= 9; index_Y(3,043)= 6; d_Y(043)=+xp( 43) 
   index_Y(1,044)= 0; index_Y(2,044)= 9; index_Y(3,044)= 7; d_Y(044)=+xp( 44) 
   index_Y(1,045)= 0; index_Y(2,045)= 9; index_Y(3,045)= 8; d_Y(045)=+xp( 45) 
   index_Y(1,046)= 1; index_Y(2,046)= 1; index_Y(3,046)= 0; d_Y(046)=+xp( 46) 
   index_Y(1,047)= 1; index_Y(2,047)= 1; index_Y(3,047)= 1; d_Y(047)=+xp( 47) 
   index_Y(1,048)= 1; index_Y(2,048)= 1; index_Y(3,048)= 2; d_Y(048)=+xp( 48) 
   index_Y(1,049)= 1; index_Y(2,049)= 1; index_Y(3,049)= 3; d_Y(049)=+xp( 49) 
   index_Y(1,050)= 1; index_Y(2,050)= 1; index_Y(3,050)= 4; d_Y(050)=+xp( 50) 
   index_Y(1,051)= 1; index_Y(2,051)= 1; index_Y(3,051)= 5; d_Y(051)=+xp( 51) 
   index_Y(1,052)= 1; index_Y(2,052)= 1; index_Y(3,052)= 6; d_Y(052)=+xp( 52) 
   index_Y(1,053)= 1; index_Y(2,053)= 1; index_Y(3,053)= 7; d_Y(053)=+xp( 53) 
   index_Y(1,054)= 1; index_Y(2,054)= 1; index_Y(3,054)= 8; d_Y(054)=+xp( 54) 
   index_Y(1,055)= 1; index_Y(2,055)= 3; index_Y(3,055)= 0; d_Y(055)=+xp( 55) 
   index_Y(1,056)= 1; index_Y(2,056)= 3; index_Y(3,056)= 1; d_Y(056)=+xp( 56) 
   index_Y(1,057)= 1; index_Y(2,057)= 3; index_Y(3,057)= 2; d_Y(057)=+xp( 57) 
   index_Y(1,058)= 1; index_Y(2,058)= 3; index_Y(3,058)= 3; d_Y(058)=+xp( 58) 
   index_Y(1,059)= 1; index_Y(2,059)= 3; index_Y(3,059)= 4; d_Y(059)=+xp( 59) 
   index_Y(1,060)= 1; index_Y(2,060)= 3; index_Y(3,060)= 5; d_Y(060)=+xp( 60) 
   index_Y(1,061)= 1; index_Y(2,061)= 3; index_Y(3,061)= 6; d_Y(061)=+xp( 61) 
   index_Y(1,062)= 1; index_Y(2,062)= 3; index_Y(3,062)= 7; d_Y(062)=+xp( 62) 
   index_Y(1,063)= 1; index_Y(2,063)= 3; index_Y(3,063)= 8; d_Y(063)=+xp( 63) 
   index_Y(1,064)= 1; index_Y(2,064)= 5; index_Y(3,064)= 0; d_Y(064)=+xp( 64) 
   index_Y(1,065)= 1; index_Y(2,065)= 5; index_Y(3,065)= 1; d_Y(065)=+xp( 65) 
   index_Y(1,066)= 1; index_Y(2,066)= 5; index_Y(3,066)= 2; d_Y(066)=+xp( 66) 
   index_Y(1,067)= 1; index_Y(2,067)= 5; index_Y(3,067)= 3; d_Y(067)=+xp( 67) 
   index_Y(1,068)= 1; index_Y(2,068)= 5; index_Y(3,068)= 4; d_Y(068)=+xp( 68) 
   index_Y(1,069)= 1; index_Y(2,069)= 5; index_Y(3,069)= 5; d_Y(069)=+xp( 69) 
   index_Y(1,070)= 1; index_Y(2,070)= 5; index_Y(3,070)= 6; d_Y(070)=+xp( 70) 
   index_Y(1,071)= 1; index_Y(2,071)= 5; index_Y(3,071)= 7; d_Y(071)=+xp( 71) 
   index_Y(1,072)= 1; index_Y(2,072)= 5; index_Y(3,072)= 8; d_Y(072)=+xp( 72) 
   index_Y(1,073)= 1; index_Y(2,073)= 7; index_Y(3,073)= 0; d_Y(073)=+xp( 73) 
   index_Y(1,074)= 1; index_Y(2,074)= 7; index_Y(3,074)= 1; d_Y(074)=+xp( 74) 
   index_Y(1,075)= 1; index_Y(2,075)= 7; index_Y(3,075)= 2; d_Y(075)=+xp( 75) 
   index_Y(1,076)= 1; index_Y(2,076)= 7; index_Y(3,076)= 3; d_Y(076)=+xp( 76) 
   index_Y(1,077)= 1; index_Y(2,077)= 7; index_Y(3,077)= 4; d_Y(077)=+xp( 77) 
   index_Y(1,078)= 1; index_Y(2,078)= 7; index_Y(3,078)= 5; d_Y(078)=+xp( 78) 
   index_Y(1,079)= 1; index_Y(2,079)= 7; index_Y(3,079)= 6; d_Y(079)=+xp( 79) 
   index_Y(1,080)= 1; index_Y(2,080)= 7; index_Y(3,080)= 7; d_Y(080)=+xp( 80) 
   index_Y(1,081)= 1; index_Y(2,081)= 7; index_Y(3,081)= 8; d_Y(081)=+xp( 81) 
   index_Y(1,082)= 2; index_Y(2,082)= 1; index_Y(3,082)= 0; d_Y(082)=+xp( 82) 
   index_Y(1,083)= 2; index_Y(2,083)= 1; index_Y(3,083)= 1; d_Y(083)=+xp( 83) 
   index_Y(1,084)= 2; index_Y(2,084)= 1; index_Y(3,084)= 2; d_Y(084)=+xp( 84) 
   index_Y(1,085)= 2; index_Y(2,085)= 1; index_Y(3,085)= 3; d_Y(085)=+xp( 85) 
   index_Y(1,086)= 2; index_Y(2,086)= 1; index_Y(3,086)= 4; d_Y(086)=+xp( 86) 
   index_Y(1,087)= 2; index_Y(2,087)= 1; index_Y(3,087)= 5; d_Y(087)=+xp( 87) 
   index_Y(1,088)= 2; index_Y(2,088)= 1; index_Y(3,088)= 6; d_Y(088)=+xp( 88) 
   index_Y(1,089)= 2; index_Y(2,089)= 1; index_Y(3,089)= 7; d_Y(089)=+xp( 89) 
   index_Y(1,090)= 2; index_Y(2,090)= 1; index_Y(3,090)= 8; d_Y(090)=+xp( 90) 
   index_Y(1,091)= 2; index_Y(2,091)= 3; index_Y(3,091)= 0; d_Y(091)=+xp( 91) 
   index_Y(1,092)= 2; index_Y(2,092)= 3; index_Y(3,092)= 1; d_Y(092)=+xp( 92) 
   index_Y(1,093)= 2; index_Y(2,093)= 3; index_Y(3,093)= 2; d_Y(093)=+xp( 93) 
   index_Y(1,094)= 2; index_Y(2,094)= 3; index_Y(3,094)= 3; d_Y(094)=+xp( 94) 
   index_Y(1,095)= 2; index_Y(2,095)= 3; index_Y(3,095)= 4; d_Y(095)=+xp( 95) 
   index_Y(1,096)= 2; index_Y(2,096)= 3; index_Y(3,096)= 5; d_Y(096)=+xp( 96) 
   index_Y(1,097)= 2; index_Y(2,097)= 3; index_Y(3,097)= 6; d_Y(097)=+xp( 97) 
   index_Y(1,098)= 2; index_Y(2,098)= 3; index_Y(3,098)= 7; d_Y(098)=+xp( 98) 
   index_Y(1,099)= 2; index_Y(2,099)= 3; index_Y(3,099)= 8; d_Y(099)=+xp( 99) 
   index_Y(1,100)= 2; index_Y(2,100)= 5; index_Y(3,100)= 0; d_Y(100)=+xp(100) 
   index_Y(1,101)= 2; index_Y(2,101)= 5; index_Y(3,101)= 1; d_Y(101)=+xp(101) 
   index_Y(1,102)= 2; index_Y(2,102)= 5; index_Y(3,102)= 2; d_Y(102)=+xp(102) 
   index_Y(1,103)= 2; index_Y(2,103)= 5; index_Y(3,103)= 3; d_Y(103)=+xp(103) 
   index_Y(1,104)= 2; index_Y(2,104)= 5; index_Y(3,104)= 4; d_Y(104)=+xp(104) 
   index_Y(1,105)= 2; index_Y(2,105)= 5; index_Y(3,105)= 5; d_Y(105)=+xp(105) 
   index_Y(1,106)= 2; index_Y(2,106)= 5; index_Y(3,106)= 6; d_Y(106)=+xp(106) 
   index_Y(1,107)= 2; index_Y(2,107)= 5; index_Y(3,107)= 7; d_Y(107)=+xp(107) 
   index_Y(1,108)= 2; index_Y(2,108)= 5; index_Y(3,108)= 8; d_Y(108)=+xp(108) 
   index_Y(1,109)= 2; index_Y(2,109)= 7; index_Y(3,109)= 0; d_Y(109)=+xp(109) 
   index_Y(1,110)= 2; index_Y(2,110)= 7; index_Y(3,110)= 1; d_Y(110)=+xp(110) 
   index_Y(1,111)= 2; index_Y(2,111)= 7; index_Y(3,111)= 2; d_Y(111)=+xp(111) 
   index_Y(1,112)= 2; index_Y(2,112)= 7; index_Y(3,112)= 3; d_Y(112)=+xp(112) 
   index_Y(1,113)= 2; index_Y(2,113)= 7; index_Y(3,113)= 4; d_Y(113)=+xp(113) 
   index_Y(1,114)= 2; index_Y(2,114)= 7; index_Y(3,114)= 5; d_Y(114)=+xp(114) 
   index_Y(1,115)= 2; index_Y(2,115)= 7; index_Y(3,115)= 6; d_Y(115)=+xp(115) 
   index_Y(1,116)= 2; index_Y(2,116)= 7; index_Y(3,116)= 7; d_Y(116)=+xp(116) 
   index_Y(1,117)= 2; index_Y(2,117)= 7; index_Y(3,117)= 8; d_Y(117)=+xp(117) 
   index_Y(1,118)= 3; index_Y(2,118)= 1; index_Y(3,118)= 0; d_Y(118)=+xp(118) 
   index_Y(1,119)= 3; index_Y(2,119)= 1; index_Y(3,119)= 1; d_Y(119)=+xp(119) 
   index_Y(1,120)= 3; index_Y(2,120)= 1; index_Y(3,120)= 2; d_Y(120)=+xp(120) 
   index_Y(1,121)= 3; index_Y(2,121)= 1; index_Y(3,121)= 3; d_Y(121)=+xp(121) 
   index_Y(1,122)= 3; index_Y(2,122)= 1; index_Y(3,122)= 4; d_Y(122)=+xp(122) 
   index_Y(1,123)= 3; index_Y(2,123)= 1; index_Y(3,123)= 5; d_Y(123)=+xp(123) 
   index_Y(1,124)= 3; index_Y(2,124)= 1; index_Y(3,124)= 6; d_Y(124)=+xp(124) 
   index_Y(1,125)= 3; index_Y(2,125)= 1; index_Y(3,125)= 7; d_Y(125)=+xp(125) 
   index_Y(1,126)= 3; index_Y(2,126)= 1; index_Y(3,126)= 8; d_Y(126)=+xp(126) 
   index_Y(1,127)= 3; index_Y(2,127)= 3; index_Y(3,127)= 0; d_Y(127)=+xp(127) 
   index_Y(1,128)= 3; index_Y(2,128)= 3; index_Y(3,128)= 1; d_Y(128)=+xp(128) 
   index_Y(1,129)= 3; index_Y(2,129)= 3; index_Y(3,129)= 2; d_Y(129)=+xp(129) 
   index_Y(1,130)= 3; index_Y(2,130)= 3; index_Y(3,130)= 3; d_Y(130)=+xp(130) 
   index_Y(1,131)= 3; index_Y(2,131)= 3; index_Y(3,131)= 4; d_Y(131)=+xp(131) 
   index_Y(1,132)= 3; index_Y(2,132)= 3; index_Y(3,132)= 5; d_Y(132)=+xp(132) 
   index_Y(1,133)= 3; index_Y(2,133)= 3; index_Y(3,133)= 6; d_Y(133)=+xp(133) 
   index_Y(1,134)= 3; index_Y(2,134)= 3; index_Y(3,134)= 7; d_Y(134)=+xp(134) 
   index_Y(1,135)= 3; index_Y(2,135)= 3; index_Y(3,135)= 8; d_Y(135)=+xp(135) 
   index_Y(1,136)= 3; index_Y(2,136)= 5; index_Y(3,136)= 0; d_Y(136)=+xp(136) 
   index_Y(1,137)= 3; index_Y(2,137)= 5; index_Y(3,137)= 1; d_Y(137)=+xp(137) 
   index_Y(1,138)= 3; index_Y(2,138)= 5; index_Y(3,138)= 2; d_Y(138)=+xp(138) 
   index_Y(1,139)= 3; index_Y(2,139)= 5; index_Y(3,139)= 3; d_Y(139)=+xp(139) 
   index_Y(1,140)= 3; index_Y(2,140)= 5; index_Y(3,140)= 4; d_Y(140)=+xp(140) 
   index_Y(1,141)= 3; index_Y(2,141)= 5; index_Y(3,141)= 5; d_Y(141)=+xp(141) 
   index_Y(1,142)= 3; index_Y(2,142)= 5; index_Y(3,142)= 6; d_Y(142)=+xp(142) 
   index_Y(1,143)= 3; index_Y(2,143)= 5; index_Y(3,143)= 7; d_Y(143)=+xp(143) 
   index_Y(1,144)= 3; index_Y(2,144)= 5; index_Y(3,144)= 8; d_Y(144)=+xp(144) 
   index_Y(1,145)= 4; index_Y(2,145)= 1; index_Y(3,145)= 0; d_Y(145)=+xp(145) 
   index_Y(1,146)= 4; index_Y(2,146)= 1; index_Y(3,146)= 1; d_Y(146)=+xp(146) 
   index_Y(1,147)= 4; index_Y(2,147)= 1; index_Y(3,147)= 2; d_Y(147)=+xp(147) 
   index_Y(1,148)= 4; index_Y(2,148)= 1; index_Y(3,148)= 3; d_Y(148)=+xp(148) 
   index_Y(1,149)= 4; index_Y(2,149)= 1; index_Y(3,149)= 4; d_Y(149)=+xp(149) 
   index_Y(1,150)= 4; index_Y(2,150)= 1; index_Y(3,150)= 5; d_Y(150)=+xp(150) 
   index_Y(1,151)= 4; index_Y(2,151)= 1; index_Y(3,151)= 6; d_Y(151)=+xp(151) 
   index_Y(1,152)= 4; index_Y(2,152)= 1; index_Y(3,152)= 7; d_Y(152)=+xp(152) 
   index_Y(1,153)= 4; index_Y(2,153)= 1; index_Y(3,153)= 8; d_Y(153)=+xp(153) 
   index_Y(1,154)= 4; index_Y(2,154)= 3; index_Y(3,154)= 0; d_Y(154)=+xp(154) 
   index_Y(1,155)= 4; index_Y(2,155)= 3; index_Y(3,155)= 1; d_Y(155)=+xp(155) 
   index_Y(1,156)= 4; index_Y(2,156)= 3; index_Y(3,156)= 2; d_Y(156)=+xp(156) 
   index_Y(1,157)= 4; index_Y(2,157)= 3; index_Y(3,157)= 3; d_Y(157)=+xp(157) 
   index_Y(1,158)= 4; index_Y(2,158)= 3; index_Y(3,158)= 4; d_Y(158)=+xp(158) 
   index_Y(1,159)= 4; index_Y(2,159)= 3; index_Y(3,159)= 5; d_Y(159)=+xp(159) 
   index_Y(1,160)= 4; index_Y(2,160)= 3; index_Y(3,160)= 6; d_Y(160)=+xp(160) 
   index_Y(1,161)= 4; index_Y(2,161)= 3; index_Y(3,161)= 7; d_Y(161)=+xp(161) 
   index_Y(1,162)= 4; index_Y(2,162)= 3; index_Y(3,162)= 8; d_Y(162)=+xp(162) 
   index_Y(1,163)= 4; index_Y(2,163)= 5; index_Y(3,163)= 0; d_Y(163)=+xp(163) 
   index_Y(1,164)= 4; index_Y(2,164)= 5; index_Y(3,164)= 1; d_Y(164)=+xp(164) 
   index_Y(1,165)= 4; index_Y(2,165)= 5; index_Y(3,165)= 2; d_Y(165)=+xp(165) 
   index_Y(1,166)= 4; index_Y(2,166)= 5; index_Y(3,166)= 3; d_Y(166)=+xp(166) 
   index_Y(1,167)= 4; index_Y(2,167)= 5; index_Y(3,167)= 4; d_Y(167)=+xp(167) 
   index_Y(1,168)= 4; index_Y(2,168)= 5; index_Y(3,168)= 5; d_Y(168)=+xp(168) 
   index_Y(1,169)= 4; index_Y(2,169)= 5; index_Y(3,169)= 6; d_Y(169)=+xp(169) 
   index_Y(1,170)= 4; index_Y(2,170)= 5; index_Y(3,170)= 7; d_Y(170)=+xp(170) 
   index_Y(1,171)= 4; index_Y(2,171)= 5; index_Y(3,171)= 8; d_Y(171)=+xp(171) 
   index_Y(1,172)= 5; index_Y(2,172)= 1; index_Y(3,172)= 0; d_Y(172)=+xp(172) 
   index_Y(1,173)= 5; index_Y(2,173)= 1; index_Y(3,173)= 1; d_Y(173)=+xp(173) 
   index_Y(1,174)= 5; index_Y(2,174)= 1; index_Y(3,174)= 2; d_Y(174)=+xp(174) 
   index_Y(1,175)= 5; index_Y(2,175)= 1; index_Y(3,175)= 3; d_Y(175)=+xp(175) 
   index_Y(1,176)= 5; index_Y(2,176)= 1; index_Y(3,176)= 4; d_Y(176)=+xp(176) 
   index_Y(1,177)= 5; index_Y(2,177)= 1; index_Y(3,177)= 5; d_Y(177)=+xp(177) 
   index_Y(1,178)= 5; index_Y(2,178)= 1; index_Y(3,178)= 6; d_Y(178)=+xp(178) 
   index_Y(1,179)= 5; index_Y(2,179)= 1; index_Y(3,179)= 7; d_Y(179)=+xp(179) 
   index_Y(1,180)= 5; index_Y(2,180)= 1; index_Y(3,180)= 8; d_Y(180)=+xp(180) 
   index_Y(1,181)= 5; index_Y(2,181)= 3; index_Y(3,181)= 0; d_Y(181)=+xp(181) 
   index_Y(1,182)= 5; index_Y(2,182)= 3; index_Y(3,182)= 1; d_Y(182)=+xp(182) 
   index_Y(1,183)= 5; index_Y(2,183)= 3; index_Y(3,183)= 2; d_Y(183)=+xp(183) 
   index_Y(1,184)= 5; index_Y(2,184)= 3; index_Y(3,184)= 3; d_Y(184)=+xp(184) 
   index_Y(1,185)= 5; index_Y(2,185)= 3; index_Y(3,185)= 4; d_Y(185)=+xp(185) 
   index_Y(1,186)= 5; index_Y(2,186)= 3; index_Y(3,186)= 5; d_Y(186)=+xp(186) 
   index_Y(1,187)= 5; index_Y(2,187)= 3; index_Y(3,187)= 6; d_Y(187)=+xp(187) 
   index_Y(1,188)= 5; index_Y(2,188)= 3; index_Y(3,188)= 7; d_Y(188)=+xp(188) 
   index_Y(1,189)= 5; index_Y(2,189)= 3; index_Y(3,189)= 8; d_Y(189)=+xp(189) 
   index_Y(1,190)= 6; index_Y(2,190)= 1; index_Y(3,190)= 0; d_Y(190)=+xp(190) 
   index_Y(1,191)= 6; index_Y(2,191)= 1; index_Y(3,191)= 1; d_Y(191)=+xp(191) 
   index_Y(1,192)= 6; index_Y(2,192)= 1; index_Y(3,192)= 2; d_Y(192)=+xp(192) 
   index_Y(1,193)= 6; index_Y(2,193)= 1; index_Y(3,193)= 3; d_Y(193)=+xp(193) 
   index_Y(1,194)= 6; index_Y(2,194)= 1; index_Y(3,194)= 4; d_Y(194)=+xp(194) 
   index_Y(1,195)= 6; index_Y(2,195)= 1; index_Y(3,195)= 5; d_Y(195)=+xp(195) 
   index_Y(1,196)= 6; index_Y(2,196)= 1; index_Y(3,196)= 6; d_Y(196)=+xp(196) 
   index_Y(1,197)= 6; index_Y(2,197)= 1; index_Y(3,197)= 7; d_Y(197)=+xp(197) 
   index_Y(1,198)= 6; index_Y(2,198)= 1; index_Y(3,198)= 8; d_Y(198)=+xp(198) 
   index_Y(1,199)= 6; index_Y(2,199)= 3; index_Y(3,199)= 0; d_Y(199)=+xp(199) 
   index_Y(1,200)= 6; index_Y(2,200)= 3; index_Y(3,200)= 1; d_Y(200)=+xp(200) 
   index_Y(1,201)= 6; index_Y(2,201)= 3; index_Y(3,201)= 2; d_Y(201)=+xp(201)
   index_Y(1,202)= 6; index_Y(2,202)= 3; index_Y(3,202)= 3; d_Y(202)=+xp(202)
   index_Y(1,203)= 6; index_Y(2,203)= 3; index_Y(3,203)= 4; d_Y(203)=+xp(203)
   index_Y(1,204)= 6; index_Y(2,204)= 3; index_Y(3,204)= 5; d_Y(204)=+xp(204)
   index_Y(1,205)= 6; index_Y(2,205)= 3; index_Y(3,205)= 6; d_Y(205)=+xp(205)
   index_Y(1,206)= 6; index_Y(2,206)= 3; index_Y(3,206)= 7; d_Y(206)=+xp(206)
   index_Y(1,207)= 6; index_Y(2,207)= 3; index_Y(3,207)= 8; d_Y(207)=+xp(207)
   index_Y(1,208)= 7; index_Y(2,208)= 1; index_Y(3,208)= 0; d_Y(208)=+xp(208)
   index_Y(1,209)= 7; index_Y(2,209)= 1; index_Y(3,209)= 1; d_Y(209)=+xp(209)
   index_Y(1,210)= 7; index_Y(2,210)= 1; index_Y(3,210)= 2; d_Y(210)=+xp(210)
   index_Y(1,211)= 7; index_Y(2,211)= 1; index_Y(3,211)= 3; d_Y(211)=+xp(211)
   index_Y(1,212)= 7; index_Y(2,212)= 1; index_Y(3,212)= 4; d_Y(212)=+xp(212)
   index_Y(1,213)= 7; index_Y(2,213)= 1; index_Y(3,213)= 5; d_Y(213)=+xp(213)
   index_Y(1,214)= 7; index_Y(2,214)= 1; index_Y(3,214)= 6; d_Y(214)=+xp(214)
   index_Y(1,215)= 7; index_Y(2,215)= 1; index_Y(3,215)= 7; d_Y(215)=+xp(215)
   index_Y(1,216)= 7; index_Y(2,216)= 1; index_Y(3,216)= 8; d_Y(216)=+xp(216)
   index_Y(1,217)= 8; index_Y(2,217)= 1; index_Y(3,217)= 0; d_Y(217)=+xp(217)
   index_Y(1,218)= 8; index_Y(2,218)= 1; index_Y(3,218)= 1; d_Y(218)=+xp(218)
   index_Y(1,219)= 8; index_Y(2,219)= 1; index_Y(3,219)= 2; d_Y(219)=+xp(219)
   index_Y(1,220)= 8; index_Y(2,220)= 1; index_Y(3,220)= 3; d_Y(220)=+xp(220)
   index_Y(1,221)= 8; index_Y(2,221)= 1; index_Y(3,221)= 4; d_Y(221)=+xp(221)
   index_Y(1,222)= 8; index_Y(2,222)= 1; index_Y(3,222)= 5; d_Y(222)=+xp(222)
   index_Y(1,223)= 8; index_Y(2,223)= 1; index_Y(3,223)= 6; d_Y(223)=+xp(223)
   index_Y(1,224)= 8; index_Y(2,224)= 1; index_Y(3,224)= 7; d_Y(224)=+xp(224)
   index_Y(1,225)= 8; index_Y(2,225)= 1; index_Y(3,225)= 8; d_Y(225)=+xp(225)

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
!theta = cos_theta!acos(cos_theta)
!cos_theta is in radians
theta = dcos(cos_theta)!cosine of radians
!***************
x1 = r1 + r2
x2 = r2 - r1
x3 = PI - cos_theta
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
muY = 0._dp
do ii = 1, ncoeffs_Y
  muY = muY + d_Y(ii) * powers_x1(index_Y(1,ii)) * powers_x2(index_Y(2,ii)) * powers_x3(index_Y(3,ii))
enddo
!***************

!add damping/asymptotic terms
damp1 = damp(r1)
damp2 = damp(r2)
value1 = muOH(r1)*( 1._dp - damp2 )
value2 = muOH(r2)*( 1._dp - damp1 )

!print *, cos_theta,1-cos_theta

cos_theta_half = sqrt( (1._dp + theta)/2._dp  )
sin_theta_half = sqrt( (1._dp - theta)/2._dp  )

muY = muY*damp1*damp2 + sin_theta_half*( value1 - value2  )


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


