C/*-----------------------------------------------------
C// Neutrino flux table at 1000 km
C// Off-axis angle is 2.3 deg. for anti-nu_mu focusing beam
C// Enu(GeV), nu_mu, bar{nu_mu}, nu_e, bar{nu_e}
C-------------------------------------------------------*
      real*8 oab23a(300,5)
      data (oab23a(1,oi),oi=1,5) / 0.2500000E-01, 0.2272494E+04, 0.1427512E+04, 0.4067150E+01, 0.1567849E+02 /
      data (oab23a(2,oi),oi=1,5) / 0.7500000E-01, 0.6661004E+04, 0.7379964E+04, 0.1633067E+02, 0.1042185E+03 /
      data (oab23a(3,oi),oi=1,5) / 0.1250000E+00, 0.4453189E+04, 0.1530403E+05, 0.2795502E+02, 0.2380326E+03 /
      data (oab23a(4,oi),oi=1,5) / 0.1750000E+00, 0.3632902E+04, 0.2289867E+05, 0.3901501E+02, 0.3227491E+03 /
      data (oab23a(5,oi),oi=1,5) / 0.2250000E+00, 0.3501984E+04, 0.3260335E+05, 0.4973846E+02, 0.3988708E+03 /
      data (oab23a(6,oi),oi=1,5) / 0.2750000E+00, 0.3624739E+04, 0.4440387E+05, 0.5805918E+02, 0.4937844E+03 /
      data (oab23a(7,oi),oi=1,5) / 0.3250000E+00, 0.3864882E+04, 0.5692340E+05, 0.6072941E+02, 0.5052600E+03 /
      data (oab23a(8,oi),oi=1,5) / 0.3750000E+00, 0.3988429E+04, 0.6923054E+05, 0.6673694E+02, 0.5677686E+03 /
      data (oab23a(9,oi),oi=1,5) / 0.4250000E+00, 0.4135977E+04, 0.8140297E+05, 0.7076378E+02, 0.5560311E+03 /
      data (oab23a(10,oi),oi=1,5) / 0.4750000E+00, 0.4153602E+04, 0.9441521E+05, 0.6882989E+02, 0.5385109E+03 /
      data (oab23a(11,oi),oi=1,5) / 0.5250000E+00, 0.4391741E+04, 0.1092786E+06, 0.7073100E+02, 0.5781679E+03 /
      data (oab23a(12,oi),oi=1,5) / 0.5750000E+00, 0.4438520E+04, 0.1276725E+06, 0.7672573E+02, 0.5685948E+03 /
      data (oab23a(13,oi),oi=1,5) / 0.6250000E+00, 0.4435842E+04, 0.1388850E+06, 0.7926210E+02, 0.5107936E+03 /
      data (oab23a(14,oi),oi=1,5) / 0.6750000E+00, 0.4394840E+04, 0.1325403E+06, 0.6989678E+02, 0.4853602E+03 /
      data (oab23a(15,oi),oi=1,5) / 0.7250000E+00, 0.4546237E+04, 0.1137439E+06, 0.6546587E+02, 0.4831289E+03 /
      data (oab23a(16,oi),oi=1,5) / 0.7750000E+00, 0.4232780E+04, 0.9384086E+05, 0.7267843E+02, 0.4435017E+03 /
      data (oab23a(17,oi),oi=1,5) / 0.8250000E+00, 0.4007377E+04, 0.7884541E+05, 0.6475566E+02, 0.3500085E+03 /
      data (oab23a(18,oi),oi=1,5) / 0.8750000E+00, 0.3784522E+04, 0.6367958E+05, 0.6902974E+02, 0.3843553E+03 /
      data (oab23a(19,oi),oi=1,5) / 0.9250000E+00, 0.3774056E+04, 0.4878728E+05, 0.5919418E+02, 0.3391929E+03 /
      data (oab23a(20,oi),oi=1,5) / 0.9750000E+00, 0.3532816E+04, 0.3658025E+05, 0.6732345E+02, 0.3546200E+03 /
      data (oab23a(21,oi),oi=1,5) / 0.1025000E+01, 0.3487181E+04, 0.2712346E+05, 0.5960320E+02, 0.3207710E+03 /
      data (oab23a(22,oi),oi=1,5) / 0.1075000E+01, 0.3133769E+04, 0.2035451E+05, 0.5799726E+02, 0.2421412E+03 /
      data (oab23a(23,oi),oi=1,5) / 0.1125000E+01, 0.3168382E+04, 0.1559549E+05, 0.6115972E+02, 0.2409011E+03 /
      data (oab23a(24,oi),oi=1,5) / 0.1175000E+01, 0.2972234E+04, 0.1218721E+05, 0.7458992E+02, 0.2532210E+03 /
      data (oab23a(25,oi),oi=1,5) / 0.1225000E+01, 0.2778857E+04, 0.9617888E+04, 0.5327583E+02, 0.1776878E+03 /
      data (oab23a(26,oi),oi=1,5) / 0.1275000E+01, 0.2568050E+04, 0.7716423E+04, 0.5058818E+02, 0.1448740E+03 /
      data (oab23a(27,oi),oi=1,5) / 0.1325000E+01, 0.2549859E+04, 0.6545165E+04, 0.5534229E+02, 0.1644611E+03 /
      data (oab23a(28,oi),oi=1,5) / 0.1375000E+01, 0.2631745E+04, 0.5705115E+04, 0.5000788E+02, 0.1402886E+03 /
      data (oab23a(29,oi),oi=1,5) / 0.1425000E+01, 0.2290783E+04, 0.4872247E+04, 0.4709473E+02, 0.1318131E+03 /
      data (oab23a(30,oi),oi=1,5) / 0.1475000E+01, 0.2234116E+04, 0.4180376E+04, 0.5115573E+02, 0.1498778E+03 /
      data (oab23a(31,oi),oi=1,5) / 0.1525000E+01, 0.2319600E+04, 0.3708072E+04, 0.5589729E+02, 0.1231315E+03 /
      data (oab23a(32,oi),oi=1,5) / 0.1575000E+01, 0.2201238E+04, 0.3462650E+04, 0.5208450E+02, 0.1184803E+03 /
      data (oab23a(33,oi),oi=1,5) / 0.1625000E+01, 0.1799368E+04, 0.3218143E+04, 0.4962585E+02, 0.9493902E+02 /
      data (oab23a(34,oi),oi=1,5) / 0.1675000E+01, 0.1519030E+04, 0.2630474E+04, 0.5279267E+02, 0.1085450E+03 /
      data (oab23a(35,oi),oi=1,5) / 0.1725000E+01, 0.1503290E+04, 0.2659590E+04, 0.5946003E+02, 0.9642793E+02 /
      data (oab23a(36,oi),oi=1,5) / 0.1775000E+01, 0.1678727E+04, 0.2566742E+04, 0.4655968E+02, 0.8373260E+02 /
      data (oab23a(37,oi),oi=1,5) / 0.1825000E+01, 0.1571368E+04, 0.2145501E+04, 0.4263981E+02, 0.7131949E+02 /
      data (oab23a(38,oi),oi=1,5) / 0.1875000E+01, 0.1598228E+04, 0.1904668E+04, 0.4329712E+02, 0.8238255E+02 /
      data (oab23a(39,oi),oi=1,5) / 0.1925000E+01, 0.1428777E+04, 0.1689963E+04, 0.4146578E+02, 0.7246562E+02 /
      data (oab23a(40,oi),oi=1,5) / 0.1975000E+01, 0.1436516E+04, 0.1765334E+04, 0.5389365E+02, 0.7362229E+02 /
      data (oab23a(41,oi),oi=1,5) / 0.2025000E+01, 0.1313794E+04, 0.1517564E+04, 0.4608866E+02, 0.9137145E+02 /
      data (oab23a(42,oi),oi=1,5) / 0.2075000E+01, 0.1613520E+04, 0.1537797E+04, 0.4561788E+02, 0.7075756E+02 /
      data (oab23a(43,oi),oi=1,5) / 0.2125000E+01, 0.1174191E+04, 0.1359686E+04, 0.5102536E+02, 0.7936016E+02 /
      data (oab23a(44,oi),oi=1,5) / 0.2175000E+01, 0.1231461E+04, 0.1281513E+04, 0.3601154E+02, 0.6494575E+02 /
      data (oab23a(45,oi),oi=1,5) / 0.2225000E+01, 0.1142929E+04, 0.1197134E+04, 0.3383991E+02, 0.6577424E+02 /
      data (oab23a(46,oi),oi=1,5) / 0.2275000E+01, 0.8729816E+03, 0.1349805E+04, 0.3867523E+02, 0.6631066E+02 /
      data (oab23a(47,oi),oi=1,5) / 0.2325000E+01, 0.9620885E+03, 0.1118354E+04, 0.4421553E+02, 0.6173801E+02 /
      data (oab23a(48,oi),oi=1,5) / 0.2375000E+01, 0.1228129E+04, 0.9165787E+03, 0.4064490E+02, 0.6618898E+02 /
      data (oab23a(49,oi),oi=1,5) / 0.2425000E+01, 0.9323180E+03, 0.1015139E+04, 0.3343372E+02, 0.5506008E+02 /
      data (oab23a(50,oi),oi=1,5) / 0.2475000E+01, 0.7538932E+03, 0.1008119E+04, 0.3998028E+02, 0.5764668E+02 /
      data (oab23a(51,oi),oi=1,5) / 0.2525000E+01, 0.7217544E+03, 0.9892656E+03, 0.3597905E+02, 0.5660242E+02 /
      data (oab23a(52,oi),oi=1,5) / 0.2575000E+01, 0.1041687E+04, 0.1001128E+04, 0.5058412E+02, 0.6870188E+02 /
      data (oab23a(53,oi),oi=1,5) / 0.2625000E+01, 0.8524692E+03, 0.7722435E+03, 0.3106560E+02, 0.5041685E+02 /
      data (oab23a(54,oi),oi=1,5) / 0.2675000E+01, 0.6191868E+03, 0.8052379E+03, 0.3167125E+02, 0.4809278E+02 /
      data (oab23a(55,oi),oi=1,5) / 0.2725000E+01, 0.9215569E+03, 0.7112730E+03, 0.3950115E+02, 0.5630985E+02 /
      data (oab23a(56,oi),oi=1,5) / 0.2775000E+01, 0.8136641E+03, 0.7128776E+03, 0.3381687E+02, 0.4447126E+02 /
      data (oab23a(57,oi),oi=1,5) / 0.2825000E+01, 0.7825109E+03, 0.7303279E+03, 0.2961143E+02, 0.4115511E+02 /
      data (oab23a(58,oi),oi=1,5) / 0.2875000E+01, 0.6858547E+03, 0.5494213E+03, 0.3186280E+02, 0.4449167E+02 /
      data (oab23a(59,oi),oi=1,5) / 0.2925000E+01, 0.5255838E+03, 0.7016171E+03, 0.2953781E+02, 0.4423482E+02 /
      data (oab23a(60,oi),oi=1,5) / 0.2975000E+01, 0.5640757E+03, 0.7194802E+03, 0.2813533E+02, 0.4102601E+02 /
      data (oab23a(61,oi),oi=1,5) / 0.3025000E+01, 0.6617826E+03, 0.6702088E+03, 0.3328701E+02, 0.4616931E+02 /
      data (oab23a(62,oi),oi=1,5) / 0.3075000E+01, 0.5860769E+03, 0.9538921E+03, 0.2854997E+02, 0.3891004E+02 /
      data (oab23a(63,oi),oi=1,5) / 0.3125000E+01, 0.8441940E+03, 0.4649981E+03, 0.2931628E+02, 0.4318315E+02 /
      data (oab23a(64,oi),oi=1,5) / 0.3175000E+01, 0.6287080E+03, 0.7319430E+03, 0.4879466E+02, 0.4623818E+02 /
      data (oab23a(65,oi),oi=1,5) / 0.3225000E+01, 0.4497620E+03, 0.6296363E+03, 0.3014691E+02, 0.4132180E+02 /
      data (oab23a(66,oi),oi=1,5) / 0.3275000E+01, 0.7462478E+03, 0.7026593E+03, 0.3788932E+02, 0.4435075E+02 /
      data (oab23a(67,oi),oi=1,5) / 0.3325000E+01, 0.6102305E+03, 0.5138965E+03, 0.2250935E+02, 0.3165540E+02 /
      data (oab23a(68,oi),oi=1,5) / 0.3375000E+01, 0.4349159E+03, 0.6267437E+03, 0.2726116E+02, 0.3471105E+02 /
      data (oab23a(69,oi),oi=1,5) / 0.3425000E+01, 0.3630156E+03, 0.6664876E+03, 0.2906470E+02, 0.4417700E+02 /
      data (oab23a(70,oi),oi=1,5) / 0.3475000E+01, 0.3599201E+03, 0.5466796E+03, 0.2333711E+02, 0.3468612E+02 /
      data (oab23a(71,oi),oi=1,5) / 0.3525000E+01, 0.7116003E+03, 0.6049808E+03, 0.2317715E+02, 0.3123705E+02 /
      data (oab23a(72,oi),oi=1,5) / 0.3575000E+01, 0.2962486E+03, 0.5753127E+03, 0.4080868E+02, 0.3291350E+02 /
      data (oab23a(73,oi),oi=1,5) / 0.3625000E+01, 0.4007325E+03, 0.4123209E+03, 0.2766419E+02, 0.3258060E+02 /
      data (oab23a(74,oi),oi=1,5) / 0.3675000E+01, 0.3207969E+03, 0.5140878E+03, 0.1612474E+02, 0.2324579E+02 /
      data (oab23a(75,oi),oi=1,5) / 0.3725000E+01, 0.3574339E+03, 0.6040859E+03, 0.1841630E+02, 0.2289914E+02 /
      data (oab23a(76,oi),oi=1,5) / 0.3775000E+01, 0.4501667E+03, 0.4921981E+03, 0.1807370E+02, 0.2189534E+02 /
      data (oab23a(77,oi),oi=1,5) / 0.3825000E+01, 0.2828122E+03, 0.4503603E+03, 0.2952050E+02, 0.3242660E+02 /
      data (oab23a(78,oi),oi=1,5) / 0.3875000E+01, 0.5775602E+03, 0.5831017E+03, 0.2406061E+02, 0.2490874E+02 /
      data (oab23a(79,oi),oi=1,5) / 0.3925000E+01, 0.3543093E+03, 0.3966169E+03, 0.2323373E+02, 0.2321240E+02 /
      data (oab23a(80,oi),oi=1,5) / 0.3975000E+01, 0.1509533E+03, 0.3994618E+03, 0.2027602E+02, 0.2310367E+02 /
      data (oab23a(81,oi),oi=1,5) / 0.4025000E+01, 0.3702803E+03, 0.4853477E+03, 0.2397093E+02, 0.2146048E+02 /
      data (oab23a(82,oi),oi=1,5) / 0.4075000E+01, 0.2608558E+03, 0.4393660E+03, 0.2221367E+02, 0.2393253E+02 /
      data (oab23a(83,oi),oi=1,5) / 0.4125000E+01, 0.4680160E+03, 0.5138404E+03, 0.2116061E+02, 0.2153994E+02 /
      data (oab23a(84,oi),oi=1,5) / 0.4175000E+01, 0.3539931E+03, 0.5249521E+03, 0.1528872E+02, 0.1868442E+02 /
      data (oab23a(85,oi),oi=1,5) / 0.4225000E+01, 0.2729001E+03, 0.4232255E+03, 0.2272669E+02, 0.2692778E+02 /
      data (oab23a(86,oi),oi=1,5) / 0.4275000E+01, 0.4151220E+03, 0.3638153E+03, 0.1947203E+02, 0.1787956E+02 /
      data (oab23a(87,oi),oi=1,5) / 0.4325000E+01, 0.6752990E+03, 0.4474738E+03, 0.1274586E+02, 0.1544061E+02 /
      data (oab23a(88,oi),oi=1,5) / 0.4375000E+01, 0.3694874E+03, 0.4366998E+03, 0.2094656E+02, 0.2088498E+02 /
      data (oab23a(89,oi),oi=1,5) / 0.4425000E+01, 0.4322847E+03, 0.3545191E+03, 0.1437929E+02, 0.1759935E+02 /
      data (oab23a(90,oi),oi=1,5) / 0.4475000E+01, 0.3273586E+03, 0.4481287E+03, 0.1155180E+02, 0.1311659E+02 /
      data (oab23a(91,oi),oi=1,5) / 0.4525000E+01, 0.4872922E+03, 0.3458086E+03, 0.1745872E+02, 0.1573449E+02 /
      data (oab23a(92,oi),oi=1,5) / 0.4575000E+01, 0.8798606E+02, 0.3597303E+03, 0.9752135E+01, 0.9569557E+01 /
      data (oab23a(93,oi),oi=1,5) / 0.4625000E+01, 0.4178306E+03, 0.3409008E+03, 0.1659740E+02, 0.1146720E+02 /
      data (oab23a(94,oi),oi=1,5) / 0.4675000E+01, 0.2195709E+03, 0.3714457E+03, 0.1407429E+02, 0.1291900E+02 /
      data (oab23a(95,oi),oi=1,5) / 0.4725000E+01, 0.2679369E+03, 0.3005158E+03, 0.1755051E+02, 0.1569733E+02 /
      data (oab23a(96,oi),oi=1,5) / 0.4775000E+01, 0.8442319E+02, 0.2557901E+03, 0.1103412E+02, 0.9512788E+01 /
      data (oab23a(97,oi),oi=1,5) / 0.4825000E+01, 0.1183399E+03, 0.3940980E+03, 0.1191302E+02, 0.1199079E+02 /
      data (oab23a(98,oi),oi=1,5) / 0.4875000E+01, 0.2166915E+03, 0.2623988E+03, 0.1306047E+02, 0.1190099E+02 /
      data (oab23a(99,oi),oi=1,5) / 0.4925000E+01, 0.1154491E+03, 0.3037491E+03, 0.1783684E+02, 0.1145226E+02 /
      data (oab23a(100,oi),oi=1,5) / 0.4975000E+01, 0.2701507E+03, 0.2819786E+03, 0.1606937E+02, 0.1183866E+02 /
      data (oab23a(101,oi),oi=1,5) / 0.5025000E+01, 0.2163878E+03, 0.2616739E+03, 0.1455809E+02, 0.1203936E+02 /
      data (oab23a(102,oi),oi=1,5) / 0.5075000E+01, 0.1052044E+03, 0.3196707E+03, 0.9129355E+01, 0.7995713E+01 /
      data (oab23a(103,oi),oi=1,5) / 0.5125000E+01, 0.1272405E+03, 0.2476335E+03, 0.6509817E+01, 0.7444425E+01 /
      data (oab23a(104,oi),oi=1,5) / 0.5175000E+01, 0.8895573E+02, 0.3570842E+03, 0.1317297E+02, 0.1161152E+02 /
      data (oab23a(105,oi),oi=1,5) / 0.5225000E+01, 0.1546794E+03, 0.1900282E+03, 0.9616251E+01, 0.7334017E+01 /
      data (oab23a(106,oi),oi=1,5) / 0.5275000E+01, 0.8827828E+02, 0.2100859E+03, 0.1408910E+02, 0.1306667E+02 /
      data (oab23a(107,oi),oi=1,5) / 0.5325000E+01, 0.4334795E+03, 0.1884312E+03, 0.8969620E+01, 0.5321748E+01 /
      data (oab23a(108,oi),oi=1,5) / 0.5375000E+01, 0.9170897E+02, 0.1826933E+03, 0.1733751E+02, 0.1675937E+02 /
      data (oab23a(109,oi),oi=1,5) / 0.5425000E+01, 0.2408284E+03, 0.1977820E+03, 0.7780325E+01, 0.5353585E+01 /
      data (oab23a(110,oi),oi=1,5) / 0.5475000E+01, 0.9831573E+02, 0.1794235E+03, 0.1121111E+02, 0.8691399E+01 /
      data (oab23a(111,oi),oi=1,5) / 0.5525000E+01, 0.2641590E+03, 0.2008966E+03, 0.7656675E+01, 0.5001439E+01 /
      data (oab23a(112,oi),oi=1,5) / 0.5575000E+01, 0.1898943E+03, 0.1892457E+03, 0.6587331E+01, 0.4088040E+01 /
      data (oab23a(113,oi),oi=1,5) / 0.5625000E+01, 0.1124800E+03, 0.1975519E+03, 0.9428531E+01, 0.7462829E+01 /
      data (oab23a(114,oi),oi=1,5) / 0.5675000E+01, 0.1948796E+03, 0.1485228E+03, 0.1819024E+02, 0.1493971E+02 /
      data (oab23a(115,oi),oi=1,5) / 0.5725000E+01, 0.7495311E+02, 0.1593937E+03, 0.7151326E+01, 0.6133707E+01 /
      data (oab23a(116,oi),oi=1,5) / 0.5775000E+01, 0.9028636E+02, 0.1538588E+03, 0.1063994E+02, 0.7836826E+01 /
      data (oab23a(117,oi),oi=1,5) / 0.5825000E+01, 0.7116173E+02, 0.1275726E+03, 0.1378111E+02, 0.8060252E+01 /
      data (oab23a(118,oi),oi=1,5) / 0.5875000E+01, 0.8181016E+02, 0.1379935E+03, 0.1470101E+02, 0.1501155E+02 /
      data (oab23a(119,oi),oi=1,5) / 0.5925000E+01, 0.1941379E+03, 0.1225687E+03, 0.3486899E+01, 0.3203714E+01 /
      data (oab23a(120,oi),oi=1,5) / 0.5975000E+01, 0.8232356E+02, 0.1322325E+03, 0.4418652E+01, 0.2899295E+01 /
      data (oab23a(121,oi),oi=1,5) / 0.6025000E+01, 0.9298823E+02, 0.1204503E+03, 0.5572413E+01, 0.3853077E+01 /
      data (oab23a(122,oi),oi=1,5) / 0.6075000E+01, 0.7020931E+02, 0.1201960E+03, 0.1603940E+02, 0.1277621E+02 /
      data (oab23a(123,oi),oi=1,5) / 0.6125000E+01, 0.1819167E+03, 0.1269364E+03, 0.7712647E+01, 0.3360583E+01 /
      data (oab23a(124,oi),oi=1,5) / 0.6175000E+01, 0.7100435E+02, 0.1258972E+03, 0.5436451E+01, 0.3204702E+01 /
      data (oab23a(125,oi),oi=1,5) / 0.6225000E+01, 0.8776857E+02, 0.7709408E+02, 0.9764536E+01, 0.7119602E+01 /
      data (oab23a(126,oi),oi=1,5) / 0.6275000E+01, 0.8821363E+02, 0.2663282E+03, 0.2513268E+01, 0.1204408E+01 /
      data (oab23a(127,oi),oi=1,5) / 0.6325000E+01, 0.9898589E+02, 0.9971867E+02, 0.6038327E+01, 0.3532992E+01 /
      data (oab23a(128,oi),oi=1,5) / 0.6375000E+01, 0.7745529E+02, 0.9714801E+02, 0.6935168E+01, 0.4916983E+01 /
      data (oab23a(129,oi),oi=1,5) / 0.6425000E+01, 0.7686380E+02, 0.7600601E+02, 0.1287120E+02, 0.1174577E+02 /
      data (oab23a(130,oi),oi=1,5) / 0.6475000E+01, 0.8775706E+02, 0.1053467E+03, 0.4974392E+01, 0.5650676E+01 /
      data (oab23a(131,oi),oi=1,5) / 0.6525000E+01, 0.2757439E+03, 0.6212659E+02, 0.3741884E+01, 0.1677410E+01 /
      data (oab23a(132,oi),oi=1,5) / 0.6575000E+01, 0.9150942E+02, 0.7590613E+02, 0.1172576E+02, 0.9106334E+01 /
      data (oab23a(133,oi),oi=1,5) / 0.6625000E+01, 0.1097189E+03, 0.6675583E+02, 0.9780799E+01, 0.5845293E+01 /
      data (oab23a(134,oi),oi=1,5) / 0.6675000E+01, 0.8277703E+02, 0.9514115E+02, 0.6359079E+01, 0.4452257E+01 /
      data (oab23a(135,oi),oi=1,5) / 0.6725000E+01, 0.6795301E+02, 0.6309396E+02, 0.6244654E+01, 0.3700520E+01 /
      data (oab23a(136,oi),oi=1,5) / 0.6775000E+01, 0.6979656E+02, 0.8445041E+02, 0.1335089E+02, 0.9341256E+01 /
      data (oab23a(137,oi),oi=1,5) / 0.6825000E+01, 0.6493889E+02, 0.6798453E+02, 0.3827916E+01, 0.2131400E+01 /
      data (oab23a(138,oi),oi=1,5) / 0.6875000E+01, 0.2090107E+03, 0.5992272E+02, 0.6095472E+01, 0.3357322E+01 /
      data (oab23a(139,oi),oi=1,5) / 0.6925000E+01, 0.7308408E+02, 0.3847201E+02, 0.9804083E+01, 0.6562249E+01 /
      data (oab23a(140,oi),oi=1,5) / 0.6975000E+01, 0.4943960E+02, 0.4426333E+02, 0.6399619E+01, 0.5647045E+01 /
      data (oab23a(141,oi),oi=1,5) / 0.7025000E+01, 0.8882552E+02, 0.6682397E+02, 0.4557494E+01, 0.2483726E+01 /
      data (oab23a(142,oi),oi=1,5) / 0.7075000E+01, 0.3145869E+03, 0.5083967E+02, 0.8545288E+01, 0.5654335E+01 /
      data (oab23a(143,oi),oi=1,5) / 0.7125000E+01, 0.4403603E+02, 0.3657516E+02, 0.5530160E+01, 0.3027945E+01 /
      data (oab23a(144,oi),oi=1,5) / 0.7175000E+01, 0.8267647E+02, 0.5334288E+02, 0.2850128E+01, 0.1554858E+01 /
      data (oab23a(145,oi),oi=1,5) / 0.7225000E+01, 0.2106733E+03, 0.6835564E+02, 0.3618555E+01, 0.0000000E+00 /
      data (oab23a(146,oi),oi=1,5) / 0.7275000E+01, 0.7276029E+02, 0.3587982E+02, 0.5803178E+01, 0.3942118E+01 /
      data (oab23a(147,oi),oi=1,5) / 0.7325000E+01, 0.7838400E+02, 0.2091218E+03, 0.3047331E+01, 0.8491431E+00 /
      data (oab23a(148,oi),oi=1,5) / 0.7375000E+01, 0.6444854E+02, 0.5757228E+02, 0.6701171E+01, 0.4290043E+01 /
      data (oab23a(149,oi),oi=1,5) / 0.7425000E+01, 0.6644323E+02, 0.6301353E+02, 0.6936690E+01, 0.3590938E+01 /
      data (oab23a(150,oi),oi=1,5) / 0.7475000E+01, 0.7335781E+02, 0.4250699E+02, 0.8421788E+01, 0.6578301E+01 /
      data (oab23a(151,oi),oi=1,5) / 0.7525000E+01, 0.6106562E+02, 0.4557328E+02, 0.2307741E+01, 0.2499164E+01 /
      data (oab23a(152,oi),oi=1,5) / 0.7575000E+01, 0.7517302E+02, 0.5949384E+02, 0.4030155E+01, 0.1343981E+01 /
      data (oab23a(153,oi),oi=1,5) / 0.7625000E+01, 0.5908269E+02, 0.4305918E+02, 0.4748396E+01, 0.2510498E+01 /
      data (oab23a(154,oi),oi=1,5) / 0.7675000E+01, 0.9970403E+02, 0.5230828E+02, 0.1991964E+01, 0.3578282E+00 /
      data (oab23a(155,oi),oi=1,5) / 0.7725000E+01, 0.7184025E+02, 0.2525968E+02, 0.2046957E+01, 0.0000000E+00 /
      data (oab23a(156,oi),oi=1,5) / 0.7775000E+01, 0.5757905E+02, 0.4725501E+02, 0.3462844E+01, 0.2874560E+00 /
      data (oab23a(157,oi),oi=1,5) / 0.7825000E+01, 0.6478725E+02, 0.5055019E+02, 0.1587032E+01, 0.0000000E+00 /
      data (oab23a(158,oi),oi=1,5) / 0.7875000E+01, 0.6161862E+02, 0.3540814E+02, 0.3704289E+01, 0.1919761E+01 /
      data (oab23a(159,oi),oi=1,5) / 0.7925000E+01, 0.5047771E+02, 0.1993154E+02, 0.3346521E+01, 0.2469032E+01 /
      data (oab23a(160,oi),oi=1,5) / 0.7975000E+01, 0.3466691E+03, 0.2820184E+02, 0.1995901E+01, 0.2708953E+01 /
      data (oab23a(161,oi),oi=1,5) / 0.8025000E+01, 0.6113799E+02, 0.2998793E+02, 0.3168330E+01, 0.2689223E+01 /
      data (oab23a(162,oi),oi=1,5) / 0.8075000E+01, 0.7177876E+02, 0.1793160E+02, 0.3443328E+01, 0.1739374E+01 /
      data (oab23a(163,oi),oi=1,5) / 0.8125000E+01, 0.5867872E+02, 0.1399817E+02, 0.4811559E+01, 0.2092236E+01 /
      data (oab23a(164,oi),oi=1,5) / 0.8175000E+01, 0.5657533E+02, 0.3818298E+02, 0.2539792E+01, 0.0000000E+00 /
      data (oab23a(165,oi),oi=1,5) / 0.8225000E+01, 0.6871333E+02, 0.3142724E+02, 0.1313719E+01, 0.6783638E+00 /
      data (oab23a(166,oi),oi=1,5) / 0.8275000E+01, 0.4048043E+02, 0.2754707E+02, 0.4907990E+01, 0.3405932E+01 /
      data (oab23a(167,oi),oi=1,5) / 0.8325000E+01, 0.4981522E+02, 0.2786968E+02, 0.1171162E+01, 0.3848761E+00 /
      data (oab23a(168,oi),oi=1,5) / 0.8375000E+01, 0.6229785E+02, 0.2969240E+02, 0.1275996E+01, 0.0000000E+00 /
      data (oab23a(169,oi),oi=1,5) / 0.8425000E+01, 0.7195865E+02, 0.3901539E+02, 0.2668705E+01, 0.2675977E+01 /
      data (oab23a(170,oi),oi=1,5) / 0.8475000E+01, 0.5772918E+02, 0.3033123E+02, 0.2549020E+01, 0.1295469E+01 /
      data (oab23a(171,oi),oi=1,5) / 0.8525000E+01, 0.5838503E+02, 0.3071318E+02, 0.2804942E+01, 0.3091548E+01 /
      data (oab23a(172,oi),oi=1,5) / 0.8575000E+01, 0.2646819E+02, 0.2487476E+02, 0.1893119E+01, 0.2843538E+00 /
      data (oab23a(173,oi),oi=1,5) / 0.8625000E+01, 0.6468881E+02, 0.2671958E+02, 0.4036763E+01, 0.2044250E+01 /
      data (oab23a(174,oi),oi=1,5) / 0.8675000E+01, 0.7639541E+02, 0.7937341E+01, 0.2696654E+01, 0.2028442E+01 /
      data (oab23a(175,oi),oi=1,5) / 0.8725000E+01, 0.6769695E+02, 0.3054096E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(176,oi),oi=1,5) / 0.8775000E+01, 0.2436985E+02, 0.2279754E+02, 0.1831826E+01, 0.1540953E+01 /
      data (oab23a(177,oi),oi=1,5) / 0.8825000E+01, 0.4449837E+02, 0.2467838E+02, 0.8950277E+00, 0.3258652E+00 /
      data (oab23a(178,oi),oi=1,5) / 0.8875000E+01, 0.4835330E+02, 0.3349822E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(179,oi),oi=1,5) / 0.8925000E+01, 0.4548551E+02, 0.1011011E+02, 0.1010435E+01, 0.1035533E+01 /
      data (oab23a(180,oi),oi=1,5) / 0.8975000E+01, 0.6819138E+02, 0.3579974E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(181,oi),oi=1,5) / 0.9025000E+01, 0.3098565E+02, 0.2068594E+02, 0.4240367E+01, 0.2352072E+01 /
      data (oab23a(182,oi),oi=1,5) / 0.9075000E+01, 0.3999629E+02, 0.2091621E+02, 0.6459433E+01, 0.5913202E+01 /
      data (oab23a(183,oi),oi=1,5) / 0.9125000E+01, 0.5284180E+02, 0.1585347E+02, 0.1138640E+01, 0.4942367E+00 /
      data (oab23a(184,oi),oi=1,5) / 0.9175000E+01, 0.6053506E+02, 0.1069093E+02, 0.9707994E+00, 0.0000000E+00 /
      data (oab23a(185,oi),oi=1,5) / 0.9225000E+01, 0.4317542E+02, 0.1800713E+02, 0.3522060E+01, 0.2201978E+01 /
      data (oab23a(186,oi),oi=1,5) / 0.9275000E+01, 0.4542897E+02, 0.1274658E+02, 0.8434738E+00, 0.0000000E+00 /
      data (oab23a(187,oi),oi=1,5) / 0.9325000E+01, 0.2022318E+02, 0.2022359E+02, 0.4067725E+00, 0.0000000E+00 /
      data (oab23a(188,oi),oi=1,5) / 0.9375000E+01, 0.6700178E+02, 0.3703815E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(189,oi),oi=1,5) / 0.9425000E+01, 0.2249875E+02, 0.2255221E+02, 0.5091725E+00, 0.0000000E+00 /
      data (oab23a(190,oi),oi=1,5) / 0.9475000E+01, 0.5320565E+02, 0.1137925E+02, 0.2595191E+01, 0.0000000E+00 /
      data (oab23a(191,oi),oi=1,5) / 0.9525000E+01, 0.4414329E+02, 0.1338060E+02, 0.8649122E+00, 0.4316453E+00 /
      data (oab23a(192,oi),oi=1,5) / 0.9575000E+01, 0.3101458E+02, 0.1938704E+02, 0.2556669E+01, 0.2556669E+01 /
      data (oab23a(193,oi),oi=1,5) / 0.9625000E+01, 0.5094925E+02, 0.3938888E+01, 0.0000000E+00, 0.5388589E+00 /
      data (oab23a(194,oi),oi=1,5) / 0.9675000E+01, 0.5142639E+02, 0.2371533E+02, 0.1048385E+01, 0.0000000E+00 /
      data (oab23a(195,oi),oi=1,5) / 0.9725000E+01, 0.3600660E+02, 0.3005981E+02, 0.8052153E+01, 0.6075395E+01 /
      data (oab23a(196,oi),oi=1,5) / 0.9775000E+01, 0.4452067E+02, 0.1817342E+02, 0.3373637E+01, 0.3373637E+01 /
      data (oab23a(197,oi),oi=1,5) / 0.9825000E+01, 0.4292239E+02, 0.1021365E+02, 0.7360679E+00, 0.0000000E+00 /
      data (oab23a(198,oi),oi=1,5) / 0.9875000E+01, 0.2267519E+02, 0.4117169E+01, 0.2109138E+01, 0.0000000E+00 /
      data (oab23a(199,oi),oi=1,5) / 0.9925000E+01, 0.3954922E+02, 0.8321030E+01, 0.1017735E+01, 0.0000000E+00 /
      data (oab23a(200,oi),oi=1,5) / 0.9975000E+01, 0.2950400E+02, 0.0000000E+00, 0.1872035E+01, 0.0000000E+00 /
      data (oab23a(201,oi),oi=1,5) / 0.1002500E+02, 0.5940182E+02, 0.2120831E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(202,oi),oi=1,5) / 0.1007500E+02, 0.4300545E+02, 0.1290362E+02, 0.3680413E+01, 0.3680413E+01 /
      data (oab23a(203,oi),oi=1,5) / 0.1012500E+02, 0.3039248E+02, 0.1521363E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(204,oi),oi=1,5) / 0.1017500E+02, 0.2849575E+02, 0.4360761E+01, 0.1758412E+01, 0.0000000E+00 /
      data (oab23a(205,oi),oi=1,5) / 0.1022500E+02, 0.4201151E+02, 0.1547535E+02, 0.2835635E+01, 0.0000000E+00 /
      data (oab23a(206,oi),oi=1,5) / 0.1027500E+02, 0.5359323E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(207,oi),oi=1,5) / 0.1032500E+02, 0.4513476E+02, 0.0000000E+00, 0.7270710E+00, 0.0000000E+00 /
      data (oab23a(208,oi),oi=1,5) / 0.1037500E+02, 0.2733011E+02, 0.0000000E+00, 0.7223342E+00, 0.0000000E+00 /
      data (oab23a(209,oi),oi=1,5) / 0.1042500E+02, 0.3455541E+02, 0.0000000E+00, 0.1262428E+01, 0.0000000E+00 /
      data (oab23a(210,oi),oi=1,5) / 0.1047500E+02, 0.1624580E+02, 0.0000000E+00, 0.7319242E+00, 0.6175871E+00 /
      data (oab23a(211,oi),oi=1,5) / 0.1052500E+02, 0.2814220E+02, 0.4700901E+01, 0.7263591E+01, 0.7263591E+01 /
      data (oab23a(212,oi),oi=1,5) / 0.1057500E+02, 0.2839955E+02, 0.4736185E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(213,oi),oi=1,5) / 0.1062500E+02, 0.8357262E+02, 0.4790153E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(214,oi),oi=1,5) / 0.1067500E+02, 0.1204384E+02, 0.0000000E+00, 0.4067189E+01, 0.2735793E+01 /
      data (oab23a(215,oi),oi=1,5) / 0.1072500E+02, 0.1952436E+02, 0.7270121E+01, 0.7860950E+00, 0.0000000E+00 /
      data (oab23a(216,oi),oi=1,5) / 0.1077500E+02, 0.4908047E+02, 0.0000000E+00, 0.5228421E+00, 0.0000000E+00 /
      data (oab23a(217,oi),oi=1,5) / 0.1082500E+02, 0.7415141E+01, 0.2230989E+02, 0.5683377E+00, 0.0000000E+00 /
      data (oab23a(218,oi),oi=1,5) / 0.1087500E+02, 0.1504612E+02, 0.5017205E+01, 0.1479589E+01, 0.0000000E+00 /
      data (oab23a(219,oi),oi=1,5) / 0.1092500E+02, 0.3285802E+02, 0.1013292E+02, 0.5654345E+00, 0.0000000E+00 /
      data (oab23a(220,oi),oi=1,5) / 0.1097500E+02, 0.2287233E+02, 0.5118617E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(221,oi),oi=1,5) / 0.1102500E+02, 0.1797073E+02, 0.1284934E+02, 0.2361906E+01, 0.2361906E+01 /
      data (oab23a(222,oi),oi=1,5) / 0.1107500E+02, 0.3107554E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(223,oi),oi=1,5) / 0.1112500E+02, 0.4183852E+02, 0.1307341E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(224,oi),oi=1,5) / 0.1117500E+02, 0.5292485E+01, 0.0000000E+00, 0.6267873E+00, 0.0000000E+00 /
      data (oab23a(225,oi),oi=1,5) / 0.1122500E+02, 0.2137240E+02, 0.1336493E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(226,oi),oi=1,5) / 0.1127500E+02, 0.1882129E+02, 0.5382865E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(227,oi),oi=1,5) / 0.1132500E+02, 0.2981202E+02, 0.8137540E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(228,oi),oi=1,5) / 0.1137500E+02, 0.3287078E+02, 0.5473425E+01, 0.5220169E+00, 0.0000000E+00 /
      data (oab23a(229,oi),oi=1,5) / 0.1142500E+02, 0.4963699E+02, 0.1104171E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(230,oi),oi=1,5) / 0.1147500E+02, 0.5013131E+02, 0.1389206E+02, 0.1043507E+01, 0.0000000E+00 /
      data (oab23a(231,oi),oi=1,5) / 0.1152500E+02, 0.1965292E+02, 0.1405826E+02, 0.8120261E+00, 0.6378054E+00 /
      data (oab23a(232,oi),oi=1,5) / 0.1157500E+02, 0.2271291E+02, 0.1134701E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(233,oi),oi=1,5) / 0.1162500E+02, 0.2000329E+02, 0.5703881E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(234,oi),oi=1,5) / 0.1167500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(235,oi),oi=1,5) / 0.1172500E+02, 0.1163401E+02, 0.5835070E+01, 0.9102843E+00, 0.0000000E+00 /
      data (oab23a(236,oi),oi=1,5) / 0.1177500E+02, 0.2929294E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(237,oi),oi=1,5) / 0.1182500E+02, 0.8902565E+01, 0.8908361E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(238,oi),oi=1,5) / 0.1187500E+02, 0.5959449E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(239,oi),oi=1,5) / 0.1192500E+02, 0.2101449E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(240,oi),oi=1,5) / 0.1197500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(241,oi),oi=1,5) / 0.1202500E+02, 0.3975686E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(242,oi),oi=1,5) / 0.1207500E+02, 0.1849574E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(243,oi),oi=1,5) / 0.1212500E+02, 0.1553133E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(244,oi),oi=1,5) / 0.1217500E+02, 0.3759887E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(245,oi),oi=1,5) / 0.1222500E+02, 0.2534212E+02, 0.9452159E+01, 0.9794071E+00, 0.0000000E+00 /
      data (oab23a(246,oi),oi=1,5) / 0.1227500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(247,oi),oi=1,5) / 0.1232500E+02, 0.1934293E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(248,oi),oi=1,5) / 0.1237500E+02, 0.2910363E+02, 0.0000000E+00, 0.6443614E+00, 0.0000000E+00 /
      data (oab23a(249,oi),oi=1,5) / 0.1242500E+02, 0.2616217E+02, 0.6555685E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(250,oi),oi=1,5) / 0.1247500E+02, 0.9901900E+01, 0.0000000E+00, 0.1053069E+01, 0.0000000E+00 /
      data (oab23a(251,oi),oi=1,5) / 0.1252500E+02, 0.1990929E+02, 0.1326029E+02, 0.6672470E+00, 0.0000000E+00 /
      data (oab23a(252,oi),oi=1,5) / 0.1257500E+02, 0.3016586E+02, 0.0000000E+00, 0.8541138E+00, 0.0000000E+00 /
      data (oab23a(253,oi),oi=1,5) / 0.1262500E+02, 0.1351456E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(254,oi),oi=1,5) / 0.1267500E+02, 0.1020958E+02, 0.6777994E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(255,oi),oi=1,5) / 0.1272500E+02, 0.1367867E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(256,oi),oi=1,5) / 0.1277500E+02, 0.3104764E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(257,oi),oi=1,5) / 0.1282500E+02, 0.1738261E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(258,oi),oi=1,5) / 0.1287500E+02, 0.4902654E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(259,oi),oi=1,5) / 0.1292500E+02, 0.7062350E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(260,oi),oi=1,5) / 0.1297500E+02, 0.5349853E+02, 0.7100194E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(261,oi),oi=1,5) / 0.1302500E+02, 0.2514126E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(262,oi),oi=1,5) / 0.1307500E+02, 0.7224550E+01, 0.1087013E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(263,oi),oi=1,5) / 0.1312500E+02, 0.3276988E+02, 0.0000000E+00, 0.9525578E+00, 0.0000000E+00 /
      data (oab23a(264,oi),oi=1,5) / 0.1317500E+02, 0.0000000E+00, 0.7337834E+01, 0.7400570E+00, 0.0000000E+00 /
      data (oab23a(265,oi),oi=1,5) / 0.1322500E+02, 0.1851223E+02, 0.7406698E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(266,oi),oi=1,5) / 0.1327500E+02, 0.2613479E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(267,oi),oi=1,5) / 0.1332500E+02, 0.4884746E+02, 0.0000000E+00, 0.0000000E+00, 0.1177906E+01 /
      data (oab23a(268,oi),oi=1,5) / 0.1337500E+02, 0.3029012E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(269,oi),oi=1,5) / 0.1342500E+02, 0.7620738E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(270,oi),oi=1,5) / 0.1347500E+02, 0.2306891E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(271,oi),oi=1,5) / 0.1352500E+02, 0.1936221E+02, 0.7724310E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(272,oi),oi=1,5) / 0.1357500E+02, 0.7770818E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(273,oi),oi=1,5) / 0.1362500E+02, 0.1570068E+02, 0.7869074E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(274,oi),oi=1,5) / 0.1367500E+02, 0.0000000E+00, 0.7937074E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(275,oi),oi=1,5) / 0.1372500E+02, 0.7948822E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(276,oi),oi=1,5) / 0.1377500E+02, 0.2006997E+02, 0.0000000E+00, 0.1233475E+01, 0.0000000E+00 /
      data (oab23a(277,oi),oi=1,5) / 0.1382500E+02, 0.8069506E+01, 0.1209944E+02, 0.1000647E+01, 0.0000000E+00 /
      data (oab23a(278,oi),oi=1,5) / 0.1387500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(279,oi),oi=1,5) / 0.1392500E+02, 0.0000000E+00, 0.0000000E+00, 0.1488552E+01, 0.0000000E+00 /
      data (oab23a(280,oi),oi=1,5) / 0.1397500E+02, 0.5357038E+02, 0.8258934E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(281,oi),oi=1,5) / 0.1402500E+02, 0.1244651E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(282,oi),oi=1,5) / 0.1407500E+02, 0.0000000E+00, 0.8403946E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(283,oi),oi=1,5) / 0.1412500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(284,oi),oi=1,5) / 0.1417500E+02, 0.8490330E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(285,oi),oi=1,5) / 0.1422500E+02, 0.2144874E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(286,oi),oi=1,5) / 0.1427500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(287,oi),oi=1,5) / 0.1432500E+02, 0.8652650E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(288,oi),oi=1,5) / 0.1437500E+02, 0.8722642E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(289,oi),oi=1,5) / 0.1442500E+02, 0.3518510E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(290,oi),oi=1,5) / 0.1447500E+02, 0.8880150E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(291,oi),oi=1,5) / 0.1452500E+02, 0.2235653E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(292,oi),oi=1,5) / 0.1457500E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(293,oi),oi=1,5) / 0.1462500E+02, 0.9067535E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(294,oi),oi=1,5) / 0.1467500E+02, 0.9100482E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(295,oi),oi=1,5) / 0.1472500E+02, 0.1374961E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(296,oi),oi=1,5) / 0.1477500E+02, 0.2310459E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(297,oi),oi=1,5) / 0.1482500E+02, 0.4656715E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(298,oi),oi=1,5) / 0.1487500E+02, 0.0000000E+00, 0.0000000E+00, 0.1345789E+01, 0.0000000E+00 /
      data (oab23a(299,oi),oi=1,5) / 0.1492500E+02, 0.1417536E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab23a(300,oi),oi=1,5) / 0.1497500E+02, 0.9517726E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /