C/*-----------------------------------------------------
C// Neutrino flux table at 1000 km
C// Off-axis angle is 1.3 deg. for anti-nu_mu focusing beam
C// Enu(GeV), nu_mu, bar{nu_mu}, nu_e, bar{nu_e}
C-------------------------------------------------------*
      real*8 oab13a(300,5)
      data (oab13a(1,oi),oi=1,5) / 0.2500000E-01, 0.2242214E+04, 0.1425167E+04, 0.4103605E+01, 0.1779729E+02 /
      data (oab13a(2,oi),oi=1,5) / 0.7500000E-01, 0.6973362E+04, 0.7324554E+04, 0.1717964E+02, 0.9491532E+02 /
      data (oab13a(3,oi),oi=1,5) / 0.1250000E+00, 0.4058831E+04, 0.1518469E+05, 0.2817810E+02, 0.2269295E+03 /
      data (oab13a(4,oi),oi=1,5) / 0.1750000E+00, 0.3304797E+04, 0.2207593E+05, 0.4035328E+02, 0.3184018E+03 /
      data (oab13a(5,oi),oi=1,5) / 0.2250000E+00, 0.3221549E+04, 0.3055133E+05, 0.5074343E+02, 0.4092247E+03 /
      data (oab13a(6,oi),oi=1,5) / 0.2750000E+00, 0.3342443E+04, 0.4047488E+05, 0.5777635E+02, 0.4788209E+03 /
      data (oab13a(7,oi),oi=1,5) / 0.3250000E+00, 0.3383431E+04, 0.5045359E+05, 0.6389758E+02, 0.4928611E+03 /
      data (oab13a(8,oi),oi=1,5) / 0.3750000E+00, 0.3622759E+04, 0.5957358E+05, 0.6565826E+02, 0.5518737E+03 /
      data (oab13a(9,oi),oi=1,5) / 0.4250000E+00, 0.3891236E+04, 0.6781047E+05, 0.6729898E+02, 0.5780750E+03 /
      data (oab13a(10,oi),oi=1,5) / 0.4750000E+00, 0.4075671E+04, 0.7503148E+05, 0.7271097E+02, 0.5319582E+03 /
      data (oab13a(11,oi),oi=1,5) / 0.5250000E+00, 0.4385902E+04, 0.8170005E+05, 0.6762715E+02, 0.5165025E+03 /
      data (oab13a(12,oi),oi=1,5) / 0.5750000E+00, 0.4430155E+04, 0.8800203E+05, 0.7264653E+02, 0.5298910E+03 /
      data (oab13a(13,oi),oi=1,5) / 0.6250000E+00, 0.4557453E+04, 0.9330663E+05, 0.6752398E+02, 0.5581161E+03 /
      data (oab13a(14,oi),oi=1,5) / 0.6750000E+00, 0.4865802E+04, 0.1003835E+06, 0.6578859E+02, 0.5625415E+03 /
      data (oab13a(15,oi),oi=1,5) / 0.7250000E+00, 0.4916653E+04, 0.1080668E+06, 0.6889132E+02, 0.5337137E+03 /
      data (oab13a(16,oi),oi=1,5) / 0.7750000E+00, 0.4819363E+04, 0.1154766E+06, 0.6499781E+02, 0.5194378E+03 /
      data (oab13a(17,oi),oi=1,5) / 0.8250000E+00, 0.4692517E+04, 0.1168653E+06, 0.7933822E+02, 0.5062839E+03 /
      data (oab13a(18,oi),oi=1,5) / 0.8750000E+00, 0.4959300E+04, 0.1148828E+06, 0.6571728E+02, 0.4615289E+03 /
      data (oab13a(19,oi),oi=1,5) / 0.9250000E+00, 0.5052053E+04, 0.1127864E+06, 0.6538092E+02, 0.4641575E+03 /
      data (oab13a(20,oi),oi=1,5) / 0.9750000E+00, 0.5134177E+04, 0.1107276E+06, 0.7698386E+02, 0.4782470E+03 /
      data (oab13a(21,oi),oi=1,5) / 0.1025000E+01, 0.5034690E+04, 0.1071830E+06, 0.6750224E+02, 0.4071020E+03 /
      data (oab13a(22,oi),oi=1,5) / 0.1075000E+01, 0.5127857E+04, 0.1018629E+06, 0.5776080E+02, 0.3989594E+03 /
      data (oab13a(23,oi),oi=1,5) / 0.1125000E+01, 0.4953127E+04, 0.9393115E+05, 0.6312409E+02, 0.3235289E+03 /
      data (oab13a(24,oi),oi=1,5) / 0.1175000E+01, 0.5111818E+04, 0.8625695E+05, 0.6142782E+02, 0.3538949E+03 /
      data (oab13a(25,oi),oi=1,5) / 0.1225000E+01, 0.4779138E+04, 0.7770749E+05, 0.6242895E+02, 0.3295844E+03 /
      data (oab13a(26,oi),oi=1,5) / 0.1275000E+01, 0.4500130E+04, 0.6965713E+05, 0.6070906E+02, 0.2755357E+03 /
      data (oab13a(27,oi),oi=1,5) / 0.1325000E+01, 0.4666040E+04, 0.6126198E+05, 0.6080239E+02, 0.3304609E+03 /
      data (oab13a(28,oi),oi=1,5) / 0.1375000E+01, 0.4245363E+04, 0.5477893E+05, 0.5034427E+02, 0.3070416E+03 /
      data (oab13a(29,oi),oi=1,5) / 0.1425000E+01, 0.4362183E+04, 0.4904873E+05, 0.5963334E+02, 0.2670541E+03 /
      data (oab13a(30,oi),oi=1,5) / 0.1475000E+01, 0.4110504E+04, 0.4375441E+05, 0.5244037E+02, 0.2788709E+03 /
      data (oab13a(31,oi),oi=1,5) / 0.1525000E+01, 0.4313080E+04, 0.3978424E+05, 0.6069466E+02, 0.2063112E+03 /
      data (oab13a(32,oi),oi=1,5) / 0.1575000E+01, 0.4063857E+04, 0.3474154E+05, 0.5531745E+02, 0.2697520E+03 /
      data (oab13a(33,oi),oi=1,5) / 0.1625000E+01, 0.4196835E+04, 0.3099942E+05, 0.4669062E+02, 0.2379954E+03 /
      data (oab13a(34,oi),oi=1,5) / 0.1675000E+01, 0.3488142E+04, 0.2820306E+05, 0.5648587E+02, 0.1925800E+03 /
      data (oab13a(35,oi),oi=1,5) / 0.1725000E+01, 0.3504354E+04, 0.2429009E+05, 0.5811551E+02, 0.2476507E+03 /
      data (oab13a(36,oi),oi=1,5) / 0.1775000E+01, 0.3277553E+04, 0.2251241E+05, 0.5705496E+02, 0.1470485E+03 /
      data (oab13a(37,oi),oi=1,5) / 0.1825000E+01, 0.3537216E+04, 0.1993689E+05, 0.4485220E+02, 0.1499981E+03 /
      data (oab13a(38,oi),oi=1,5) / 0.1875000E+01, 0.3746147E+04, 0.1764145E+05, 0.5079119E+02, 0.1329572E+03 /
      data (oab13a(39,oi),oi=1,5) / 0.1925000E+01, 0.3107423E+04, 0.1591643E+05, 0.4580638E+02, 0.1281759E+03 /
      data (oab13a(40,oi),oi=1,5) / 0.1975000E+01, 0.3272701E+04, 0.1486877E+05, 0.6441414E+02, 0.1665453E+03 /
      data (oab13a(41,oi),oi=1,5) / 0.2025000E+01, 0.3038562E+04, 0.1310822E+05, 0.4588360E+02, 0.1415672E+03 /
      data (oab13a(42,oi),oi=1,5) / 0.2075000E+01, 0.3092718E+04, 0.1202985E+05, 0.5552280E+02, 0.1272514E+03 /
      data (oab13a(43,oi),oi=1,5) / 0.2125000E+01, 0.3126578E+04, 0.1162211E+05, 0.7327701E+02, 0.1377720E+03 /
      data (oab13a(44,oi),oi=1,5) / 0.2175000E+01, 0.3054503E+04, 0.9698718E+04, 0.4458267E+02, 0.1118876E+03 /
      data (oab13a(45,oi),oi=1,5) / 0.2225000E+01, 0.2657683E+04, 0.9498566E+04, 0.3909021E+02, 0.9655803E+02 /
      data (oab13a(46,oi),oi=1,5) / 0.2275000E+01, 0.3095627E+04, 0.9141891E+04, 0.3312098E+02, 0.8015245E+02 /
      data (oab13a(47,oi),oi=1,5) / 0.2325000E+01, 0.2692323E+04, 0.8091939E+04, 0.3658580E+02, 0.8125948E+02 /
      data (oab13a(48,oi),oi=1,5) / 0.2375000E+01, 0.2693390E+04, 0.7347801E+04, 0.4570151E+02, 0.1117056E+03 /
      data (oab13a(49,oi),oi=1,5) / 0.2425000E+01, 0.2381089E+04, 0.8047973E+04, 0.4440054E+02, 0.7194501E+02 /
      data (oab13a(50,oi),oi=1,5) / 0.2475000E+01, 0.2800117E+04, 0.7206872E+04, 0.4632458E+02, 0.8909934E+02 /
      data (oab13a(51,oi),oi=1,5) / 0.2525000E+01, 0.2488947E+04, 0.6628337E+04, 0.4743925E+02, 0.1173325E+03 /
      data (oab13a(52,oi),oi=1,5) / 0.2575000E+01, 0.2363640E+04, 0.5403245E+04, 0.4162835E+02, 0.7186691E+02 /
      data (oab13a(53,oi),oi=1,5) / 0.2625000E+01, 0.2305754E+04, 0.5113275E+04, 0.3594361E+02, 0.5619931E+02 /
      data (oab13a(54,oi),oi=1,5) / 0.2675000E+01, 0.2080292E+04, 0.4574983E+04, 0.4346027E+02, 0.9191129E+02 /
      data (oab13a(55,oi),oi=1,5) / 0.2725000E+01, 0.2208178E+04, 0.5187050E+04, 0.3500060E+02, 0.8929556E+02 /
      data (oab13a(56,oi),oi=1,5) / 0.2775000E+01, 0.2237096E+04, 0.4686392E+04, 0.4683915E+02, 0.7211985E+02 /
      data (oab13a(57,oi),oi=1,5) / 0.2825000E+01, 0.1884586E+04, 0.4598677E+04, 0.3764564E+02, 0.5095240E+02 /
      data (oab13a(58,oi),oi=1,5) / 0.2875000E+01, 0.2246693E+04, 0.4155818E+04, 0.6005852E+02, 0.6338731E+02 /
      data (oab13a(59,oi),oi=1,5) / 0.2925000E+01, 0.2340911E+04, 0.3706368E+04, 0.3221250E+02, 0.5988557E+02 /
      data (oab13a(60,oi),oi=1,5) / 0.2975000E+01, 0.1902673E+04, 0.3307813E+04, 0.3711173E+02, 0.9996787E+02 /
      data (oab13a(61,oi),oi=1,5) / 0.3025000E+01, 0.1660162E+04, 0.3810251E+04, 0.4470062E+02, 0.5793211E+02 /
      data (oab13a(62,oi),oi=1,5) / 0.3075000E+01, 0.1251799E+04, 0.2947537E+04, 0.3686790E+02, 0.4630442E+02 /
      data (oab13a(63,oi),oi=1,5) / 0.3125000E+01, 0.2198709E+04, 0.2997991E+04, 0.5033475E+02, 0.6546242E+02 /
      data (oab13a(64,oi),oi=1,5) / 0.3175000E+01, 0.1633713E+04, 0.3076358E+04, 0.3709670E+02, 0.6639474E+02 /
      data (oab13a(65,oi),oi=1,5) / 0.3225000E+01, 0.1598314E+04, 0.2604531E+04, 0.3157730E+02, 0.1178945E+03 /
      data (oab13a(66,oi),oi=1,5) / 0.3275000E+01, 0.1412868E+04, 0.2835124E+04, 0.3729031E+02, 0.5841880E+02 /
      data (oab13a(67,oi),oi=1,5) / 0.3325000E+01, 0.1421621E+04, 0.2482741E+04, 0.3766940E+02, 0.5188911E+02 /
      data (oab13a(68,oi),oi=1,5) / 0.3375000E+01, 0.1866286E+04, 0.2909746E+04, 0.3811247E+02, 0.4879020E+02 /
      data (oab13a(69,oi),oi=1,5) / 0.3425000E+01, 0.1519132E+04, 0.2603158E+04, 0.3542053E+02, 0.4396412E+02 /
      data (oab13a(70,oi),oi=1,5) / 0.3475000E+01, 0.1038919E+04, 0.2426221E+04, 0.5085137E+02, 0.5169521E+02 /
      data (oab13a(71,oi),oi=1,5) / 0.3525000E+01, 0.1643925E+04, 0.2396029E+04, 0.3403336E+02, 0.4420342E+02 /
      data (oab13a(72,oi),oi=1,5) / 0.3575000E+01, 0.1214526E+04, 0.2417464E+04, 0.3769621E+02, 0.4542923E+02 /
      data (oab13a(73,oi),oi=1,5) / 0.3625000E+01, 0.1203645E+04, 0.2484040E+04, 0.3396639E+02, 0.4455866E+02 /
      data (oab13a(74,oi),oi=1,5) / 0.3675000E+01, 0.1219876E+04, 0.2379402E+04, 0.3705147E+02, 0.4690197E+02 /
      data (oab13a(75,oi),oi=1,5) / 0.3725000E+01, 0.1433289E+04, 0.1907931E+04, 0.3038201E+02, 0.3823063E+02 /
      data (oab13a(76,oi),oi=1,5) / 0.3775000E+01, 0.1608399E+04, 0.1988407E+04, 0.2879718E+02, 0.4841548E+02 /
      data (oab13a(77,oi),oi=1,5) / 0.3825000E+01, 0.1356161E+04, 0.2165051E+04, 0.2781622E+02, 0.3249537E+02 /
      data (oab13a(78,oi),oi=1,5) / 0.3875000E+01, 0.1150237E+04, 0.2157653E+04, 0.3255430E+02, 0.5576867E+02 /
      data (oab13a(79,oi),oi=1,5) / 0.3925000E+01, 0.1252688E+04, 0.2404197E+04, 0.2074528E+02, 0.3487094E+02 /
      data (oab13a(80,oi),oi=1,5) / 0.3975000E+01, 0.1606403E+04, 0.1952376E+04, 0.2898476E+02, 0.3496257E+02 /
      data (oab13a(81,oi),oi=1,5) / 0.4025000E+01, 0.1220156E+04, 0.2082819E+04, 0.2888390E+02, 0.3517829E+02 /
      data (oab13a(82,oi),oi=1,5) / 0.4075000E+01, 0.8053239E+03, 0.1488115E+04, 0.2874502E+02, 0.3568352E+02 /
      data (oab13a(83,oi),oi=1,5) / 0.4125000E+01, 0.1118268E+04, 0.1511071E+04, 0.3135945E+02, 0.3701502E+02 /
      data (oab13a(84,oi),oi=1,5) / 0.4175000E+01, 0.1132216E+04, 0.1506540E+04, 0.2475976E+02, 0.3660266E+02 /
      data (oab13a(85,oi),oi=1,5) / 0.4225000E+01, 0.1108126E+04, 0.1188623E+04, 0.2421417E+02, 0.2717459E+02 /
      data (oab13a(86,oi),oi=1,5) / 0.4275000E+01, 0.1293401E+04, 0.1605260E+04, 0.1977409E+02, 0.2432965E+02 /
      data (oab13a(87,oi),oi=1,5) / 0.4325000E+01, 0.1202605E+04, 0.1494149E+04, 0.1789905E+02, 0.2492859E+02 /
      data (oab13a(88,oi),oi=1,5) / 0.4375000E+01, 0.9299556E+03, 0.1253048E+04, 0.1987109E+02, 0.2724561E+02 /
      data (oab13a(89,oi),oi=1,5) / 0.4425000E+01, 0.7836597E+03, 0.1128292E+04, 0.1516899E+02, 0.2116659E+02 /
      data (oab13a(90,oi),oi=1,5) / 0.4475000E+01, 0.1148546E+04, 0.9753897E+03, 0.1830257E+02, 0.2099709E+02 /
      data (oab13a(91,oi),oi=1,5) / 0.4525000E+01, 0.1238302E+04, 0.1349206E+04, 0.3439491E+02, 0.3494181E+02 /
      data (oab13a(92,oi),oi=1,5) / 0.4575000E+01, 0.8271250E+03, 0.8132479E+03, 0.1657377E+02, 0.2130043E+02 /
      data (oab13a(93,oi),oi=1,5) / 0.4625000E+01, 0.1154061E+04, 0.1157174E+04, 0.3277519E+02, 0.2813333E+02 /
      data (oab13a(94,oi),oi=1,5) / 0.4675000E+01, 0.5272261E+03, 0.1279395E+04, 0.2415414E+02, 0.2490869E+02 /
      data (oab13a(95,oi),oi=1,5) / 0.4725000E+01, 0.9727507E+03, 0.1317589E+04, 0.2727894E+02, 0.2702287E+02 /
      data (oab13a(96,oi),oi=1,5) / 0.4775000E+01, 0.1412676E+04, 0.1029448E+04, 0.1956240E+02, 0.2548131E+02 /
      data (oab13a(97,oi),oi=1,5) / 0.4825000E+01, 0.8301450E+03, 0.7264343E+03, 0.2851571E+02, 0.2424641E+02 /
      data (oab13a(98,oi),oi=1,5) / 0.4875000E+01, 0.1117807E+04, 0.9055133E+03, 0.1811626E+02, 0.1764235E+02 /
      data (oab13a(99,oi),oi=1,5) / 0.4925000E+01, 0.8764760E+03, 0.7830425E+03, 0.2658304E+02, 0.2448531E+02 /
      data (oab13a(100,oi),oi=1,5) / 0.4975000E+01, 0.1071074E+04, 0.1103314E+04, 0.2273397E+02, 0.2282015E+02 /
      data (oab13a(101,oi),oi=1,5) / 0.5025000E+01, 0.6304527E+03, 0.7793923E+03, 0.1971307E+02, 0.2015302E+02 /
      data (oab13a(102,oi),oi=1,5) / 0.5075000E+01, 0.8739373E+03, 0.1078597E+04, 0.1804782E+02, 0.1964376E+02 /
      data (oab13a(103,oi),oi=1,5) / 0.5125000E+01, 0.7010463E+03, 0.7093548E+03, 0.1193495E+02, 0.1269129E+02 /
      data (oab13a(104,oi),oi=1,5) / 0.5175000E+01, 0.1061551E+04, 0.6858936E+03, 0.2861998E+02, 0.2344800E+02 /
      data (oab13a(105,oi),oi=1,5) / 0.5225000E+01, 0.8391272E+03, 0.5217153E+03, 0.2077274E+02, 0.1749122E+02 /
      data (oab13a(106,oi),oi=1,5) / 0.5275000E+01, 0.6497331E+03, 0.9562585E+03, 0.2294206E+02, 0.6057210E+02 /
      data (oab13a(107,oi),oi=1,5) / 0.5325000E+01, 0.6028345E+03, 0.9606273E+03, 0.1872959E+02, 0.1527030E+02 /
      data (oab13a(108,oi),oi=1,5) / 0.5375000E+01, 0.8123988E+03, 0.8118585E+03, 0.1952602E+02, 0.1897865E+02 /
      data (oab13a(109,oi),oi=1,5) / 0.5425000E+01, 0.5867706E+03, 0.9667067E+03, 0.2328529E+02, 0.2091301E+02 /
      data (oab13a(110,oi),oi=1,5) / 0.5475000E+01, 0.3656605E+03, 0.4726156E+03, 0.1566257E+02, 0.1658309E+02 /
      data (oab13a(111,oi),oi=1,5) / 0.5525000E+01, 0.5590761E+03, 0.8349283E+03, 0.2000694E+02, 0.1655919E+02 /
      data (oab13a(112,oi),oi=1,5) / 0.5575000E+01, 0.4734026E+03, 0.2834768E+03, 0.1782422E+02, 0.1690953E+02 /
      data (oab13a(113,oi),oi=1,5) / 0.5625000E+01, 0.1186466E+04, 0.6971380E+03, 0.1494825E+02, 0.1243729E+02 /
      data (oab13a(114,oi),oi=1,5) / 0.5675000E+01, 0.5251662E+03, 0.4700270E+03, 0.3115843E+02, 0.2908438E+02 /
      data (oab13a(115,oi),oi=1,5) / 0.5725000E+01, 0.7409878E+03, 0.8063817E+03, 0.2016863E+02, 0.1458298E+02 /
      data (oab13a(116,oi),oi=1,5) / 0.5775000E+01, 0.2148016E+03, 0.1044140E+04, 0.1329299E+02, 0.9975746E+01 /
      data (oab13a(117,oi),oi=1,5) / 0.5825000E+01, 0.2925370E+03, 0.5587238E+03, 0.1243389E+02, 0.9366533E+01 /
      data (oab13a(118,oi),oi=1,5) / 0.5875000E+01, 0.3303644E+03, 0.5564350E+03, 0.1369616E+02, 0.1152523E+02 /
      data (oab13a(119,oi),oi=1,5) / 0.5925000E+01, 0.4753815E+03, 0.8061350E+03, 0.2153099E+02, 0.1789831E+02 /
      data (oab13a(120,oi),oi=1,5) / 0.5975000E+01, 0.1310327E+04, 0.3153117E+03, 0.1216664E+02, 0.1088693E+02 /
      data (oab13a(121,oi),oi=1,5) / 0.6025000E+01, 0.4117168E+03, 0.6533462E+03, 0.1740741E+02, 0.1208461E+02 /
      data (oab13a(122,oi),oi=1,5) / 0.6075000E+01, 0.2521710E+03, 0.3360802E+03, 0.1264593E+02, 0.8564751E+01 /
      data (oab13a(123,oi),oi=1,5) / 0.6125000E+01, 0.4337896E+03, 0.8300868E+03, 0.1251568E+02, 0.9657792E+01 /
      data (oab13a(124,oi),oi=1,5) / 0.6175000E+01, 0.3987112E+03, 0.4181184E+03, 0.1561263E+02, 0.1292867E+02 /
      data (oab13a(125,oi),oi=1,5) / 0.6225000E+01, 0.6409974E+03, 0.3669729E+03, 0.1572302E+02, 0.1419975E+02 /
      data (oab13a(126,oi),oi=1,5) / 0.6275000E+01, 0.3712352E+03, 0.7312496E+03, 0.1917262E+02, 0.1274757E+02 /
      data (oab13a(127,oi),oi=1,5) / 0.6325000E+01, 0.6914606E+03, 0.7967653E+03, 0.1165680E+02, 0.1093306E+02 /
      data (oab13a(128,oi),oi=1,5) / 0.6375000E+01, 0.6035706E+03, 0.4688978E+03, 0.2015304E+02, 0.1684405E+02 /
      data (oab13a(129,oi),oi=1,5) / 0.6425000E+01, 0.3355263E+03, 0.1973751E+03, 0.1153725E+02, 0.7910561E+01 /
      data (oab13a(130,oi),oi=1,5) / 0.6475000E+01, 0.2771866E+03, 0.5714342E+03, 0.9331743E+01, 0.9270815E+01 /
      data (oab13a(131,oi),oi=1,5) / 0.6525000E+01, 0.3778221E+03, 0.3479647E+03, 0.1335978E+02, 0.1144568E+02 /
      data (oab13a(132,oi),oi=1,5) / 0.6575000E+01, 0.4958033E+03, 0.3575840E+03, 0.1283133E+02, 0.9451298E+01 /
      data (oab13a(133,oi),oi=1,5) / 0.6625000E+01, 0.4809019E+03, 0.2827010E+03, 0.1158040E+02, 0.5925372E+01 /
      data (oab13a(134,oi),oi=1,5) / 0.6675000E+01, 0.4290791E+03, 0.4348261E+03, 0.1239938E+02, 0.9545209E+01 /
      data (oab13a(135,oi),oi=1,5) / 0.6725000E+01, 0.3882048E+03, 0.3750219E+03, 0.1299380E+02, 0.1031285E+02 /
      data (oab13a(136,oi),oi=1,5) / 0.6775000E+01, 0.3851914E+03, 0.2214499E+03, 0.1158724E+02, 0.4683333E+01 /
      data (oab13a(137,oi),oi=1,5) / 0.6825000E+01, 0.1181858E+03, 0.3836406E+03, 0.9777253E+01, 0.9910574E+01 /
      data (oab13a(138,oi),oi=1,5) / 0.6875000E+01, 0.4519331E+03, 0.4698085E+03, 0.1556514E+02, 0.1143249E+02 /
      data (oab13a(139,oi),oi=1,5) / 0.6925000E+01, 0.2508742E+03, 0.3072509E+03, 0.9969158E+01, 0.6617005E+01 /
      data (oab13a(140,oi),oi=1,5) / 0.6975000E+01, 0.6227973E+03, 0.3913519E+03, 0.9212806E+01, 0.8635554E+01 /
      data (oab13a(141,oi),oi=1,5) / 0.7025000E+01, 0.1032836E+03, 0.4208042E+03, 0.1901974E+02, 0.1737284E+02 /
      data (oab13a(142,oi),oi=1,5) / 0.7075000E+01, 0.4077283E+03, 0.3015941E+03, 0.7202165E+01, 0.2331826E+01 /
      data (oab13a(143,oi),oi=1,5) / 0.7125000E+01, 0.2137271E+03, 0.1298355E+03, 0.1099685E+02, 0.9138451E+01 /
      data (oab13a(144,oi),oi=1,5) / 0.7175000E+01, 0.2757730E+03, 0.3109384E+03, 0.9075270E+01, 0.5679317E+01 /
      data (oab13a(145,oi),oi=1,5) / 0.7225000E+01, 0.4320600E+03, 0.1524620E+03, 0.1328103E+02, 0.1038978E+02 /
      data (oab13a(146,oi),oi=1,5) / 0.7275000E+01, 0.9743583E+02, 0.1712396E+03, 0.1655274E+02, 0.1331155E+02 /
      data (oab13a(147,oi),oi=1,5) / 0.7325000E+01, 0.2762388E+03, 0.1633365E+03, 0.6618530E+01, 0.7221057E+01 /
      data (oab13a(148,oi),oi=1,5) / 0.7375000E+01, 0.9097085E+02, 0.1540473E+03, 0.6738301E+01, 0.3776255E+01 /
      data (oab13a(149,oi),oi=1,5) / 0.7425000E+01, 0.2603537E+03, 0.1131612E+03, 0.1352847E+02, 0.6535026E+01 /
      data (oab13a(150,oi),oi=1,5) / 0.7475000E+01, 0.4262902E+03, 0.1453169E+03, 0.8101768E+01, 0.2508000E+01 /
      data (oab13a(151,oi),oi=1,5) / 0.7525000E+01, 0.8750464E+02, 0.1389289E+03, 0.1111929E+02, 0.8578288E+01 /
      data (oab13a(152,oi),oi=1,5) / 0.7575000E+01, 0.1093174E+03, 0.2827649E+03, 0.1001352E+02, 0.5839381E+01 /
      data (oab13a(153,oi),oi=1,5) / 0.7625000E+01, 0.2701566E+03, 0.9722003E+02, 0.1558677E+02, 0.1099844E+02 /
      data (oab13a(154,oi),oi=1,5) / 0.7675000E+01, 0.2618789E+03, 0.7724592E+02, 0.8238507E+01, 0.5272577E+01 /
      data (oab13a(155,oi),oi=1,5) / 0.7725000E+01, 0.1036623E+03, 0.2915451E+03, 0.1098462E+02, 0.4332010E+01 /
      data (oab13a(156,oi),oi=1,5) / 0.7775000E+01, 0.9466081E+02, 0.1098701E+03, 0.9768231E+01, 0.4695020E+01 /
      data (oab13a(157,oi),oi=1,5) / 0.7825000E+01, 0.8811759E+02, 0.1256079E+03, 0.1317320E+02, 0.8006061E+01 /
      data (oab13a(158,oi),oi=1,5) / 0.7875000E+01, 0.8278643E+02, 0.9314294E+02, 0.6798954E+01, 0.4260711E+01 /
      data (oab13a(159,oi),oi=1,5) / 0.7925000E+01, 0.7182678E+02, 0.3891091E+03, 0.1030276E+02, 0.4712577E+01 /
      data (oab13a(160,oi),oi=1,5) / 0.7975000E+01, 0.1143995E+03, 0.1077217E+03, 0.6831966E+01, 0.3114939E+01 /
      data (oab13a(161,oi),oi=1,5) / 0.8025000E+01, 0.1117618E+03, 0.1104474E+03, 0.1206087E+02, 0.8561650E+01 /
      data (oab13a(162,oi),oi=1,5) / 0.8075000E+01, 0.8694540E+02, 0.1213141E+03, 0.7167905E+01, 0.5667746E+01 /
      data (oab13a(163,oi),oi=1,5) / 0.8125000E+01, 0.2979024E+03, 0.5437191E+02, 0.1517362E+01, 0.1077318E+01 /
      data (oab13a(164,oi),oi=1,5) / 0.8175000E+01, 0.9895139E+02, 0.9887023E+02, 0.8357230E+01, 0.7763618E+01 /
      data (oab13a(165,oi),oi=1,5) / 0.8225000E+01, 0.3104568E+03, 0.9292537E+02, 0.5619139E+01, 0.2986156E+01 /
      data (oab13a(166,oi),oi=1,5) / 0.8275000E+01, 0.4888410E+03, 0.8108567E+02, 0.5188307E+01, 0.2741232E+01 /
      data (oab13a(167,oi),oi=1,5) / 0.8325000E+01, 0.2995183E+03, 0.7482749E+02, 0.4760756E+01, 0.3378124E+01 /
      data (oab13a(168,oi),oi=1,5) / 0.8375000E+01, 0.5275003E+03, 0.5048572E+02, 0.1242289E+02, 0.1194869E+02 /
      data (oab13a(169,oi),oi=1,5) / 0.8425000E+01, 0.1186614E+03, 0.5559123E+02, 0.5355233E+01, 0.3083461E+01 /
      data (oab13a(170,oi),oi=1,5) / 0.8475000E+01, 0.6380125E+02, 0.6678207E+02, 0.2746364E+01, 0.9174029E+00 /
      data (oab13a(171,oi),oi=1,5) / 0.8525000E+01, 0.6769640E+02, 0.8460519E+02, 0.7172127E+01, 0.3122361E+01 /
      data (oab13a(172,oi),oi=1,5) / 0.8575000E+01, 0.5443161E+02, 0.2813347E+03, 0.3645166E+01, 0.0000000E+00 /
      data (oab13a(173,oi),oi=1,5) / 0.8625000E+01, 0.8492422E+02, 0.7081939E+02, 0.6881779E+01, 0.4364363E+01 /
      data (oab13a(174,oi),oi=1,5) / 0.8675000E+01, 0.2954621E+03, 0.2799605E+03, 0.2905975E+01, 0.1418225E+01 /
      data (oab13a(175,oi),oi=1,5) / 0.8725000E+01, 0.9342713E+02, 0.5807039E+02, 0.3584899E+01, 0.1249960E+01 /
      data (oab13a(176,oi),oi=1,5) / 0.8775000E+01, 0.8313027E+02, 0.7808549E+02, 0.5640407E+01, 0.2825152E+01 /
      data (oab13a(177,oi),oi=1,5) / 0.8825000E+01, 0.1104156E+03, 0.6766986E+02, 0.1149266E+02, 0.1009539E+02 /
      data (oab13a(178,oi),oi=1,5) / 0.8875000E+01, 0.1000480E+03, 0.4336922E+02, 0.2747644E+01, 0.0000000E+00 /
      data (oab13a(179,oi),oi=1,5) / 0.8925000E+01, 0.9102753E+02, 0.4714184E+02, 0.5090677E+01, 0.2990907E+01 /
      data (oab13a(180,oi),oi=1,5) / 0.8975000E+01, 0.8858679E+02, 0.6304286E+02, 0.7597477E+01, 0.6048689E+01 /
      data (oab13a(181,oi),oi=1,5) / 0.9025000E+01, 0.3168945E+03, 0.6369623E+02, 0.1271350E+02, 0.8541464E+01 /
      data (oab13a(182,oi),oi=1,5) / 0.9075000E+01, 0.7664831E+02, 0.2991211E+03, 0.3365330E+01, 0.2244437E+01 /
      data (oab13a(183,oi),oi=1,5) / 0.9125000E+01, 0.8451616E+02, 0.2643030E+02, 0.1334894E+02, 0.1209937E+02 /
      data (oab13a(184,oi),oi=1,5) / 0.9175000E+01, 0.8720714E+02, 0.2920733E+03, 0.8669024E+01, 0.7321453E+01 /
      data (oab13a(185,oi),oi=1,5) / 0.9225000E+01, 0.7022388E+02, 0.3242382E+02, 0.2953010E+01, 0.1567170E+01 /
      data (oab13a(186,oi),oi=1,5) / 0.9275000E+01, 0.5826112E+02, 0.5271228E+02, 0.1001146E+02, 0.4085705E+01 /
      data (oab13a(187,oi),oi=1,5) / 0.9325000E+01, 0.9199079E+02, 0.5885469E+02, 0.9281814E+00, 0.2052352E+01 /
      data (oab13a(188,oi),oi=1,5) / 0.9375000E+01, 0.7627248E+02, 0.3530806E+02, 0.1322176E+02, 0.8046574E+01 /
      data (oab13a(189,oi),oi=1,5) / 0.9425000E+01, 0.9766804E+02, 0.4331885E+02, 0.4766645E+01, 0.2007607E+01 /
      data (oab13a(190,oi),oi=1,5) / 0.9475000E+01, 0.5135394E+02, 0.2904505E+03, 0.1664886E+01, 0.7942067E+00 /
      data (oab13a(191,oi),oi=1,5) / 0.9525000E+01, 0.1055945E+03, 0.4792269E+02, 0.1575685E+01, 0.3627578E+00 /
      data (oab13a(192,oi),oi=1,5) / 0.9575000E+01, 0.3333454E+03, 0.2994036E+03, 0.9096301E+01, 0.6253794E+01 /
      data (oab13a(193,oi),oi=1,5) / 0.9625000E+01, 0.3651104E+03, 0.3135681E+02, 0.2024385E+01, 0.0000000E+00 /
      data (oab13a(194,oi),oi=1,5) / 0.9675000E+01, 0.6138203E+02, 0.1981065E+02, 0.1206379E+01, 0.1052413E+01 /
      data (oab13a(195,oi),oi=1,5) / 0.9725000E+01, 0.6402494E+02, 0.4564510E+03, 0.7450344E+01, 0.4064093E+01 /
      data (oab13a(196,oi),oi=1,5) / 0.9775000E+01, 0.9098187E+02, 0.3191450E+03, 0.6819413E+01, 0.6819413E+01 /
      data (oab13a(197,oi),oi=1,5) / 0.9825000E+01, 0.3854941E+03, 0.3464432E+02, 0.2083802E+01, 0.9951539E+00 /
      data (oab13a(198,oi),oi=1,5) / 0.9875000E+01, 0.7419286E+02, 0.3718696E+02, 0.8719138E+01, 0.6422497E+01 /
      data (oab13a(199,oi),oi=1,5) / 0.9925000E+01, 0.3846823E+03, 0.7496444E+02, 0.3429014E+01, 0.3049645E+01 /
      data (oab13a(200,oi),oi=1,5) / 0.9975000E+01, 0.1010590E+03, 0.2528482E+02, 0.2757865E+01, 0.1494609E+01 /
      data (oab13a(201,oi),oi=1,5) / 0.1002500E+02, 0.7233872E+02, 0.4696601E+03, 0.1510086E+01, 0.6971945E+00 /
      data (oab13a(202,oi),oi=1,5) / 0.1007500E+02, 0.3225075E+02, 0.4518730E+02, 0.9361116E+01, 0.5597854E+01 /
      data (oab13a(203,oi),oi=1,5) / 0.1012500E+02, 0.4337982E+02, 0.3687742E+02, 0.5621790E+01, 0.2792565E+01 /
      data (oab13a(204,oi),oi=1,5) / 0.1017500E+02, 0.3719756E+03, 0.3501135E+02, 0.2041314E+01, 0.2041314E+01 /
      data (oab13a(205,oi),oi=1,5) / 0.1022500E+02, 0.9064625E+02, 0.2656980E+02, 0.2257873E+01, 0.4487493E+00 /
      data (oab13a(206,oi),oi=1,5) / 0.1027500E+02, 0.9150703E+02, 0.2005522E+02, 0.1003793E+01, 0.8286030E+00 /
      data (oab13a(207,oi),oi=1,5) / 0.1032500E+02, 0.6531615E+02, 0.2707983E+02, 0.2200438E+01, 0.8628407E+00 /
      data (oab13a(208,oi),oi=1,5) / 0.1037500E+02, 0.6829206E+02, 0.4784758E+02, 0.1741038E+01, 0.0000000E+00 /
      data (oab13a(209,oi),oi=1,5) / 0.1042500E+02, 0.5065985E+02, 0.1375175E+02, 0.6222422E+00, 0.7876589E+00 /
      data (oab13a(210,oi),oi=1,5) / 0.1047500E+02, 0.6963349E+02, 0.0000000E+00, 0.1779605E+01, 0.0000000E+00 /
      data (oab13a(211,oi),oi=1,5) / 0.1052500E+02, 0.4674481E+02, 0.3983896E+02, 0.3165642E+01, 0.0000000E+00 /
      data (oab13a(212,oi),oi=1,5) / 0.1057500E+02, 0.5919702E+02, 0.9475250E+01, 0.2081971E+01, 0.0000000E+00 /
      data (oab13a(213,oi),oi=1,5) / 0.1062500E+02, 0.5493182E+02, 0.3579039E+02, 0.3358143E+01, 0.0000000E+00 /
      data (oab13a(214,oi),oi=1,5) / 0.1067500E+02, 0.1156602E+03, 0.2650258E+02, 0.5273661E+00, 0.0000000E+00 /
      data (oab13a(215,oi),oi=1,5) / 0.1072500E+02, 0.6565518E+02, 0.3651222E+02, 0.6587098E+00, 0.1190464E+01 /
      data (oab13a(216,oi),oi=1,5) / 0.1077500E+02, 0.8601173E+02, 0.1962865E+02, 0.6215404E+01, 0.4326273E+01 /
      data (oab13a(217,oi),oi=1,5) / 0.1082500E+02, 0.5700601E+02, 0.3645669E+03, 0.3048406E+01, 0.2379679E+01 /
      data (oab13a(218,oi),oi=1,5) / 0.1087500E+02, 0.7760083E+02, 0.2750887E+02, 0.6291493E+00, 0.5409489E+00 /
      data (oab13a(219,oi),oi=1,5) / 0.1092500E+02, 0.7326226E+02, 0.1517695E+02, 0.6537645E+00, 0.0000000E+00 /
      data (oab13a(220,oi),oi=1,5) / 0.1097500E+02, 0.6120007E+02, 0.5086729E+01, 0.4766784E+01, 0.2315268E+01 /
      data (oab13a(221,oi),oi=1,5) / 0.1102500E+02, 0.5399522E+02, 0.4889794E+02, 0.3595325E+01, 0.1502168E+01 /
      data (oab13a(222,oi),oi=1,5) / 0.1107500E+02, 0.5443010E+02, 0.7766897E+01, 0.2564861E+01, 0.1337683E+01 /
      data (oab13a(223,oi),oi=1,5) / 0.1112500E+02, 0.3657354E+02, 0.5500793E+02, 0.0000000E+00, 0.5228801E+00 /
      data (oab13a(224,oi),oi=1,5) / 0.1117500E+02, 0.2381987E+02, 0.2115542E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(225,oi),oi=1,5) / 0.1122500E+02, 0.7199230E+02, 0.7962515E+01, 0.2851049E+01, 0.0000000E+00 /
      data (oab13a(226,oi),oi=1,5) / 0.1127500E+02, 0.6724822E+02, 0.2423921E+02, 0.6603398E+00, 0.0000000E+00 /
      data (oab13a(227,oi),oi=1,5) / 0.1132500E+02, 0.5418978E+02, 0.1627940E+02, 0.3253791E+01, 0.0000000E+00 /
      data (oab13a(228,oi),oi=1,5) / 0.1137500E+02, 0.4655435E+02, 0.1368961E+02, 0.7738722E+00, 0.0000000E+00 /
      data (oab13a(229,oi),oi=1,5) / 0.1142500E+02, 0.7738073E+02, 0.1657510E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(230,oi),oi=1,5) / 0.1147500E+02, 0.3060169E+02, 0.0000000E+00, 0.3498952E+01, 0.0000000E+00 /
      data (oab13a(231,oi),oi=1,5) / 0.1152500E+02, 0.4219669E+02, 0.1408184E+02, 0.2047774E+01, 0.0000000E+00 /
      data (oab13a(232,oi),oi=1,5) / 0.1157500E+02, 0.5384741E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(233,oi),oi=1,5) / 0.1162500E+02, 0.1998723E+02, 0.2289231E+02, 0.1421662E+01, 0.0000000E+00 /
      data (oab13a(234,oi),oi=1,5) / 0.1167500E+02, 0.5485058E+02, 0.0000000E+00, 0.6330927E+01, 0.6330927E+01 /
      data (oab13a(235,oi),oi=1,5) / 0.1172500E+02, 0.5518588E+02, 0.3492476E+02, 0.1530700E+01, 0.0000000E+00 /
      data (oab13a(236,oi),oi=1,5) / 0.1177500E+02, 0.6746321E+02, 0.1468879E+02, 0.2428080E+01, 0.0000000E+00 /
      data (oab13a(237,oi),oi=1,5) / 0.1182500E+02, 0.2074127E+02, 0.5912878E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(238,oi),oi=1,5) / 0.1187500E+02, 0.3874477E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(239,oi),oi=1,5) / 0.1192500E+02, 0.4507796E+02, 0.1808927E+02, 0.2252085E+01, 0.0000000E+00 /
      data (oab13a(240,oi),oi=1,5) / 0.1197500E+02, 0.7892348E+02, 0.2121272E+02, 0.6245739E+01, 0.5307585E+01 /
      data (oab13a(241,oi),oi=1,5) / 0.1202500E+02, 0.8263485E+02, 0.9186478E+01, 0.1373543E+01, 0.0000000E+00 /
      data (oab13a(242,oi),oi=1,5) / 0.1207500E+02, 0.4925459E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(243,oi),oi=1,5) / 0.1212500E+02, 0.4968662E+02, 0.1244059E+02, 0.1629587E+01, 0.0000000E+00 /
      data (oab13a(244,oi),oi=1,5) / 0.1217500E+02, 0.7211679E+02, 0.6267466E+01, 0.8708537E+00, 0.0000000E+00 /
      data (oab13a(245,oi),oi=1,5) / 0.1222500E+02, 0.6641749E+02, 0.9467908E+01, 0.6824621E+00, 0.0000000E+00 /
      data (oab13a(246,oi),oi=1,5) / 0.1227500E+02, 0.6376957E+02, 0.6369989E+01, 0.8421706E+01, 0.8421706E+01 /
      data (oab13a(247,oi),oi=1,5) / 0.1232500E+02, 0.3215984E+02, 0.6408334E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(248,oi),oi=1,5) / 0.1237500E+02, 0.7451665E+02, 0.6487257E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(249,oi),oi=1,5) / 0.1242500E+02, 0.4576664E+02, 0.2605321E+02, 0.6419561E+01, 0.5080829E+01 /
      data (oab13a(250,oi),oi=1,5) / 0.1247500E+02, 0.5926176E+02, 0.2960354E+02, 0.3515425E+01, 0.3515425E+01 /
      data (oab13a(251,oi),oi=1,5) / 0.1252500E+02, 0.7305977E+02, 0.2984024E+02, 0.6728522E+00, 0.0000000E+00 /
      data (oab13a(252,oi),oi=1,5) / 0.1257500E+02, 0.6697570E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(253,oi),oi=1,5) / 0.1262500E+02, 0.4384695E+02, 0.1012535E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(254,oi),oi=1,5) / 0.1267500E+02, 0.8164723E+02, 0.0000000E+00, 0.4406445E+01, 0.2768093E+01 /
      data (oab13a(255,oi),oi=1,5) / 0.1272500E+02, 0.6171548E+02, 0.2748135E+02, 0.9476798E+00, 0.0000000E+00 /
      data (oab13a(256,oi),oi=1,5) / 0.1277500E+02, 0.4484337E+02, 0.2415991E+02, 0.3133295E+01, 0.0000000E+00 /
      data (oab13a(257,oi),oi=1,5) / 0.1282500E+02, 0.4864475E+02, 0.1744950E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(258,oi),oi=1,5) / 0.1287500E+02, 0.5620807E+02, 0.1047868E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(259,oi),oi=1,5) / 0.1292500E+02, 0.5655268E+02, 0.0000000E+00, 0.1457759E+01, 0.0000000E+00 /
      data (oab13a(260,oi),oi=1,5) / 0.1297500E+02, 0.2492987E+02, 0.1424678E+02, 0.8918766E+00, 0.0000000E+00 /
      data (oab13a(261,oi),oi=1,5) / 0.1302500E+02, 0.6810260E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(262,oi),oi=1,5) / 0.1307500E+02, 0.2524219E+02, 0.1805942E+02, 0.1051234E+01, 0.0000000E+00 /
      data (oab13a(263,oi),oi=1,5) / 0.1312500E+02, 0.1818278E+02, 0.0000000E+00, 0.7746802E+00, 0.0000000E+00 /
      data (oab13a(264,oi),oi=1,5) / 0.1317500E+02, 0.5871736E+02, 0.2206884E+02, 0.9083562E+00, 0.0000000E+00 /
      data (oab13a(265,oi),oi=1,5) / 0.1322500E+02, 0.5547953E+02, 0.0000000E+00, 0.3445293E+01, 0.3445293E+01 /
      data (oab13a(266,oi),oi=1,5) / 0.1327500E+02, 0.4110425E+02, 0.1115897E+02, 0.1215118E+01, 0.0000000E+00 /
      data (oab13a(267,oi),oi=1,5) / 0.1332500E+02, 0.2626710E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(268,oi),oi=1,5) / 0.1337500E+02, 0.3026880E+02, 0.1136909E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(269,oi),oi=1,5) / 0.1342500E+02, 0.5338152E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(270,oi),oi=1,5) / 0.1347500E+02, 0.8456583E+02, 0.7688158E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(271,oi),oi=1,5) / 0.1352500E+02, 0.7354125E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(272,oi),oi=1,5) / 0.1357500E+02, 0.6240967E+02, 0.0000000E+00, 0.1132989E+02, 0.1058770E+02 /
      data (oab13a(273,oi),oi=1,5) / 0.1362500E+02, 0.2356545E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(274,oi),oi=1,5) / 0.1367500E+02, 0.1977738E+02, 0.7936710E+01, 0.9500106E+00, 0.0000000E+00 /
      data (oab13a(275,oi),oi=1,5) / 0.1372500E+02, 0.1993599E+02, 0.0000000E+00, 0.9727378E+00, 0.0000000E+00 /
      data (oab13a(276,oi),oi=1,5) / 0.1377500E+02, 0.6837186E+02, 0.0000000E+00, 0.1886865E+01, 0.0000000E+00 /
      data (oab13a(277,oi),oi=1,5) / 0.1382500E+02, 0.8491045E+02, 0.8083918E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(278,oi),oi=1,5) / 0.1387500E+02, 0.4072248E+02, 0.8157866E+01, 0.6148209E+01, 0.6148209E+01 /
      data (oab13a(279,oi),oi=1,5) / 0.1392500E+02, 0.3682977E+02, 0.8179570E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(280,oi),oi=1,5) / 0.1397500E+02, 0.3721707E+02, 0.0000000E+00, 0.7946882E+00, 0.0000000E+00 /
      data (oab13a(281,oi),oi=1,5) / 0.1402500E+02, 0.2913017E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(282,oi),oi=1,5) / 0.1407500E+02, 0.8370006E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(283,oi),oi=1,5) / 0.1412500E+02, 0.6748445E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(284,oi),oi=1,5) / 0.1417500E+02, 0.3398480E+02, 0.8498802E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(285,oi),oi=1,5) / 0.1422500E+02, 0.4705011E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(286,oi),oi=1,5) / 0.1427500E+02, 0.6033557E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(287,oi),oi=1,5) / 0.1432500E+02, 0.4336177E+02, 0.8687794E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(288,oi),oi=1,5) / 0.1437500E+02, 0.5246232E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(289,oi),oi=1,5) / 0.1442500E+02, 0.1760347E+02, 0.1762346E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(290,oi),oi=1,5) / 0.1447500E+02, 0.4876277E+02, 0.0000000E+00, 0.2039261E+01, 0.0000000E+00 /
      data (oab13a(291,oi),oi=1,5) / 0.1452500E+02, 0.4468342E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(292,oi),oi=1,5) / 0.1457500E+02, 0.7189564E+02, 0.0000000E+00, 0.8687226E+00, 0.0000000E+00 /
      data (oab13a(293,oi),oi=1,5) / 0.1462500E+02, 0.1354590E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(294,oi),oi=1,5) / 0.1467500E+02, 0.1821360E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(295,oi),oi=1,5) / 0.1472500E+02, 0.3207491E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(296,oi),oi=1,5) / 0.1477500E+02, 0.0000000E+00, 0.1386696E+02, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(297,oi),oi=1,5) / 0.1482500E+02, 0.2322306E+02, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(298,oi),oi=1,5) / 0.1487500E+02, 0.9362142E+01, 0.9384598E+01, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(299,oi),oi=1,5) / 0.1492500E+02, 0.9443782E+01, 0.0000000E+00, 0.0000000E+00, 0.0000000E+00 /
      data (oab13a(300,oi),oi=1,5) / 0.1497500E+02, 0.2846616E+02, 0.0000000E+00, 0.1035417E+01, 0.0000000E+00 /