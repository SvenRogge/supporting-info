%chk=gaussian.chk
%nproc=20
%mem=90GB
#p opt b3lyp/gen int=grid=ultrafine pseudo=read scf=(tight,maxcycle=5000)

Title Card Required

0 1
Zr  9.171771025 -1.692510032  0.582031015
 O  8.853050027 -0.343765032 -0.976597983
 O  9.960886005 -3.027187034 -1.054891984
 C  9.682806028 -3.149381036 -2.277843981
 O  7.377803032 -2.352269025 -0.679058985
 H  7.374336021 -3.272642024 -0.953378983
 O  7.384587031 -1.036934027  1.459400016
 O  11.329661010 -1.084632044  0.328330017
 C  11.885687030 -0.016967045 -0.040649983
 O  8.497712029 -3.705651028  1.323105018
 C  7.363431002 -4.220711026  1.502085016
 O  9.872246018 -1.793960034  2.704186018
 C  9.568842001 -1.153059034  3.744600015
 O  9.387508000  0.461663965  1.339471019
 H  10.170821028  0.641947962  1.865203019
Zr  9.179346038  1.689316966 -0.576369983
 O  9.956298038  1.687338960 -2.678587981
 O  11.330206010  1.053511958 -0.403974983
 O  9.896282040  3.081167962  1.043982017
 O  8.529974058  3.745165965 -1.238836982
Zr  5.588560009 -1.672660018  0.588283016
 O  4.793679009 -3.014256012 -1.035098983
 O  6.236022029 -3.686579021  1.335485018
 O  4.877887015 -1.759824014  2.715857016
 O  3.461862019 -1.055264011  0.305323017
Zr  7.380555018 -0.845320024 -2.377660979
 O  7.392789013  1.424951971 -1.978200980
 H  7.394141008  1.987295975 -2.756866980
 O  8.821579011 -2.504917033 -2.933291978
 O  5.944983031  0.101677979 -3.808005981
 C  5.096004028  1.024370984 -3.695593983
 O  5.918926028 -2.497093021 -2.923771980
 C  5.062029010 -3.138121013 -2.259964983
 O  8.819869028  0.089527968 -3.808940979
 C  9.680565015  1.001222963 -3.696554979
 O  5.921040038 -0.331324020 -0.975138983
Zr  5.608137026  1.702927981 -0.574394983
 O  7.395212009  1.738209972  0.510479017
 O  6.270339012  3.754711978 -1.234406980
 C  7.402061044  4.283381970 -1.396089980
 O  4.896550036  3.101745980  1.046174017
 C  5.197589034  3.232572980  2.262500017
 O  4.835451030  1.719492985 -2.679836982
 O  5.392256008  0.481377982  1.342160017
 H  4.602804027  0.666578986  1.857084017
 O  3.473317028  1.070585988 -0.426948982
 C  2.883955025 -0.002046008 -0.097651982
 C  1.391965022 -0.015417002 -0.156506983
 C  0.691610020 -1.165211996  0.223736017
 H  1.247097015 -2.042515998  0.527688017
 C  0.696621027  1.129350000 -0.561076984
 H  1.255935031  2.007031999 -0.857087983
Zr  7.390999018  0.835201975  2.395013017
 O  6.019922038 -0.179179021  3.860302014
 C  5.185373020 -1.115336017  3.752910013
 O  8.772857031  2.534346964  2.927153015
 C  9.602751041  3.209272963  2.261759015
 O  6.026374031  2.554286980  2.925206018
 O  8.753272023 -0.201018032  3.855268014
Zr -5.581521000 -1.668410974  0.591379015
 O -5.381000991  0.478585026  1.336338015
 H -4.590127968  0.657173025  1.851413015
 O -3.459557976 -1.051233983  0.290960018
 C -2.883596973  0.001671016 -0.115703983
 C -1.391091975 -0.012715990 -0.167120982
 C -0.695246981 -1.163689992  0.218643017
 H -1.254609982 -2.039695991  0.518974017
 C -0.690995968  1.130670006 -0.566616983
 H -1.246404967  2.009360008 -0.867134984
 O -4.797539986 -3.014794974 -1.032088981
 C -5.069095986 -3.133745972 -2.256178979
 O -7.372795004 -2.347170963 -0.665975983
 H -7.366002010 -3.268661966 -0.936082980
 O -5.914430985 -0.329696971 -0.976870980
 O -4.852458981 -1.767776976  2.715743016
 C -5.162238975 -1.137086973  3.762078014
 O -7.377540983 -1.036994967  1.467541016
 O -6.246495981 -3.670066967  1.347591016
 C -7.379126980 -4.186432964  1.523668016
Zr -5.613459966  1.703483024 -0.592256984
 O -6.308725950  3.760108030 -1.262260984
 O -4.850393967  1.732054022 -2.713281982
 O -4.910855961  3.097930025  1.018793017
 O -3.473995968  1.071904020 -0.451470982
Zr -7.374057992 -0.847394964 -2.372693979
 O -8.853643975 -0.348693959 -0.971551982
 O -5.953502996  0.099909029 -3.823890981
 C -5.125003965  1.043044023 -3.730492982
 O -8.807796961  0.065026040 -3.805599981
 C -9.673417948  0.972463045 -3.697465979
 O -7.396982954  1.417931033 -1.984812980
 H -7.402492959  1.975551035 -2.766641979
 O -5.922802992 -2.488543972 -2.919253981
 O -8.806118993 -2.500056958 -2.909030981
 C -9.665731966 -3.137247955 -2.247160979
Zr -9.176624003 -1.678562958  0.597939019
 O -9.386651950  0.461340043  1.345512016
 H -10.173074952  0.641364046  1.866500016
 O -9.888549967 -1.764534956  2.726758014
 C -9.587769981 -1.110742955  3.759408013
 O -11.346245953 -1.079529951  0.336367017
 C -11.902261971 -0.018092947 -0.050087983
 O -8.511301988 -3.666646958  1.347155016
 O -9.953210977 -3.008992953 -1.027990980
Zr -7.374949973  0.825513032  2.396820015
 O -8.741765965  2.522428041  2.920184016
 C -9.580512951  3.188925041  2.258975014
 O -6.012960977  2.533650027  2.912038015
 C -5.193987961  3.213215021  2.240612015
 O -7.391976990  1.724860032  0.503896017
 O -8.741205989 -0.185075960  3.867178015
 O -5.992408951 -0.197215971  3.871779015
Zr -9.184604948  1.676714040 -0.576814984
 O -9.883564959  3.065434041  1.044150014
 O -9.962353942  1.656530046 -2.683307983
 O -8.568137948  3.716782037 -1.257285983
 C -7.447878953  4.264389033 -1.440375984
 O -11.342011954  1.054150051 -0.402969982
 H -7.482513972  5.295021019 -1.820691982
 H -4.646527956  3.984372023  2.791673018
 H -4.597588965  1.307508022 -4.653343980
 H -4.507845987 -3.905333973 -2.800374983
 H -7.385741007 -5.225245962  1.846821017
 H -10.217099004 -3.911411954 -2.790368983
 H -10.223787962  1.202440046 -4.608656981
 H -13.003738986 -0.016105942 -0.064286983
 H -4.629894979 -1.436120976  4.678008015
 H -10.146454958 -1.364178954  4.667198014
 H -10.128091979  3.955040043  2.814722014
 H  7.404711058  5.331356972 -1.731166983
 H  4.658992037  4.018548981  2.812967014
 H  12.984276007 -0.008220050 -0.023639983
 H  7.353497023 -5.268456022  1.817952019
 H  10.093949995 -1.451759036  4.657714016
 H  4.647089016 -1.391814016  4.666568013
 H  10.244354011  1.237926961 -4.601565980
 H  4.526943030  1.264787985 -4.596500982
 H  10.245162012 -3.921495034 -2.828938979
 H  4.490915005 -3.911020014 -2.798897981
 H  10.149708022  3.988133960  2.812185017

H C O 0
6-311g(d,p)
****
Zr 0
LanL2DZ
****

Zr 0
LanL2DZ ecp












