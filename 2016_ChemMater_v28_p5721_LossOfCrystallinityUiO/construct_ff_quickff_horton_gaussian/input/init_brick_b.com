%chk=gaussian.chk
%nproc=24
%mem=62GB
#p opt b3lyp/gen gfinput int=grid=ultrafine iop(6/7=3) pseudo=read scf=(tight,maxcycle=5000)

Title Card Required

0 1
 C    -1     -6.500102    -6.124534     0.312768
 H    -1     -3.059748     5.764947     0.345002
 H    -1      5.552769    -2.782691     0.293150
 H    -1      4.414686     7.640618     0.312793
 O    -1     -3.874788    -1.845456     0.321544
 C    -1     -3.450735    -0.098574    -2.922054
 O    -1      1.681253     0.024652    -3.312232
 H    -1     -3.120299    -0.047169    -5.603622
 H    -1      5.557004    -0.011273    -2.543048
 O    -1      0.869441    -1.036932     1.283326
 O    -1     -0.412708    -3.591843    -1.655734
 C    -1      2.894020     0.021543    -2.920105
 H    -1      1.744174     1.963290     2.230588
 H    -1     -4.868502    -7.526494     0.307642
 C    -1     -0.318105     3.778865     5.881642
 H    -1      7.341711    -0.084026    -4.317182
 C    -1     -0.384337    -3.207493    -2.869335
 H    -1     -0.316945     5.747787     3.102905
 O    -1      3.341783     1.992657     0.157913
 C    -1      5.369447     3.874242     0.241786
 O    -1     -0.344149    -1.996310    -3.275959
 O    -1     -2.237554    -0.068432    -3.326870
 C    -1     -0.361299    -4.257977    -3.933090
 C    -1     -0.342999    -5.508015     4.143991
 C    -1      3.946461    -4.212826     0.297251
 O    -1     -0.314639    -1.872204     3.874170
 C    -1     -0.249382     5.563047    -3.602280
 O    -1      3.370030    -0.032108     2.185865
 H    -1     -0.283307    -7.645429    -4.337133
 C    -1      2.963708    -0.068170     3.398153
 O    -1      1.747812    -0.089742     3.772964
 C    -1      5.946316     0.031116     6.483056
 H    -1     -6.117047    -0.080407    -2.531863
 H    -1     -0.369388    -5.890747    -2.534235
 C    -1      6.361497     4.844295     0.293697
 H    -1     -0.232452     2.777621    -5.565006
 C    -1     -0.265570    -6.204173    -5.933318
 C    -1      4.662735     6.580080     0.293486
 C    -1     -0.334147    -4.157749     4.524411
 C    -1     -0.281351     4.812546    -6.297703
 C    -1     -4.541082    -4.126591     0.306705
 C    -1     -0.359750    -4.816676     6.851488
 H    -1     -0.315863    -2.765984     6.163885
 O    -1      1.099145    -1.435630    -1.135689
 C    -1     -0.236042     6.475279     5.135950
 C    -1     -5.141657    -6.471704     0.309190
 H    -1     -0.293634     4.546637    -7.354215
 H    -1     -6.123865     2.757298     0.342535
 H    -1     -4.814661     7.557954     0.220102
 C    -1      5.300952    -0.023171    -3.601695
 O    -1      0.735947     0.990099    -0.757774
 C    -1      5.923111    -0.078300    -5.933554
 C    -1      6.284746    -0.063755    -4.579758
 C    -1     -6.872236    -4.772273     0.303474
 C    -1     -0.302321     5.489587     4.160843
 H    -1     -0.340731    -2.817895    -5.530361
 C    -1      4.019769     4.252967     0.216605
 C    -1      5.933920    -6.170985     0.204502
 C    -1     -6.519339     0.036869    -5.920237
 C    -1     -0.384204     3.091446     3.457245
 H    -1      5.651978    -0.031253     3.081514
 H    -1     -0.327969    -5.754388     3.082814
 H    -1     -6.172749    -2.725122     0.286371
 O    -1     -0.474603     1.862025     3.812178
 O    -1     -0.236168     3.546412    -1.685442
 C    -1     -0.299223     6.167771    -5.941214
 C    -1     -0.337059    -3.818816     5.884407
 C    -1      3.622740    -0.040828     5.817376
 C    -1     -6.445221     6.155793     0.214258
 C    -1     -6.863527     0.000987    -4.561731
 H    -1     -0.236543     5.829546    -2.546017
 C    -1     -5.868801    -0.052054    -3.592416
 O    -1      1.743615     3.603941     0.204786
 H    -1      4.303370     0.008798     7.870601
 C    -1      2.889422    -3.155945     0.317414
Zr    -1      2.303177     0.027628     0.221529
Zr    -1     -0.312834    -2.431635     0.311667
 H    -1     -4.920433     0.047521    -7.357473
 C    -1     -5.168691     0.019500    -6.297013
 O    -1      1.194697     1.424593     1.643578
 C    -1     -0.303470    -3.089034     3.479842
 C    -1      6.332253     0.019612     5.135000
 H    -1     -0.223277     4.509370     7.916818
 C    -1     -0.281666    -4.845978    -6.277928
 C    -1      4.582038    -6.542922     0.226388
 C    -1      4.571624    -0.052247    -6.303409
 H    -1     -0.203177     7.531047     4.868863
 H    -1     -2.270890     1.918526    -1.693378
 H    -1      5.618565     2.814000     0.221887
 C    -1     -5.899210    -3.779651     0.300022
 O    -1     -1.294061     1.158727     1.316117
 C    -1     -0.283763     6.538938    -4.589598
 C    -1      3.949416     0.001471    -3.971405
 H    -1     -0.343504     2.723148     6.150487
 C    -1      6.012019     6.201068     0.318472
Zr    -1     -0.273420    -0.026764    -2.263773
 C    -1     -5.088055     6.503706     0.249167
 H    -1     -7.879255     4.553165     0.209907
 H    -1      2.536676    -5.836055     0.290015
 H    -1      4.315895    -0.063364    -7.362207
 H    -1     -7.917545     0.015021    -4.285735
 C    -1     -6.819710     4.804621     0.244693
 O    -1     -1.687914     1.373061    -1.146324
 C    -1      3.587692    -0.011967    -5.325178
 C    -1     -4.118911     5.512476     0.321632
 O    -1     -2.228625     3.531883     0.316189
 C    -1     -0.378259    -6.166462     6.471151
 O    -1     -0.241273    -3.470407     2.267815
 C    -1     -0.214292     3.156552    -2.899406
Zr    -1     -2.616767    -0.023754     0.236595
 H    -1      2.560486    -0.064470     6.058419
 O    -1      1.673471    -3.533162     0.338639
 C    -1     -0.211749     6.119332     6.491258
 H    -1      7.394643     0.043060     4.894063
 C    -1     -4.174308    -0.033289    -5.328018
 C    -1      5.369161    -0.021880     4.133756
Zr    -1     -0.222964     2.537768     0.251444
 C    -1      3.592784    -5.568344     0.273708
 C    -1     -5.848837     3.811494     0.321153
 C    -1      2.955571     3.207681     0.183587
 H    -1     -7.932045    -4.519757     0.296797
 H    -1      7.343682    -4.547287     0.213230
 O    -1     -1.669446    -1.213870     1.677734
 O    -1     -3.827238    -0.124478    -1.707068
 H    -1     -0.297424     7.598239    -4.335814
 O    -1     -0.281886     3.519378     2.261050
 C    -1     -3.427891     3.110459     0.410392
 C    -1     -4.492064     4.160079     0.364002
 O    -1      3.283531     0.014539    -1.706464
 O    -1     -3.789263     1.882904     0.500846
 C    -1      6.287729    -4.815363     0.229491
 C    -1     -0.248486     4.766336     6.858403
 O    -1      3.303741    -1.947244     0.282006
 H    -1      1.636858    -2.032151    -1.675747
 C    -1      3.670627     5.609926     0.242967
 O    -1     -2.279304    -3.453840     0.233134
 C    -1     -0.328962    -3.878530    -5.282962
 C    -1      4.008225    -0.051167     4.469199
 C    -1     -0.348149     4.135347     4.526727
 C    -1     -0.366390    -6.506083     5.110776
 H    -1      2.530563     0.008335    -5.587477
 H    -1      4.329109    -7.602359     0.206592
 O    -1     -0.200111     1.943306    -3.292237
 C    -1     -4.169086    -5.479159     0.304747
 C    -1     -0.345448    -5.616042    -3.588406
 C    -1     -3.488859    -3.064833     0.295639
 C    -1      4.585568     0.000395     6.818307
 H    -1     -0.373080    -7.559345     4.831669
 H    -1      2.615925     5.882118     0.223290
 C    -1     -0.296695    -6.584403    -4.583833
Zr    -1     -0.220235    -0.039480     2.630772
 O    -1     -1.342514    -1.092139    -0.836033
 C    -1     -4.518230    -0.068515    -3.969605
 C    -1     -0.246534     3.836715    -5.310408
 C    -1     -0.230751     4.207728    -3.958401
 H    -1     -3.109357    -5.731857     0.296748
 H    -1     -0.257691    -4.570685    -7.331871
 H    -1      7.416385     4.572754     0.314330
 C    -1      5.298857    -3.841682     0.275158
 H    -1     -0.360714    -4.570576     7.912801
 H     0     -7.256842     6.592889     0.369949
 H     0      6.817754     6.927698     0.456312
 H     0      6.799707    -6.847859     0.322024
 H     0     -7.047910    -6.806295     0.369425
 H     0     -0.278192     6.884940     7.249030
 H     0     -0.152114    -6.926720    -6.583336
 H     0     -0.457319     6.758733    -6.560121
 H     0     -0.468903    -7.018895     7.400488
 H     0      6.618375     0.022286    -6.685727
 H     0      6.761820     0.207729     7.289460
 H     0     -7.208513     0.110372    -6.616637

H C O 0
6-31g(d)
****
Zr 0
LanL2DZ
****

Zr 0
LanL2DZ ecp












