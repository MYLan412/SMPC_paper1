clear all;
close all;
%% [e(k) Δx(k)]  20*20 这是例1中预见长度为 0 的程序
I_N = eye(5);
I_m = eye(1);
I_mN = eye(5);
I_s = eye(20);

A1 =  [  0.0571    0.0200    0.7000
        -0.0571    0.0400    0.0100
        0.0013    -0.0083    0.1000 ] ;
A2 = [ 0.0571     0.0200    0.7000
       -0.0571    0.0400    0.0100
       0.0026    -0.0165    0.1000] ;
B = [  2;
       0;
       0 ];
  
A1_d = [    0.0063    0.0040    0.0040
           -0.0063    0.0500    0.0800
            0.0023    0.0075   -0.2000 ];
     
A2_d = [    0.0063    0.0020    0.0060
           -0.0063    0.0200    0.0700
            0.0047    0.0052   -0.2000 ];
     
C = [-0.02  -0.01   -0.01 ];
     
L = [ 4  -1  -1  -1  -1 ;
      0  3  -1  -1  -1;
     -1  0  3   -1   -1;
      -1  -1  -1   4  -1;
     -1  -1  -1  -1   4];
M = [ 1 0 0 0 0;
      0 1 0 0 0;
      0 0 1 0 0;
      0 0 0 1 0;
      0 0 0 0 1 ];
H = L + M;

A_R0_1 = [    I_N  kron(-H,C*A1);
            zeros(15,5)  kron(I_N,A1)  ];
A_R0_1_d = [ zeros(5)  kron(-H,C*A1_d);
             zeros(15,5)  kron(I_N,A1_d) ];

A_R0_2 = [    I_N  kron(-H,C*A2);
            zeros(15,5)  kron(I_N,A2)   ];
A_R0_2_d = [ zeros(5)  kron(-H,C*A2_d);
             zeros(15,5)  kron(I_N,A2_d) ];

B_R0 = [    kron(-H,C*B);
            kron(I_N,B) ];
G_hat = eye(20) - B_R0*inv(B_R0'*B_R0)*B_R0';

F = [ kron(-H,C);
      eye(15)];
D = [ eye(5)  zeros(5,15) ];

%%
setlmis([]);

%%
% 
gama = lmivar(2,[1,1]);
% 表示矩阵
P0 = lmivar(1,[20,1]); % P
Q0 = lmivar(1,[20,1]); % Q
R0 = lmivar(1,[20,1]);  % R
Y01 = lmivar(2,[5,20]); % Y1
Y02 = lmivar(2,[5,20]); % Y2


%% 第一个大矩阵 < 0
% (1,1)
lmiterm([1 1 1 Q0],1,1);
% lmiterm([1 1 1 R0'],-(A_R0_1-eye(35)),1);
lmiterm([1 1 1 Y01],-B_R0,1,'s');
lmiterm([1 1 1 R0],-1,(A_R0_1-I_s)','s');
% lmiterm([1 1 1 Y01'],-1,B_R0');

% (1,2)
lmiterm([1 1 2 P0],1,1);
lmiterm([1 1 2 R0'],1,1);
lmiterm([1 1 2 -R0],1,(A_R0_1-I_s)');
lmiterm([1 1 2 -Y01'],1,B_R0');

% (1,3)
lmiterm([1 1 3 R0'],-G_hat*A_R0_1_d,1);

% (1,4)
lmiterm([1 1 4 0],-G_hat*F);

% (1,5)
lmiterm([1 1 5 R0],1,D');

% (2,2)
lmiterm([1 2 2 P0],1,1);
lmiterm([1 2 2 R0],1,1,'s');
% lmiterm([1 2 2 R0'],1,1);

% (2,3)
lmiterm([1 2 3 R0'],-G_hat*A_R0_1_d,1);

% (2,4)
lmiterm([1 2 4 0],-G_hat*F);

% (2,5)
lmiterm([1 2 5 0],0);

% (3,3)
lmiterm([1 3 3 Q0],-1,1);

% (3,4)
lmiterm([1 3 4 0],0);

% (3,5)
lmiterm([1 3 5 0],0);

% (4,4)
lmiterm([1 4 4 gama^2],-1,eye(15));

% (4,5)
lmiterm([1 4 5 0],0);

% (5,5)
lmiterm([1 5 5 0],-eye(5));

%% 第二个大矩阵 < 0
% (1,1)
lmiterm([2 1 1 Q0],1,1);
% lmiterm([1 1 1 R0'],-(A_R0_1-eye(35)),1);
lmiterm([2 1 1 Y02],-B_R0,1,'s');
lmiterm([2 1 1 R0],-1,(A_R0_2-I_s)','s');
% lmiterm([1 1 1 Y01'],-1,B_R0');

% (1,2)
lmiterm([2 1 2 P0],1,1);
lmiterm([2 1 2 R0'],1,1);
lmiterm([2 1 2 -R0],1,(A_R0_2-I_s)');
lmiterm([2 1 2 -Y02'],1,B_R0');

% (1,3)
lmiterm([2 1 3 R0'],-G_hat*A_R0_2_d,1);

% (1,4)
lmiterm([2 1 4 0],-G_hat*F);

% (1,5)
lmiterm([2 1 5 R0],1,D');

% (2,2)
lmiterm([2 2 2 P0],1,1);
lmiterm([2 2 2 R0],1,1,'s');
% lmiterm([1 2 2 R0'],1,1);

% (2,3)
lmiterm([2 2 3 R0'],-G_hat*A_R0_2_d,1);

% (2,4)
lmiterm([2 2 4 0],-G_hat*F);

% (2,5)
lmiterm([2 2 5 0],0);

% (3,3)
lmiterm([2 3 3 Q0],-1,1);

% (3,4)
lmiterm([2 3 4 0],0);

% (3,5)
lmiterm([2 3 5 0],0);

% (4,4)
lmiterm([2 4 4 gama^2],-1,eye(15));

% (4,5)
lmiterm([2 4 5 0],0);

% (5,5)
lmiterm([2 5 5 0],-eye(5));

%% 第三个大矩阵 < 0
% (1,1)
lmiterm([3 1 1 Q0],2,1);
lmiterm([3 1 1 R0'],-(3/2)*(A_R0_1-I_s),1,'s');
lmiterm([3 1 1 Y02],-(1/2)*B_R0,1,'s');
lmiterm([3 1 1 R0'],-(1/2)*(A_R0_2-I_s),1,'s');
lmiterm([3 1 1 Y01],-(3/2)*B_R0,1,'s');

% (1,2)
lmiterm([3 1 2 P0],2,1);
lmiterm([3 1 2 R0'],2,1);
lmiterm([3 1 2 -R0],3/2,(A_R0_1-I_s)');
lmiterm([3 1 2 -Y02'],1/2,B_R0');
lmiterm([3 1 2 -R0],1/2,(A_R0_2-I_s)');
lmiterm([3 1 2 -Y01'],3/2,B_R0');

% (1,3)
lmiterm([3 1 3 R0'],-(3/2)*G_hat*A_R0_1_d,1);
lmiterm([3 1 3 R0'],-(1/2)*G_hat*A_R0_2_d,1);

% (1,4)
lmiterm([3 1 4 0],-2*G_hat*F);

% (1,5)
lmiterm([3 1 5 R0],2,D');

% (2,2)
lmiterm([3 2 2 P0],2,1);
lmiterm([3 2 2 R0],2,1,'s');

% (2,3)
lmiterm([3 2 3 R0'],-(3/2)*G_hat*A_R0_1_d,1);
lmiterm([3 2 3 R0'],-(1/2)*G_hat*A_R0_2_d,1);

% (2,4)
lmiterm([3 2 4 0],-2*G_hat*F);

% (2,5)
lmiterm([3 2 5 0],0);

% (3,3)
lmiterm([3 3 3 Q0],-2,1);

% (3,4)
lmiterm([3 3 4 0],0);

% (3,5)
lmiterm([3 3 5 0],0);

% (4,4)
lmiterm([3 4 4 gama^2],-2,eye(15));

% (4,5)
lmiterm([3 4 5 0],0);

% (5,5)
lmiterm([3 5 5 0],-2*eye(5));

%% 第四个大矩阵 < 0
% (1,1)
lmiterm([4 1 1 Q0],2,1);
lmiterm([4 1 1 R0'],-(1/2)*(A_R0_1-I_s),1,'s');
lmiterm([4 1 1 Y02],-(3/2)*B_R0,1,'s');
lmiterm([4 1 1 R0'],-(3/2)*(A_R0_2-I_s),1,'s');
lmiterm([4 1 1 Y01],-(1/2)*B_R0,1,'s');

% (1,2)
lmiterm([4 1 2 P0],2,1);
lmiterm([4 1 2 R0'],2,1);
lmiterm([4 1 2 -R0],1/2,(A_R0_1-I_s)');
lmiterm([4 1 2 -Y02'],3/2,B_R0');
lmiterm([4 1 2 -R0],3/2,(A_R0_2-I_s)');
lmiterm([4 1 2 -Y01'],1/2,B_R0');

% (1,3)
lmiterm([4 1 3 R0'],-(1/2)*G_hat*A_R0_1_d,1);
lmiterm([4 1 3 R0'],-(3/2)*G_hat*A_R0_2_d,1);

% (1,4)
lmiterm([4 1 4 0],-2*G_hat*F);

% (1,5)
lmiterm([4 1 5 R0],2,D');

% (2,2)
lmiterm([4 2 2 P0],2,1);
lmiterm([4 2 2 R0],2,1,'s');

% (2,3)
lmiterm([4 2 3 R0'],-(1/2)*G_hat*A_R0_1_d,1);
lmiterm([4 2 3 R0'],-(3/2)*G_hat*A_R0_2_d,1);

% (2,4)
lmiterm([4 2 4 0],-2*G_hat*F);

% (2,5)
lmiterm([4 2 5 0],0);

% (3,3)
lmiterm([4 3 3 Q0],-2,1);

% (3,4)
lmiterm([4 3 4 0],0);

% (3,5)
lmiterm([4 3 5 0],0);

% (4,4)
lmiterm([4 4 4 gama^2],-2,eye(15));

% (4,5)
lmiterm([4 4 5 0],0);

% (5,5)
lmiterm([4 5 5 0],-2*eye(5));


%% 约束条件
lmiterm([-5 1 1 P0],1,1);
lmiterm([5 1 1 0],0);

lmiterm([-6 1 1 gama],1,1);
lmiterm([6 1 1 0],0);

lmiterm([-7 1 1 Q0],1,1);
lmiterm([7 1 1 0],0);

%%
lmisys=getlmis;
options=[0,0,0,0,0];
[tmin,xfeas]=feasp(lmisys,options); 

%%
gama = dec2mat(lmisys,xfeas,gama);
P0 = dec2mat(lmisys,xfeas,P0);
Q0 = dec2mat(lmisys,xfeas,Q0);
R0 = dec2mat(lmisys,xfeas,R0);
Y01 = dec2mat(lmisys,xfeas,Y01);
Y02 = dec2mat(lmisys,xfeas,Y02);

K1_20 = Y01 * inv(R0')
K2_20 = Y02 * inv(R0')

%K1_20 = [12.4301150967106,4.28357658543951,4.28358839425510,4.28358048300257,4.28358523003061,53.4182748063869,-9.29568381072098e-07,1.12728303442368e-10,6.59056850770293,-6.89834762503621e-07,8.54948369749455e-11,6.59053814009376,-6.90025108634331e-07,8.55216066544410e-11,6.59052338089769,-6.89881299415090e-07,8.54659665137181e-11,6.59053657665800,-6.89708596400645e-07,8.54905390119279e-11,-1.69167421211824e-06,-2.13118033335455e-07,-2.13118766096860e-07,-2.13142922589421e-07,-2.13130667603088e-07,6.22932226387772,3.24924107111279,3.24923988737876,3.24924006410111,3.24924017376708,2.27672244078404,2.11586439407449,2.11586211982326,2.11586176697022,2.11586189087538,1.45260696786980,1.38937609618120,1.38937507719869,1.38937509712393,1.38937514111805;4.28357945689669,12.4301033999621,4.28358934417837,4.28357984499292,4.28358537716934,6.59057272980418,-6.88529463066458e-07,8.53409907058262e-11,53.4182654096822,-9.29510767534888e-07,1.12679894032322e-10,6.59053119425992,-6.89345351375015e-07,8.54196052134332e-11,6.59052277433066,-6.89199027789711e-07,8.53729143760759e-11,6.59053635943333,-6.89094075482999e-07,8.53934033943873e-11,-2.13147188209560e-07,-1.69167604702140e-06,-2.13096498939083e-07,-2.13138143788974e-07,-2.13121300775494e-07,3.24924083701686,6.22932104953379,3.24923942162442,3.24923976615444,3.24924010299668,2.11586368593385,2.27672148485006,2.11586169296881,2.11586157367154,2.11586177482457,1.38937609306314,1.45260650470371,1.38937472802283,1.38937508576433,1.38937506883177;4.28358882766332,4.28358689624743,12.4301129305282,4.28357742311747,4.28358152927641,6.59054484531414,-6.89015539756117e-07,8.54257215312888e-11,6.59053380863721,-6.89433908529479e-07,8.54450732913074e-11,53.4182478153317,-9.30134678329994e-07,1.12777021900394e-10,6.59054545391934,-6.89541407775720e-07,8.54187934069818e-11,6.59056336188920,-6.89370441162247e-07,8.54565778323052e-11,-2.13172795781895e-07,-2.13125740928016e-07,-1.69164360962823e-06,-2.13115406207597e-07,-2.13108979142515e-07,3.24924016018317,3.24923988906747,6.22932179136269,3.24924194894044,3.24924235619388,2.11586240894051,2.11586253663200,2.27672214292810,2.11586495095526,2.11586504693435,1.38937648891222,1.38937609773549,1.45260680503634,1.38937660196321,1.38937642557713;4.28358505163521,4.28358186885840,4.28358165236671,12.4301029311937,4.28358095243782,6.59052901040815,-6.88658029434043e-07,8.53454112153725e-11,6.59052225211157,-6.88823800287632e-07,8.53449897610403e-11,6.59054407896424,-6.89376995336629e-07,8.54165336650264e-11,53.4182517175543,-9.29494930061417e-07,1.12647152991763e-10,6.59057145935502,-6.89053666835858e-07,8.53945012403303e-11,-2.13176735651280e-07,-2.13132072242791e-07,-2.13099963144468e-07,-1.69167226066541e-06,-2.13101069977072e-07,3.24923919448758,3.24923917596355,3.24924082833401,6.22932146472318,3.24924159321806,2.11586031047887,2.11586076463133,2.11586304523615,2.27672158099176,2.11586437931993,1.38937422471597,1.38937398127440,1.38937428790627,1.45260635684575,1.38937558734411;4.28358708554607,4.28358447719547,4.28358312065060,4.28357829972426,12.4301088544589,6.59053626722253,-6.88836391707049e-07,8.53890404590291e-11,6.59053085162927,-6.89332467702269e-07,8.54209425490781e-11,6.59055505265351,-6.89420770240509e-07,8.54297941010044e-11,6.59056449266092,-6.89243186378354e-07,8.53565105315825e-11,53.4182437532048,-9.29578249312249e-07,1.12704669654467e-10,-2.13174399178546e-07,-2.13134510220516e-07,-2.13093077422796e-07,-2.13107026985101e-07,-1.69166566493284e-06,3.24923979113982,3.24923992690329,3.24924168748662,3.24924200205939,6.22932122043319,2.11586138413119,2.11586190195639,2.11586429809370,2.11586535495744,2.27672111785942,1.38937518180169,1.38937507679286,1.38937522028136,1.38937664517758,1.45260652113109];


%K2_20 = [12.4301061321777,4.28357788573038,4.28358342750278,4.28358285214037,4.28358420676128,53.4182669944473,-9.40979119855426e-07,1.14492834870123e-10,6.59056928724113,-7.03499131397056e-07,8.73341210750688e-11,6.59053758535373,-7.03586006417759e-07,8.73439254330241e-11,6.59052386379371,-7.03607318324428e-07,8.73461698699374e-11,6.59053509605017,-7.03556955929288e-07,8.73414913420893e-11,-1.69167810346176e-06,-2.13120290760996e-07,-2.13124367016720e-07,-2.13144000243294e-07,-2.13132243913392e-07,6.22932174134406,3.24924085218087,3.24923956025282,3.24923984302743,3.24923984608680,2.27672182810505,2.11586386007915,2.11586147666480,2.11586125382871,2.11586132250180,1.45260677779044,1.38937584579802,1.38937477797824,1.38937485354168,1.38937486975180;4.28357909083183,12.4301043695020,4.28358451829635,4.28358303729239,4.28358451559940,6.59056731426397,-7.03482855696907e-07,8.73321041731883e-11,53.4182671241091,-9.40880911221422e-07,1.14482410667543e-10,6.59053169020461,-7.03700151863679e-07,8.73571108420950e-11,6.59052447315522,-7.03709629520349e-07,8.73582442433668e-11,6.59053593236521,-7.03561018840415e-07,8.73419190585605e-11,-2.13145225728781e-07,-1.69167750816245e-06,-2.13101045722947e-07,-2.13137786932784e-07,-2.13121803522099e-07,3.24924080724641,6.22932117994223,3.24923945154210,3.24923991309882,3.24924013253176,2.11586374244328,2.27672146037386,2.11586156662848,2.11586158811286,2.11586172372435,1.38937629415479,1.45260657043919,1.38937474713979,1.38937516395557,1.38937511592573;4.28358398331187,4.28358384174649,12.4301082124045,4.28358044168222,4.28358127985762,6.59053741902265,-7.03556153245474e-07,8.73403280797437e-11,6.59053390471481,-7.03459827214914e-07,8.73298214129605e-11,53.4182474530660,-9.40941782299458e-07,1.14488439043016e-10,6.59054622850375,-7.03624775237149e-07,8.73477502756049e-11,6.59056220109865,-7.03447335303449e-07,8.73279442007297e-11,-2.13173627267358e-07,-2.13130331833511e-07,-1.69164892586186e-06,-2.13115980989961e-07,-2.13109978200203e-07,3.24923971321948,3.24923963123043,6.22932149902217,3.24924176913499,3.24924207175928,2.11586194182217,2.11586200103673,2.27672157718318,2.11586452143432,2.11586456418959,1.38937636750779,1.38937585652781,1.45260654785536,1.38937640239296,1.38937619875856;4.28358463272980,4.28358395030953,4.28358174417978,12.4301047821507,4.28357977762629,6.59052428413427,-7.03502627592678e-07,8.73337436917212e-11,6.59052472830607,-7.03528723789183e-07,8.73373409485650e-11,6.59054679481724,-7.03422054563815e-07,8.73246589581908e-11,53.4182536143820,-9.40772520090640e-07,1.14468176374215e-10,6.59057148299082,-7.03438070966278e-07,8.73259306293749e-11,-2.13174281515596e-07,-2.13132709320772e-07,-2.13101818553307e-07,-1.69167219509279e-06,-2.13101301849776e-07,3.24923935212332,3.24923952543169,3.24924114704948,6.22932178161776,3.24924180941985,2.11586059680169,2.11586100188282,2.11586323114653,2.27672181198620,2.11586456075833,1.38937455658109,1.38937420870243,1.38937447312739,1.45260657153685,1.38937577603300;4.28358476851053,4.28358395156808,4.28358143308733,4.28357858105986,12.4301076836062,6.59053026991772,-7.03658871873140e-07,8.73528966147467e-11,6.59053273040002,-7.03620861414287e-07,8.73486304001904e-11,6.59055663989912,-7.03548915956141e-07,8.73401581913625e-11,6.59056589410910,-7.03620293609980e-07,8.73473121363620e-11,53.4182429950444,-9.40914604645789e-07,1.14486292758409e-10,-2.13173589247532e-07,-2.13138032055886e-07,-2.13096762677869e-07,-2.13109485564440e-07,-1.69166668540061e-06,3.24923966224555,3.24924000830761,3.24924173103520,3.24924207636520,6.22932118878736,2.11586121739621,2.11586166996360,2.11586403549400,2.11586514244880,2.27672088041047,1.38937524609323,1.38937503077643,1.38937514480813,1.38937660915089,1.45260646004639];

%% [e(k) Δx(k) ] 20*20
A_R0_1_3 = [    I_N  kron(-H,C*A1);
            zeros(15,5)  kron(I_N,A1) ];
A_R0_1_d_3 = [ zeros(5)  kron(-H,C*A1_d);
             zeros(15,5)  kron(I_N,A1_d) ];

A_R0_2_3 = [    I_N  kron(-H,C*A2);
            zeros(15,5)  kron(I_N,A2) ];
        
A_R0_2_d_3 = [ zeros(5)  kron(-H,C*A2_d);
             zeros(15,5)  kron(I_N,A2_d) ];

B_R0_3 = [    kron(-H,C*B);
            kron(I_N,B) ];
G_hat_3 = eye(20) - B_R0*inv(B_R0'*B_R0)*B_R0';

F_3 = [ kron(-H,C);
      eye(15) ];
D_3 = [ eye(5)  zeros(5,15) ];


Y = zeros(15,1);
  
X= [ kron(-H,C)*Y;
     Y ]; % [e x1 x2 x3 ]
 
%%  R(k)
for k=1:1:500
    if k<25
        R(:,k) =zeros(5,1);
    elseif k>=25 & k<=50
        R(:,k) = ( (k-25)/10 )*100*ones(5,1);
    else
        R(:,k) = 2.5*100*ones(5,1);
    end
end

%% 滑模面、控制律、输出
for k=1:1:201
    time(k)=(k-1)/20;
%%
      %预见长度(1,1)
      Y1(:,k)=[Y(1,k) Y(4,k) Y(7,k) Y(10,k) Y(13,k)];
      Y2(:,k)=[Y(2,k) Y(5,k) Y(8,k) Y(11,k) Y(14,k)];
      Y3(:,k)=[Y(3,k) Y(6,k) Y(9,k) Y(12,k) Y(15,k)];
      Tr1(:,k)=kron(I_N,C)*Y(:,k);
      
      if k>1
          % 误差 e
          x11(:,k)= -kron(H,I_m) * (kron(I_N,C)*Y(:,k) - R(:,k));
          % Δx 
          x21(:,k)=Y1(:,k)-Y1(:,k-1);
          x31(:,k)=Y2(:,k)-Y2(:,k-1);
          x41(:,k)=Y3(:,k)-Y3(:,k-1);
      else
          x11(:,k)= -kron(H,I_m) * (kron(I_N,C)*Y(:,k) - R(:,k));
          
          x21(:,k)=Y1(:,k);
          x31(:,k)=Y2(:,k);
          x41(:,k)=Y3(:,k);
      end
      
      %h1(:,k) =  ( 0.02*(   sin( (Y(2,k)+Y(4,k)+Y(6,k)+Y(8,k)+Y(10,k))/5  )   )^2  ) / (exp(2)+exp(5));
      h1(:,k) =  ( 1-1/(2+cos(-0.05*(k-0.5*pi) ) )  ) * ( 1-1/(2+cos(-0.05*(k+0.5*pi) ) )  );
      h2(:,k) = 1-h1(:,k);
      
      if k>1
         
          S3(:,k) = B_R0_3'*[x11(:,k);x21(:,k);x31(:,k);x41(:,k)] - B_R0_3'* ...
                  (   h1(:,k-1)*h1(:,k-1)*( A_R0_1_3 + B_R0_3*K1_20 )*[x11(:,k-1);x21(:,k-1);x31(:,k-1);x41(:,k-1)] + ...
                      h1(:,k-1)*h2(:,k-1)*( A_R0_1_3 + B_R0_3*K2_20 )*[x11(:,k-1);x21(:,k-1);x31(:,k-1);x41(:,k-1)] + ...
                      h2(:,k-1)*h1(:,k-1)*( A_R0_2_3 + B_R0_3*K1_20 )*[x11(:,k-1);x21(:,k-1);x31(:,k-1);x41(:,k-1)] + ...
                      h2(:,k-1)*h2(:,k-1)*( A_R0_2_3 + B_R0_3*K2_20 )*[x11(:,k-1);x21(:,k-1);x31(:,k-1);x41(:,k-1)]...
                  );
      else
          S3(:,k) = B_R0_3'*[x11(:,k);x21(:,k);x31(:,k);x41(:,k)];
      end
      
      e1=0;
      e2=0;
      
      % Δu(k) 的线性项
      for j=1:1:k
          e1=e1+x11(:,j);
      end

      % Δu(k) 的非线性项
      for j=1:1:k
         for jj=1:1:5
               sat(jj,j)=sign(S3(jj,j));
         end
          %W(:,j) = kron(ones(5,1),(  [ exp(-5*j)*cos(j) ;0]  )   );
          
          if j>1
             if j>4
                W(:,j) = ( h1(:,j)-h1(:,j-1)+h2(:,j)-h2(:,j-1) ) * (  kron(I_N,A1)*[Y1(:,j-1);Y2(:,j-1);Y3(:,j-1)] + kron(I_N,A1_d)*[Y1(:,j-3-1);Y2(:,j-3-1);Y3(:,j-3-1)] ) + ...
                          1*[ 0;0; -0.01*cos(Y(3,j))+0.04*Y(3,j);  0;0; -0.01*cos(Y(6,j))+0.04*Y(6,j);  0;0; -0.01*cos(Y(9,j))+0.04*Y(9,j); 0;0; -0.01*cos(Y(12,j))+0.04*Y(12,j); 0;0; -0.01*cos(Y(15,j))+0.04*Y(15,j) ]; 
                
                %e2=e2 + inv(B_R0_3'*B_R0_3)* (  (0.54-1)*S3(:,j) + 0.5*sat + norm(inv(B_R0_3'*B_R0_3)*B_R0_3'*F_3) * norm( W(:,j)) ) +  0.01*norm( [x21(:,j) zeros(5,1) zeros(5,1) ] );
                
                gbgf=inv(B_R0_3'*B_R0_3)*B_R0_3'*F_3;
                e2=e2 + inv(B_R0_3'*B_R0_3)* (  (0.54-1)*S3(:,j) + 0.5*sat(:,j) + ...
                   [ norm(gbgf(:,1:5)) ,norm(gbgf(:,6:10)) ,norm(gbgf(:,11:15)) ]* [ norm( W(1:5,j)); norm( W(5:10,j)); norm( W(11:15,j)) ] )...
                   +  0.01*norm( [x21(:,j) zeros(5,1) zeros(5,1) ] ) ;
             else
                 W(:,j) = ( h1(:,j)-h1(:,j-1)+h2(:,j)-h2(:,j-1) ) * (  kron(I_N,A1)*[Y1(:,j-1);Y2(:,j-1);Y3(:,j-1)] ) + ...
                          1*[ 0;0; -0.01*cos(Y(3,j))+0.04*Y(3,j);  0;0; -0.01*cos(Y(6,j))+0.04*Y(6,j);  0;0; -0.01*cos(Y(9,j))+0.04*Y(9,j); 0;0; -0.01*cos(Y(12,j))+0.04*Y(12,j); 0;0; -0.01*cos(Y(15,j))+0.04*Y(15,j) ] ;
                
                %e2=e2 + inv(B_R0_3'*B_R0_3)* (  (0.54-1)*S3(:,j) + 0.5*sat + norm(inv(B_R0_3'*B_R0_3)*B_R0_3'*F_3) * norm( W(:,j)) ) +  0.01*norm( [x21(:,j) zeros(5,1) zeros(5,1) ] );
                
                gbgf=inv(B_R0_3'*B_R0_3)*B_R0_3'*F_3;
                e2=e2 + inv(B_R0_3'*B_R0_3)* (  (0.54-1)*S3(:,j) + 0.5*sat(:,j) + ...
                   [ norm(gbgf(:,1:5)) ,norm(gbgf(:,6:10)) ,norm(gbgf(:,11:15)) ]* [ norm( W(1:5,j)); norm( W(5:10,j)); norm( W(11:15,j)) ] )...
                   +  0.01*norm( [x21(:,j) zeros(5,1) zeros(5,1) ] ) ;
             end
             
          else
               W(:,j) = 1*[ 0;0; -0.01*cos(Y(3,j))+0.04*Y(3,j);  0;0; -0.01*cos(Y(6,j))+0.04*Y(6,j);  0;0; -0.01*cos(Y(9,j))+0.04*Y(9,j); 0;0; -0.01*cos(Y(12,j))+0.04*Y(12,j); 0;0; -0.01*cos(Y(15,j))+0.04*Y(15,j) ];
               
              %e2=e2 + inv(B_R0_3'*B_R0_3)* (  (0.54-1)*S3(:,j) + 0.5*sat + norm(inv(B_R0_3'*B_R0_3)*B_R0_3'*F_3) * norm( W(:,j)) ) + 0.01*norm([x21(:,j) zeros(5,1) zeros(5,1) ]);
              
              gbgf=inv(B_R0_3'*B_R0_3)*B_R0_3'*F_3;
              e2=e2 + inv(B_R0_3'*B_R0_3)* (  (0.54-1)*S3(:,j) + 0.5*sat(:,j) + ...
                 [ norm(gbgf(:,1:5)) ,norm(gbgf(:,6:10)) ,norm(gbgf(:,11:15)) ]* [ norm( W(1:5,j)); norm( W(5:10,j)); norm( W(11:15,j)) ] )...
                 +  0.01*norm( [x21(:,j) zeros(5,1) zeros(5,1) ] ) ;
          end
      end
      
      if k>3
         U3(:,k)=( h1(:,k)*K1_20+h2(:,k)*K2_20 ) * [e1;Y1(:,k);Y2(:,k);Y3(:,k)] - ...
                 ( h1(:,k)*inv(B_R0_3'*B_R0_3)*B_R0_3'*A_R0_1_d +  h2(:,k)*inv(B_R0_3'*B_R0_3)*B_R0_3'*A_R0_2_d ) *...
                 [ e1;Y1(:,k-3);Y2(:,k-3);Y3(:,k-3) ] -e2 ;
             
         Y(:,k+1) =  h1(:,k) * (  kron(I_N,A1)*Y(:,k) + kron(I_N,A1_d)*Y(:,k-3) + kron(I_N,B)*( U3(:,k)+0.01*[sin(Y(1,k));sin(Y(4,k)); sin(Y(7,k)); sin(Y(10,k)); sin(Y(13,k)) ] ) ...
                      + 1*[ 0;0; -0.01*cos(Y(3,k))+0.04*Y(3,k);  0;0; -0.01*cos(Y(6,k))+0.04*Y(6,k);  0;0; -0.01*cos(Y(9,k))+0.04*Y(9,k); 0;0; -0.01*cos(Y(12,k))+0.04*Y(12,k); 0;0; -0.01*cos(Y(15,k))+0.04*Y(15,k) ]   ) +...
                      h2(:,k) * (  kron(I_N,A2)*Y(:,k) + kron(I_N,A2_d)*Y(:,k-3) + kron(I_N,B)*( U3(:,k)+0.01*[sin(Y(1,k)); sin(Y(4,k));sin(Y(7,k)); sin(Y(10,k)); sin(Y(13,k)) ] )...
                      + 1*[ 0;0; -0.01*cos(Y(3,k))+0.04*Y(3,k);  0;0; -0.01*cos(Y(6,k))+0.04*Y(6,k);  0;0; -0.01*cos(Y(9,k))+0.04*Y(9,k); 0;0; -0.01*cos(Y(12,k))+0.04*Y(12,k); 0;0; -0.01*cos(Y(15,k))+0.04*Y(15,k) ] ) ;
      else
          U3(:,k)=( h1(:,k)*K1_20+h2(:,k)*K2_20 ) * [e1;Y1(:,k);Y2(:,k);Y3(:,k)] - e2 ;
          
          Y(:,k+1) =  h1(:,k) * (  kron(I_N,A1)*Y(:,k) + kron(I_N,B)*( U3(:,k)+0.01*[sin(Y(1,k));sin(Y(4,k)); sin(Y(7,k)); sin(Y(10,k)); sin(Y(13,k)) ] )...
                      + 1*[ 0;0; -0.01*cos(Y(3,k))+0.04*Y(3,k);  0;0; -0.01*cos(Y(6,k))+0.04*Y(6,k);  0;0; -0.01*cos(Y(9,k))+0.04*Y(9,k); 0;0; -0.01*cos(Y(12,k))+0.04*Y(12,k); 0;0; -0.01*cos(Y(15,k))+0.04*Y(15,k) ] ) + ...
                      h2(:,k) * (  kron(I_N,A2)*Y(:,k) + kron(I_N,B)*( U3(:,k)+0.01*[sin(Y(1,k)); sin(Y(4,k)); sin(Y(7,k));sin(Y(10,k)); sin(Y(13,k)) ] )...
                      + 1*[ 0;0; -0.01*cos(Y(3,k))+0.04*Y(3,k);  0;0; -0.01*cos(Y(6,k))+0.04*Y(6,k);  0;0; -0.01*cos(Y(9,k))+0.04*Y(9,k); 0;0; -0.01*cos(Y(12,k))+0.04*Y(12,k); 0;0; -0.01*cos(Y(15,k))+0.04*Y(15,k) ] ) ;
      end
     
end

%% R1(k)
for k=1:1:201
    if k<25
       R1(:,k) =zeros(5,1);
    elseif k>=25 & k<=50
        R1(:,k) = ( (k-25)/10 )*100*ones(5,1);
    else
       R1(:,k) = 2.5*100*ones(5,1);
    end
end

%% 
figure(1);
plot(time,R1(1,:)/100,'k-',time,Tr1(1,:)/100,time,Tr1(2,:)/100,time,Tr1(3,:)/100,time,Tr1(4,:)/100,time,Tr1(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('output response','FontName','Times New Roman','FontSize',12);
f1=legend('$r_{k}$','agent1','agent2','agent3','agent4','agent5','Location','NorthWest','FontName','Times New Roman','FontSize',10);
set(f1,'Interpreter','latex','fontsize',10);

axes('Position',[0.4,0.5,0.5,0.2]); % 生成子图   
plot(time,R1(1,:)/100,'k-',time,Tr1(1,:)/100,time,Tr1(2,:)/100,time,Tr1(3,:)/100,time,Tr1(4,:)/100,time,Tr1(5,:)/100,'LineWidth',1);
ylim([1.7411,1.7412]); % 放大的横坐标
xlim([2.05,2.43]); % 放大的横坐标

figure(2);
plot(time,S3(1,:)/100,time,S3(2,:)/100,time,S3(3,:)/100,time,S3(4,:)/100,time,S3(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('sliding mode surface','FontName','Times New Roman','FontSize',12);
legend('agent1','agent2','agent3','agent4','agent5','FontName','Times New Roman','FontSize',10);

figure(3);
plot(time,U3(1,:)/100,time,U3(2,:)/100,time,U3(3,:)/100,time,U3(4,:)/100,time,U3(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('control input','FontName','Times New Roman','FontSize',12);
legend('agent1','agent2','agent3','agent4','agent5','FontName','Times New Roman','FontSize',10);

axes('Position',[0.4,0.4,0.5,0.2]); % 生成子图   
plot(time,U3(1,:)/100,time,U3(2,:)/100,time,U3(3,:)/100,time,U3(4,:)/100,time,U3(5,:)/100,'LineWidth',1);
ylim([-33.852,-33.851]); % 放大的横坐标
xlim([2.09,2.19]); % 放大的横坐标

figure(4);
plot(time,x11(1,:)/100,time,x11(2,:)/100,time,x11(3,:)/100,time,x11(4,:)/100,time,x11(5,:)/100,'LineWidth',1);
xlabel('time(s)','FontName','Times New Roman','FontSize',12);
ylabel('tracking error','FontName','Times New Roman','FontSize',12);
legend('agent1','agent2','agent3','agent4','agent5','FontName','Times New Roman','FontSize',10);