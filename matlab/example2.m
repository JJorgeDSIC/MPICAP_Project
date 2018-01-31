N = 5;
%Basic problem

Atest = [
    0.7577,0.7060,0.8235,0.4387,0.4898
    0.7431,0.0318,0.6948,0.3816,0.4456
    0.3922,0.2769,0.3171,0.7655,0.6463
    0.6555,0.0462,0.9502,0.7952,0.7094
    0.1712,0.0971,0.0344,0.1869,0.7547]';

Btest = [
    0.8147    ,0.9058   ,0.1270   ,0.9134  ,0.6324
    0.0975   ,0.2785   ,0.5469  ,0.9575    ,0.9649
    0.1576   ,0.9706   ,0.9572  ,0.4854    ,0.8003
    0.1419   ,0.4218   ,0.9157  ,0.7922    ,0.9595
    0.6557   ,0.0357   ,0.8491  ,0.9340    ,0.6787]';

Ctest = [
    0.2760,0.4984,0.7513,0.9593,0.8407
    0.6797,0.9597,0.2551,0.5472,0.2543
    0.6551,0.3404,0.5060,0.1386,0.8143
    0.1626,0.5853,0.6991,0.1493,0.2435
    0.1190,0.2238,0.8909,0.2575,0.9293]';


[Xc, Tc]=solver_AX_BX_C(Atest,Btest,Ctest);

compute_norm(Atest,Btest,Ctest,Xc)