function [NS,QN,EL]=QS_Data(RAD,UNIV)

%****************
%Assign the input
%****************
ATOM=RAD.ATOM;

hbar=UNIV.hbar;
q=UNIV.q;
c=UNIV.c;        

%//////////////////////////////////////////////////////////////////////////
%**************************************************************************
%***************************NIST ENERGY LEVELS*****************************
%**************************************************************************
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Hydrogen Energy Levels - No Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'H')==1 && RAD.FS==0)
    LEVELS(1,1)=0;

    LEVELS(2,1:3)=82259.158;

    LEVELS(3,1:5)=97492.304;

    LEVELS(4,1:7)=102823.904;

    LEVELS(5,1:9)=105291.657;

    LEVELS(6,1:11)=106632.1681;

    LEVELS(7,1:13)=107440.4508;

    LEVELS(8,1:15)=107965.0568;

    LEVELS(9,1:17)=108324.7253;

    LEVELS(10,1:19)=108581.9945;

    LEVELS(11,21)=108772.3445;

    LEVELS(12,23)=108917.1209;
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Hydrogen Energy Levels - No Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Deutrium Energy Levels - No Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'D')==1 && RAD.FS==0)
    LEVELS(1,1)=0;

    LEVELS(2,1:3)=82281.493;

    LEVELS(3,1:5)=97518.836;

    LEVELS(4,1:7)=102851.878;

    LEVELS(5,1:9)=105320.308;

    LEVELS(6,1:11)=106661.1812;

    LEVELS(7,1:13)=107469.6848;

    LEVELS(8,1:15)=107994.4344;

    LEVELS(9,1:17)=108354.2009;

    LEVELS(10,1:19)=108611.5396;

    LEVELS(11,1:21)=108801.9411;

    LEVELS(12,1:23)=108946.7571;
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Deutrium Energy Levels - No Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Hydrogen Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'H')==1 && RAD.FS==1)
    LEVELS(1,1)=0;

    LEVELS(2,1)=82258.9543992821;
    LEVELS(2,2)=82258.919113;
    LEVELS(2,3)=82259.2850014;

    LEVELS(3,1)=97492.221701;
    LEVELS(3,2)=97492.2112;
    LEVELS(3,3)=97492.319611;
    LEVELS(3,4)=97492.319433;
    LEVELS(3,5)=97492.355566;

    LEVELS(4,1)=102823.8530211;
    LEVELS(4,2)=102823.8485825;
    LEVELS(4,3)=102823.8943175;
    LEVELS(4,4)=102823.89425;
    LEVELS(4,5)=102823.9094871;
    LEVELS(4,6)=102823.90949;
    LEVELS(4,7)=102823.917091;

    LEVELS(5,1)=105291.63094;
    LEVELS(5,2)=105291.62867;
    LEVELS(5,3)=105291.65209;
    LEVELS(5,4)=105291.651993;
    LEVELS(5,5)=105291.659796;
    LEVELS(5,6)=105291.65983494;
    LEVELS(5,7)=105291.6637;
    LEVELS(5,8)=105291.66373033;
    LEVELS(5,9)=105291.666072;

    LEVELS(6,1)=106632.1498416;
    LEVELS(6,2)=106632.1485242;
    LEVELS(6,3)=106632.1620756;
    LEVELS(6,4)=106632.1620536;
    LEVELS(6,5)=106632.1665697;
    LEVELS(6,6)=106632.16656761;
    LEVELS(6,7)=106632.16882614;
    LEVELS(6,8)=106632.16882188;
    LEVELS(6,9)=106632.17017701;
    LEVELS(6,10)=106632.17017434;
    LEVELS(6,11)=106632.17107779;

    LEVELS(7,1)=107440.4393319;
    LEVELS(7,2)=107440.4385022;
    LEVELS(7,3)=107440.4470360;
    LEVELS(7,4)=107440.4470215;
    LEVELS(7,5)=107440.4498661;
    LEVELS(7,6)=107440.4498609;
    LEVELS(7,7)=107440.4512832;
    LEVELS(7,8)=107440.4512805;
    LEVELS(7,9)=107440.4521339;
    LEVELS(7,10)=107440.4521322;
    LEVELS(7,11)=107440.4527011;
    LEVELS(7,12)=107440.4527000;
    LEVELS(7,13)=107440.4531064;

    LEVELS(8,1)=107965.0497146;
    LEVELS(8,2)=107965.0491587;
    LEVELS(8,3)=107965.0548757;
    LEVELS(8,4)=107965.0548659;
    LEVELS(8,5)=107965.0567716;
    LEVELS(8,6)=107965.0567681;
    LEVELS(8,7)=107965.0577209;
    LEVELS(8,8)=107965.0577191;
    LEVELS(8,9)=107965.0582908;
    LEVELS(8,10)=107965.0582897;
    LEVELS(8,11)=107965.0586708;
    LEVELS(8,12)=107965.0586701;
    LEVELS(8,13)=107965.0589423;
    LEVELS(8,14)=107965.0589417;
    LEVELS(8,15)=107965.0591459; 
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Hydrogen Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Deutrium Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'D')==1 && RAD.FS==1)
    LEVELS(1,1)=0;

    LEVELS(2,1)=82281.336361094;
    LEVELS(2,2)=82281.3010273;
    LEVELS(2,3)=82281.667009;

    LEVELS(3,1)=97518.748607;
    LEVELS(3,2)=97518.738088;
    LEVELS(3,3)=97518.846525;
    LEVELS(3,4)=97518.84627;
    LEVELS(3,5)=97518.88241;

    LEVELS(4,1)=102851.8306186;
    LEVELS(4,2)=102851.826163;
    LEVELS(4,3)=102851.871920;
    LEVELS(4,4)=102851.871846;
    LEVELS(4,5)=102851.8870958;
    LEVELS(4,6)=102851.88706860;
    LEVELS(4,7)=102851.89469315;

    LEVELS(5,1)=105320.27991;
    LEVELS(5,2)=105320.27762;
    LEVELS(5,3)=105320.30105;
    LEVELS(5,4)=105320.30125;
    LEVELS(5,5)=105320.30906;
    LEVELS(5,6)=105320.30891596;
    LEVELS(5,7)=105320.31281973;
    LEVELS(5,8)=105320.31281242;
    LEVELS(5,9)=105320.31515469;

    LEVELS(6,1)=106661.16360;
    LEVELS(6,2)=106661.16228;
    LEVELS(6,3)=106661.17584;
    LEVELS(6,4)=106661.17594;
    LEVELS(6,5)=106661.18046;
    LEVELS(6,6)=106661.18039908;
    LEVELS(6,7)=106661.18265821;
    LEVELS(6,8)=106661.18265396;
    LEVELS(6,9)=106661.18400944;
    LEVELS(6,10)=106661.18400680;
    LEVELS(6,11)=106661.18491047;
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Deutrium Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Tritium Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'T')==1 && RAD.FS==1)
    LEVELS(1,1)=0;

    LEVELS(2,1)=82288.78325;
    LEVELS(2,2)=82288.74801;
    LEVELS(2,3)=82289.11417;

    LEVELS(3,1)=97527.584;
    LEVELS(3,2)=97527.56400;
    LEVELS(3,3)=97527.67248;
    LEVELS(3,4)=97527.67231;
    LEVELS(3,5)=97527.70863;
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Tritium Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Helium II Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'HeII')==1 && RAD.FS==1)
    LEVELS(1,1)=0;

    LEVELS(2,1)=329179.76197;
    LEVELS(2,2)=329179.293589;
    LEVELS(2,3)=329185.150747;

    LEVELS(3,1)=390140.964175;
    LEVELS(3,2)=390140.824628;
    LEVELS(3,3)=390142.5600890;
    LEVELS(3,4)=390142.55723373;
    LEVELS(3,5)=390143.13566549;

    LEVELS(4,1)=411477.181412;
    LEVELS(4,2)=411477.1224141;
    LEVELS(4,3)=411477.8545580;
    LEVELS(4,4)=411477.853339;
    LEVELS(4,5)=411478.097366;
    LEVELS(4,6)=411478.096926;
    LEVELS(4,7)=411478.218936;

    LEVELS(5,1)=421352.708831;
    LEVELS(5,2)=421352.6785905;
    LEVELS(5,3)=421353.0534464;
    LEVELS(5,4)=421353.052818;
    LEVELS(5,5)=421353.177760;
    LEVELS(5,6)=421353.177533;
    LEVELS(5,7)=421353.240002;
    LEVELS(5,8)=421353.2398840;
    LEVELS(5,9)=421353.2773653;

    LEVELS(6,1)=426717.153325;
    LEVELS(6,2)=426717.135815;
    LEVELS(6,3)=426717.352745;
    LEVELS(6,4)=426717.352377;
    LEVELS(6,5)=426717.424682;
    LEVELS(6,6)=426717.424550;
    LEVELS(6,7)=426717.460701;
    LEVELS(6,8)=426717.4606322;
    LEVELS(6,9)=426717.4823229;
    LEVELS(6,10)=426717.4822804;
    LEVELS(6,11)=426717.4967408;

    LEVELS(7,1)=429951.7257273;
    LEVELS(7,2)=429951.714697;
    LEVELS(7,3)=429951.851305;
    LEVELS(7,4)=429951.851072;
    LEVELS(7,5)=429951.896605;
    LEVELS(7,6)=429951.896522;
    LEVELS(7,7)=429951.919288;
    LEVELS(7,8)=429951.9192438;
    LEVELS(7,9)=429951.9329032;
    LEVELS(7,10)=429951.9328764;
    LEVELS(7,11)=429951.9419827; 
    LEVELS(7,12)=429951.94196436;
    LEVELS(7,13)=429951.94846880;

    LEVELS(8,1)=432051.0782102;
    LEVELS(8,2)=432051.0708182;
    LEVELS(8,3)=432051.1623352;
    LEVELS(8,4)=432051.162178;
    LEVELS(8,5)=432051.192681;
    LEVELS(8,6)=432051.192626;
    LEVELS(8,7)=432051.207877;
    LEVELS(8,8)=432051.2078481;
    LEVELS(8,9)=432051.2169988;
    LEVELS(8,10)=432051.2169808;
    LEVELS(8,11)=432051.2230813; 
    LEVELS(8,12)=432051.22306901;
    LEVELS(8,13)=432051.22742651; 
    LEVELS(8,14)=432051.22741759;
    LEVELS(8,15)=432051.23068568;  

    LEVELS(9,1)=433490.3821779;
    LEVELS(9,2)=433490.3769859;
    LEVELS(9,3)=433490.4412599;
    LEVELS(9,4)=433490.4411499;
    LEVELS(9,5)=433490.4625729;
    LEVELS(9,6)=433490.4625339;
    LEVELS(9,7)=433490.4732459;
    LEVELS(9,8)=433490.4732251;
    LEVELS(9,9)=433490.4796520;
    LEVELS(9,10)=433490.4796393;
    LEVELS(9,11)=433490.4839239; 
    LEVELS(9,12)=433490.48391531;
    LEVELS(9,13)=433490.48697561; 
    LEVELS(9,14)=433490.48696940;
    LEVELS(9,15)=433490.48926469;
    LEVELS(9,16)=433490.489259934;
    LEVELS(9,17)=433490.491045154;  

    LEVELS(10,1)=434519.9048042;
    LEVELS(10,2)=434519.9010184;
    LEVELS(10,3)=434519.9478747;
    LEVELS(10,4)=434519.9477941;
    LEVELS(10,5)=434519.9634121;
    LEVELS(10,6)=434519.9633831;
    LEVELS(10,7)=434519.9711918;
    LEVELS(10,8)=434519.9711768;
    LEVELS(10,9)=434519.9758620;
    LEVELS(10,10)=434519.9758527;
    LEVELS(10,11)=434519.9789762; 
    LEVELS(10,12)=434519.97896989;
    LEVELS(10,13)=434519.98120092; 
    LEVELS(10,14)=434519.98119635;
    LEVELS(10,15)=434519.98286961;
    LEVELS(10,16)=434519.98286614;
    LEVELS(10,17)=434519.98416757; 
    LEVELS(10,18)=434519.984164842;
    LEVELS(10,19)=434519.985205983;

    LEVELS(11,1)=435281.6332465;
    LEVELS(11,2)=435281.6304019;
    LEVELS(11,3)=435281.6656057;
    LEVELS(11,4)=435281.6655450;
    LEVELS(11,5)=435281.6772787;
    LEVELS(11,6)=435281.6772572;
    LEVELS(11,7)=435281.6831240;
    LEVELS(11,8)=435281.6831127;
    LEVELS(11,9)=435281.6866328;
    LEVELS(11,10)=435281.6866258;
    LEVELS(11,11)=435281.6889725; 
    LEVELS(11,12)=435281.68896776;
    LEVELS(11,13)=435281.69064396; 
    LEVELS(11,14)=435281.69064053;
    LEVELS(11,15)=435281.69189768;
    LEVELS(11,16)=435281.69189507;
    LEVELS(11,17)=435281.69287285; 
    LEVELS(11,18)=435281.692870804;
    LEVELS(11,19)=435281.693653024;
    LEVELS(11,20)=435281.693651372;
    LEVELS(11,21)=435281.694291373; 
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Helium II Energy Levels - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Helium I Energy Levels - Singlet States (S=0) - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'HeI')==1 && RAD.SPIN==0)
    LEVELS(1,1)=0;

    LEVELS(2,1)=166277.440141;
    LEVELS(2,2)=171134.896946;

    LEVELS(3,1)=184864.82932;
    LEVELS(3,2)=186209.364940;
    LEVELS(3,3)=186104.9666893;

    LEVELS(4,1)=190940.226355;
    LEVELS(4,2)=191492.711909;
    LEVELS(4,3)=191446.4557405;
    LEVELS(4,4)=191451.89746084;

    LEVELS(5,1)=193663.512095;
    LEVELS(5,2)=193942.462294;
    LEVELS(5,3)=193918.28990114;
    LEVELS(5,4)=193921.13088112;
    LEVELS(5,5)=193921.621933101;

    LEVELS(6,1)=195114.868700;
    LEVELS(6,2)=195274.9084660;
    LEVELS(6,3)=195260.77050808;
    LEVELS(6,4)=195262.431784358;
    LEVELS(6,5)=195262.727124343;
    LEVELS(6,6)=195262.795750108;

    LEVELS(7,1)=195978.895183;
    LEVELS(7,2)=196079.0875698;
    LEVELS(7,3)=196070.12838342;
    LEVELS(7,4)=196071.181064965;
    LEVELS(7,5)=196071.371283543;
    LEVELS(7,6)=196071.416212368;
    LEVELS(7,7)=196071.430106489;

    LEVELS(8,1)=196534.5641527;
    LEVELS(8,2)=196601.4002470;
    LEVELS(8,3)=196595.37405237;
    LEVELS(8,4)=196596.082092602;
    LEVELS(8,5)=196596.211361727;
    LEVELS(8,6)=196596.242205021;
    LEVELS(8,7)=196596.251866245;
    LEVELS(8,8)=196596.255509048;
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Helium I Energy Levels - Singlet States (S=0) - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Helium I Energy Levels - Triplet States (S=1) - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'HeI')==1 && RAD.SPIN==1)
    LEVELS(2,1)=159855.9743297;
    LEVELS(2,2)=169087.8308131;
    LEVELS(2,3)=169086.8428979;
    LEVELS(2,4)=169086.7664725;

    LEVELS(3,1)=183236.79170;
    LEVELS(3,2)=185564.854540;
    LEVELS(3,3)=185564.583895;
    LEVELS(3,4)=185564.561920;
    LEVELS(3,5)=186101.5928903;
    LEVELS(3,6)=186101.5486891;
    LEVELS(3,7)=186101.5461767;

    LEVELS(4,1)=190298.113260;
    LEVELS(4,2)=191217.160290;
    LEVELS(4,3)=191217.049963;
    LEVELS(4,4)=191217.040967;
    LEVELS(4,5)=191444.5006512;
    LEVELS(4,6)=191444.4821307;
    LEVELS(4,7)=191444.4809292;
    LEVELS(4,8)=191451.88970691;
    LEVELS(4,9)=191451.87394878;
    LEVELS(4,10)=191451.88108855;

    LEVELS(5,1)=193346.991344;
    LEVELS(5,2)=193800.767563;
    LEVELS(5,3)=193800.712118;
    LEVELS(5,4)=193800.707595;
    LEVELS(5,5)=193917.16138710;
    LEVELS(5,6)=193917.15192855;
    LEVELS(5,7)=193917.15128741;
    LEVELS(5,8)=193921.12575310;
    LEVELS(5,9)=193921.11826439;
    LEVELS(5,10)=193921.12134275;
    LEVELS(5,11)=193921.620238142;
    LEVELS(5,12)=193921.614948739;
    LEVELS(5,13)=193921.617719262;

    LEVELS(6,1)=194936.119697;
    LEVELS(6,2)=195192.7772813;
    LEVELS(6,3)=195192.7455590;
    LEVELS(6,4)=195192.7429742;
    LEVELS(6,5)=195260.07720298;
    LEVELS(6,6)=195260.07173649;
    LEVELS(6,7)=195260.07135803;
    LEVELS(6,8)=195262.428373005;
    LEVELS(6,9)=195262.424216836;
    LEVELS(6,10)=195262.425821479;
    LEVELS(6,11)=195262.726141773;
    LEVELS(6,12)=195262.723082377;
    LEVELS(6,13)=195262.724684231;
    LEVELS(6,14)=195262.795042911;
    LEVELS(6,15)=195262.793050975;
    LEVELS(6,16)=195262.794095149;

    LEVELS(7,1)=195868.236975;
    LEVELS(7,2)=196027.3364649;
    LEVELS(7,3)=196027.3166408;
    LEVELS(7,4)=196027.3150267;
    LEVELS(7,5)=196069.67648966;
    LEVELS(7,6)=196069.67304998;
    LEVELS(7,7)=196069.67280898;
    LEVELS(7,8)=196071.178729980;
    LEVELS(7,9)=196071.176177974;
    LEVELS(7,10)=196071.177123492;
    LEVELS(7,11)=196071.370663981;
    LEVELS(7,12)=196071.368738129;
    LEVELS(7,13)=196071.369746173;
    LEVELS(7,14)=196071.415767010;
    LEVELS(7,15)=196071.414512625;
    LEVELS(7,16)=196071.415170183;
    LEVELS(7,17)=196071.429772608;
    LEVELS(7,18)=196071.428891148;
    LEVELS(7,19)=196071.429352355;

    LEVELS(8,1)=196461.3618105;
    LEVELS(8,2)=196566.7261046;
    LEVELS(8,3)=196566.7128971;
    LEVELS(8,4)=196566.7118220;
    LEVELS(8,5)=196595.06466814;
    LEVELS(8,6)=196595.06236505;
    LEVELS(8,7)=196595.06220247;
    LEVELS(8,8)=196596.080442094;
    LEVELS(8,9)=196596.078760123;
    LEVELS(8,10)=196596.079366013;
    LEVELS(8,11)=196596.210946275;
    LEVELS(8,12)=196596.209656481;
    LEVELS(8,13)=196596.210331441;
    LEVELS(8,14)=196596.241906662;
    LEVELS(8,15)=196596.241066328;
    LEVELS(8,16)=196596.241506841;
    LEVELS(8,17)=196596.251642570;
    LEVELS(8,18)=196596.251052062;
    LEVELS(8,19)=196596.251361035;
    LEVELS(8,20)=196596.255335763;
    LEVELS(8,21)=196596.254898015;
    LEVELS(8,22)=196596.255126492;
    
    %***************
    %Assign the QN's
    %***************
    n=RAD.PQN;
    s=RAD.SPIN;
    l=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Helium I Energy Levels - Triplet States (S=1) - Fine Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Argon II Energy - 4726.86 A - Fine Structure
%             
%                 CONFIGURATION                    TERM
%
%            LL -> 3s^2 3p^4 (3P) 4s                2P
%            UL -> 3s^2 3p^4 (3P) 4p                2Do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'ArII')==1 && strcmpi(RAD.WAVE,'4726A')==1)
    %********************
    %Lower level energies
    %********************
    LEVELS(1,1)=139258.3389;
    LEVELS(1,2)=138243.6444;

    %********************
    %Upper level energies
    %********************
    LEVELS(2,1)=159393.3850;
    LEVELS(2,2)=158730.2995;
    
    %***************
    %Assign the QN's
    %***************
    n=[4 4];
    s=0.5;
    l=[1 2];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Argon II Energy - 4726.86 A - Fine Structure
%             
%                 CONFIGURATION                    TERM
%
%            LL -> 3s^2 3p^4 (3P) 4s                2P
%            UL -> 3s^2 3p^4 (3P) 4p                2Do
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Argon II Energy - 4806.02 A - Fine Structure
%             
%                 CONFIGURATION                    TERM
%
%            LL -> 3s^2 3p^4 (3P) 4s                4P
%            UL -> 3s^2 3p^4 (3P) 4p                4Po
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmpi(ATOM,'ArII')==1 && strcmpi(RAD.WAVE,'4806A')==1)
    %********************
    %Lower level energies
    %********************
    LEVELS(1,1)=135601.7336;
    LEVELS(1,2)=135085.9958;
    LEVELS(1,3)=134241.7389;

    %********************
    %Upper level energies
    %********************
    LEVELS(2,1)=155708.1075;
    LEVELS(2,2)=155351.1206;
    LEVELS(2,3)=155043.1619;
    
    %***************
    %Assign the QN's
    %***************
    n=[4 4];
    s=1.5;
    l=[1 1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Argon II Energy - 4806.02 A - Fine Structure
%             
%                 CONFIGURATION                    TERM
%
%            LL -> 3s^2 3p^4 (3P) 4s                4P
%            UL -> 3s^2 3p^4 (3P) 4p                4Po
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%//////////////////////////////////////////////////////////////////////////
%**************************************************************************
%***************************NIST ENERGY LEVELS*****************************
%**************************************************************************
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%***************************************
%Conversion factor from wavenumber to eV
%***************************************
C=100*2*pi*hbar*c/q;

%**************************
%Calc. the number of states
%**************************
NS=Num_States(n,s,l);

QN=cell(2,max(NS));
EL(1:2,1:max(NS))=0;
if isempty(l)==0
    %*******************************************************
    %Calc. the QN's |n,s,l,j,mj> with BOTH n and l specified
    %*******************************************************
    for ii=1:2
        mm=0;
        nn=0;
        for kk=1:(abs(l(ii)+s)-abs(l(ii)-s)+1)
            mm=mm+1;
            
            %*******************************
            %Calc. total angular momentum QN
            %*******************************
            j=abs(l(ii)-s)+(kk-1);
            
            for ll=1:2*j+1
                nn=nn+1;
                
                %******************************************
                %Calc. total angular momentum projection QN
                %******************************************
                mj=-j+(ll-1);

                %***********
                %Assign QN's
                %***********
                QN{ii,nn}(1:5)=[n(ii),s,l(ii),j,mj];

                %********************
                %Assign energy levels
                %********************
                EL(ii,nn)=LEVELS(ii,mm);
            end
        end
    end 
elseif isempty(l)==1
    %*****************************************************
    %Calc. the QN's |n,s,l,j,mj> with ONLY n and specified
    %*****************************************************
    for ii=1:2
        mm=0;
        nn=0;
        for jj=1:n(ii)
            %*************************
            %Calc. angular momentum QN
            %*************************            
            l=jj-1;
            
            for kk=1:(abs(l+s)-abs(l-s)+1)
                mm=mm+1;
                
                %*******************************
                %Calc. total angular momentum QN
                %*******************************
                j=abs(l-s)+(kk-1);
                
                for ll=1:2*j+1
                    nn=nn+1;

                    %******************************************
                    %Calc. total angular momentum projection QN
                    %******************************************
                    mj=-j+(ll-1);
                    
                    %***********
                    %Assign QN's
                    %***********
                    QN{ii,nn}(1:5)=[n(ii),s,l,j,mj];
                    
                    %********************
                    %Assign energy levels
                    %********************
                    EL(ii,nn)=LEVELS(n(ii),mm);
                end
            end
        end
    end
end

%************************
%Convert wavenumber to eV
%************************
EL=C*EL;

end