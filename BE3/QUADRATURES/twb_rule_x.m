function x = twb_rule_x ( strength )

%*****************************************************************************80
%
%% twb_rule_x returns the x abscissas of a TWB rule of given strength.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 April 2019
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Mark Taylor, Beth Wingate, Len Bos, 
%    Several new quadrature formulas for polynomial integration in the triangle, 
%    http://arxiv.org/abs/math/0501496v2,
%    08 February 2007.
%
%  Parameters:
%
%    Input, integer STRENGTH, the desired strength.
%    1 <= STRENGTH.
%
%    Output, real X(N), the x abscissas of the rule, if it exists.
%
  n = twb_rule_n ( strength );

  if ( n < 1 )
    x = [];
    return
  end

  if ( strength == 1 )
    x = [ ...
    0.3333333333333 ];
  elseif ( strength == 2 )
    x = [ ...
    0.1666666666667; ...
    0.6666666666667; ...
    0.1666666666667 ];
  elseif ( strength == 4 )
    x = [ ...
    0.0915762135098; ...
    0.8168475729805; ...
    0.0915762135098; ...
    0.1081030181681; ...
    0.4459484909160; ...
    0.4459484909160 ];
  elseif ( strength == 5 )
    x = [ ...
    0.0000000000000; ...
    1.0000000000000; ...
    0.0000000000000; ...
    0.2673273531185; ...
    0.6728175529461; ...
    0.0649236350054; ...
    0.6716498539042; ...
    0.0654032456800; ...
    0.2693767069140; ...
    0.3386738503896 ];
  elseif ( strength == 7 )
    x = [ ...
    1.0000000000000; ...
    0.0000000000000; ...
    0.0000000000000; ...
    0.7839656651012; ...
    0.1738960507345; ...
    0.1738960507345; ...
    0.0421382841642; ...
    0.7839656651012; ...
    0.0421382841642; ...
    0.4743880861752; ...
    0.4743880861752; ...
    0.0512238276497; ...
    0.2385615300181; ...
    0.5228769399639; ...
    0.2385615300181 ];
  elseif ( strength == 9 )
    x = [ ...
    0.0451890097844; ...
    0.0451890097844; ...
    0.9096219804312; ...
    0.7475124727339; ...
    0.2220631655373; ...
    0.7475124727339; ...
    0.2220631655373; ...
    0.0304243617288; ...
    0.0304243617288; ...
    0.1369912012649; ...
    0.6447187277637; ...
    0.1369912012649; ...
    0.2182900709714; ...
    0.2182900709714; ...
    0.6447187277637; ...
    0.0369603304334; ...
    0.4815198347833; ...
    0.4815198347833; ...
    0.4036039798179; ...
    0.4036039798179; ...
    0.1927920403641 ];
  elseif ( strength == 11 )
    x = [ ...
    0.0000000000000; ...
    0.9451704450173; ...
    0.9289002405719; ...
    0.0685505797224; ...
    0.0243268355615; ...
    0.1279662835335; ...
    0.0277838749488; ...
    0.0287083428360; ...
    0.7498347588656; ...
    0.7228007909707; ...
    0.2497602062386; ...
    0.0865562992839; ...
    0.8325513856998; ...
    0.3061619157672; ...
    0.0303526617491; ...
    0.4868610595047; ...
    0.6657904293017; ...
    0.1765456154221; ...
    0.0293121007360; ...
    0.5295657488667; ...
    0.1444673824391; ...
    0.3299740111411; ...
    0.5361815729052; ...
    0.5511507516862; ...
    0.1437790861923; ...
    0.3348066587327; ...
    0.1529619437161; ...
    0.3430183498147 ];
  elseif ( strength == 13 )
    x = [ ...
    0.0242935351590; ...
    0.0265193427722; ...
    0.9492126023551; ...
    0.0033775763749; ...
    0.4757672298101; ...
    0.5190783193471; ...
    0.8616839745321; ...
    0.1249209759926; ...
    0.0138565453861; ...
    0.0211887064222; ...
    0.8432296787219; ...
    0.1354231797865; ...
    0.3088853510679; ...
    0.6685057595169; ...
    0.0226545012557; ...
    0.2808515408772; ...
    0.6922446749051; ...
    0.0268617447119; ...
    0.1141778485470; ...
    0.7974807922061; ...
    0.0892807293894; ...
    0.1052487892455; ...
    0.6663022280740; ...
    0.2307803737547; ...
    0.1705059157540; ...
    0.5086593973043; ...
    0.3141823862281; ...
    0.4617460817864; ...
    0.0693087496081; ...
    0.4651955259268; ...
    0.2578625857893; ...
    0.6112627766779; ...
    0.1305182135934; ...
    0.4281437991828; ...
    0.3356995783730; ...
    0.2305424298836 ];
  elseif ( strength == 14 )
    x = [ ...
    0.0000000000000; ...
    1.0000000000000; ...
    0.0000000000000; ...
    0.0573330873026; ...
    0.0573330873026; ...
    0.9275286857160; ...
    0.0151382269814; ...
    0.9275286857160; ...
    0.0151382269814; ...
    0.8159625040711; ...
    0.8159625040711; ...
    0.1659719969565; ...
    0.0180654989724; ...
    0.1659719969565; ...
    0.0180654989724; ...
    0.3165475556378; ...
    0.6647637544849; ...
    0.0186886898773; ...
    0.0186886898773; ...
    0.3165475556378; ...
    0.6647637544849; ...
    0.0192662192492; ...
    0.4903668903754; ...
    0.4903668903754; ...
    0.0875134669581; ...
    0.0875134669581; ...
    0.8249730660837; ...
    0.0935526036219; ...
    0.0935526036219; ...
    0.2079865423167; ...
    0.6984608540613; ...
    0.6984608540613; ...
    0.2079865423167; ...
    0.0974892983467; ...
    0.3645018421383; ...
    0.5380088595149; ...
    0.5380088595149; ...
    0.3645018421383; ...
    0.0974892983467; ...
    0.2217145894873; ...
    0.5565708210253; ...
    0.2217145894873; ...
    0.3860471669296; ...
    0.2279056661408; ...
    0.3860471669296 ];
  elseif ( strength == 16 )
    x = [ ...
    1.0000000000000; ...
    0.0000000000000; ...
    0.0000000000000; ...
    0.9398863583577; ...
    0.0543806683058; ...
    0.0093940049164; ...
    0.0164345086362; ...
    0.9469487269862; ...
    0.0426604005768; ...
    0.0122269495439; ...
    0.8673696521047; ...
    0.8456744021389; ...
    0.1395759632103; ...
    0.1317821743231; ...
    0.0157955126300; ...
    0.7365462884436; ...
    0.0139688430330; ...
    0.2547895186039; ...
    0.7316386522555; ...
    0.0157253728951; ...
    0.2662302843647; ...
    0.8673504065214; ...
    0.0741493666957; ...
    0.0159285948360; ...
    0.0156061028068; ...
    0.5910094817484; ...
    0.4034771496889; ...
    0.5694745628526; ...
    0.0678493700650; ...
    0.4265968590272; ...
    0.0670982507890; ...
    0.7528310231480; ...
    0.7753727783557; ...
    0.1689073157787; ...
    0.1687335832919; ...
    0.0821244708436; ...
    0.6288705363345; ...
    0.0811413015266; ...
    0.2969112065080; ...
    0.0767542314171; ...
    0.6223022333845; ...
    0.3103786288051; ...
    0.0819218215187; ...
    0.4717022665013; ...
    0.4546603415250; ...
    0.1701091339237; ...
    0.6406004329487; ...
    0.1912267583717; ...
    0.1885315767070; ...
    0.4772929957691; ...
    0.3126974621760; ...
    0.4961225945946; ...
    0.1928805312867; ...
    0.3360041453816; ...
    0.3337280550848 ];
  elseif ( strength == 18 )
    x = [ ...
    0.0116731059668; ...
    0.9810030858388; ...
    0.0106966317092; ...
    0.9382476983551; ...
    0.0126627518417; ...
    0.0598109409984; ...
    0.0137363297927; ...
    0.9229527959405; ...
    0.0633107354993; ...
    0.0117265100335; ...
    0.1554720587323; ...
    0.8343293888982; ...
    0.8501638031957; ...
    0.0128816350522; ...
    0.1510801608959; ...
    0.0101917879217; ...
    0.2813372399303; ...
    0.7124374628501; ...
    0.2763025250863; ...
    0.0109658368561; ...
    0.4289110517884; ...
    0.4215420555115; ...
    0.5711258590444; ...
    0.5826868270511; ...
    0.0130567806713; ...
    0.0130760400964; ...
    0.7263437062407; ...
    0.0687230068637; ...
    0.8652302101529; ...
    0.0648599071037; ...
    0.1483494943362; ...
    0.0624359898396; ...
    0.7871369011735; ...
    0.0519104921610; ...
    0.1543129927444; ...
    0.2617842745603; ...
    0.7667257872813; ...
    0.2582103676627; ...
    0.0679065925147; ...
    0.5293578274804; ...
    0.0666036150484; ...
    0.0585675461899; ...
    0.0644535360411; ...
    0.6748138429151; ...
    0.3914602310369; ...
    0.6487701492307; ...
    0.3946498220408; ...
    0.5390137151933; ...
    0.1627895082785; ...
    0.6812436322641; ...
    0.1542832878020; ...
    0.2522727750445; ...
    0.2547981532407; ...
    0.1485580549194; ...
    0.2930239606436; ...
    0.2808991272310; ...
    0.4820989592971; ...
    0.5641878245444; ...
    0.1307699644344; ...
    0.1479692221948; ...
    0.5638684222946; ...
    0.4361157428790; ...
    0.3603263935285; ...
    0.4224188334674; ...
    0.3719001833052; ...
    0.2413645006928 ];
  elseif ( strength == 20 )
    x = [ ...
    0.0089411337112; ...
    0.9792622629807; ...
    0.0105475382112; ...
    0.0023777061947; ...
    0.0630425115795; ...
    0.9308422496730; ...
    0.0629076555490; ...
    0.9315962246381; ...
    0.0061951689415; ...
    0.0287125819237; ...
    0.9293844478305; ...
    0.0375457566621; ...
    0.0086895739064; ...
    0.1547597053965; ...
    0.8331025294185; ...
    0.8374231073526; ...
    0.1559362505234; ...
    0.0098599642095; ...
    0.4055873733289; ...
    0.5964727898618; ...
    0.0080747800416; ...
    0.0075073977721; ...
    0.3936764519237; ...
    0.5846530726212; ...
    0.4870804112120; ...
    0.2683512811785; ...
    0.7223956288748; ...
    0.2716826742357; ...
    0.0112580842046; ...
    0.0115034734370; ...
    0.7140525900564; ...
    0.4902871053112; ...
    0.0201423425209; ...
    0.0361107464859; ...
    0.8607998819851; ...
    0.1005891526001; ...
    0.0918740717058; ...
    0.8604888296191; ...
    0.0439842178673; ...
    0.2011017606735; ...
    0.7449993726263; ...
    0.0532186641310; ...
    0.7453984647401; ...
    0.1957289932876; ...
    0.1092532057988; ...
    0.0567625702001; ...
    0.0483837933475; ...
    0.1080612809760; ...
    0.6185605900991; ...
    0.7721296013497; ...
    0.6115734801133; ...
    0.3381326103376; ...
    0.1173084128254; ...
    0.2674551260596; ...
    0.6542100160026; ...
    0.0538297481158; ...
    0.1848840324117; ...
    0.3376267104744; ...
    0.6067102034499; ...
    0.4612614085496; ...
    0.1525465365671; ...
    0.0700582543543; ...
    0.4704201379032; ...
    0.1216461693746; ...
    0.6371404052702; ...
    0.2379904515119; ...
    0.1483929857177; ...
    0.3598069571550; ...
    0.4941441055095; ...
    0.1440630687981; ...
    0.5019764440004; ...
    0.3555423834298; ...
    0.2443439540771; ...
    0.2437064989342; ...
    0.5122200807321; ...
    0.2526038315178; ...
    0.3759895652851; ...
    0.3729077987144 ];
  elseif ( strength == 21 )
    x = [ ...
    0.0035524391922; ...
    0.0035524391922; ...
    0.9928951216156; ...
    0.9553548273730; ...
    0.0358552797177; ...
    0.9553548273730; ...
    0.0087898929093; ...
    0.0087898929093; ...
    0.0358552797177; ...
    0.8865264879047; ...
    0.8865264879047; ...
    0.0052405375935; ...
    0.0052405375935; ...
    0.1082329745017; ...
    0.1082329745017; ...
    0.0466397432150; ...
    0.0466397432150; ...
    0.9067205135700; ...
    0.2075720456946; ...
    0.2075720456946; ...
    0.7841520301770; ...
    0.0082759241284; ...
    0.0082759241284; ...
    0.7841520301770; ...
    0.0858119489725; ...
    0.8827043562574; ...
    0.0314836947701; ...
    0.0858119489725; ...
    0.8827043562574; ...
    0.0314836947701; ...
    0.6688778233826; ...
    0.0095150760625; ...
    0.0095150760625; ...
    0.6688778233826; ...
    0.3216071005550; ...
    0.3216071005550; ...
    0.4379999543113; ...
    0.0099859785681; ...
    0.4379999543113; ...
    0.0099859785681; ...
    0.5520140671206; ...
    0.5520140671206; ...
    0.7974931072148; ...
    0.0405093994119; ...
    0.0405093994119; ...
    0.1619974933734; ...
    0.7974931072148; ...
    0.1619974933734; ...
    0.3864215551955; ...
    0.3864215551955; ...
    0.2271568896090; ...
    0.8090129379329; ...
    0.0954935310336; ...
    0.0954935310336; ...
    0.2745425238718; ...
    0.0479840480721; ...
    0.6774734280561; ...
    0.6774734280561; ...
    0.2745425238718; ...
    0.0479840480721; ...
    0.4053472446667; ...
    0.0516677930989; ...
    0.4053472446667; ...
    0.5429849622344; ...
    0.0516677930989; ...
    0.5429849622344; ...
    0.1877738615539; ...
    0.7054113116872; ...
    0.7054113116872; ...
    0.1068148267588; ...
    0.1877738615539; ...
    0.1068148267588; ...
    0.1195059712009; ...
    0.1195059712009; ...
    0.5747817297348; ...
    0.5747817297348; ...
    0.3057122990643; ...
    0.3057122990643; ...
    0.5981245743363; ...
    0.2009377128319; ...
    0.2009377128319; ...
    0.2160775200005; ...
    0.3121360256673; ...
    0.2160775200005; ...
    0.3121360256673; ...
    0.4717864543321; ...
    0.4717864543321; ...
    0.4376579903849; ...
    0.4376579903849; ...
    0.1246840192303; ...
    0.3333333333333 ];
  elseif ( strength == 23 )
    x = [ ...
    0.0087809303836; ...
    0.9903675314220; ...
    0.0027029276450; ...
    0.0335909214524; ...
    0.0091675068606; ...
    0.9675568182558; ...
    0.0084737200688; ...
    0.0078781948792; ...
    0.0676785477700; ...
    0.9470266955047; ...
    0.0442974755680; ...
    0.9144243214882; ...
    0.0081735424459; ...
    0.2497452292741; ...
    0.3833232646055; ...
    0.8876850353557; ...
    0.1035329228297; ...
    0.0077255923618; ...
    0.1403192425107; ...
    0.8104591009652; ...
    0.1809643003717; ...
    0.8330767948684; ...
    0.0083010907126; ...
    0.0348407706147; ...
    0.2740287679608; ...
    0.7173982224778; ...
    0.2394976858234; ...
    0.0081859185845; ...
    0.0068836152075; ...
    0.4843741485699; ...
    0.4960767772741; ...
    0.6112936776245; ...
    0.3804323980345; ...
    0.7303890713524; ...
    0.0083987168639; ...
    0.6128525675612; ...
    0.0075475961037; ...
    0.0079525316513; ...
    0.3559774870460; ...
    0.9110236977966; ...
    0.0437233605166; ...
    0.0388480061835; ...
    0.0967032117936; ...
    0.0873226911312; ...
    0.0421445202084; ...
    0.8485617974961; ...
    0.8477921333864; ...
    0.1067435889398; ...
    0.1833966521991; ...
    0.0416340541167; ...
    0.7611632251560; ...
    0.1941599254144; ...
    0.7579378747173; ...
    0.0439826512395; ...
    0.0369760535918; ...
    0.5363187134342; ...
    0.1001256948921; ...
    0.7912266693524; ...
    0.0379866714177; ...
    0.4157414028965; ...
    0.6507106491463; ...
    0.0420141133438; ...
    0.0425548444254; ...
    0.2920627107240; ...
    0.5389729538180; ...
    0.4193031828489; ...
    0.6549472009700; ...
    0.3007352790917; ...
    0.3752400771585; ...
    0.3453980282786; ...
    0.0994532168761; ...
    0.1598309359585; ...
    0.1797326661667; ...
    0.7124584461943; ...
    0.1066065678636; ...
    0.7001701904096; ...
    0.0993303629801; ...
    0.6065648052521; ...
    0.1023223542704; ...
    0.2533382324938; ...
    0.6166226715217; ...
    0.2769500693109; ...
    0.0904184571873; ...
    0.4981522767248; ...
    0.0928231860168; ...
    0.3738418699229; ...
    0.2521678840407; ...
    0.5087500218708; ...
    0.3905579116731; ...
    0.1706141469096; ...
    0.5266737761312; ...
    0.3487581527629; ...
    0.2588053596017; ...
    0.1696614558053; ...
    0.3013521806875; ...
    0.2580202409759; ...
    0.4584740860198; ...
    0.1848898683498; ...
    0.6130740338465; ...
    0.1921611750994; ...
    0.4180541160599; ...
    0.1650612642036; ...
    0.5159205739625; ...
    0.2982718935750; ...
    0.4098894602340 ];
  elseif ( strength == 25 )
    x = [ ...
    0.0082881595033; ...
    0.4618422030241; ...
    0.0071066441239; ...
    0.9847613141699; ...
    0.5374447869049; ...
    0.0000000000000; ...
    0.4914131929361; ...
    0.0070345937020; ...
    0.9564734714228; ...
    0.0370198792045; ...
    0.1024124542747; ...
    0.5928065811509; ...
    0.0050948422371; ...
    0.0081562023689; ...
    0.0424936107568; ...
    0.9495543500844; ...
    0.8932787471239; ...
    0.0069317612927; ...
    0.9035839030665; ...
    0.0905665738209; ...
    0.0083929332787; ...
    0.6261245686071; ...
    0.0062801592979; ...
    0.8272539257367; ...
    0.0062005875353; ...
    0.1676900311185; ...
    0.7199353069567; ...
    0.2749740090237; ...
    0.0079257582005; ...
    0.0069981220752; ...
    0.8125248773263; ...
    0.0073536969970; ...
    0.7283665935411; ...
    0.1800642304565; ...
    0.2658102467762; ...
    0.0070892364520; ...
    0.3774054302043; ...
    0.0369649608668; ...
    0.9203194109805; ...
    0.0425477806431; ...
    0.6191278394983; ...
    0.3762697209178; ...
    0.0956111149690; ...
    0.0302473410377; ...
    0.8739905691754; ...
    0.8604133734958; ...
    0.0347307852352; ...
    0.1043606608343; ...
    0.7797622824754; ...
    0.0185865164256; ...
    0.0324585286618; ...
    0.8371293901157; ...
    0.0836602075315; ...
    0.0784070242501; ...
    0.4929238648458; ...
    0.1870637584073; ...
    0.4892636967025; ...
    0.0401982618372; ...
    0.7894259278865; ...
    0.1686260456429; ...
    0.3750901913174; ...
    0.0356362876880; ...
    0.5887548164804; ...
    0.0373308082182; ...
    0.2820769993374; ...
    0.6819277603320; ...
    0.0374938324382; ...
    0.6984079204127; ...
    0.2654390894079; ...
    0.1429848440800; ...
    0.7623554007647; ...
    0.0934222022749; ...
    0.5759004479923; ...
    0.3822427332525; ...
    0.0411414081675; ...
    0.0802462538379; ...
    0.7625229819410; ...
    0.1524941445131; ...
    0.0622159195833; ...
    0.1109539036076; ...
    0.4575627212057; ...
    0.4322865136374; ...
    0.5865002850241; ...
    0.0869359250818; ...
    0.0929594906936; ...
    0.6661932141454; ...
    0.4780306362227; ...
    0.4372215294577; ...
    0.6779224504669; ...
    0.2423431255660; ...
    0.2288925420305; ...
    0.3315065049959; ...
    0.3424200526607; ...
    0.0862630046475; ...
    0.5113188946635; ...
    0.1538977841001; ...
    0.6779951348472; ...
    0.1664600469411; ...
    0.0950910318888; ...
    0.3436048136712; ...
    0.5560417025366; ...
    0.1452404029513; ...
    0.1619685156238; ...
    0.5800164844262; ...
    0.2450201223288; ...
    0.2557621891794; ...
    0.2205239985511; ...
    0.4940183111285; ...
    0.2531570689798; ...
    0.5846891116357; ...
    0.1660333602278; ...
    0.2505426292461; ...
    0.3519336802182; ...
    0.3502668835419; ...
    0.4400892485512; ...
    0.4680855471546; ...
    0.1770237763947; ...
    0.3900920779501; ...
    0.2805847774120; ...
    0.3361523347440 ];
  end

  return
end