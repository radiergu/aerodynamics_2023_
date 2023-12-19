clear
close all

%%Thin airfoil theory%%

alpha = deg2rad(5);

xc = linspace(0,1,100);

deltaCp = 4*alpha*(2-2*xc./sqrt(1-(1-2*xc).^2));


%%Data extracted from paper%%

x1_paper = [0 0.02392344498 0.04625199362 0.0701754386 0.09728867624 0.1228070175 0.1467304625 0.1722488038 0.1977671451 0.2216905901 0.2472089314 0.2727272727 0.298245614 0.3237639553 0.3476874003 0.3732057416 0.3987240829 0.4226475279 0.4481658692 0.4720893142 0.4976076555 0.5231259968 0.5486443381 0.5725677831 0.5980861244 0.6236044657 0.6475279107 0.6746411483 0.6985645933 0.7256778309 0.7496012759 0.7735247209 0.7990430622 0.8245614035 0.8484848485 0.8740031898 0.8995215311 0.9234449761 0.9489633174 0.9728867624 0.9984051037];

delta_cp_paper1 = [3.896202532,1.074141049,0.9764918626,0.8983725136,0.8397830018,0.8007233273,0.7909584087,0.7616636528,0.7421338156,0.703074141,0.7226039783,0.703074141,0.6737793852,0.6835443038,0.6737793852,0.6640144665,0.6542495479,0.6640144665,0.6347197107,0.6054249548,0.6054249548,0.5858951175,0.5858951175,0.5566003617,0.5370705244,0.5175406872,0.4882459313,0.4589511754,0.4296564195,0.4003616637,0.3613019892,0.3417721519,0.3027124774,0.2734177215,0.2148282098,0.1952983725,0.1464737794,0.1171790235,0.06835443038,0,0];


x3_paper = [0.00159744408945681 0.0239616613418530 0.0479233226837061 0.0734824281150160 0.0990415335463258 0.123003194888179 0.148562300319489 0.170926517571885 0.198083067092652 0.223642172523962 0.249201277955272 0.274760383386581 0.298722044728434 0.324281150159744 0.349840255591054 0.373801916932907 0.397763578274760 0.424920127795527 0.448881789137380 0.474440894568690 0.498402555910543 0.522364217252396 0.549520766773163 0.573482428115016 0.599041533546326 0.624600638977636 0.650159744408946 0.675718849840256 0.699680511182109 0.723642172523962 0.747603833865815 0.774760383386582 0.800319488817891 0.825878594249201 0.849840255591054 0.873801916932907 0.899361022364217 0.923322683706070 0.95047923322683 0.971246006389776 1];

delta_cp_paper3 = [0.009782608696 0.7630434783 0.7239130435 0.7043478261 0.6847826087 0.675 0.6652173913 0.6456521739 0.6456521739 0.6260869565 0.6163043478 0.6065217391 0.6065217391 0.5967391304 0.577173913 0.5869565217 0.5673913043 0.5576086957 0.5576086957 0.547826087 0.5380434783 0.5282608696 0.5086956522 0.4793478261 0.4695652174 0.4597826087 0.45 0.4304347826 0.4010869565 0.3815217391 0.3423913043 0.3032608696 0.2739130435 0.2445652174 0.1956521739 0.1565217391 0.1173913043 0.06847826087 0.02934782609 -0.04891304348 -0.01956521739];


x5_paper = [0.00149713386534784 0.0254209284028323 0.0493447229403170 0.0748634371136338 0.100382151286951 0.124305945824435 0.149824659997752 0.173748454535237 0.200862088344386 0.224785882881870 0.250304597055187 0.274228391592672 0.299747105765988 0.326860739575138 0.349189614476790 0.374708328650107 0.400227042823424 0.425745756996740 0.449669551534225 0.473593346071709 0.499112060245026 0.523035854782511 0.550149488591660 0.574073283129144 0.599591997302461 0.623515791839946 0.649034506013263 0.672958300550747 0.700071934359897 0.723995728897381 0.751109362706530 0.773438237608182 0.800551871417331 0.826070585590648 0.849994380128133 0.873918174665617 0.899436888838934 0.924955603012251 0.948879397549736 0.972803192087220 0.999916825896369];

delta_cp_paper5 = [0 0.5471840981 0.5569552427 0.5471840981 0.5569552427 0.5569552427 0.5667263873 0.5667263873 0.5569552427 0.5569552427 0.5569552427 0.5667263873 0.5667263873 0.5667263873 0.5569552427 0.5374129535 0.5569552427 0.5374129535 0.5374129535 0.5276418089 0.5178706643 0.5080995196 0.4885572304 0.4787860858 0.4592437966 0.4397015074 0.4299303628 0.4103880736 0.3908457843 0.3615323505 0.3224477721 0.3029054829 0.2638209044 0.2442786152 0.2149651814 0.1856517476 0.1367960245 0.07816915687 0.09771144609 -0.03908457843 -0.03908457843];

%%Plot%%

plot(xc,deltaCp)
hold on
plot(x1_paper, delta_cp_paper1)
plot(x3_paper, delta_cp_paper3)
plot(x5_paper, delta_cp_paper5)

legend({'Thin airfoil theory','t/c = 1%', 't/c = 3%', 't/c = 5%'},'Location','northeast')

xlabel('x/c') 
ylabel('\DeltaCp')