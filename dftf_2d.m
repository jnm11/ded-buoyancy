
x=linspace(-pi,pi,101);
y=linspace(-pi,pi,100)';

mesh(x,y,f);



mesh(x,y,f);

f3=(-(10.*(90518.*cos(y)-840.*cos(2.*y)-12321.*cos(3.*y)+354.*cos(4.*y)+907.*cos(5.*y)+79590)).*cos(x)+(60360.*cos(3.*y)-8240.*cos(4.*y)-7320.*cos(5.*y)+8400.*cos(y)+125120.*cos(2.*y)-55440).*cos(2.*x)+(-12975.*cos(3.*y)-12930.*cos(4.*y)+1125.*cos(5.*y)+123210.*cos(y)+60360.*cos(2.*y)+63930).*cos(3.*x)+(-12930.*cos(3.*y)-3700.*cos(4.*y)+1110.*cos(5.*y)-3540.*cos(y)-8240.*cos(2.*y)-3420).*cos(4.*x)+(1125.*cos(3.*y)+1110.*cos(4.*y)-503.*cos(5.*y)-9070.*cos(y)-7320.*cos(2.*y)-2238).*cos(5.*x)-795900.*cos(y)-55440.*cos(2.*y)+63930.*cos(3.*y)-3420.*cos(4.*y)-2238.*cos(5.*y)+2077164)/(2*1376256);

f4=-1355.*cos(4.*y).*(1/442368)+13.*cos(4.*x).*cos(3.*y).*(1/55296)+13.*cos(3.*x).*cos(4.*y).*(1/55296)-67.*cos(x).*cos(y).*(1/256)+19.*cos(3.*x).*cos(y).*(1/2304)-277.*cos(4.*x).*cos(y).*(1/55296)+1141.*cos(2.*x).*cos(y).*(1/13824)-277.*cos(4.*y).*cos(x).*(1/55296)-61.*cos(2.*x).*cos(4.*y).*(1/36864)+1141.*cos(2.*y).*cos(x).*(1/13824)-109.*cos(2.*y).*cos(3.*x).*(1/13824)-61.*cos(2.*y).*cos(4.*x).*(1/36864)+1255.*cos(2.*x).*cos(2.*y).*(1/27648)+19.*cos(3.*y).*cos(x).*(1/2304)-11.*cos(3.*x).*cos(3.*y).*(1/2304)-109.*cos(2.*x).*cos(3.*y).*(1/13824)-1355.*cos(4.*x).*(1/442368)+205.*cos(3.*x).*(1/18432)+3419.*cos(2.*x).*(1/110592)-6101.*cos(x).*(1/18432)+205.*cos(3.*y).*(1/18432)+3419.*cos(2.*y).*(1/110592)-25.*cos(4.*x).*cos(4.*y).*(1/442368)-6101.*cos(y).*(1/18432)+96181/147456;


f5=(-(120.*(-20529.*cos(2.*y)-18180.*cos(3.*y)+5490.*cos(4.*y)+940.*cos(5.*y)-415.*cos(6.*y)+130904.*cos(y)+129118)).*cos(x)+(3428745.*cos(2.*y)+333900.*cos(3.*y)-537810.*cos(4.*y)+13500.*cos(5.*y)+19335.*cos(6.*y)+2463480.*cos(y)-99390).*cos(2.*x)+(333900.*cos(2.*y)-640720.*cos(3.*y)-204120.*cos(4.*y)+77040.*cos(5.*y)+7540.*cos(6.*y)+2181600.*cos(y)+1480600).*cos(3.*x)+(-537810.*cos(2.*y)-204120.*cos(3.*y)+29220.*cos(4.*y)+33480.*cos(5.*y)+690.*cos(6.*y)-658800.*cos(y)-321540).*cos(4.*x)+(13500.*cos(2.*y)+77040.*cos(3.*y)+33480.*cos(4.*y)-1104.*cos(5.*y)-1020.*cos(6.*y)-112800.*cos(y)-82824).*cos(5.*x)+(19335.*cos(2.*y)+7540.*cos(3.*y)+690.*cos(4.*y)-1020.*cos(5.*y)+105.*cos(6.*y)+49800.*cos(y)+36190).*cos(6.*x)-99390.*cos(2.*y)+1480600.*cos(3.*y)-321540.*cos(4.*y)-82824.*cos(5.*y)+36190.*cos(6.*y)-15494160.*cos(y)+34522852)/48234496;







f2=((-cos(y)-1).*cos(x)-cos(y)+3)/4;

f3=((cos(y)+1).*cos(2.*x)+(cos(x)+1).*cos(2.*y)+(-4.*cos(y)-5).*cos(x)-5.*cos(y)+10)/16;

f4=((6.*cos(y)-cos(2.*y)+7).*cos(2.*x)+(6.*cos(x)+7).*cos(2.*y)+(-cos(y)-1).*cos(3.*x)+(-cos(x)-1).*cos(3.*y)+(-14.*cos(y)-21).*cos(x)-21.*cos(y)+35)/4^3;

f5=((27.*cos(y)-8.*cos(2.*y)+cos(3.*y)+36).*cos(2.*x)+(27.*cos(x)+cos(3.*x)+36).*cos(2.*y)+(-8.*cos(y)-9).*cos(3.*x)+(-8.*cos(x)-9).*cos(3.*y)+(cos(y)+1).*cos(4.*x)+(cos(x)+1).*cos(4.*y)+(-48.*cos(y)-84).*cos(x)-84.*cos(y)+126)/4^4;

f6=((10.*cos(3.*y)-cos(4.*y)+110.*cos(y)-44.*cos(2.*y)+165).*cos(2.*x)+(110.*cos(x)+10.*cos(3.*x)-cos(4.*x)+165).*cos(2.*y)+(-cos(3.*y)-44.*cos(y)-55).*cos(3.*x)+(-44.*cos(x)-55).*cos(3.*y)+(10.*cos(y)+11).*cos(4.*x)+(10.*cos(x)+11).*cos(4.*y)+(-cos(y)-1).*cos(5.*x)+(-cos(x)-1).*cos(5.*y)+(-165.*cos(y)-330).*cos(x)-330.*cos(y)+462)/4^5;

f7=((-208.*cos(2.*y)+65.*cos(3.*y)-12.*cos(4.*y)+cos(5.*y)+429.*cos(y)+715).*cos(2.*x)+(429.*cos(x)+65.*cos(3.*x)-12.*cos(4.*x)+cos(5.*x)+715).*cos(2.*y)+(-12.*cos(3.*y)+cos(4.*y)-208.*cos(y)-286).*cos(3.*x)+(-208.*cos(x)+cos(4.*x)-286).*cos(3.*y)+(65.*cos(y)+78).*cos(4.*x)+(65.*cos(x)+78).*cos(4.*y)+(-12.*cos(y)-13).*cos(5.*x)+(-12.*cos(x)-13).*cos(5.*y)+(cos(y)+1).*cos(6.*x)+(cos(x)+1).*cos(6.*y)+(-572.*cos(y)-1287).*cos(x)-1287.*cos(y)+1716)/4^6;

f8=((-910.*cos(2.*y)+350.*cos(3.*y)-90.*cos(4.*y)+14.*cos(5.*y)-cos(6.*y)+1638.*cos(y)+3003).*cos(2.*x)+(1638.*cos(x)+350.*cos(3.*x)-90.*cos(4.*x)+14.*cos(5.*x)-cos(6.*x)+3003).*cos(2.*y)+(-90.*cos(3.*y)+14.*cos(4.*y)-cos(5.*y)-910.*cos(y)-1365).*cos(3.*x)+(-910.*cos(x)+14.*cos(4.*x)-cos(5.*x)-1365).*cos(3.*y)+(-cos(4.*y)+350.*cos(y)+455).*cos(4.*x)+(350.*cos(x)+455).*cos(4.*y)+(-90.*cos(y)-105).*cos(5.*x)+(-90.*cos(x)-105).*cos(5.*y)+(14.*cos(y)+15).*cos(6.*x)+(14.*cos(x)+15).*cos(6.*y)+(-cos(y)-1).*cos(7.*x)+(-cos(x)-1).*cos(7.*y)+(-2002.*cos(y)-5005).*cos(x)-5005.*cos(y)+6435)/4^7;




contour(x,y,f9,linspace(0.01,0.99,10));


f4=((-4.*cos(y)-cos(2.*y)-3).*cos(2.*x)+(-4.*cos(x)-3).*cos(2.*y)+(-16.*cos(y)-12).*cos(x)-12.*cos(y)+55)/4^3;

f5 =( (-2.*cos(y)+4.*cos(2.*y)+2.*cos(3.*y)-4).*cos(2.*x)+(-2.*cos(x)+2.*cos(3.*x)-4).*cos(2.*y)+(8.*cos(y)+6).*cos(3.*x)+(8.*cos(x)+6).*cos(3.*y)+(-80.*cos(y)-70).*cos(x)-70.*cos(y)+196)/4^4;

f6 =( (-3.*cos(4.*y)+48.*cos(y)+24.*cos(2.*y)+27).*cos(2.*x)+(48.*cos(x)-3.*cos(4.*x)+27).*cos(2.*y)+(36.*cos(y)-4.*cos(3.*y)+32).*cos(3.*x)+(-12.*cos(y)-9).*cos(4.*x)+(36.*cos(x)+32).*cos(3.*y)+(-12.*cos(x)-9).*cos(4.*y)+(-324.*cos(y)-324).*cos(x)-324.*cos(y)+714)/4^5;

f7 =( (-36.*cos(3.*y)-8.*cos(4.*y)+4.*cos(5.*y)+352.*cos(y)+64.*cos(2.*y)+264).*cos(2.*x)+(352.*cos(x)-36.*cos(3.*x)-8.*cos(4.*x)+4.*cos(5.*x)+264).*cos(2.*y)+(-16.*cos(3.*y)+6.*cos(4.*y)+96.*cos(y)+110).*cos(3.*x)+(96.*cos(x)+6.*cos(4.*x)+110).*cos(3.*y)+(-86.*cos(y)-72).*cos(4.*x)+(-86.*cos(x)-72).*cos(4.*y)+(16.*cos(y)+12).*cos(5.*x)+(16.*cos(x)+12).*cos(5.*y)+(-1232.*cos(y)-1386).*cos(x)-1386.*cos(y)+2640)/4^6;

f8 =( (26.*cos(2.*y)-220.*cos(3.*y)+26.*cos(4.*y)+20.*cos(5.*y)-5.*cos(6.*y)+1820.*cos(y)+1573).*cos(2.*x)+(1820.*cos(x)-220.*cos(3.*x)+26.*cos(4.*x)+20.*cos(5.*x)-5.*cos(6.*x)+1573).*cos(2.*y)+(-16.*cos(3.*y)+40.*cos(4.*y)-8.*cos(5.*y)+104.*cos(y)+260).*cos(3.*x)+(104.*cos(x)+40.*cos(4.*x)-8.*cos(5.*x)+260).*cos(3.*y)+(-9.*cos(4.*y)-400.*cos(y)-377).*cos(4.*x)+(-400.*cos(x)-377).*cos(4.*y)+(152.*cos(y)+124).*cos(5.*x)+(152.*cos(x)+124).*cos(5.*y)+(-20.*cos(y)-15).*cos(6.*x)+(-20.*cos(x)-15).*cos(6.*y)+(-4576.*cos(y)-5720).*cos(x)-5720.*cos(y)+9867)/4^7;



mesh(f8);




% 2nd and 4th derivatives vanish at boundary
f6=((-90.*cos(y)-36.*cos(2.*y)-6.*cos(3.*y)-60).*cos(2.*x)+(-15.*cos(y)-6.*cos(2.*y)-cos(3.*y)-10).*cos(3.*x)+(-90.*cos(x)-60).*cos(2.*y)+(-15.*cos(x)-10).*cos(3.*y)+(-225.*cos(y)-150).*cos(x)-150.*cos(y)+924)/4^5;

f7=((3.*cos(4.*y)+108.*cos(y)+60.*cos(2.*y)+20.*cos(3.*y)+65).*cos(3.*x)+(108.*cos(x)+60.*cos(2.*x)+3.*cos(4.*x)+65).*cos(3.*y)+(18.*cos(4.*y)-252.*cos(y)-210).*cos(2.*x)+(45.*cos(y)+18.*cos(2.*y)+30).*cos(4.*x)+(-252.*cos(x)-210).*cos(2.*y)+(45.*cos(x)+30).*cos(4.*y)+(-1260.*cos(y)-945).*cos(x)-945.*cos(y)+3396)/4^7;

f8=((252.*cos(3.*y)-36.*cos(4.*y)-36.*cos(5.*y)-216.*cos(y)+432.*cos(2.*y)-396).*cos(2.*x)+(-28.*cos(3.*y)-36.*cos(4.*y)-6.*cos(5.*y)+834.*cos(y)+252.*cos(2.*y)+584).*cos(3.*x)+(-216.*cos(x)-36.*cos(4.*x)-36.*cos(5.*x)-396).*cos(2.*y)+(834.*cos(x)-36.*cos(4.*x)-6.*cos(5.*x)+584).*cos(3.*y)+(-9.*cos(4.*y)+36.*cos(y)+45).*cos(4.*x)+(36.*cos(x)+45).*cos(4.*y)+(-90.*cos(y)-60).*cos(5.*x)+(-90.*cos(x)-60).*cos(5.*y)+(-5544.*cos(y)-4620).*cos(x)-4620.*cos(y)+12639)/4^8;

mesh(f8);



% These one are maximally round and then they can be widened by taking a larger power



f=((-90.*cos(y)-36.*cos(2.*y)-6.*cos(3.*y)-60).*cos(2.*x)+(-15.*cos(y)-6.*cos(2.*y)-cos(3.*y)-10).*cos(3.*x)+(-90.*cos(x)-60).*cos(2.*y)+(-15.*cos(x)-10).*cos(3.*y)+(-225.*cos(y)-150).*cos(x)-150.*cos(y)+924)/1024;

% This is larger than 1 % Don't use;
f=((-18.*cos(4.*y)-828.*cos(y)-432.*cos(2.*y)-132.*cos(3.*y)-510).*cos(2.*x)+(-3.*cos(4.*y)-288.*cos(y)-132.*cos(2.*y)-32.*cos(3.*y)-185).*cos(3.*x)+(-828.*cos(x)-18.*cos(4.*x)-510).*cos(2.*y)+(-288.*cos(x)-3.*cos(4.*x)-185).*cos(3.*y)+(-45.*cos(y)-30).*cos(4.*x)+(-45.*cos(x)-30).*cos(4.*y)+(-1440.*cos(y)-855).*cos(x)-855.*cos(y)+7692)/8192;


% This is larger than 1 % Don't use;
f=((-41088.*cos(3.*y)+2124.*cos(4.*y)+864.*cos(5.*y)-601056.*cos(y)-245808.*cos(2.*y)-397596).*cos(2.*x)+(32.*cos(3.*y)+2424.*cos(4.*y)+144.*cos(5.*y)-129456.*cos(y)-41088.*cos(2.*y)-90616).*cos(3.*x)+(-601056.*cos(x)+2124.*cos(4.*x)+864.*cos(5.*x)-397596).*cos(2.*y)+(-129456.*cos(x)+2424.*cos(4.*x)+144.*cos(5.*x)-90616).*cos(3.*y)+(621.*cos(4.*y)-3384.*cos(y)-3705).*cos(4.*x)+(-3384.*cos(x)-3705).*cos(4.*y)+(2160.*cos(y)+1440).*cos(5.*x)+(2160.*cos(x)+1440).*cos(5.*y)+(-1347264.*cos(y)-870120).*cos(x)-870120.*cos(y)+5848149)/(6256.*1024);


