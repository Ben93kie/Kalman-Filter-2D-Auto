%Im Folgenden Skript wird ein einfacher extended 3-D Kalman-Filter
%implementiert, der Vorhersagen über die Trajektorie bzw die Position von
%Autos basierend auf einem vereinfachten mathematischen Modell und
%Positionsmessungen macht. Das Modell berücksichtigt dabei x-Position,
%y-Position und Gierrate. Die nächste Position ergibt sich dabei aus der
%alten geschätzten Position plus die Geschwindigkeit in x- bzw. y-Richtung
%in Richtung der geschätzten Gierrate. Wichtig: Es ist linear, sprich, in
%jedem Schritt wird die nächste Position anhand des Gradienten bestimmt.
%Geplottet wird Auto1, Auto2, die Vorhersage von Auto 2 (rot) in Form von
%Rechtecken und weiterhin: grüner Kreis bzw. grünes Quadrat ist die
%korrigierte Kalman-Filter-Prädiktion in jedem Zeitschritt mit Messfehlern
%im mathematischen Modell und in den Positionsmessungen. Der gelbe Kreis
%bzw das gelbe Quadrat ist die Prädiktion, wo das jeweilige Auto in einer
%bestimmten Zeit sein wird basierend auf exakten mathematischen Modell. Der
%blaue Kreis bzw das blaue Quadrat ist die geschätzte Position zu jedem
%Zeitpunkt basierend auf dem mathematischen Modell mit der letzten 
%Prädiktion und mit Messfehlern. Der rote Kreis bzw. das rote Quadrat ist
%die korrigierte Positionsschätzung zu jedem Zeitpunkt, nachdem die Messung
%und das mathematische Modell fusioniert wurden mit Messfehlern. 

clc; clear all ;
%Zeit
tfor = linspace(0,4*pi) ;
tplot=zeros(length(tfor));

syms t
%Definiere Trajektorie von Auto1
%HINWEIS: zum linearen Term wird ein nichtlinearer hinzugefügt nur aus dem
%Grund, dass die symbolische Ableitung nicht konstant ist. Denn wenn diese
%konstant ist, ergibt sich nach Auswerten eines Vektor nur ein Skalar (da
%die Aleitung an jedem Punkt konstant ist. Kann leicht geändert werden. Der
%Einfachkeit halber wurde ein winiziger nichtlinearer Teil addiert.
x=t+(1/(t+2))+4;
y=1/(t+8)+2*t;
%Trajektorie Auto 2
w=-t+(1/(t+6))+10;
z=4+1/(t+100)+t/8;%+sin(t); %uncomment für eine nichtlineare Beispieltraj
%ektorie von Auto 2

%Symbolische Ableitungen der Trajektorien für Gierrate
dw=diff(w,t);
dz=diff(z,t);
dx=diff(x,t);
dxjac=diff(x,t);
dy=diff(y,t);
dyjac=diff(y,t);
t=linspace(0,4*pi);
xnew=double(subs(x));
ynew=double(subs(y));
dxnew=double(subs(dx));
dynew=double(subs(dy));
wnew=double(subs(w));
znew=double(subs(z));
dwnew=double(subs(dw));
dznew=double(subs(dz));
%Definiere Auto 1
  xi=[-0.5 0.5 0.5 -0.5];
  yi=[-0.25 -0.25 0.25 0.25];
 %Auto 2 
  wi=[-0.5 0.5 0.5 -0.5];
  zi=[-0.25 -0.25 0.25 0.25];
%Initialisierung
 mittelpunktaltx=0;
 mittelpunktalty=0;
 firstmittelpunktx=mittelpunktaltx;
 firstmittelpunkty=mittelpunktalty;
 sum=0;
 %Aus Plotgründen extra Variable
 xiplot=xi;
 yiplot=yi;
 wiplot=wi;
 ziplot=zi;
 %Gierrate basierend auf Ableitungen der Trajektorie
    if dxnew(1)<0
         gierrate=atan(dynew(1)/dxnew(1))+pi;
     else
         gierrate=atan(dynew(1)/dxnew(1));
    end
    
     if dwnew(1)<0
         gierrate2=atan(dznew(1)/dwnew(1))+pi;
     else
         gierrate2=atan(dznew(1)/dwnew(1));
     end
    
    %kx ist der state-Vektor von Auto 1
 kx=[xnew(1) ynew(1) gierrate]';
 %kxx ist der state-Vektor von Auto 2
 kkx=[wnew(1) znew(1) gierrate2]';
 
 %kxjustmathpredict ist die Trajektorie nur basierend auf dem
 %mathematischen Modell für Auto 1
 
 kxjustmathpredict=[xnew(1) ynew(1) gierrate]';
 
 %entsprechend für Auto 2
 kkxjustmathpredict=[wnew(1) znew(1) gierrate2]';
 
 %kP bzw kkP ist die Varianz (=Fehler des mathematischen Modells am Anfang, sprich
 %Initialfehler)
 kP=zeros(3,3)+0.000001*eye(3);
 kkP=zeros(3,3)+0.000001*eye(3);
 t=gierrate;
 tt=gierrate2;
 %kF bzw kkF sind die linearisierten Übergangsmatrizen des extended
 %Kalman-Filters. Es wurde also per Hand schon abgeleitet. Die
 %Übergangsfunktion ergibt sich aus dem neuen kx bzw kkx unten.
 kF=[[1 0 -0.35*sin(gierrate)];[0 1 0.35*cos(gierrate)];[0 0 1]];
 kkF=[[1 0 -0.35*sin(gierrate2)];[0 1 0.35*cos(gierrate2)];[0 0 1]];
 %"Fehler" der Messung
 kQ=0.1*eye(3);
 kkQ=0.1*eye(3);
 %Zufallsvektor der auf die Messung addiert wird
 kw=mvnrnd([0 0 0],kQ)';
 kkw=mvnrnd([0 0 0],kkQ)';
 %%%%Wird nicht benutzt, kann aber hilfreich zum Verstehen sein
 drehungvector=zeros(length(tfor),1);
 drehungvector2=zeros(length(tfor),1);
 approxdreh=zeros(length(tfor),1);
 approxdreh2=zeros(length(tfor),1);
 %%%%
 
 %Für die plots am Ende
 fehlermessung=zeros(length(tfor),1);
 fehlermessung2=zeros(length(tfor),1);
 fehlercorrected=zeros(length(tfor),1);
 fehlercorrected2=zeros(length(tfor),1);
 fehlerpredict=zeros(length(tfor),1);
 %Matrix, die bestimmt, was gemessen wird
 H=[[1 0 0];[0 1 0];[0 0 1]];
 for i = 1:length(tfor)
     
%in jedem Schritt, wo der Mittelpunkt des Autos ist
     mittelpunktx=xnew(i);
     mittelpunkty=ynew(i);
     mittelpunktw=wnew(i);
     mittelpunktz=znew(i);
     ableitungsrichtungx=dxnew(i);
     ableitungsrichtungy=dynew(i);
     ableitungsrichtungw=dwnew(i);
     ableitungsrichtungz=dznew(i);
     mittelpunktdiffx=mittelpunktx-mittelpunktaltx;
     mittelpunktdiffy=mittelpunkty-mittelpunktalty;
     %Update der Übergangsfunktion
     kF=[[1 0 -0.25*sin(kx(3))];[0 1 0.25*cos(kx(3))];[0 0 1]];
     kkF=[[1 0 -0.13*sin(kkx(3))];[0 1 0.13*cos(kkx(3))];[0 0 1]];
     mittelpunktaltx=mittelpunktx;
     mittelpunktalty=mittelpunkty;
     newmatrix=vertcat(xi,yi);
     newmatrix2=vertcat(wi,zi);
    %Update der Fahrtrichtung/Gierrate
     if ableitungsrichtungx<0
         drehung=atan(ableitungsrichtungy/ableitungsrichtungx)+pi;
     else
         drehung=atan(ableitungsrichtungy/ableitungsrichtungx);
     end
     if ableitungsrichtungw<0
         drehung2=atan(ableitungsrichtungz/ableitungsrichtungw)+pi;
     else
         drehung2=atan(ableitungsrichtungz/ableitungsrichtungw);
     end
     %Muss Auto in jedem Zeitpunkt entsprechend der Gierrate(=drehung)
     %drehen
     rotationmatrix=[cos(drehung) -sin(drehung); sin(drehung) cos(drehung)];
     rotationmatrix2=[cos(drehung2) -sin(drehung2); sin(drehung2) cos(drehung2)];
     
     newmatrixrotated=rotationmatrix*newmatrix;
     newmatrixrotated2=rotationmatrix2*newmatrix2;
     
     newmatrixrotatedx=newmatrixrotated(1,:);
     newmatrixrotatedy=newmatrixrotated(2,:);
     
     newmatrixrotatedw=newmatrixrotated2(1,:);
     newmatrixrotatedz=newmatrixrotated2(2,:);
     
     drehungvector(i)=drehung;
     drehungvector2(i)=drehung2;
   %Plot der Autos zu jedem Zeitpunkt
     xiplot=newmatrixrotatedx+mittelpunktx;
     yiplot=newmatrixrotatedy+mittelpunkty;
     wiplot=newmatrixrotatedw+mittelpunktw;
     ziplot=newmatrixrotatedz+mittelpunktz;
     
     %Prädiktion anhand des mathematischen Modells und der vorherigen
     %korrigierten Position
     kx=[kx(1)+0.25*cos(kx(3)),kx(2)+0.25*sin(kx(3)),kx(3)]';
     kkx=[kkx(1)+0.13*cos(kkx(3)),kkx(2)+0.13*sin(kkx(3)),kkx(3)]';
     %Kalman-Schritte im FOlgenden
     kP=kF*kP*kF'+kQ;
     kkP=kkF*kkP*kkF'+kkQ;
     
     kR=0.005*eye(3);
     kkR=0.005*eye(3);
     
     
     kv=mvnrnd([0,0,0],kR);
     kkv=mvnrnd([0,0,0],kkR);
     
     kz=[mittelpunktx,mittelpunkty,drehung]'+kv';
     kkz=[mittelpunktw,mittelpunktz,drehung2]'+kkv';
     
     kxjustmathpredict=[mittelpunktx+cos(drehung)*2,mittelpunkty+sin(drehung)*2,kxjustmathpredict(3)];
     kkxjustmathpredict=[mittelpunktw+cos(drehung2)*2,mittelpunktz+sin(drehung2)*2,kkxjustmathpredict(3)];
     
     ky=kz-H*kx;
     kky=kkz-H*kkx;
     
     kS=H*kP*H'+kR;
     kkS=H*kkP*H'+kkR;
     
     Kalman=kP*H'*inv(kS);
     kKalman=kkP*H'*inv(kkS);
     %Korrigierte Position nach Anwendung des Kalman-Filters
     kxcorrected=kx+Kalman*ky;
     kkxcorrected=kkx+kKalman*kky;
     %Geschätzte Position in der Zukunft basierend auf korrigierter
     %Position (hier 2 units in die Zukunft)
     predictionbycar1=[kkxcorrected(1)+cos(kkxcorrected(3))*2,kkxcorrected(2)+sin(kkxcorrected(3))*2,kkxcorrected(3)]';
     predictionbycar2=[kxcorrected(1)+cos(kxcorrected(3))*2,kxcorrected(2)+sin(kxcorrected(3))*2,kxcorrected(3)]';
     %fürs Plotten der Prediction
     rotationmatrixpred=[cos(kkxcorrected(3)) -sin(kkxcorrected(3)); sin(kkxcorrected(3)) cos(kkxcorrected(3))];
     newmatrixrotatedpred=rotationmatrixpred*newmatrix2;
     newmatrixrotatedpredx=newmatrixrotatedpred(1,:);
     newmatrixrotatedpredy=newmatrixrotatedpred(2,:);
     
     predplotx=newmatrixrotatedpredx+kkxcorrected(1)+cos(kkxcorrected(3))*2;
     predploty=newmatrixrotatedpredy+kkxcorrected(2)+sin(kkxcorrected(3))*2;
     
      approxdreh(i)=kxcorrected(3);
      %Neuer Zufallsvektor für nächste Schleife
     kw=mvnrnd([0 0 0],kQ)';
     kkw=mvnrnd([0 0 0],kkQ)';
     %Update des mathematischen Fehlers
     kP=kP-Kalman*H*kP;
     kkP=kkP-Kalman*H*kkP;
     
    
     %Plotten. c ist für die Farben der Autos da
     c = [0; 6 ; 4 ; 3 ];
     %Erklärung zu plots am Anfang
     plot(xnew,ynew)
     hold on
     plot(kx(1), kx(2), 'bo', 'MarkerSize', 20);
     hold on
     plot(kkx(1),kkx(2), 'bs', 'MarkerSize', 20);
     hold on
     plot(kxcorrected(1),kxcorrected(2), 'ro','MarkerSize',20)
     hold on
     plot(kkxcorrected(1),kkxcorrected(2), 'rs', 'MarkerSize',20)
     hold on
     plot(kz(1),kz(2),'go','MarkerSize',15)
     hold on
     plot(kkz(1),kkz(2),'gs','MarkerSize',15)
     axis([0 14 0 14])
     hold on
     plot(kxjustmathpredict(1),kxjustmathpredict(2),'yo','MarkerSize',20)
     hold on
     plot(kkxjustmathpredict(1),kkxjustmathpredict(2),'ys','MarkerSize',20)
     hold on
     %uncomment, wenn die Prädiktion der beiden Autos in Form von Kreisen
     %bzw. Quadraten angezeigt werden soll.
%      plot(predictionbycar1(1),predictionbycar1(2),'cs','MarkerSize',20)
%      hold on
%      plot(predictionbycar2(1),predictionbycar2(2),'co','MarkerSize',20)
%      hold on
     plot(wnew,znew)
     hold on
     patch(xiplot,yiplot,c)  
     hold on
     patch(wiplot,ziplot,c)
     hold on
    
     patch(predplotx,predploty,'red')
     
     hold off
     pause(0.5)
     
     
     fehlermessung(i)=norm([mittelpunktx,mittelpunkty]-[kz(1) kz(2)]);
     fehlercorrected(i)=norm([mittelpunktx,mittelpunkty]-[kxcorrected(1),kxcorrected(2)]);
     fehlerpredict(i)=norm([mittelpunktx,mittelpunkty]-[kx(1),kx(2)]);
         
     kx=kxcorrected;
     kkx=kkxcorrected;
 end
 %fehlermessung gibt an, wie schlecht die Messung der Position in jedem 
 %Zeitpunkt ist
plot(tfor,fehlermessung,'r')
hold on
%fehlercorrected gibt an, wie gut die korrgierte Messung (nach
%Kalman-Filter) ist
plot(tfor, fehlercorrected,'b')

