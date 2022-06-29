%Console and data clear
clear all; clc;

%Liczba probek 801
N_probek=801;

%DOWLOAND DATA
Data=load('heating_system.dat');

%experiment
Data_x=Data(:,1);
Data_voltage=Data(:,2);
Data_results=Data(:,3);

%CODE

%first plot
figure(1)
Plot1=plot(Data_x,Data_results);
ylabel('*C');
xlabel('time s');
title('Data results in case of samples');

%2nd plot
figure(2)
Plot2=plot(Data_x,Data_voltage);
ylabel('V');
xlabel('time s');
title('voltage in case of samples');

%metoda nieparametryczna
f=0.5*(0:N_probek-1)/N_probek; %wektor czêstotliwoœci w [Hz]
[Px,wx] = periodogram(Data_voltage,[],N_probek,0.5);
[Py,wy] = periodogram(Data_results,[],N_probek,0.5);
figure(3)
Widmowa=Py/Px;
plot(Widmowa(:,1));
title('OdpowiedŸ impulsowa na podstawie widmowej gêstoœci mocy');
xlim([0 50])

%opóŸnieie transportowe wynosi oko³o 2s
%obiekt jest nieliniowy

%METODA LS

%pierwsza próbka
fi_1 = [0 0];
fi_2 = [0 0 0 0];
fi_3 = [0 0 0 0 0 0];
%druga próbka
fi_1 = [fi_1; Data_voltage(1) Data_results(1)]
fi_2 = [fi_2; Data_voltage(1) 0 Data_results(1) 0]
fi_3 = [fi_3; Data_voltage(1) 0 0 Data_results(1) 0 0];
%trzecia próbka
fi_1 = [fi_1; Data_voltage(2) Data_results(2)];
fi_2 = [fi_2; Data_voltage(2) Data_voltage(1) Data_results(2) Data_results(1)];
fi_3 = [fi_3; Data_voltage(2) Data_voltage(1) 0 Data_results(2) Data_results(1) 0];


for i=4:N_probek
    fi_1 = [fi_1; Data_voltage(i-1) Data_results(i-1)];
    fi_2 = [fi_2; Data_voltage(i-1) Data_voltage(i-2) Data_results(i-1) Data_results(i-2)];
    fi_3 = [fi_3; Data_voltage(i-1) Data_voltage(i-2) Data_voltage(i-3) Data_results(i-1) Data_results(i-2) Data_results(i-3)];
end 

pLS_1=pinv(fi_1)*Data_results;
pLS_2=pinv(fi_2)*Data_results;
pLS_3=pinv(fi_3)*Data_results;

%pierwsza próbka
y_1 = [0];
y_2 = [0];
y_3 = [0];
%druga próbka
y_1 = [y_1; pLS_1(1)*Data_voltage(1) + pLS_1(2)*Data_results(1)];
y_2 = [y_2; pLS_2(1)*Data_voltage(1)+ 0 + pLS_2(3)*Data_results(1) + 0];
y_3 = [y_3; pLS_3(1)*Data_voltage(1)+ pLS_3(4)*Data_results(1)];
%trzecia próbka
y_1 = [y_1; pLS_1(1)*Data_voltage(2)+ pLS_1(2)*Data_results(2)];
y_2 = [y_2; pLS_2(1)*Data_voltage(2)+ pLS_2(2)*Data_voltage(1) + pLS_2(3)*Data_results(2) + pLS_2(4)*Data_results(1)];
y_3 = [y_3; pLS_3(1)*Data_voltage(2)+ pLS_3(2)*Data_voltage(1) + pLS_3(4)*Data_results(2) + pLS_3(5)*Data_results(1)];

for i=4:N_probek
    y_1 = [y_1; pLS_1(1)*Data_voltage(i-1)+ pLS_1(2)*Data_results(i-1)];
    y_2 = [y_2; pLS_2(1)*Data_voltage(i-1)+ pLS_2(2)*Data_voltage(i-2) + pLS_2(3)*Data_results(i-1) + pLS_2(4)*Data_results(i-2)];
    y_3 = [y_3; pLS_3(1)*Data_voltage(i-1)+ pLS_3(2)*Data_voltage(i-2) + pLS_3(3)*Data_voltage(i-3) + pLS_3(4)*Data_results(i-1) + pLS_3(5)*Data_results(i-2) + pLS_3(6)*Data_results(i-3)];
end 

%pierwsza próbka
ys_1 = [0];
ys_2 = [0];
ys_3 = [0];
%druga próbka
ys_1 = [ys_1; pLS_1(1)*Data_voltage(1) + pLS_1(2)*ys_1(1)];
ys_2 = [ys_2; pLS_2(1)*Data_voltage(1)+ 0 + pLS_2(3)*ys_2(1) + 0];
ys_3 = [ys_3; pLS_3(1)*Data_voltage(1)+ pLS_3(4)*ys_3(1)];
%trzecia próbka
ys_1 = [ys_1; pLS_1(1)*Data_voltage(2)+ pLS_1(2)*ys_1(2)];
ys_2 = [ys_2; pLS_2(1)*Data_voltage(2)+ pLS_2(2)*Data_voltage(1) + pLS_2(3)*ys_2(2) + pLS_2(4)*ys_2(1)];
ys_3 = [ys_3; pLS_3(1)*Data_voltage(2)+ pLS_3(2)*Data_voltage(1) + pLS_3(4)*ys_3(2) + pLS_3(5)*ys_3(1)];

for i=4:N_probek
    ys_1 = [ys_1; pLS_1(1)*Data_voltage(i-1)+ pLS_1(2)*ys_1(i-1)];
    ys_2 = [ys_2; pLS_2(1)*Data_voltage(i-1)+ pLS_2(2)*Data_voltage(i-2) + pLS_2(3)*ys_2(i-1) + pLS_2(4)*ys_2(i-2)];
    ys_3 = [ys_3; pLS_3(1)*Data_voltage(i-1)+ pLS_3(2)*Data_voltage(i-2) + pLS_3(3)*Data_voltage(i-3) + pLS_3(4)*ys_3(i-1) + pLS_3(5)*ys_3(i-2) + pLS_3(6)*ys_3(i-3)];
end

figure(4)
hold on;
plot(y_1,'r-')
plot(y_2,'g-')
plot(y_3,'b-')
title('Predykcja jednokrokowa wraz z odpowiedzi¹ obiektu metoda LS')
plot(Data_results,'k')
legend('y_1','Y_2','y_3','obiekt')
hold off

figure(5)
hold on;
plot(ys_1,'r-')
plot(ys_2,'g-')
plot(ys_3,'b-')
title('OdpowiedŸ modelu wraz z odpowiedzi¹ obiektu metoda LS')
plot(Data_results,'k')
legend('ys_1','ys_2','ys_3','obiekt')
hold off

%METODA IV
Gp=tf([1],[1 1 1 0.1 0.1 0.1],2)
Yp = lsim(Gp, Data_voltage, Data_x);

%pierwsza próbka
z_1 = [0 0];
z_2 = [0 0 0 0];
z_3 = [0 0 0 0 0 0];
%druga próbka
z_1 = [z_1; Data_voltage(1) Yp(1)]
z_2 = [z_2; Data_voltage(1) 0 Yp(1) 0]
z_3 = [z_3; Data_voltage(1) 0 0 Yp(1) 0 0];
%trzecia próbka
z_1 = [z_1; Data_voltage(2) Yp(2)];
z_2 = [z_2; Data_voltage(2) Data_voltage(1) Yp(2) Yp(1)];
z_3 = [z_3; Data_voltage(2) Data_voltage(1) 0 Yp(2) Yp(1) 0];

for i=4:N_probek
    z_1 = [z_1; Data_voltage(i-1) Yp(i-1)];
    z_2 = [z_2; Data_voltage(i-1) Data_voltage(i-2) Yp(i-1) Yp(i-2)];
    z_3 = [z_3; Data_voltage(i-1) Data_voltage(i-2) Data_voltage(i-3) Yp(i-1) Yp(i-2) Yp(i-3)];
end 

piV_1=inv(z_1'*fi_1)*z_1'*Data_results
piV_2=inv(z_2'*fi_2)*z_2'*Data_results
piV_3=inv(z_3'*fi_3)*z_3'*Data_results

%pierwsza próbka
yiv_1 = [0];
yiv_2 = [0];
yiv_3 = [0];
%druga próbka
yiv_1 = [yiv_1; piV_1(1)*Data_voltage(1) + piV_1(2)*Data_results(1)];
yiv_2 = [yiv_2; piV_2(1)*Data_voltage(1)+ 0 + piV_2(3)*Data_results(1) + 0];
yiv_3 = [yiv_3; piV_3(1)*Data_voltage(1)+ piV_3(4)*Data_results(1)];
%trzecia próbka
yiv_1 = [yiv_1; piV_1(1)*Data_voltage(2)+ piV_1(2)*Data_results(2)];
yiv_2 = [yiv_2; piV_2(1)*Data_voltage(2)+ piV_2(2)*Data_voltage(1) + piV_2(3)*Data_results(2) + piV_2(4)*Data_results(1)];
yiv_3 = [yiv_3; piV_3(1)*Data_voltage(2)+ piV_3(2)*Data_voltage(1) + piV_3(4)*Data_results(2) + piV_3(5)*Data_results(1)];

for i=4:N_probek
    yiv_1 = [yiv_1; piV_1(1)*Data_voltage(i-1)+ piV_1(2)*Data_results(i-1)];
    yiv_2 = [yiv_2; piV_2(1)*Data_voltage(i-1)+ piV_2(2)*Data_voltage(i-2) + piV_2(3)*Data_results(i-1) + piV_2(4)*Data_results(i-2)];
    yiv_3 = [yiv_3; piV_3(1)*Data_voltage(i-1)+ piV_3(2)*Data_voltage(i-2) + piV_3(3)*Data_voltage(i-3) + piV_3(4)*Data_results(i-1) + piV_3(5)*Data_results(i-2) + piV_3(6)*Data_results(i-3)];
end 

%pierwsza próbka
yivs_1 = [0];
yivs_2 = [0];
yivs_3 = [0];
%druga próbka
yivs_1 = [yivs_1; piV_1(1)*Data_voltage(1) + piV_1(2)*yivs_1(1)];
yivs_2 = [yivs_2; piV_2(1)*Data_voltage(1)+ 0 + piV_2(3)*yivs_2(1) + 0];
yivs_3 = [yivs_3; piV_3(1)*Data_voltage(1)+ piV_3(4)*yivs_3(1)];
%trzecia próbka
yivs_1 = [yivs_1; piV_1(1)*Data_voltage(2)+ piV_1(2)*yivs_1(2)];
yivs_2 = [yivs_2; piV_2(1)*Data_voltage(2)+ piV_2(2)*Data_voltage(1) + piV_2(3)*yivs_2(2) + piV_2(4)*yivs_2(1)];
yivs_3 = [yivs_3; piV_3(1)*Data_voltage(2)+ piV_3(2)*Data_voltage(1) + piV_3(4)*yivs_3(2) + piV_3(5)*yivs_3(1)];

for i=4:N_probek
    yivs_1 = [yivs_1; piV_1(1)*Data_voltage(i-1)+ piV_1(2)*yivs_1(i-1)];
    yivs_2 = [yivs_2; piV_2(1)*Data_voltage(i-1)+ piV_2(2)*Data_voltage(i-2) + piV_2(3)*yivs_2(i-1) + piV_2(4)*yivs_2(i-2)];
    yivs_3 = [yivs_3; piV_3(1)*Data_voltage(i-1)+ piV_3(2)*Data_voltage(i-2) + piV_3(3)*Data_voltage(i-3) + piV_3(4)*yivs_3(i-1) + piV_3(5)*yivs_3(i-2) + piV_3(6)*yivs_3(i-3)];
end 

figure(6)
hold on;
plot(yiv_1,'r-')
plot(yiv_2,'g-')
plot(yiv_3,'b-')
title('Predykcja jednokrokowa wraz z odpowiedzi¹ obiektu metoda IV')
plot(Data_results,'k')
legend('y_1','y_2','y_3','obiekt')
hold off


figure(7)
hold on;
plot(yivs_1,'r-')
plot(yivs_2,'g-')
plot(yivs_3,'b-')
title('OdpowiedŸ modelu wraz z odpowiedzi¹ obiektu metoda IV')
plot(Data_results,'k')
legend('ys_1','ys_2','ys_3','obiekt')
hold off

ynLS_d1 = fi_1*pLS_1;
ynLS_d2 = fi_2*pLS_2;
ynLS_d3 = fi_3*pLS_3;
ynIV_d1 = z_1*piV_1;
ynIV_d2 = z_2*piV_2;
ynIV_d3 = z_3*piV_3;

enLS_1 = y_1 - ynLS_d1;
enLS_2 = y_2 - ynLS_d2;
enLS_3 = y_3 - ynLS_d3;
enIV_1 = yiv_1 - ynIV_d1;
enIV_2 = yiv_2 - ynIV_d2;
enIV_3 = yiv_3 - ynIV_d3;

Vlsn_1 = (enLS_1'*enLS_1)/N_probek ;
Vlsn_2 = (enLS_2'*enLS_2)/N_probek;
Vlsn_3 = (enLS_3'*enLS_3)/N_probek;
Vivn_1 = (enIV_1'*enIV_1)/N_probek
Vivn_2 = (enIV_2'*enIV_2)/N_probek
Vivn_3 = (enIV_3'*enIV_3)/N_probek

AICLS_1 = N_probek * log(Vlsn_1) + 2*1
AICLS_2 = N_probek * log(Vlsn_2) + 2*2
AICLS_3 = N_probek * log(Vlsn_3) + 2*3
AICIV_1 = N_probek * log(Vivn_1) + 2*1
AICIV_2 = N_probek * log(Vivn_2) + 2*2
AICIV_3 = N_probek * log(Vivn_3) + 2*3

figure(8)
hold on;
plot(yivs_1,'r-')
plot(ys_1,'g-')
plot(Data_results,'k')
legend('yIV_1','yLS_1','obiekt')
hold off



