#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

//Propietats Fluid Intern (Water)
double pw(double T) {
    return 847.2+1.298*T-2.657e-3*T*T; 
}
double yw(double T) {
    return -0.722+7.168e-3*T-9.137e-6*T*T; 
}
double cpw(double T) {
    return 5648.79-9.14*T+14.21e-3*T*T; 
}
double uw(double T) {
    return exp(7.867-0.077*T+9.04e-5*T*T); 
}
double Bw(double T, double p) {
    return -(1.298-2.657e-3*T*2)/p; 
}

//Propietats Fluid Extern (Air)

double pa(double T, double P) {
    return P/(287*T);
}
double ya(double T) {
    return 2.728e-3+7.776e-5*T;
}
double cpa(double T) {
    return 1031.5-0.21*T+4.143e-4*T*T;
}
double ua(double T) {
    return (2.5393e-5*sqrt(T/273.15))/(1+(122/T));
}
double Ba(double T) {
    return 1/T;
}

int main () {
    const double pi=M_PI;

    //Dades Físiques
    const double Ri=0.02/2, Rext=0.024/2, L=20, y=100000, Text=-40+273, Pext=1e5, vin=0.1, Pin=2e5, Tin=95+273;
    
    //Dades Numèriques
    const int n=50; const double d=1e-7;
    
    //Càlculs Geometria
    const double dx=L/n, S=pi*(Rext*Rext-Ri*Ri), Ai=2*pi*Ri*dx, Aext=2*pi*Rext*dx, Aint=pi*Ri*Ri;
    
    //Definició Vectors
    const int size=n+4;
    vector<double> ap(size,0), aw(size,0), ae(size,0), bp(size,0), P(size,0), R(size,0);
    vector<double> af(size,0), an(size,0), T(size,Tin), Tint(size, 0), Pint(size, 0), pint(size, 0), vint(size, 0);

    //Càlculs Previs
    const double pin=pw(Tin), min=pin*vin*Aint;

    //Valors inicial
    Tint[2]=Tin, Pint[2]=Pin, pint[2]=pin, vint[2]=vin;

    double em;

    do {
        em=0;
        vector<double> Tb=T;

        //Càlcul fluid interior
        double difT=0, difP=0, difv=0;
        for (int i=2; i<=n+1; ++i){
            Tint[i+1]=Tint[i];
            Pint[i+1]=Pint[i];
            vint[i+1]=vint[i];
            do{
                double Tg=Tint[i+1], Pg=Pint[i+1], vg=vint[i+1];
                double Tm=(Tint[i+1]+Tint[i])/2, Pm=(Pint[i+1]+Pint[i])/2, vm=(vint[i+1]+vint[i])/2;

                double Re=pw(Tm)*vm*2*Ri/uw(Tm), Pr=uw(Tm)*cpw(Tm)/yw(Tm), Gz=Re*Pr*2*Ri/L;

                if (Re<2000 && Gz>10){
                    af[i]=(1.86*pow(Re*Pr*2*Ri/L,1.0/3.0)*pow(uw(Tm)/uw(T[i]),14.0/100.0))*yw(Tm)/(2*Ri);
                }
                else if (Re<2000 && Gz<10){
                    af[i]=(3.66)*yw(Tm)/(2*Ri);
                }
                else if (Re>2000 && Pr>0.6){
                    af[i]=(0.027*pow(Re,8.0/10.0)*pow(Pr,3.3/10.0)*pow(uw(Tm)/uw(T[i]),14.0/100.0))*yw(Tm)/(2*Ri);
                }
                else{
                    af[i]=(0.023*pow(Re,8.0/10.0)*pow(Pr,4.0/10.0))*yw(Tm)/(2*Ri);
                }

                Pint[i+1]=Pint[i]+(min*(vint[i]-vint[i+1])-0.005*pw(Tm)*vm*vm/2*pi*2*Ri*dx)/Aint;
                Tint[i+1]=(af[i]*(T[i]-Tm)*pi*2*Ri*dx+min*(vint[i]*vint[i]-vint[i+1]*vint[i+1])/2)/(min*cpw(Tm))+Tint[i];
                pint[i+1]=pw(Tint[i+1]);
                vint[i+1]=min/(Aint*pint[i+1]);

                difT=Tint[i+1]-Tg, difP=Pint[i+1]-Pg, difv=vint[i+1]-vg;
                if (difT<0) difT=-difT; if (difP<0) difP=-difP; if (difv<0) difv=-difv;
        } while(difT>d || difP>d || difv>d);
        };

        //Càlcul fluid exterior
        for (int i=2; i<=n+1; ++i){
            double Tm=(Text+T[i])/2;
            double Ra=(ua(Tm)*cpa(Tm)/ya(Tm))*(9.81*Ba(Tm)*pa(Tm,Pext)*pa(Tm,Pext)*(T[i]-Text)*8*Rext*Rext*Rext/(ua(Tm)*ua(Tm)));
            if (Ra>=1e9){
                an[i]=(0.1*pow(Ra,1.0/3.0))*ya(Tm)/(2*Rext);
            }
            else{
                an[i]=(0.47*pow(Ra,1.0/4.0))*ya(Tm)/(2*Rext);
            }
        };


        //Càlcul Tub
        for (int i=1; i<=n+2; ++i){
            if(i==1){
                ae[i]=1;
                ap[i]=1;
            }
            else if(i==n+2){
                aw[i]=1;
                ap[i]=1;
            }
            else{
                double xvc=(i-1)*dx;
                double xvcp=(i-2)*dx;
                double xp=(xvc+xvcp)/2;
                if(i == n+1) {
                    ae[i] = y*S/(dx/2);
                    aw[i] = y*S/dx; 
                }
                else if(i == 2) {
                    aw[i] = y*S/(dx/2);
                    ae[i] = y*S/dx; 
                }
                else {
                    aw[i] = y*S/dx;
                    ae[i] = y*S/dx;
                }
            bp[i]=af[i]*((Tint[i]+Tint[i+1])/2*Ai)+an[i]*(Text*Aext);
            ap[i]=ae[i]+aw[i]+af[i]*Ai+an[i]*Aext;
            }
            P[i]=ae[i]/(ap[i]-aw[i]*P[i-1]);
            R[i]=(bp[i]+aw[i]*R[i-1])/(ap[i]-aw[i]*P[i-1]);
        };

        for (int i=n+2; i>=1; i=i-1){
            T[i]=P[i]*T[i+1]+R[i];

            double dif=T[i]-Tb[i];
            if (dif<0) dif=-dif;
            if (dif>em) em=dif;
        };

    } while (em>d);

    cout<<"Las Temperaturas son: \n";
    for (int i=1; i<=n+2; ++i){
        cout<<fixed<<setprecision(3)<<T[i]<<", ";
    }
    cout<<"\n";
    cout<<"Temperatura sortida: "<<Tint[n+2];
    cout<<" Pressió sortida: "<<Pint[n+2];
    cout<<" Velocitat sortida: "<<vint[n+2];
    cout<<"\n";

    double Qconva=0;
    for (int i=2; i<=n+1; i++){
        Qconva=Qconva+an[i]*(T[i]-Text)*pi*2*Rext*dx;
    }
    double Qconvw=0;
    for (int i=2; i<=n+1; i++){
        Qconvw=Qconvw+cpw((Tint[i]+Tint[i+1])/2)*(Tint[i+1]-Tint[i])*min;
    }


    cout<<"Qsortida: "<<Qconva;
    cout<<"\n";
    cout<<"Qentrada: "<<Qconvw;
    cout<<"\n";

    return 0;
}
