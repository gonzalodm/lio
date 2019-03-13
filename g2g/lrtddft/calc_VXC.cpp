#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

void calc_VXC(double* dens, double* diff,
              double* vrho,double* vsigma,
              double* v2rho2,double* v2rhosigma,double *v2sigma2,
              double* v3rho3,double* v3rho2sigma,double* v3rhosigma2,double* v3sigma3,
              double* DDUM,double* PDUM)
{

// DIFFERENCE RELAXED DENSITY
   double C[6];
   double DUMNV[2], DXV[2], DYV[2],DZV[2];
   DUMNV[0]=DUMNV[1]=diff[0];
   DXV[0]=DXV[1]=diff[1];
   DYV[0]=DYV[1]=diff[2];
   DZV[0]=DZV[1]=diff[3];
//     -- DUMGRV
//     1:A*A 2:B*B 3:A*B 4:B*A
   double DUMGRV[4], DUMXX[4];
   DUMGRV[0]=DXV[0]*dens[1]+DYV[0]*dens[2]+DZV[0]*dens[3];
   DUMGRV[1]=DXV[1]*dens[1]+DYV[1]*dens[2]+DZV[1]*dens[3];
   DUMGRV[2]=DXV[0]*dens[1]+DYV[0]*dens[2]+DZV[0]*dens[3];
   DUMGRV[3]=DXV[1]*dens[1]+DYV[1]*dens[2]+DZV[1]*dens[3];

   DUMXX[0]=DXV[0]*DXV[0]+DYV[0]*DYV[0]+DZV[0]*DZV[0];
   DUMXX[1]=DXV[1]*DXV[1]+DYV[1]*DYV[1]+DZV[1]*DZV[1];
   DUMXX[2]=DXV[0]*DXV[1]+DYV[0]*DYV[1]+DZV[0]*DZV[1];
   DUMXX[3]=DUMXX[2];

/*
   DUMGRV y DUMXX todo bien
   cout << DUMGRV[0] << " " << DUMGRV[1] << " " << DUMGRV[2] << " ";
   cout << DUMGRV[3] << endl;
   exit(-1);
*/

// GROUND STATE
   double GDUMA,GDUMAG1,GDUMAG2;
   //GDUMA=EX(IIPT,KRA)+EC(IIPT,IRA)
   GDUMA=vrho[0]+vrho[2];
   //GDUMAG1=TWO*EX(IIPT,KGA)+TWO*EC(IIPT,IGA)
   GDUMAG1=2.0f*(vsigma[0]+vsigma[2]);
   //GDUMAG2=EC(IIPT,IGC)
   GDUMAG2=vsigma[4];

   //double GDUMB,GDUMBG1,GDUMBG2;
   //double GDUMBG1,GDUMBG2; no use
   //GDUMB=EX(IIPT,KRB)+EC(IIPT,IRB)
   //GDUMB=vrho[1]+vrho[3];
   //GDUMBG1=TWO*EX(IIPT,KGB)+TWO*EC(IIPT,IGB)
   //GDUMBG1=2.0f*(vsigma[1]+vsigma[3]); no use
   //GDUMBG2=EC(IIPT,IGC)
   //GDUMBG2=vsigma[4]; no uso

// CONTRACTION
   double GDUMAX,GDUMAY,GDUMAZ;
   //double GDUMBX,GDUMBY,GDUMBZ;
   GDUMAX=GDUMAG1*dens[1]+GDUMAG2*dens[1];
   GDUMAY=GDUMAG1*dens[2]+GDUMAG2*dens[2];
   GDUMAZ=GDUMAG1*dens[3]+GDUMAG2*dens[3];
   //GDUMBX=GDUMBG1*dens[1]+GDUMBG2*dens[1];
   //GDUMBY=GDUMBG1*dens[2]+GDUMBG2*dens[2];
   //GDUMBZ=GDUMBG1*dens[3]+GDUMBG2*dens[3];

/*
   cout << "GDUMAX " << GDUMAX << endl;
   cout << "GDUMAY " << GDUMAY << endl;
   cout << "GDUMAZ " << GDUMAZ << endl;
   cout << "GDUMBX " << GDUMBX << endl;
   cout << "GDUMBY " << GDUMBY << endl;
   cout << "GDUMBZ " << GDUMBZ << endl;
   exit(-1);
   // todo bien hasta aqui
*/

// V NON CORE CONTRIBUTION
   double DUMNV1; //, DUMNV2;no uso
   //double DUMGRV1,DUMGRV2,DUMGRV3,DUMGRV4; no uso DUMGRV2
   double DUMGRV1,DUMGRV3,DUMGRV4;
   //DUMNV1=EX(IIPT,KRA)+EC(IIPT,IRA)
   DUMNV1=vrho[0]+vrho[2];
   //DUMGRV1=TWO*(EX(IIPT,KGA)+EC(IIPT,IGA))
   DUMGRV1=2.0f*(vsigma[0]+vsigma[2]);
   //DUMGRV3=EC(IIPT,IGC)
   DUMGRV3=vsigma[4];
   //DUMNV2=EX(IIPT,KRB)+EC(IIPT,IRB)
   //DUMNV2=vrho[1]+vrho[3]; no uso
   //DUMGRV2=TWO*(EX(IIPT,KGB)+EC(IIPT,IGB))
   //DUMGRV2=2.0f*(vsigma[1]+vsigma[3]); no uso
   //DUMGRV4=EC(IIPT,IGC)
   DUMGRV4=vsigma[4];

   double VNCDOMA,VNCDOMAX,VNCDOMAY,VNCDOMAZ;
   double VNCDUMAX,VNCDUMAY,VNCDUMAZ;
   VNCDOMA=DUMNV1;
   VNCDOMAX=DUMGRV1*dens[1]+DUMGRV3*dens[1];
   VNCDOMAY=DUMGRV1*dens[2]+DUMGRV3*dens[2];
   VNCDOMAZ=DUMGRV1*dens[3]+DUMGRV3*dens[3];
   VNCDUMAX=DUMGRV1*DXV[0]+DUMGRV4*DXV[1];
   VNCDUMAY=DUMGRV1*DYV[0]+DUMGRV4*DYV[1];
   VNCDUMAZ=DUMGRV1*DZV[0]+DUMGRV4*DZV[1];

/*
   cout << "VNCDOMA " << VNCDOMA << endl;
   cout << "VNCDOMAX " << VNCDOMAX << endl;
   cout << "VNCDOMAY " << VNCDOMAY << endl;
   cout << "VNCDOMAZ " << VNCDOMAZ << endl;
   cout << "VNCDUMAX " << VNCDUMAX << endl;
   cout << "VNCDUMAY " << VNCDUMAY << endl;
   cout << "VNCDUMAZ " << VNCDUMAZ << endl;
   exit(-1);
   // todo bien
*/

/*
   no las uso nunca
   double VNCDOMB,VNCDOMBX,VNCDOMBY,VNCDOMBZ;
   double VNCDUMBX,VNCDUMBY,VNCDUMBZ;
   VNCDOMB=DUMNV2;
   VNCDOMBX=DUMGRV2*dens[1]+DUMGRV4*dens[1];
   VNCDOMBY=DUMGRV2*dens[2]+DUMGRV4*dens[2];
   VNCDOMBZ=DUMGRV2*dens[3]+DUMGRV4*dens[3];
   VNCDUMBX=DUMGRV2*DXV[1]+DUMGRV3*DXV[0];
   VNCDUMBY=DUMGRV2*DYV[1]+DUMGRV3*DYV[0];
   VNCDUMBZ=DUMGRV2*DZV[1]+DUMGRV3*DZV[0];
*/

/*
   cout << "VNCDOMB " << VNCDOMB << endl;
   cout << "VNCDOMBX " << VNCDOMBX << endl;
   cout << "VNCDOMBY " << VNCDOMBY << endl;
   cout << "VNCDOMBZ " << VNCDOMBZ << endl;
   cout << "VNCDUMBX " << VNCDUMBX << endl;
   cout << "VNCDUMBY " << VNCDUMBY << endl;
   cout << "VNCDUMBZ " << VNCDUMBZ << endl;
   exit(-1);
   // todo bien
*/
// END V NON CORE

// V CORE CONTRIBUTION
   C[0]=DUMNV[0];
   C[1]=2.0f*DUMGRV[0];
   C[2]=DUMGRV[2];
   C[3]=DUMNV[1];
   C[4]=2.0f*DUMGRV[1];
   C[5]=DUMGRV[3];
   double DUMA,DUMAG1,DUMAG2;
// IRA C1
   //DUMA=C1*(EX(IIPT,KRARA)+EC(IIPT,IRARA))
   DUMA=C[0]*(v2rho2[0]+v2rho2[2]);
   //DUMAG1=C1*TWO*(EX(IIPT,KRAGA)+EC(IIPT,IRAGA))
   DUMAG1=2.0f*C[0]*(v2rhosigma[0]+v2rhosigma[2]);
   //DUMAG2=C1*EC(IIPT,IRAGC)
   DUMAG2=C[0]*v2rhosigma[4];
// IGA C2
   //DUMA=DUMA+C2*(EX(IIPT,KRAGA)+EC(IIPT,IRAGA))
   DUMA=DUMA+C[1]*(v2rhosigma[0]+v2rhosigma[2]);
   //DUMAG1=DUMAG1+C2*TWO*(EX(IIPT,KGAGA)+EC(IIPT,IGAGA))
   DUMAG1=DUMAG1+2.0f*C[1]*(v2sigma2[0]+v2sigma2[2]);
   //DUMAG2=DUMAG2+C2*EC(IIPT,IGAGC)
   DUMAG2=DUMAG2+C[1]*v2sigma2[4];
// IGC C3
   //DUMA=DUMA+C3*EC(IIPT,IRAGC)
   DUMA=DUMA+C[2]*v2rhosigma[4];
   //DUMAG1=DUMAG1+C3*TWO*EC(IIPT,IGAGC)
   DUMAG1=DUMAG1+2.0f*C[2]*v2sigma2[4];
   //DUMAG2=DUMAG2+C3*EC(IIPT,IGCGC)
   DUMAG2=DUMAG2+C[2]*v2sigma2[7];
// IRB C4
   //DUMA=DUMA+C4*EC(IIPT,IRARB)
   DUMA=DUMA+C[3]*v2rho2[3];
   //DUMAG1=DUMAG1+C4*TWO*EC(IIPT,IRBGA)
   DUMAG1=DUMAG1+2.0f*C[3]*v2rhosigma[5];
   //DUMAG2=DUMAG2+C4*EC(IIPT,IRBGC)
   DUMAG2=DUMAG2+C[3]*v2rhosigma[7];
// IGB C5
   //DUMA=DUMA+C5*EC(IIPT,IRAGB)
   DUMA=DUMA+C[4]*v2rhosigma[3];
   //DUMAG1=DUMAG1+C5*TWO*EC(IIPT,IGAGB)
   DUMAG1=DUMAG1+2.0f*C[4]*v2sigma2[3];
   //DUMAG2=DUMAG2+C5*EC(IIPT,IGBGC)
   DUMAG2=DUMAG2+C[4]*v2sigma2[6];
// IGC C6
   //DUMA=DUMA+C6*EC(IIPT,IRAGC)
   DUMA=DUMA+C[5]*v2rhosigma[4];
   //DUMAG1=DUMAG1+C6*TWO*EC(IIPT,IGAGC)
   DUMAG1=DUMAG1+2.0f*C[5]*v2sigma2[4];
   //DUMAG2=DUMAG2+C6*EC(IIPT,IGCGC)
   DUMAG2=DUMAG2+C[5]*v2sigma2[7];

/*
   cout << "DUMA " << DUMA << endl;
   cout << "DUMAG1 " << DUMAG1 << endl;
   cout << "DUMAG2 " << DUMAG2 << endl;
   exit(-1);
   // todo bien
*/
   double DUMB,DUMBG1,DUMBG2;
// IRA C1
   //DUMB=C1*EC(IIPT,IRARB)
   DUMB=C[0]*v2rho2[3];
   //DUMBG1=C1*TWO*EC(IIPT,IRAGB)
   DUMBG1=2.0f*C[0]*v2rhosigma[3];
   //DUMBG2=C1*EC(IIPT,IRAGC)
   DUMBG2=C[0]*v2rhosigma[4];
// IGA C2
   //DUMB=DUMB+C2*EC(IIPT,IRBGA)
   DUMB=DUMB+C[1]*v2rhosigma[5];
   //DUMBG1=DUMBG1+C2*TWO*EC(IIPT,IGAGB)
   DUMBG1=DUMBG1+2.0f*C[1]*v2sigma2[3];
   //DUMBG2=DUMBG2+C2*EC(IIPT,IGAGC)
   DUMBG2=DUMBG2+C[1]*v2sigma2[4];
// IGC C3
   //DUMB=DUMB+C3*EC(IIPT,IRBGC)
   DUMB=DUMB+C[2]*v2rhosigma[7];
   //DUMBG1=DUMBG1+C3*TWO*EC(IIPT,IGBGC)
   DUMBG1=DUMBG1+2.0f*C[2]*v2sigma2[6];
   //DUMBG2=DUMBG2+C3*EC(IIPT,IGCGC)
   DUMBG2=DUMBG2+C[2]*v2sigma2[7];
// IRB C4
   //DUMB=DUMB+C4*(EX(IIPT,KRBRB)+EC(IIPT,IRBRB))
   DUMB=DUMB+C[3]*(v2rho2[1]+v2rho2[4]);
   //DUMBG1=DUMBG1+C4*TWO*(EX(IIPT,KRBGB)+EC(IIPT,IRBGB))
   DUMBG1=DUMBG1+2.0f*C[3]*(v2rhosigma[1]+v2rhosigma[6]);
   //DUMBG2=DUMBG2+C4*EC(IIPT,IRBGC)
   DUMBG2=DUMBG2+C[3]*v2rhosigma[7];
// IGB C5
   //DUMB=DUMB+C5*(EX(IIPT,KRBGB)+EC(IIPT,IRBGB))
   DUMB=DUMB+C[4]*(v2rhosigma[1]+v2rhosigma[6]);
   //DUMBG1=DUMBG1+C5*TWO*(EX(IIPT,KGBGB)+EC(IIPT,IGBGB))
   DUMBG1=DUMBG1+2.0f*C[4]*(v2sigma2[1]+v2sigma2[5]);
   //DUMBG2=DUMBG2+C5*EC(IIPT,IGBGC)
   DUMBG2=DUMBG2+C[4]*v2sigma2[6];
// IGC C6
   //DUMB=DUMB+C6*EC(IIPT,IRBGC)
   DUMB=DUMB+C[5]*v2rhosigma[7];
   //DUMBG1=DUMBG1+C6*TWO*EC(IIPT,IGBGC)
   DUMBG1=DUMBG1+2.0f*C[5]*v2sigma2[6];
   //DUMBG2=DUMBG2+C6*EC(IIPT,IGCGC)
   DUMBG2=DUMBG2+C[5]*v2sigma2[7];

/*
   cout << "DUMB " << DUMB << endl;
   cout << "DUMBG1 " << DUMBG1 << endl;
   cout << "DUMBG2 " << DUMBG2 << endl;
   exit(-1);
   // todo bien
*/

// CONTRACTION OF VC
   double VCDUMA,VCDUMAX,VCDUMAY,VCDUMAZ;
   //double VCDUMB,VCDUMBX,VCDUMBY,VCDUMBZ; no lo uso
   VCDUMA=DUMA;
   VCDUMAX=DUMAG1*dens[1]+DUMAG2*dens[1];
   VCDUMAY=DUMAG1*dens[2]+DUMAG2*dens[2];
   VCDUMAZ=DUMAG1*dens[3]+DUMAG2*dens[3];
/*
   no uso 
   VCDUMB=DUMB;
   VCDUMBX=DUMBG1*dens[1]+DUMBG2*dens[1];
   VCDUMBY=DUMBG1*dens[2]+DUMBG2*dens[2];
   VCDUMBZ=DUMBG1*dens[3]+DUMBG2*dens[3];
*/

/*
   cout << "VCDUMA  " << VCDUMA << endl;
   cout << "VCDUMAX " << VCDUMAX << endl;
   cout << "VCDUMAY " << VCDUMAY << endl;
   cout << "VCDUMAZ " << VCDUMAZ << endl;
   cout << "VCDUMB  " << VCDUMB << endl;
   cout << "VCDUMBX " << VCDUMBX << endl;
   cout << "VCDUMBY " << VCDUMBY << endl;
   cout << "VCDUMBZ " << VCDUMBZ << endl;
   exit(-1);
   // todo bien
*/
// END V CORE

   DDUM[0]=GDUMA+VCDUMA;
   DDUM[1]=GDUMAX+VNCDUMAX+VCDUMAX;
   DDUM[2]=GDUMAY+VNCDUMAY+VCDUMAY;
   DDUM[3]=GDUMAZ+VNCDUMAZ+VCDUMAZ;

// PA
   PDUM[0]=VNCDOMA;
   PDUM[1]=VNCDOMAX;
   PDUM[2]=VNCDOMAY;
   PDUM[3]=VNCDOMAZ;
  
/*
   cout << DDUM[0] << endl;
   cout << DDUM[1] << endl;
   cout << DDUM[2] << endl;
   cout << DDUM[3] << endl;
   cout << PDUM[0] << endl;
   cout << PDUM[1] << endl;
   cout << PDUM[2] << endl;
   cout << PDUM[3] << endl;
   exit(-1);
   // todo bien
*/

// FOR GRDFUN
// GROUND STATE
   double DUMEXC = 0.0f;
   DUMEXC = -0.24536672503123141;
   //DUMEXC=EX0(IIPT)+EC0(IIPT)
   // no me coinciden las energias de correlacion de games con libxc
   // games = libxc / 2.4830838553290406 no se que es para el punto 1
   // al menos, los otros no probe

// FORM VXC
   double DUM1,DUM2,DUM3,DUM4,DUM5,DUM6;
   double DUMVA,DUMVB,DUMV;
   //DUM1=EX(IIPT,KRA)+EC(IIPT,IRA)
   DUM1=vrho[0]+vrho[2];
   //DUM2=TWO*(EX(IIPT,KGA)+EC(IIPT,IGA))
   DUM2=2.0f*(vsigma[0]+vsigma[2]);
   //DUM3=EC(IIPT,IGC)
   DUM3=vsigma[4];
   //DUM4=EX(IIPT,KRB)+EC(IIPT,IRB)
   DUM4=vrho[1]+vrho[3];
   //DUM5=TWO*(EX(IIPT,KGB)+EC(IIPT,IGB))
   DUM5=2.0f*(vsigma[1]+vsigma[3]);
   //DUM6=EC(IIPT,IGC)
   DUM6=vsigma[4];
   DUMVA=DUM1*DUMNV[0]+DUM2*DUMGRV[0]+DUM3*DUMGRV[2];
   DUMVB=DUM4*DUMNV[1]+DUM5*DUMGRV[1]+DUM6*DUMGRV[3];
   DUMV=DUMVA+DUMVB;
   double GRDFUN = 0.0f;
   GRDFUN=GRDFUN+DUMEXC+DUMV;
/*
   cout << GRDFUN << endl;
   exit(-1);
*/

// ASI ESTABA ANTES Y LO CORRIGIERON
//     DUMVA=DUM1*DUMNV(1)+DUM2*DUMGRV(1)
//     DUMVB=DUM4*DUMNV(2)+DUM2*DUMGRV(1)
//     DUMV=DUMVA+DUMVB
//     GRDFUN(IPT)=GRDFUN(IPT)+DUMEXC+DUMV
}
