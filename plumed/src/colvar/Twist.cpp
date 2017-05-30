/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Pbc.h"
#include "../tools/Matrix.h"
#include <string>
#include <cmath>
#include <math.h> 
#include <vector>
#define PI 3.14159265

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR TWIST
/*
 *
 * DNA twisting code works analogously to any collective variable (colvar) implemented in PLUMED. TWIST colvar monitors or controls the value of total twist between any chosen base pair levels i and j. As in imput to TWIST colvar, provide 48 atom numbers from bases that are restrained: (3×4 bases ×3 atoms: C1',N1/N9 and C6/C8 depending whether purine or pyrimidine), plus 6×2 auxiliary atoms, force constant, desired value of total twist,and how many turns the restricted DNA fragment has. The energy penalty will be added to the potential energy functional: E_tw=0.5*k*(tw0-tw)^2.
 *
 * The numbers that represent DNA bases (desriptors) are taken from bases i-1, i, i+1 and j-1, j, j+1 from the Watson DNA strand (5'->3'), and from the Crick DNA strand (3'->5'). The auxiliary atoms are used as a support for a correct implementation of the periodic boundary condition, and have to be selected to be evenly destributed along the restricted DNA fragment.
 *
 * Example of plumed.dat input file:
 *
 * tw: TWIST ATOMS=44,42,45,74,72,75,106,104,107,167,169,261,263,358,360,423,421,424,455,453,456,487,485,488,579,577,580,612,610,613,644,642,645,704,706,801,803,896,898,960,958,961,993,991,994,1025,1023,1026 N_TURNS=1.0
 * tw_r: RESTRAINT ARG=tw KAPPA=0.25 AT=358.8
 * PRINT STRIDE=250 ARG=tw,tw_r.bias FILE=Twist_358.8
 * 
 * First 9 atom numbers represent bases i-1, i, i+1 from the Watson strand, followed by 6 auxiliary atoms, follwed by 9 atoms from bases j-1, j, j+1 from the Watson strand, follwed by 9 atoms from bases j-1, j, j+1 from the Crick strand, followed by 6 auxiliary atoms, and finally the last 9 atoms represent bases i-1, i, i+1 from the Crick strand.
 *
 */
//+ENDPLUMEDOC
   
class Twist : public Colvar {
  bool pbc;
  std::vector<double> n_turns;
  Vector rotVector6Atoms(int arr[6]);
  Vector rotVectorBasis(Matrix<double> R1,  Matrix<double> R2);
  double rotAngleBasis(Matrix<double> P1,  Matrix<double> P2);
  double traceBasis( Matrix<double> P1, Matrix<double> P2);
  Matrix<double> basis(Vector a, Vector b, Vector c);
  Matrix<double> basis2Vectors(Vector e1, Vector diff);
  Matrix<double> basis2Vectors_2(Vector e1, Vector diff);
  Matrix<double> dCrosProduct(Matrix<double> da, Vector a, Matrix<double> db, Vector b);
  Matrix<double> dUnit(Vector a, Matrix<double> da);

public:
  explicit Twist(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Twist,"TWIST")

void Twist::registerKeywords(Keywords& keys){
  Colvar::registerKeywords(keys);
  keys.add("atoms","ATOMS","the keyword with which you specify what atoms to use");
  keys.add("compulsory","N_TURNS","0.0","specifies number of DNA turns between restrained base pairs");

}

Twist::Twist(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
n_turns(1,0.0)
{
  parseVector("N_TURNS",n_turns);
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=48)
  {
	  error("Number of specified atoms should be 48");
	  log.printf("number or atoms %d \n ",atoms.size());
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  checkRead();
}

// calculator
void Twist::calculate(){

  if(pbc) makeWhole();

  const double N=n_turns[0];

// defining rotation vectors for each of the nucleotides pair
  int atomNmbs[10][6]=
  {
		  { 0, 1, 2, 3, 4, 5}, //Basis 1 - strand 1 begins
		  { 3, 4, 5, 6, 7, 8}, //Basis 2
		  { 9,10,11,12,13,14}, //Middle atoms strand #1
		  {15,16,17,18,19,20}, //Basis 3
		  {18,19,20,21,22,23}, //Basis 4
		  {24,25,26,27,28,29}, //Basis 5 - strand 2 begins
		  {27,28,29,30,31,32}, //Basis 6
		  {33,34,35,36,37,38}, //Middle atoms strand #2
		  {39,40,41,42,43,44}, //Basis 7
		  {42,43,44,45,46,47}  //Basis 8
  };

  double atomDerivatives[48][3];
  for (int i=0;i<48;i++){
	  for(int j=0; j<3;j++){
		  atomDerivatives[i][j] = 0.0;
	  }
  }

  int Direction[10]={1,1,0, 1, 1,-1,-1,0,-1,-1};
  int Top_Btm[10]  ={1,1,0,-1,-1,-1,-1,0, 1, 1};

  double U_Derivatives[10][4][3][3];
  for(int i=0;i<10;i++){
	for(int j=0;j<4;j++){
		for(int k=0;k<3;k++){
			for(int l=0;l<3;l++){
				U_Derivatives[i][j][l][l] = 0.0 ;}}}} 
  

  Vector U; Vector U_avTop; Vector U_avBtm;

  Matrix<double> dU(3,3);
  Matrix<double> dU_final(3,3);
  Matrix<double> dU_avTop(3,3);
  Matrix<double> dU_avBtm(3,3);
  Vector dQ; Vector dTrace; Vector dTwist; Vector drotAng;
  Matrix<double> dScrew(3,3);
  Matrix<double> dU_final_norm(3,3);

for(int n_basis=0;n_basis<10;n_basis++){
if (Direction[n_basis]==1 || Direction[n_basis]==-1) {
  U = rotVector6Atoms(atomNmbs[n_basis]);

  // calculating dU

  Matrix<double> Basis1 = basis(getPosition(atomNmbs[n_basis][0]),getPosition(atomNmbs[n_basis][1]),getPosition(atomNmbs[n_basis][2]));
  Matrix<double> Basis2 = basis(getPosition(atomNmbs[n_basis][3]),getPosition(atomNmbs[n_basis][4]),getPosition(atomNmbs[n_basis][5]));

  Vector a1 = delta(getPosition(atomNmbs[n_basis][0]),getPosition(atomNmbs[n_basis][1]));
  Vector b1 = delta(getPosition(atomNmbs[n_basis][0]),getPosition(atomNmbs[n_basis][2]));
  Vector e2_1 = crossProduct(a1,b1);
  Vector e3_1 = crossProduct(a1,e2_1);

  Vector a2 = delta(getPosition(atomNmbs[n_basis][3]),getPosition(atomNmbs[n_basis][4]));
  Vector b2 = delta(getPosition(atomNmbs[n_basis][3]),getPosition(atomNmbs[n_basis][5]));
  Vector e2_2 = crossProduct(a2,b2);
  Vector e3_2 = crossProduct(a2,e2_2);

  Matrix<double> da1(3,3);
  Matrix<double> da2(3,3);
  Matrix<double> db1(3,3);
  Matrix<double> db2(3,3);
  da1 = 0.0; da2 = 0.0; db1 = 0.0; db2 = 0.0;

  Matrix<double> dR11(3,3); Matrix<double> dR12(3,3); Matrix<double> dR13(3,3);
  Matrix<double> dR21(3,3); Matrix<double> dR22(3,3); Matrix<double> dR23(3,3);

  Matrix<double> de2_1(3,3); Matrix<double> de3_1(3,3);
  Matrix<double> de2_2(3,3); Matrix<double> de3_2(3,3);

  for(int k=0;k<4;k++){
	  if(k==0)		  { da1 = 0.0; da2 = 0.0; db1 = 0.0; db2 = 0.0; da1(0,0) = 1.0; da1(1,1) = 1.0; da1(2,2) = 1.0;}
	  else if(k==1)	{ da1 = 0.0; da2 = 0.0; db1 = 0.0; db2 = 0.0; db1(0,0) = 1.0; db1(1,1) = 1.0; db1(2,2) = 1.0;}
	  else if(k==2)	{ da1 = 0.0; da2 = 0.0; db1 = 0.0; db2 = 0.0; da2(0,0) = 1.0; da2(1,1) = 1.0; da2(2,2) = 1.0;}
	  else if(k==3)	{ da1 = 0.0; da2 = 0.0; db1 = 0.0; db2 = 0.0; db2(0,0) = 1.0; db2(1,1) = 1.0; db2(2,2) = 1.0;}

	  de2_1 = dCrosProduct(da1,a1,db1,b1);
	  de3_1 = dCrosProduct(da1,a1,de2_1,e2_1);
	  dR11  = dUnit(a1,da1);
	  dR12  = dUnit(e2_1,de2_1);
	  dR13  = dUnit(e3_1,de3_1);

	  de2_2 = dCrosProduct(da2,a2,db2,b2);
	  de3_2 = dCrosProduct(da2,a2,de2_2,e2_2);
	  dR21  = dUnit(a2,da2);
	  dR22  = dUnit(e2_2,de2_2);
	  dR23  = dUnit(e3_2,de3_2);

	  for(int i=0;i<3;i++) {dQ(i) = 0.0; }
	  for(int i=0;i<3;i++)
	  {
		  dQ(i) += dR11(1,i)*Basis2(0,2)+Basis1(0,1)*dR21(2,i)
		  	  -dR11(2,i)*Basis2(0,1)-Basis1(0,2)*dR21(1,i)
			  +dR12(1,i)*Basis2(1,2)+Basis1(1,1)*dR22(2,i)
		  	  -dR12(2,i)*Basis2(1,1)-Basis1(1,2)*dR22(1,i)
			  +dR13(1,i)*Basis2(2,2)+Basis1(2,1)*dR23(2,i)
			  -dR13(2,i)*Basis2(2,1)-Basis1(2,2)*dR23(1,i);
	  }
	  for(int i=0;i<3;i++)
	  {
		 dU(0,i)=dQ(i);
	  }
  
	  for(int i=0;i<3;i++) {dQ(i) = 0.0; }
	  for(int i=0;i<3;i++)
	  {
		  dQ(i) += dR11(2,i)*Basis2(0,0)+Basis1(0,2)*dR21(0,i)
	  		  -dR11(0,i)*Basis2(0,2)-Basis1(0,0)*dR21(2,i)
	  		  +dR12(2,i)*Basis2(1,0)+Basis1(1,2)*dR22(0,i)
	  		  -dR12(0,i)*Basis2(1,2)-Basis1(1,0)*dR22(2,i)
	  		  +dR13(2,i)*Basis2(2,0)+Basis1(2,2)*dR23(0,i)
	  		  -dR13(0,i)*Basis2(2,2)-Basis1(2,0)*dR23(2,i);
	  }
	  for(int i=0;i<3;i++)
	  {
		  dU(1,i)=dQ(i);
	  }
 
	  for(int i=0;i<3;i++) {dQ(i) = 0.0; }
	  for(int i=0;i<3;i++)
	  {
		  dQ(i) += dR11(0,i)*Basis2(0,1)+Basis1(0,0)*dR21(1,i)
		  	  -dR11(1,i)*Basis2(0,0)-Basis1(0,1)*dR21(0,i)
		  	  +dR12(0,i)*Basis2(1,1)+Basis1(1,0)*dR22(1,i)
		  	  -dR12(1,i)*Basis2(1,0)-Basis1(1,1)*dR22(0,i)
		  	  +dR13(0,i)*Basis2(2,1)+Basis1(2,0)*dR23(1,i)
		  	  -dR13(1,i)*Basis2(2,0)-Basis1(2,1)*dR23(0,i);
	  }
	  for(int i=0;i<3;i++)
	  {
		  dU(2,i)=dQ(i);
 	  }

	  Matrix<double> dU_norm = dUnit(U,dU);

	  for(int i=0;i<3;i++) {
		  for(int j=0;j<3;j++) {
			  U_Derivatives[n_basis][k][i][j] = Direction[n_basis]*dU_norm(i,j); } }

  } // finish loop k
  // finish calculating dU for each n_basis

  U *= 1.0/U.modulo();	

  // calculating averaged rotation vectors for top and bottom base pairs
  if     (Top_Btm[n_basis] == 1){
	  U_avTop += Direction[n_basis]*U;}
  else if(Top_Btm[n_basis] == -1){
	  U_avBtm += Direction[n_basis]*U;}

  } // finish if loop
  } // finish loop n_basis

  Vector diffTop = delta(getPosition(3),getPosition(42));
  Vector diffBtm = delta(getPosition(18),getPosition(27));

  Matrix<double> Top = basis2Vectors_2(U_avTop,diffTop);
  Matrix<double> Btm = basis2Vectors_2(U_avBtm,diffBtm);

  // calculating resulting rotation vector, and twist angle
  Vector U_final = rotVectorBasis(Top,Btm);
  U_final *= 1.0/U_final.modulo();
  double rotAng  = rotAngleBasis(Top, Btm)*180/PI; //+ 360*N;
  Vector Screw   = U_final*rotAng; //Euler's axis
  double Trace = (traceBasis(Top, Btm) - 1)*0.5;
  if (Trace*Trace > 1.0) fprintf(stderr,"trace is too big %g \n", Trace);
  double twst = (Screw(0)*Top(2,0) + Screw(1)*Top(2,1) + Screw(2)*Top(2,2))+ 360*N;

  // Matrix describing derivatives of which atoms are not 0
  // for 1 basis, you have 6 atoms, which provide 4 distances (k=1,2,3,4)
  // for the calculation of R1 and R2 (basis associated with restricted bases)

  Matrix<int> dCoord(6,4);
  dCoord = 0;
  dCoord(0,0)=-1; dCoord(0,1)=-1; dCoord(1,0)=1; dCoord(2,1)=1;
  dCoord(3,2)=-1; dCoord(3,3)=-1; dCoord(4,2)=1; dCoord(5,3)=1;

  // calculating dU_final derivatives

  Vector e1_1 = crossProduct(U_avTop,diffTop);
  Vector e2_1 = crossProduct(U_avTop,e1_1);

  Vector e1_2 = crossProduct(U_avBtm,diffBtm);
  Vector e2_2 = crossProduct(U_avBtm,e1_2);

  for(int n_basis=0; n_basis<10; n_basis++){
  if (Direction[n_basis]==1 || Direction[n_basis]==-1) {

	  Matrix<double> dR11(3,3);
	  Matrix<double> dR12(3,3);
	  Matrix<double> dR13(3,3);

	  Matrix<double> dR21(3,3);
	  Matrix<double> dR22(3,3);
	  Matrix<double> dR23(3,3);

	  Matrix<double> de1_1(3,3);
	  Matrix<double> de1_2(3,3);
	  Matrix<double> de2_1(3,3);
	  Matrix<double> de2_2(3,3);

	  for(int k=0;k<4;k++){

		  if (Top_Btm[n_basis] == 1){
			  for(int i=0; i<3; i++){
				  for(int j=0; j<3; j++){
					  dU_avTop(i,j) = U_Derivatives[n_basis][k][i][j];
				  }
			  }

			  for(int i=0;i<3;i++){
				  de1_1(0,i) =  dU_avTop(1,i)*diffTop(2)-dU_avTop(2,i)*diffTop(1);
				  de1_1(1,i) = -dU_avTop(0,i)*diffTop(2)+dU_avTop(2,i)*diffTop(0);
				  de1_1(2,i) =  dU_avTop(0,i)*diffTop(1)-dU_avTop(1,i)*diffTop(0);
			  }
			  de2_1 = dCrosProduct(dU_avTop,U_avTop,de1_1,e1_1);

			  dR11 = dUnit(e1_1,de1_1);
			  dR12 = dUnit(e2_1,de2_1);
			  dR13 = dUnit(U_avTop,dU_avTop);
		  } 
		  else if (Top_Btm[n_basis] == -1){
			  for(int i=0; i<3; i++){
				  for(int j=0; j<3; j++){
					  dU_avBtm(i,j) = U_Derivatives[n_basis][k][i][j];
			  	  }
			  }
			  for(int i=0;i<3;i++){
			  		de1_2(0,i) =  dU_avBtm(1,i)*diffBtm(2)-dU_avBtm(2,i)*diffBtm(1);
			  		de1_2(1,i) = -dU_avBtm(0,i)*diffBtm(2)+dU_avBtm(2,i)*diffBtm(0);
			  		de1_2(2,i) =  dU_avBtm(0,i)*diffBtm(1)-dU_avBtm(1,i)*diffBtm(0);
			  }
			  de2_2 = dCrosProduct(dU_avBtm,U_avBtm,de1_2,e1_2);

			  dR21  = dUnit(e1_2,de1_2);
			  dR22  = dUnit(e2_2,de2_2);
			  dR23  = dUnit(U_avBtm,dU_avBtm);
		  }

		  for(int i=0;i<3;i++) {dQ(i) = 0.0; }
		  for(int i=0;i<3;i++)
		  {
			  dQ(i) += dR11(1,i)*Btm(0,2)+Top(0,1)*dR21(2,i)
		  		  -dR11(2,i)*Btm(0,1)-Top(0,2)*dR21(1,i)
		  		  +dR12(1,i)*Btm(1,2)+Top(1,1)*dR22(2,i)
		  		  -dR12(2,i)*Btm(1,1)-Top(1,2)*dR22(1,i)
		  		  +dR13(1,i)*Btm(2,2)+Top(2,1)*dR23(2,i)
		  		  -dR13(2,i)*Btm(2,1)-Top(2,2)*dR23(1,i);
		  }
		  for(int i=0;i<3;i++)
		  {
			  dU_final(0,i)=dQ(i);
		  }

		  for(int i=0;i<3;i++) {dQ(i) = 0.0; }
		  for(int i=0;i<3;i++)
		  {
			  dQ(i) += dR11(2,i)*Btm(0,0)+Top(0,2)*dR21(0,i)
		  	  	  -dR11(0,i)*Btm(0,2)-Top(0,0)*dR21(2,i)
		  	  	  +dR12(2,i)*Btm(1,0)+Top(1,2)*dR22(0,i)
		  	  	  -dR12(0,i)*Btm(1,2)-Top(1,0)*dR22(2,i)
		  	  	  +dR13(2,i)*Btm(2,0)+Top(2,2)*dR23(0,i)
		  	  	  -dR13(0,i)*Btm(2,2)-Top(2,0)*dR23(2,i);
		  }
		  for(int i=0;i<3;i++)
		  {
			  dU_final(1,i)=dQ(i);
		  }

		  for(int i=0;i<3;i++) {dQ(i) = 0.0; }
		  for(int i=0;i<3;i++)
		  {
			  dQ(i) += dR11(0,i)*Btm(0,1)+Top(0,0)*dR21(1,i)
		  		  -dR11(1,i)*Btm(0,0)-Top(0,1)*dR21(0,i)
		  		  +dR12(0,i)*Btm(1,1)+Top(1,0)*dR22(1,i)
		  		  -dR12(1,i)*Btm(1,0)-Top(1,1)*dR22(0,i)
		  		  +dR13(0,i)*Btm(2,1)+Top(2,0)*dR23(1,i)
		  		  -dR13(1,i)*Btm(2,0)-Top(2,1)*dR23(0,i);
		  }
		  for(int i=0;i<3;i++)
		  {
			  dU_final(2,i)=dQ(i);
		  }
		  U_final = rotVectorBasis(Top,Btm);
		  dU_final_norm = dUnit(U_final,dU_final);
		  U_final *= 1.0/U_final.modulo();

		  //-------------------------------
		  // calculating dScrew and drotAng
		  //-------------------------------
		  for(int i=0;i<3;i++) {dTrace(i) = 0.0;}

		  for (int i=0; i<3; i++){
			  for (int j=0; j<3; j++){
				  dQ(j) = 0.0;
				  dQ(j)+=  Top(0,i)*dR21(i,j) + Top(1,i)*dR22(i,j) + Top(2,i)*dR23(i,j)
					 + Btm(0,i)*dR11(i,j) + Btm(1,i)*dR12(i,j) + Btm(2,i)*dR13(i,j);
			  }
			  for (int i=0; i<3; i++){
				  dTrace(i) += dQ(i);
			  }
		  }

		  for (int i=0; i<3; i++){
			  drotAng(i) = -0.5*dTrace(i)/sqrt(1-Trace*Trace);
		  }

		  for (int i=0; i<3; i++){
			  for (int j=0; j<3; j++){
				  dScrew(i,j) = dU_final_norm(i,j)*rotAng + (U_final(i)*drotAng(j))*180/PI;
			  }
		  }

		  for (int j=0; j<3; j++){
			  dTwist(j) = 0.0;
			  for (int i=0; i<3; i++){
				  dTwist(j) += dScrew(i,j)*Top(2,i) + Screw(i)*dR13(i,j);
			  }
		  }
	// derivatives
		  for (int i=0; i<6; i++){
			  for (int j=0; j<3; j++){
				  atomDerivatives[atomNmbs[n_basis][i]][j] -= dTwist(j)*dCoord(i,k);
			  }
		  }
	  } // finish k
  } // close loop if loop
  } // finish n_basis

  //-----------------------------------------------------------------
  // additional derivatives with respect to diffTop and diffBtm atoms
  //-----------------------------------------------------------------
  Matrix<double> d_diffTop(3,3);
  Matrix<double> R(3,3);

  Vector e1; Vector e2; Vector e3;
  Matrix<double> de2(3,3); Matrix<double> de3(3,3);
  Matrix<double> de1_norm(3,3); Matrix<double> de2_norm(3,3); Matrix<double> de3_norm(3,3);

  for (int k=0; k<2; k++){
	  d_diffTop = 0.0;
	  for (int i=0; i<3; i++){
			  d_diffTop(i,i) = -1.0;
	  }
	  if (k == 0){
		  R = Btm;
		  e1 = U_avTop;
		  e2 = e1_1;
		  e3 = e2_1;
	  }
	  else if (k == 1){
		  R = Top;
		  e1 = U_avBtm;
		  e2 = e1_2;
		  e3 = e2_2;
	  }
	  for (int i=0; i<3; i++){
		  de2(0,i) =  e1(1)*d_diffTop(2,i) - e1(2)*d_diffTop(1,i);
		  de2(1,i) = -e1(0)*d_diffTop(2,i) + e1(2)*d_diffTop(0,i);
		  de2(2,i) =  e1(0)*d_diffTop(1,i) - e1(1)*d_diffTop(0,i);
	  }
	  for (int i=0; i<3; i++){
	  	  de3(0,i) =  e1(1)*de2(2,i) - e1(2)*de2(1,i);
	  	  de3(1,i) = -e1(0)*de2(2,i) + e1(2)*de2(0,i);
	  	  de3(2,i) =  e1(0)*de2(1,i) - e1(1)*de2(0,i);
	  }

	  de2_norm = dUnit(e2,de2);
	  de3_norm = dUnit(e3,de3);

	  for (int i=0;i<3;i++) {dQ(i) = 0.0; }
	  for (int i=0; i<3; i++){
		  dQ(i) += de2_norm(1,i)*R(0,2)+de3_norm(1,i)*R(1,2)-de2_norm(2,i)*R(0,1)-de3_norm(2,i)*R(1,1);
	  }
	  for (int i=0; i<3; i++){
		  dU_final(0,i) = dQ(i);
	  }

	  for (int i=0;i<3;i++) {dQ(i) = 0.0; }
	  for (int i=0; i<3; i++){
		  dQ(i) += de2_norm(2,i)*R(0,0)+de3_norm(2,i)*R(1,0)-de2_norm(0,i)*R(0,2)-de3_norm(0,i)*R(1,2);
	  }
	  for (int i=0; i<3; i++){
		  dU_final(1,i) = dQ(i);
	  }

	  for (int i=0;i<3;i++) {dQ(i) = 0.0; }
	  for (int i=0; i<3; i++){
		  dQ(i) += de2_norm(0,i)*R(0,1)+de3_norm(0,i)*R(1,1)-de2_norm(1,i)*R(0,0)-de3_norm(1,i)*R(1,0);
	  }
	  for (int i=0; i<3; i++){
		  dU_final(2,i) = dQ(i);
	  }

	  if (k == 1) {
		  for (int i=0; i<3; i++){
			  for (int j=0; j<3; j++){
				  dU_final(i,j) = -dU_final(i,j);
			  }
		  }
	  }
	  U_final = rotVectorBasis(Top,Btm);
	  dU_final_norm = dUnit(U_final,dU_final);
	  U_final *= 1.0/U_final.modulo();

	  for (int i=0;i<3;i++) {dTrace(i) = 0.0; }
	  for (int i=0;i<3;i++) {dQ(i) = 0.0; }

	  for (int n=0; n<3; n++){
		  for (int i=0; i<3; i++){
			  dQ(i) = de2_norm(n,i)*R(0,n) + de3_norm(n,i)*R(1,n);
			  dTrace(i) += dQ(i);
			  drotAng(i) = -0.5*dTrace(i)/sqrt(1-Trace*Trace);
		  }
	  }

	  for (int i=0; i<3; i++){
		  for (int j=0; j<3; j++){
			  dScrew(i,j) = dU_final_norm(i,j)*rotAng + (U_final(i)*drotAng(j))*180/PI;
	  	  }
	  }

	  for (int j=0; j<3; j++){
		  dTwist(j) = 0.0;
	  	  for (int i=0; i<3; i++){
	  		  dTwist(j) += dScrew(i,j)*Top(2,i);
	  	  }
	  }

	  if (k == 0){
		  for (int i=0; i<3; i++){
			  atomDerivatives[3][i]  -= dTwist(i);
			  atomDerivatives[42][i] += dTwist(i); 
		  }
	  }
	  else if (k == 1){
		  for (int i=0; i<3; i++){
			  atomDerivatives[18][i] -= dTwist(i);
			  atomDerivatives[27][i] += dTwist(i);
		  }
	  }
  }
/*
atomDerivatives[3][i]  -= dTwist(i); 
atomDerivatives[18][i] -= dTwist(i);
atomDerivatives[27][i] += dTwist(i);
atomDerivatives[42][i] += dTwist(i);
*/
  for(int atom=0; atom < 48; atom++){
	  Vector derivatives;
	  for (int i=0; i<3; i++){
		  derivatives(i) -= atomDerivatives[atom][i];
	  }
	  setAtomsDerivatives(atom,derivatives);
   }

  setValue(twst);
  setBoxDerivativesNoPbc();

}

 Matrix<double> Twist::basis(Vector a, Vector b, Vector c){
   return basis2Vectors(delta(a,b),delta(a,c));
  // Vector e1 = b - a ;
  // Vector e2 = crossProduct(e1,(c - a));
  // Vector e3 = crossProduct(e1,e2);
  // Matrix<double> M(3,3);
  // for (int i=0;i<3;i++) 
  // {
  //   M(0,i) = e1(i)/e1.modulo();
  //   M(1,i) = e2(i)/e2.modulo();
  //   M(2,i) = e3(i)/e3.modulo();
  // }
  // return M;
}

Matrix<double> Twist::basis2Vectors(Vector e1, Vector diff){
  Vector e2 = crossProduct(e1,diff);
  Vector e3 = crossProduct(e1,e2);
  Matrix<double> M(3,3);
  for (int i=0;i<3;i++) 
  {
    M(0,i) = e1(i)/e1.modulo();
    M(1,i) = e2(i)/e2.modulo();
    M(2,i) = e3(i)/e3.modulo();
  }
  return M;
}

// basis2Vectors_2 is a different version of the basis function from 2 vectors
Matrix<double> Twist::basis2Vectors_2(Vector RotVect, Vector N1N9){
  Vector e1 = crossProduct(RotVect,N1N9);
  Vector e2 = crossProduct(RotVect,e1);
    Matrix<double> M(3,3);
    for (int i=0;i<3;i++)
    {
      M(0,i) = e1(i)/e1.modulo();
      M(1,i) = e2(i)/e2.modulo();
      M(2,i) = RotVect(i)/RotVect.modulo();
    }
    return M;
}

Vector Twist::rotVector6Atoms(int arr[6]){
  int i0=arr[0]; int i1=arr[1]; int i2=arr[2];
  int i3=arr[3]; int i4=arr[4]; int i5=arr[5];
  Matrix<double> B1 = basis(getPosition(i0),getPosition(i1),getPosition(i2));
  Matrix<double> B2 = basis(getPosition(i3),getPosition(i4),getPosition(i5));
  return rotVectorBasis(B1,B2);
}

Vector Twist::rotVectorBasis( Matrix<double> R1,  Matrix<double> R2){
  Matrix<double> Q(3,3);
  Matrix<double> R1T(3,3);
  transpose(R1,R1T);
  mult(R1T,R2,Q);
  Vector u;
  u(0) = Q(1,2) - Q(2,1);
  u(1) = Q(2,0) - Q(0,2);
  u(2) = Q(0,1) - Q(1,0);
  //u *= 1/u.modulo();
  return u;
} 

double Twist::rotAngleBasis( Matrix<double> P1, Matrix<double> P2){
  double trace = traceBasis(P1,P2);
  double tr = (trace-1)*0.5;
  if(tr > 1.0) tr = 1.0;
  else if (tr < -1.0) tr =-1.0;
  double theta=acos(tr);//*180.0/PI;
  if(theta < 0.0) { fprintf(stderr,"negative theta %g \n", theta);}
  return theta;
}

double Twist:: traceBasis( Matrix<double> P1, Matrix<double> P2){
  Matrix<double> Q(3,3);
  Matrix<double> P1T(3,3);
  transpose(P1,P1T);
  mult(P1T,P2,Q);
  double trace=Q(0,0)+Q(1,1)+Q(2,2);
  return trace;
}

Matrix<double> Twist::dCrosProduct(Matrix<double> da, Vector a, Matrix<double> db, Vector b){
	Matrix<double> dc(3,3);
	for(int j=0;j<3;j++)
	{
		dc(0,j)= da(1,j)*b(2)-da(2,j)*b(1)+a(1)*db(2,j)-a(2)*db(1,j);
	  dc(1,j)=-da(0,j)*b(2)+da(2,j)*b(0)-a(0)*db(2,j)+a(2)*db(0,j);
	  dc(2,j)= da(0,j)*b(1)-da(1,j)*b(0)+a(0)*db(1,j)-a(1)*db(0,j);
	}
	return dc;
}

Matrix<double> Twist::dUnit(Vector a, Matrix<double> da){
	double a_norm=a.modulo();
	Matrix<double> da_unit(3,3);
	Vector dk;
	for(int i=0;i<3;i++){
		dk(i)=0.0;
		for(int k=0;k<3;k++){
			dk(i) += da(k,i)*a(k);
		}
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			da_unit(i,j) = da(i,j)/a_norm-1.0/(a_norm*a_norm*a_norm)*a(i)*dk(j);
		}
	}
	return da_unit;
}

}
}
