///All structures used in BD solver
//Uses attribute aligned and padding to ensure that
//	memory alignment matches alignment for structs on host
//structures defined in StdAfx.h for host side

typedef struct FD_coeff
{
	double4 cA;			///finite difference coefficients
	double4 cB;
	double4 cC;
	double Alpha[5];
	int Type;
	int indc;
	int inde;		//for outlet inde is i-1, indw is i-2
	int indw;
	int indn;
	int inds;
} FDcoeff __attribute__((aligned(64)));

typedef struct Particle
{
	double2 pos;		///position vector
	uint Num_rep;		///number of particles represented by this particle
	short type;			///type of particle (cooresponds to Param element
	short Dep_Flag; 	//-2 if waiting for re-release, -1 for not deposited, > 0 signifies the BL location it has deposited at.
	ushort Dep_timer;	//timer set to specified value once particle deposits and decrements at 0, particle is deposited
	ushort timer;		///time for use in re-releasing
	int loc;			//Node number particle is located within
} par __attribute__((aligned(32)));

typedef struct Particle_Param
{
	double Q_A_prime[2];	///Used in deposition/rebound calculation
	double Q_A[2];			///used in dep/reb
	double tau_crit[2];	///critical shear stress
	double Dp;			///particle diameter
	double Mp;			///particle mass
	double Kth;			///Thermophoretic Coeff
	double D_dist;		///% of particles distributed in bin
	double L_coeff;		///Lift coefficient (Lift force = L_coeff*Tau^1.5)
	double D_coeff;		///Drag coefficient (Drag force = D_coeff*Tau)
} Pparam __attribute__ ((aligned(32)));

typedef struct NodeCoeff
{
	double4 CoeffT00; //coefficients used to calculate NodeVar variables in update_Nodes kernel
	double4 CoeffT10;
	double4 CoeffT01;
	double4 CoeffT11;
	double4 CoeffU00;
	double4 CoeffU10;
	double4 CoeffU01;
	double4 CoeffU11;
	int4 neigh;		//Indicies of 4 lattice points enclosing square (U and T array are of size vlb.nX*(vlb.Channel_Height+1))
	char Pad[16];
} nodeC __attribute__ ((aligned(32)));

typedef struct NodeInfo
{
	short BLind[MAX_BL_PER_NODE];	//Inidices of three closest boundary links (-1 used to represent when less than three BL's near)
	short Wall_Flag;  //0 is no walls nearby, 1 if bottom wall nearby and 2 if top wall nearby, -1 if inactive
	char Pad[14 - 2 * MAX_BL_PER_NODE];
} nodeI __attribute__ ((aligned(16)));

typedef struct Neighbors
{
	int2 ii00;
	int2 ii10;
	int2 ii01;
	int2 ii11;
} Nact __attribute__((aligned(32)));

typedef struct NodeVar
{
	double4 Temps;
	double2 U00;
	double2 U10;
	double2 U01;
	double2 U11;
} nodeV __attribute__ ((aligned(32)));

typedef struct BL_bounds
{//Defines range of BL used in particle section of domain
	int MIN_BL_BOT;
	int	MAX_BL_BOT;
	int MIN_BL_TOP;
	int MAX_BL_TOP;
} BLbound __attribute__ ((aligned(16)));

typedef struct BL_Links
{
	double2 vP0;		//Location of left node
	double2 vP1;		//Location of right node
	double2 vTvec;	//tangential vector (points to the right)
	double2 vNvec;   //normal vector pointing into domain
	double Tau;		//Shear stress at location
	double blLen;	//length of BL
	int Node_loc;	//points to location BL is located in (or majority of BL when across two)
	int Color_ind;
	short P0ind;		//index of right node in vls.C array (maybe unnecessary)
	short P1ind;		//index of left node in vls.C array
	short dir;		//direction of shear stress
	short int_type;	//Designates if interface between particle and surface is soot-soot or soot-wall
} bLinks __attribute__ ((aligned(32)));

typedef struct Tr_Param
{//Used for re-releasing particles
	double Top_location;	//Y value of uppermost node
	double Bottom_location;	//Y value of lowermost node
	double umax_val;		//Max velocity at inlet
	double bval;			//spacing between upper and lower wall at particle inlet
	double offset_y;		//Location of wall at particle inlet
	double X_release;	//X location of release
	uint BL_rel_bot;	//index of bottom BL at particle inlet 
	uint BL_rel_top;	//index of top BL at particle inlet
	uint Uvals_start;	//location of starting point for Uvals used in re-distribution
	char Pad[4];		//padding
} Trparam __attribute__((aligned(64)));

typedef struct FoulInfo
{//used for updating fouling layer
	double4 WeightsL;	//Weights applied to BL deposits to left
	double4 WeightsR;	//Weights applied to BL deposits to right
	double2 vN;			//Normal vector of C node
	double disp;		//Current displacement distance	
	uint BL_ind;		//index of BL to left of node
	uint C_ind;			//index of wall node this struct corresponds to
} foulI __attribute__((aligned(32)));

typedef struct Ramp_info
{
	double Ybegin;
	double Coeff;
	uint IOind; //Index of LS nodes at beginning and end of TR area
	uint Cind; //direction traversed from IO_end when creating ramp
	char pad[8];

} rampI __attribute__((aligned(32)));
