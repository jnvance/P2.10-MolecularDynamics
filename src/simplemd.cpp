#include "Vector.h"
#include "Random.h"
#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <mpi.h>

using namespace std;
using namespace PLMD;

#define INDEX(i0,i1,i2,n) (i0 * n[1] + i1) * n[2] + i2

inline int modulo(int a, int b){
  int r = a % b;
  return r + b * (r < 0);
}

class SimpleMD
{

  int iv[32];
  int iy;
  int iset;
  double gset;
  bool write_positions_first;
  bool write_statistics_first;
  int write_statistics_last_time_reopened;
  FILE* write_statistics_fp;

  // For linked cells
  int ndomains[3];
  int totdomains;
  double dcell[3];
  vector<vector<int> > neighbors;

  public:
  SimpleMD(){
    for(int i=0;i<32;i++) iv[i]=0.0;
    iy=0;
    iset=0;
    gset=0.0;
    write_positions_first=true;
    write_statistics_first=true;
    write_statistics_last_time_reopened=0;
    write_statistics_fp=NULL;
  }

  // For MPI
  int world_rank, world_size;
  int myrank, nprocs;
  int nrep, irep;
  MPI_Comm comm;
  int myrank_col, nprocs_col;
  MPI_Comm comm_col;


  private:

  void
  read_input(FILE*   fp,
             double& temperature,
             double& tstep,
             double& friction,
             double& forcecutoff,
             double& listcutoff,
             int&    nstep,
             int&    nconfig,
             int&    nstat,
             bool&   wrapatoms,
             string& inputfile,
             string& outputfile,
             string& trajfile,
             string& statfile,
             int&    maxneighbours,
             int&    idum,
             int&    exchangestride)
  {
    temperature=1.0;
    tstep=0.005;
    friction=0.0;
    forcecutoff=2.5;
    listcutoff=3.0;
    nstep=1;
    nconfig=10;
    nstat=1;
    maxneighbours=1000;
    idum=0;
    exchangestride=50;
    wrapatoms=false;
    statfile="";
    trajfile="";
    outputfile="";
    inputfile="";

    string line;

    line.resize(256);
    char buffer[256];
    char buffer1[256];

    while(fgets(buffer,256,fp)){
      line=buffer;
      for(int i=0;i<line.length();++i) if(line[i]=='#' || line[i]=='\n') line.erase(i);
      for(int i=line.length()-1;i>=0;--i){
        if(line[i]!=' ')break;
        line.erase(i);
      }
      buffer[0]=0;
      sscanf(line.c_str(),"%s",buffer);
      if(strlen(buffer)==0) continue;
      string keyword=buffer;
      if(keyword=="temperature")
        sscanf(line.c_str(),"%s %lf",buffer,&temperature);
      else if(keyword=="tstep")
        sscanf(line.c_str(),"%s %lf",buffer,&tstep);
      else if(keyword=="friction")
        sscanf(line.c_str(),"%s %lf",buffer,&friction);
      else if(keyword=="forcecutoff")
        sscanf(line.c_str(),"%s %lf",buffer,&forcecutoff);
      else if(keyword=="listcutoff")
        sscanf(line.c_str(),"%s %lf",buffer,&listcutoff);
      else if(keyword=="nstep")
        sscanf(line.c_str(),"%s %d",buffer,&nstep);
      else if(keyword=="nconfig")
      {
        sscanf(line.c_str(),"%s %d %s",buffer,&nconfig,buffer1);
        trajfile=buffer1;
      }
      else if(keyword=="nstat")
      {
        sscanf(line.c_str(),"%s %d %s",buffer,&nstat,buffer1);
        statfile=buffer1;
      }
      else if(keyword=="wrapatoms")
      {
        sscanf(line.c_str(),"%s %s",buffer,buffer1);
        if(buffer1[0]=='T' || buffer1[0]=='t') wrapatoms=true;
      }
      else if(keyword=="maxneighbours")
        sscanf(line.c_str(),"%s %d",buffer,&maxneighbours);
      else if(keyword=="inputfile")
      {
        sscanf(line.c_str(),"%s %s",buffer,buffer1);
        inputfile=buffer1;
      }
      else if(keyword=="outputfile")
      {
        sscanf(line.c_str(),"%s %s",buffer,buffer1);
        outputfile=buffer1;
      }
      else if(keyword=="idum")
        sscanf(line.c_str(),"%s %d",buffer,&idum);
      else if(keyword=="exchangestride")
        sscanf(line.c_str(),"%s %d",buffer,&exchangestride);
      else{
        fprintf(stderr,"Unknown keywords :%s\n",keyword.c_str());
        exit(1);
      }
    }

    if(inputfile.length()==0){
        fprintf(stderr,"Specify input file\n");
        exit(1);
    }
    if(outputfile.length()==0){
        fprintf(stderr,"Specify output file\n");
        exit(1);
    }
    if(trajfile.length()==0){
        fprintf(stderr,"Specify traj file\n");
        exit(1);
    }
    if(statfile.length()==0){
        fprintf(stderr,"Specify stat file\n");
        exit(1);
    }
  }

  void read_natoms(const string & inputfile,int & natoms){
    // read the number of atoms in file "input.xyz"
    FILE* fp=fopen(inputfile.c_str(),"r");
    fscanf(fp,"%d",&natoms);
    fclose(fp);
  }

  void read_positions(const string& inputfile,int natoms,vector<Vector>& positions,double cell[3]){
    // read positions and cell from a file called inputfile
    // natoms (input variable) and number of atoms in the file should be consistent
    FILE* fp=fopen(inputfile.c_str(),"r");
    char buffer[256];
    char atomname[256];
    fgets(buffer,256,fp);
    fscanf(fp,"%lf %lf %lf",&cell[0],&cell[1],&cell[2]);
    for(int i=0;i<natoms;i++){
      fscanf(fp,"%s %lf %lf %lf",atomname,&positions[i][0],&positions[i][1],&positions[i][2]);
    // note: atomname is read but not used
    }
    fclose(fp);
  }

  void randomize_velocities(const int natoms,const double temperature,const vector<double>&masses,vector<Vector>& velocities,Random&random){
    // randomize the velocities according to the temperature
    for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<3;i++)
        velocities[iatom][i]=sqrt(temperature/masses[iatom])*random.Gaussian();
  }

  void pbc(const double cell[3],const Vector & vin,Vector & vout){
    // apply periodic boundary condition to a vector
    for(int i=0;i<3;i++){
      vout[i]=vin[i]-floor(vin[i]/cell[i]+0.5)*cell[i];
    }
  }

  void check_list(const int natoms,const vector<Vector>& positions,const vector<Vector>&positions0,const double listcutoff,
                  const double forcecutoff,bool & recompute)
  {
    // check if the neighbour list have to be recomputed
    Vector displacement;  // displacement from positions0 to positions
    double delta2;        // square of the 'skin' thickness
    recompute=false;
    delta2=(0.5*(listcutoff-forcecutoff))*(0.5*(listcutoff-forcecutoff));
    // if ANY atom moved more than half of the skin thickness, recompute is set to .true.
    for(int iatom=0;iatom<natoms;iatom++){
      for(int k=0;k<3;k++) displacement[k]=positions[iatom][k]-positions0[iatom][k];
      double s=0.0;
      for(int k=0;k<3;k++) s+=displacement[k]*displacement[k];
      if(s>delta2) recompute=true;
    }
  }

  void compute_cells(const int natoms, const double cell[3], const double listcutoff){

    totdomains = 1;
    for(int k=0;k<3;k++) {
      ndomains[k] = cell[k]/listcutoff + 1;
      totdomains *= ndomains[k];
      dcell[k] = cell[k]/ndomains[k];
    }

    // Works efficiently only for ndomains > 3 on each side
    bool do_cell_list = (totdomains >= 64);

    // Calculate index of neighboring domains
    // vector<vector<int> > neighbors(totdomains);
    neighbors.resize(totdomains);
    for (int i0 = 0; i0 < ndomains[0]; ++i0)
      for (int i1 = 0; i1 < ndomains[1]; ++i1)
        for (int i2 = 0; i2 < ndomains[2]; ++i2){

          int mydomain = INDEX(i0,i1,i2,ndomains);

          for (int d0=-1; d0<=+1; ++d0)
            for (int d1=-1; d1<=+1; ++d1)
              for (int d2=-1; d2<=+1; ++d2){
                int neighbor = INDEX(
                                      modulo(i0+d0,ndomains[0]),
                                      modulo(i1+d1,ndomains[1]),
                                      modulo(i2+d2,ndomains[2]),
                                      ndomains);
                assert(neighbor >= 0 && neighbor < totdomains);
                neighbors[mydomain].push_back(neighbor);
              }
        }
  }

  void compute_list(const int natoms,const vector<Vector>& positions,const double cell[3],const double listcutoff,vector<vector<int> >& list){
    Vector distance;     // distance of the two atoms
    Vector distance_pbc; // minimum-image distance of the two atoms
    double listcutoff2;  // squared list cutoff
    listcutoff2=listcutoff*listcutoff;


    #ifdef _CELL_LIST

    list.assign(natoms,vector<int>());
    // Recalculate cell list
    vector<vector<int> > domain(totdomains);
    int i[3];
    for (int iatom = 0; iatom < natoms; ++iatom){
      for(int k=0;k<3;k++) i[k] = modulo((int)floor(positions[iatom][k]/dcell[k]),ndomains[k]);
      domain[INDEX(i[0],i[1],i[2],ndomains)].push_back(iatom);
    }

    for(int dom = 0; dom < totdomains; dom++){
      for(int i = 0; i < domain[dom].size(); i++ ){
        int iatom = domain[dom][i];
        for(int jnei = 0; jnei < neighbors[dom].size(); jnei++ ){
          for(int j = 0; j < domain[ neighbors[dom][jnei] ].size(); j++ ){
            int jatom = domain[ neighbors[dom][jnei] ][j];
            if(jatom<=iatom)continue;
            for(int k=0;k<3;k++) distance[k]=positions[iatom][k]-positions[jatom][k];
            pbc(cell,distance,distance_pbc);
            // if the interparticle distance is larger than the cutoff, skip
            double d2=0; for(int k=0;k<3;k++) d2+=distance_pbc[k]*distance_pbc[k];
            if(d2>listcutoff2)continue;
            list[iatom].push_back(jatom);
          }
        }
      }
    }

    #else

    list.assign(natoms/nprocs+1,vector<int>());
    for(int iatom=myrank;iatom<natoms-1;iatom+=nprocs){
      for(int jatom=iatom+1;jatom<natoms;jatom++){
        for(int k=0;k<3;k++) distance[k]=positions[iatom][k]-positions[jatom][k];
        pbc(cell,distance,distance_pbc);
        // if the interparticle distance is larger than the cutoff, skip
        double d2=0; for(int k=0;k<3;k++) d2+=distance_pbc[k]*distance_pbc[k];
        if(d2>listcutoff2)continue;
        list[iatom/nprocs].push_back(jatom);
      }
    }

    #endif

  }

  void compute_forces(const int natoms,const vector<Vector>& positions,const double cell[3],
                      double forcecutoff,const vector<vector<int> >& list,vector<Vector>& forces,double & engconf)
  {
    Vector distance;        // distance of the two atoms
    Vector distance_pbc;    // minimum-image distance of the two atoms
    double distance_pbc2;   // squared minimum-image distance
    double forcecutoff2;    // squared force cutoff
    Vector f;               // force
    double engcorrection;   // energy necessary shift the potential avoiding discontinuities

    forcecutoff2=forcecutoff*forcecutoff;
    engconf=0.0;
    // vector<Vector> forces(natoms);
    for(int i=0;i<natoms;i++)for(int k=0;k<3;k++) forces[i][k]=0.0;
    engcorrection=4.0*(1.0/pow(forcecutoff2,6.0)-1.0/pow(forcecutoff2,3));

    for(int iatom=myrank;iatom<natoms-1;iatom+=nprocs){

      #ifdef _CELL_LIST

      for(int jlist=0;jlist<list[iatom].size();jlist++){
        int jatom=list[iatom][jlist];

      #else

      for(int jlist=0;jlist<list[iatom/nprocs].size();jlist++){
        int jatom=list[iatom/nprocs][jlist];

      #endif

        for(int k=0;k<3;k++) distance[k]=positions[iatom][k]-positions[jatom][k];
        pbc(cell,distance,distance_pbc);
        distance_pbc2=0.0; for(int k=0;k<3;k++) distance_pbc2+=distance_pbc[k]*distance_pbc[k];
        // if the interparticle distance is larger than the cutoff, skip
        if(distance_pbc2>forcecutoff2) continue;
        double distance_pbc6=distance_pbc2*distance_pbc2*distance_pbc2;
        double distance_pbc8=distance_pbc6*distance_pbc2;
        double distance_pbc12=distance_pbc6*distance_pbc6;
        double distance_pbc14=distance_pbc12*distance_pbc2;
        engconf+=4.0*(1.0/distance_pbc12 - 1.0/distance_pbc6) - engcorrection;
        for(int k=0;k<3;k++) f[k]=2.0*distance_pbc[k]*4.0*(6.0/distance_pbc14-3.0/distance_pbc8);
        // same force on the two atoms, with opposite sign:
        for(int k=0;k<3;k++) forces[iatom][k]+=f[k];
        for(int k=0;k<3;k++) forces[jatom][k]-=f[k];
      }
    }

    MPI_Allreduce(MPI_IN_PLACE, &forces[0][0], 3*natoms, MPI_DOUBLE, MPI_SUM, comm);

    // printf(" >> MPI rank: %3.3d-%3.3d\tengconf: %f\n",irep,myrank,engconf);

    MPI_Allreduce(MPI_IN_PLACE, &engconf, 1, MPI_DOUBLE, MPI_SUM, comm);

    // printf(" >> MPI rank: %3.3d-%3.3d\tengconf_reduced: %f\n",irep,myrank,engconf);
    assert(not(isnan(engconf)));

  }

  void compute_engkin(const int natoms,const vector<double>& masses,const vector<Vector>& velocities,double & engkin)
  {
    // calculate the kinetic energy from the velocities
    engkin=0.0;
    for(int iatom=0;iatom<natoms;iatom++)for(int k=0;k<3;k++){
      engkin+=0.5*masses[iatom]*velocities[iatom][k]*velocities[iatom][k];
    }
  }


  void thermostat(const int natoms,const vector<double>& masses,const double dt,const double friction,
                  const double temperature,vector<Vector>& velocities,double & engint,Random & random){
    // Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
    // it is a linear combination of old velocities and new, randomly chosen, velocity,
    // with proper coefficients
    double c1,c2;
    c1=exp(-friction*dt);
    for(int iatom=0;iatom<natoms;iatom++){
      c2=sqrt((1.0-c1*c1)*temperature/masses[iatom]);
      for(int i=0;i<3;i++){
        engint+=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
        velocities[iatom][i]=c1*velocities[iatom][i]+c2*random.Gaussian();
        engint-=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
      }
    }
  }

  void write_positions(const string& trajfile,int natoms,const vector<Vector>& positions,const double cell[3],const bool wrapatoms)
  {
    // write positions on file trajfile
    // positions are appended at the end of the file
    Vector pos;
    FILE*fp;
    if(write_positions_first){
      if (myrank==0)
        fp=fopen(trajfile.c_str(),"w");
      else
        fp=fopen("/dev/null","w");
      write_positions_first=false;
    } else {
      if (myrank==0)
        fp=fopen(trajfile.c_str(),"a");
      else
        fp=fopen("/dev/null","a");
    }
    fprintf(fp,"%d\n",natoms);
    fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
    for(int iatom=0;iatom<natoms;iatom++){
    // usually, it is better not to apply pbc here, so that diffusion
    // is more easily calculated from a trajectory file:
      if(wrapatoms) pbc(cell,positions[iatom],pos);
      else for(int k=0;k<3;k++) pos[k]=positions[iatom][k];
      fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
    }
    fclose(fp);
  }

  void write_final_positions(const string& outputfile,int natoms,const vector<Vector>& positions,const double cell[3],const bool wrapatoms)
  {
    // write positions on file outputfile
    Vector pos;
    FILE*fp;
    if (myrank==0)
      fp=fopen(outputfile.c_str(),"w");
    else
      fp=fopen("/dev/null","w");
    fprintf(fp,"%d\n",natoms);
    fprintf(fp,"%f %f %f\n",cell[0],cell[1],cell[2]);
    for(int iatom=0;iatom<natoms;iatom++){
    // usually, it is better not to apply pbc here, so that diffusion
    // is more easily calculated from a trajectory file:
      if(wrapatoms) pbc(cell,positions[iatom],pos);
      else for(int k=0;k<3;k++) pos[k]=positions[iatom][k];
      fprintf(fp,"Ar %10.7f %10.7f %10.7f\n",pos[0],pos[1],pos[2]);
    }
    fclose(fp);
  }


  void write_statistics(const string & statfile,const int istep,const double tstep,
                        const int natoms,const double engkin,const double engconf,const double engint){
    // write statistics on file statfile
    if(write_statistics_first){
    // first time this routine is called, open the file
      if(myrank==0) write_statistics_fp=fopen(statfile.c_str(),"w");
      else          write_statistics_fp=fopen("/dev/null","w");
      write_statistics_first=false;
    }
    if(istep-write_statistics_last_time_reopened>100){
    // every 100 steps, reopen the file to flush the buffer
      fclose(write_statistics_fp);
      if(myrank==0) write_statistics_fp=fopen(statfile.c_str(),"a");
      else          write_statistics_fp=fopen("/dev/null","a");
      write_statistics_last_time_reopened=istep;
    }
    fprintf(write_statistics_fp,"%d %f %f %f %f %f\n",istep,istep*tstep,2.0*engkin/(3.0*natoms),engconf,engkin+engconf,engkin+engconf+engint);
  }


  int parallel_tempering(const int istep, const int nstep, const int exchangestride, const int partner,
                          Random& random, const double engconf, const double temperature,
                          vector<Vector>& positions, vector<Vector>& velocities, const int natoms)
  {
    int swap = 0;         // store result of metropolis check (0 or 1)
    double ET_buf[2];     // buffer for energy and temperature exchanges

    // perform parallel tempering on master processes only
    if (myrank == 0)
    {

      if ( irep > partner ) // process of higher rank sends data to lower rank
      {
        ET_buf[0] = engconf;
        ET_buf[1] = temperature;

        MPI_Send(ET_buf, 2, MPI_DOUBLE, partner, istep+nstep, comm_col);

        // metropolis check performed in lower-ranked process

        MPI_Recv(&swap, 1, MPI_INT, partner, istep+2*nstep, comm_col, MPI_STATUS_IGNORE);
      }
      else // lower-ranked partner performs the metropolis check
      {

        MPI_Recv(ET_buf, 2, MPI_DOUBLE, partner, istep+nstep, comm_col, MPI_STATUS_IGNORE);

        // check metropolis criterion
        double acc = (1.0/temperature - 1.0/ET_buf[1])*(engconf-ET_buf[0]);
        if (acc > 0)                        swap = 1;
        else if (exp(acc) > random.U01())   swap = 1;

        MPI_Send(&swap, 1, MPI_INT, partner, istep+2*nstep, comm_col);
      }
    }

    vector<Vector> buffer(natoms);

    // Broadcast the result of metropolis step
    MPI_Bcast(&swap, 1, MPI_INT, 0, comm);

    // swap particle positions and velocities
    if(swap==1){
      // put positions in a buffer
      buffer = positions;

      // send positions buffer and receive positions
      MPI_Sendrecv( &buffer[0][0],    3*natoms, MPI_DOUBLE, partner, 20*nstep+istep,
                    &positions[0][0], 3*natoms, MPI_DOUBLE, partner, 20*nstep+istep,
                    comm_col, MPI_STATUS_IGNORE);

      // put velocities in a buffer
      buffer.clear();
      buffer = velocities;

      // send velocities buffer and receive velocities
      MPI_Sendrecv( &buffer[0][0],     3*natoms, MPI_DOUBLE, partner, 40*nstep+istep,
                    &velocities[0][0], 3*natoms, MPI_DOUBLE, partner, 40*nstep+istep,
                    comm_col, MPI_STATUS_IGNORE);

      // send and receive temperatures for rescaling
      MPI_Sendrecv(&temperature, 1, MPI_DOUBLE, partner, istep+50*nstep,
                   &ET_buf[1],   1, MPI_DOUBLE, partner, istep+50*nstep, comm_col, MPI_STATUS_IGNORE);

      // rescale velocities accordingly
      double factor = sqrt( temperature / ET_buf[1] );
      for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<3;i++) velocities[iatom][i] *= factor;
    }

    return swap;
  }




  public:
  int main(FILE*in,FILE*out){
    int            natoms;       // number of atoms
    vector<Vector> positions;    // atomic positions
    vector<Vector> velocities;   // velocities
    vector<double> masses;       // masses
    vector<Vector> forces;       // forces
    double         cell[3];      // cell size
    double         cell9[3][3];  // cell size

    // neighbour list variables
    // see Allen and Tildesey book for details
    vector< vector<int> >  list;         // neighbour list
    vector<Vector> positions0;   // reference atomic positions, i.e. positions when the neighbour list

    // input parameters
    // all of them have a reasonable default value, set in read_input()
    double      tstep;             // simulation timestep
    double      temperature;       // temperature
    double      friction;          // friction for Langevin dynamics (for NVE, use 0)
    double      listcutoff;        // cutoff for neighbour list
    double      forcecutoff;       // cutoff for forces
    int         nstep;             // number of steps
    int         nconfig;           // stride for output of configurations
    int         nstat;             // stride for output of statistics
    int         maxneighbour;      // maximum average number of neighbours per atom
    int         idum;              // seed
    int         exchangestride;    // stride to perform parallel tempering
    bool        wrapatoms;         // if true, atomic coordinates are written wrapped in minimal cell
    string      inputfile;         // name of file with starting configuration (xyz)
    string      outputfile;        // name of file with final configuration (xyz)
    string      trajfile;          // name of the trajectory file (xyz)
    string      statfile;          // name of the file with statistics
    string      string;            // a string for parsing


    double engkin;                 // kinetic energy
    double engconf;                // configurational energy
    double engint;                 // integral for conserved energy in Langevin dynamics

    bool recompute_list;           // control if the neighbour list have to be recomputed

    Random random;                 // random numbers stream

    read_input(in,temperature,tstep,friction,forcecutoff,
               listcutoff,nstep,nconfig,nstat,
               wrapatoms,inputfile,outputfile,trajfile,statfile,
               maxneighbour,idum,exchangestride);

    // number of atoms is read from file inputfile
    read_natoms(inputfile,natoms);

    // write the parameters in output so they can be checked
    if (myrank==0){
      fprintf(stdout,"%s %s\n","Starting configuration           :",inputfile.c_str());
      fprintf(stdout,"%s %s\n","Final configuration              :",outputfile.c_str());
      fprintf(stdout,"%s %d\n","Number of atoms                  :",natoms);
      fprintf(stdout,"%s %f\n","Temperature                      :",temperature);
      fprintf(stdout,"%s %f\n","Time step                        :",tstep);
      fprintf(stdout,"%s %f\n","Friction                         :",friction);
      fprintf(stdout,"%s %f\n","Cutoff for forces                :",forcecutoff);
      fprintf(stdout,"%s %f\n","Cutoff for neighbour list        :",listcutoff);
      fprintf(stdout,"%s %d\n","Number of steps                  :",nstep);
      fprintf(stdout,"%s %d\n","Stride for trajectory            :",nconfig);
      fprintf(stdout,"%s %s\n","Trajectory file                  :",trajfile.c_str());
      fprintf(stdout,"%s %d\n","Stride for statistics            :",nstat);
      fprintf(stdout,"%s %s\n","Statistics file                  :",statfile.c_str());
      fprintf(stdout,"%s %d\n","Max average number of neighbours :",maxneighbour);
      fprintf(stdout,"%s %d\n","Seed                             :",idum);
      fprintf(stdout,"%s %d\n","Stride for parallel tempering    :",exchangestride);
      fprintf(stdout,"%s %s\n","Are atoms wrapped on output?     :",(wrapatoms?"T":"F"));
    }

    // Setting the seed
    random.setSeed(idum);

    // allocation of dynamical arrays
    positions.resize(natoms);
    positions0.resize(natoms);
    velocities.resize(natoms);
    forces.resize(natoms);
    masses.resize(natoms);
    list.resize(natoms);

    // masses are hard-coded to 1
    for(unsigned i=0;i<natoms;++i) masses[i]=1.0;

    // energy integral initialized to 0
    engint=0.0;

    // positions are read from file inputfile
    read_positions(inputfile,natoms,positions,cell);

    // velocities are randomized according to temperature
    randomize_velocities(natoms,temperature,masses,velocities,random);

  #ifdef _CELL_LIST

    // compute cell list domains and neighbors
    compute_cells(natoms,cell,listcutoff);

  #endif

    // neighbour list are computed, and reference positions are saved
    compute_list(natoms,positions,cell,listcutoff,list);

    int list_size=0;
    for(int i=0;i<list.size();i++) list_size+=list[i].size();
    if (myrank==0) fprintf(stdout,"List size: %d\n",list_size);
    for(int iatom=0;iatom<natoms;++iatom) for(int k=0;k<3;++k) positions0[iatom][k]=positions[iatom][k];

    // forces are computed before starting md
    compute_forces(natoms,positions,cell,forcecutoff,list,forces,engconf);

    // count swaps with partner of higher rank
    int swap_attempts = 0; // attempted swaps
    int swap_count = 0;    // successful swaps

    // here is the main md loop
    // Langevin thermostat is applied before and after a velocity-Verlet integrator
    // the overall structure is:
    //   thermostat
    //   update velocities
    //   update positions
    //   (eventually recompute neighbour list)
    //   compute forces
    //   update velocities
    //   thermostat
    //   (eventually dump output informations)
    //   perform parallel tempering
    for(int istep=0;istep<nstep;istep++){
      thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);

      for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
        velocities[iatom][k]+=forces[iatom][k]*0.5*tstep/masses[iatom];

      for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
        positions[iatom][k]+=velocities[iatom][k]*tstep;

      // a check is performed to decide whether to recalculate the neighbour list
      check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute_list);
      if(recompute_list){
        compute_list(natoms,positions,cell,listcutoff,list);
        for(int iatom=0;iatom<natoms;++iatom) for(int k=0;k<3;++k) positions0[iatom][k]=positions[iatom][k];
        fprintf(stdout,"Neighbour list recomputed at step %d\n",istep);
        int list_size=0;
        for(int i=0;i<list.size();i++) list_size+=list[i].size();
        fprintf(stdout,"List size: %d\n",list_size);
      }

      compute_forces(natoms,positions,cell,forcecutoff,list,forces,engconf);

      for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
        velocities[iatom][k]+=forces[iatom][k]*0.5*tstep/masses[iatom];

      thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);

      // kinetic energy is calculated
      compute_engkin(natoms,masses,velocities,engkin);

      // eventually, write positions and statistics
      if((istep+1)%nconfig==0) write_positions(trajfile,natoms,positions,cell,wrapatoms);
      if((istep+1)%nstat==0)   write_statistics(statfile,istep+1,tstep,natoms,engkin,engconf,engint);

      // perform parallel tempering
      if ((nrep>1)&&(istep%exchangestride==0)){
        // determine partner to exchange with
        int partner=irep+(((istep+1)/exchangestride+irep)%2)*2-1;
        if(partner<0) partner=0;
        if(partner>=nrep) partner=nrep-1;

        int swap = 0;
        if (partner != irep) // to reduce communication, do nothing if process is partnered with self
          swap = parallel_tempering(istep,nstep,exchangestride,partner,random,engconf,temperature,positions,velocities,natoms);

        if (partner > irep) {
          swap_attempts++;
          swap_count += swap;
          if(myrank==0 && swap==1)
            printf("\nParticles swapped between simulations %4d and %-4d\n\n",irep,partner);
        }
      }
    }

    write_final_positions(outputfile,natoms,positions,cell,wrapatoms);

    // close the statistic file if it was open:
    if(write_statistics_fp) fclose(write_statistics_fp);

    // gather data regarding swaps across master processes
    if(nrep>1 && myrank==0){
      vector<int> swap_attempts_data(nrep);
      vector<int> swap_count_data(nrep);
      vector<double> temperature_data(nrep);

      MPI_Gather(&swap_attempts,1,MPI_INT,&swap_attempts_data[0],1,MPI_INT,0,comm_col);
      MPI_Gather(&swap_count,1,MPI_INT,&swap_count_data[0],1,MPI_INT,0,comm_col);
      MPI_Gather(&temperature,1,MPI_DOUBLE,&temperature_data[0],1,MPI_DOUBLE,0,comm_col);

      if(myrank_col==0)
        for (int i = 0; i < nrep-1; ++i)
          printf("Swaps between %3d (T = %5.3f) and %-3d (T = %5.3f) :\t%5d/%-5d\n",
            i,temperature_data[i],i+1,temperature_data[i+1],swap_count_data[i],swap_attempts_data[i]);

    }

    return 0;
  }


};

int main(int argc,char*argv[]){

  // Note: disabled receiving input via stdin
  assert(argc>1);

  SimpleMD smd;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &smd.world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &smd.world_size);

  // Determine how many simulations to implement
  smd.nrep = argc - 1;
  smd.irep = smd.world_rank / (smd.world_size / smd.nrep);

  // Crash if there are extras
  assert(smd.world_size % smd.nrep == 0);

  // Split the communicator based on the irep and use the
  // original rank for ordering
  MPI_Comm_split(MPI_COMM_WORLD, smd.irep, smd.world_rank, &smd.comm);
  MPI_Comm_rank(smd.comm, &smd.myrank);
  MPI_Comm_size(smd.comm, &smd.nprocs);

  // Split the communicator transversally according to local rank
  MPI_Comm_split(MPI_COMM_WORLD, smd.myrank, smd.world_rank, &smd.comm_col);
  MPI_Comm_rank(smd.comm_col, &smd.myrank_col);
  MPI_Comm_size(smd.comm_col, &smd.nprocs_col);

  fprintf(stdout,"MPI: %3.3d\tOriginal Rank: %2d/%2d\tirep:%d \tComm Rank: %2d/%2d\tComm2 Rank: %2d/%2d\n",
          smd.world_rank, smd.world_rank, smd.world_size, smd.irep, smd.myrank, smd.nprocs, smd.myrank_col,smd.nprocs_col);

  FILE* in=stdin;
  if(argc>1) in=fopen(argv[smd.irep+1],"r");
  int r=smd.main(in,stdout);
  if(argc>1) fclose(in);

  MPI_Finalize();

  return r;
}


