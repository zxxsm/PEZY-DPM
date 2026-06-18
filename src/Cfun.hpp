/*
 * File: Cfun.hpp
 * Author: Zhuxiaoxong
 * Email: zhuxiaoxiong_cn@163.com
 */
#ifndef _CFUN_H
#define _CFUN_H

/**********************************************************************
 *
 *  Function: main simulation.
 *
 **********************************************************************/
void simulate(dTdpmsim *dpmsim, Tdpmvox *dpmvox, Tdpmsrc *dpmsrc,
			  dTdpmesc *dpmesc, Trdpmwck *dpmwck, Trdpmcmp *dpmcmp,
			  Tdpmrpt *dpmrpt, Trdpmlph *dpmlph, Trdpmpap *dpmpap,
			  Tdpmmat *dpmmat, Tdpmsub *dpmsub, Trdpmscsr *dpmscsr,
			  Trdpmlam *dpmlam, Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw,
			  Trdpmlbr *dpmlbr, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,
			  Tcutoff *cutoff, Txbre *xbre, Tsd *sd0, Trseed *rseed0,
			  Tmpi_para *mpi_para, double *tlaps);

/**********************************************************************
 * 
 *  数据进行重构，
 *
 **********************************************************************/			  
void reconstruct(Tdpmsim *&fdpmsim, dTdpmsim *dpmsim,
				 Tdpmsrc *&fdpmsrc, dTdpmsrc *dpmsrc,
				 Tdpmmat *&fdpmmat, dTdpmmat *dpmmat,
				 Txbre *&fxbre, dTxbre *xbre,
				 Tdpmvox *&fdpmvox, dTdpmvox *dpmvox,
				 Tcutoff *&fcutoff, dTcutoff *cutoff,
				 Tdpmsub *&fdpmsub, dTdpmsub *dpmsub,
				 Trdpmscsr *&rdpmscsr, dTdpmscsr *dpmscsr,
				 Trdpmscpw *&rdpmscpw, dTdpmscpw *dpmscpw,
				 Trxdpmbw *&rxdpmbw, dTxdpmbw *xdpmbw,
				 Trdpmlam *&rdpmlam, dTdpmlam *dpmlam,
				 Trdpmlbr *&rdpmlbr, dTdpmlbr *dpmlbr,
				 Trdpmstpw *&rdpmstpw, dTdpmstpw *dpmstpw,
				 Trdpmlph *&rdpmlph, dTdpmlph *dpmlph,
				 Trdpmwck *&rdpmwck, dTdpmwck *dpmwck,
				 Trdpmpap *&rdpmpap, dTdpmpap *dpmpap,
				 Trdpmcmp *&rdpmcmp, dTdpmcmp *dpmcmp);
/**********************************************************************
 *
 *  Function: this is an adapted version of subroutine ranecu written
 *		by f. james (comput. phys. commun. 60 (1990) 329-344),
 *		which has been modified to give a single random number at each
 *		call.
 *		the 'seeds' iseed1 and iseed2 must be initialized in the
 *		main program and transferred through the named common block
 *		/rseed/.
 *
 **********************************************************************/
UREAL Crng(Trseed *rseed);

/**********************************************************************
 *
 *  Function: Creates a new primary particle state
 *
 **********************************************************************/
void Csource(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmsrc *dpmsrc,
			 Tcgeom3 *cgeom3, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Locates the particle and sets vox in 3D
 *	Input:	particle position
 *	Output: absvox -> absolute address of the current voxel
 *
 **********************************************************************/
void Cwhere(Tdpmpart *dpmpart, Tcgeom3 *cgeom3, Tdpmvox *dpmvox);

/**********************************************************************
 *
 *  Function: Stores the inverse of vz to save time while in cango().
 *		It should be called every time {vx,vy,vz} changes.
 *
 **********************************************************************/
void Csetv(Tdpmpart *dpmpart, Tcgeom3 *cgeom3);

/**********************************************************************
 *
 *  Function: Informs about interface crossings in a 3D voxel
 *		geometry and prepares information to move to next voxel
 *	Input:
 *      {x,y,z} -> position vector (cm)
 *      {vx,vy,vz} -> normalized direction vector
 *      vox -> voxel where the particle starts
 *    Output:
 *      -> distance to intersection with nearest voxel wall (cm)
 *      index -> {1,2,3} depending on which walls are intersected
 *      dvox -> {+1,-1} depending on whether part. goes forw-backw
 *
 **********************************************************************/
UREAL Cinters(Tdpmpart *dpmpart, Tcgeom3 *cgeom3, Tcgo *cgo);

/**********************************************************************
 *
 *  Function: Restricted stopping power --linear interpolation
 *	Input:
 *      matid -> material id#
 *      e -> energy in eV  --kinetic energy--
 *    Output:
 *      StopPow in eV*cm^2/g
 *
 **********************************************************************/
UREAL Crstpip(int32 Lmatid, UREAL e, Trdpmstpw *dpmstpw);

/**********************************************************************
 *
 *  Function: Inverse 1st transport MFP --linear interpolation
 *	Input:
 *      matid -> material id#
 *      e -> energy in eV  --kinetic energy--
 *    Output:
 *      lambda_1^{-1} in cm^2/g
 *
 **********************************************************************/
UREAL Cscpwip(int32 Lmatid, UREAL e, Trdpmscpw *dpmscpw);

/**********************************************************************
 *
 *  Function: Bremsstrahlung inverse MFP --linear interpolation
 *	Input:
 *      matid -> material id number
 *      ie -> -1/energy in eV^-1  --kinetic energy--
 *    Output:
 *      inverse mean free path in cm^2/g
 *
 **********************************************************************/
UREAL Cilabip(int32 Lmatid, UREAL ie, Trdpmlbr *dpmlbr);

/**********************************************************************
 *
 *  Function: Deposites energy in the corresponding counters
 *	Input:
 *      edep -> energy being deposited (eV)
 *  Worthing Vectorizing ???
 **********************************************************************/
void Cscore(UREAL edep, Tdpmpart *dpmpart, Tdpmesc_new *dpmesc,
			Tdpmrpt *dpmrpt);

/**********************************************************************
 *
 *  Function: Transport of e- below the nominal Eabs until absorption
 *	Input:
 *      {x,y,z} -> take off location
 *		{vx,vy,vz} -> direction of flight
 *		energy -> initial kinetic energy
 *    Output:                                                      *
 *      escore -> energy deposited by CSDA
 *
 **********************************************************************/
void Csubabs(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmsub *dpmsub,
			 Tdpmesc_new *dpmesc_new, Tcgo *cgo, Tdpmrpt *dpmrpt, Tcgeom3 *cgeom3);

/**********************************************************************
 *
 *  Function: Changes voxel according to the information
 *		passed by inters()
 *	Input:
 *      vox -> voxel where the particle starts
 *		index -> {1,2,3} depending on which walls are intersected
 *		dvox -> {+1,-1} depending on whether part. goes forw-backw
 *    Output:
 *      vox -> voxel where the final position lies
 *
 **********************************************************************/
void Cchvox(Tdpmpart *dpmpart, Tcgo *cgo, Tdpmvox *dpmvox);

/**********************************************************************
 *
 *  Function: Transports the particle following a rectiliniar
 *		trajectory, taking care of interface crossings and keeping
 * 		track of energy losses in the corresponding counters
 *		(using CSDA).
 *	Input:
 *      e- initial state
 *      Fuel variables affecting the flight
 *      event -> event that stopped flight last time; in addition
 *                to the output codes,
 *                -1 new particle
 *                20 end of elastic step
 *    Output:
 *      e- final state
 *      Remaining fuel vars
 *      event -> kind of event that causes flight to stop:
 *                1 run out of energy; absorbed
 *                2 run out of elastic fuel
 *                3 run out of Moller fuel
 *                4 run out of bremsstrahlung fuel
 *               99 escaped from universe
 *
 ***********************************************************************/
void Cflight(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmmat *dpmmat,
			 Tdpmjmp *dpmjmp, Tdpmsrc *dpmsrc, Tdpmsub *dpmsub, Tcgo *cgo,
			 Tcgeom3 *cgeom3, Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw,
			 Trdpmlbr *dpmlbr, Tdpmesc_new *dpmesc_new, Tdpmrpt *dpmrpt);

/**********************************************************************
 *
 *  Function: 3spline interpolation for scattering strength as a
 *		function of kinetic energy; this quantity is related to
 *		the step length
 *	Input:
 *      e -> kinetic energy in eV
 *  Output:
 *      K = scattering strength = integ{ds/lambda1(s)} k
 *
 ***********************************************************************/
UREAL Cstepip(UREAL e, Trdpmscsr *dpmscsr);

/**********************************************************************
 *
 *  Function: Moller mean free path, 3-spline interpolation
 *	Input:
 *      ie -> -1/energy in eV^-1  --kinetic energy--
 *  Output:
 *      mean free path in g/cm^2
 *
 ***********************************************************************/
UREAL Clamoip(UREAL ie, Trdpmlam *dpmlam);

/**********************************************************************
 *
 *  Function: 3spline interpolation for bw as a function energy
 *	Input:
 *      ie -> -1/energy in eV^-1  --kinetic energy--
 *  Output:
 *      bw, broad screening parameter that gets flattest q surf
 *
 ***********************************************************************/
UREAL Cxbwip(UREAL ie, Trxdpmbw *xdpmbw);

/**********************************************************************
 *
 *  Function: Linearly interpolated q(u;energy) surface
 *	Input:
 *      u -> angular variable                                      *
 *      ie -> -1/energy in eV^-1  --kinetic energy--
 *  Output:
 *      bw, broad screening parameter that gets flattest q surf
 *
 ***********************************************************************/
UREAL Cxq2Dip(UREAL u, UREAL ie, Txdpmq2D *xdpmq2D);

/**********************************************************************
 *
 *  Function: Samples cos(theta) according to the G&S distribution.
 *		Uses interpolated data for bw and the q surface and this
 *		latter quantity to perform a rejection procedure.
 *	Input:
 *      e -> kinetic energy in eV
 *  Output:
 *      mu -> polar angle -cos(theta)-
 *
 ***********************************************************************/
void Csamsca(UREAL e, UREAL *mu, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Rotates a vector; the rotation is specified by giving
 *		the polar and azimuthal angles in the "self-frame", as
 *		determined by the vector to be rotated.
 *	Input:
 *      (u,v,w) -> input vector (=d) in the lab. frame
 *		costh -> cos(theta), angle between d before and after turn
 *      phi -> azimuthal angle (rad) turned by d in its self-frame
 *  Output:
 *       (u,v,w) -> rotated vector components in the lab. frame
 *
 ***********************************************************************/
void Crotate(UREAL *u, UREAL *v, UREAL *w, UREAL costh, UREAL phi);

/**********************************************************************
 *
 *  Function: Samples an energy loss according to Moller DCS
 *	Input:
 *      energy -> kinetic energy in eV
 *  Output:
 *      eloss/energy, energy loss fraction
 *
 ***********************************************************************/
UREAL Csammo(UREAL Lenergy, Tcutoff *cutoff, Trseed *rseed);

/**********************************************************************
 *
 *  Function: gives the outgoing angle of the delta ray after a Moller
 *    interaction
 *	Input:
 *      eloss -> energy loss --also kinetic energy of the recoil e-
 *		energy -> kinetic energy in eV
 *	Output:
 *		cos(theta') in the lab frame
 *
 ***********************************************************************/
UREAL Csecmo(UREAL eloss, UREAL Lenergy);

/**********************************************************************
 *
 *  Function: Creates a new secondary electron from a Moller
 *		interaction and stores its state in the secondary stack
 *	Input:
 *      elost -> kinetic energy of the secondary electron (eV)
 *
 ***********************************************************************/
void Cputmol(UREAL elost, Tdpmpart *dpmpart, Trdpmstck *dpmstck, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Canibalized from PENELOPE too.
 *		screening functions f1(b) and f2(b) in the bethe-heitler
 *		differential cross section for bremsstrahlung emission.

 *	Input:
 *      elost -> kinetic energy of the secondary electron (eV)
 *
 ***********************************************************************/
void Cschiff(UREAL *b, UREAL *f1, UREAL *f2);

/**********************************************************************
 *
 *  Function: New routine to replace sambre in dpm  for more accurate
 *		sampling of Bremsstrahlung energies.
 *		Canibalized from PENELOPE.
 *	Input:
 *      elost -> kinetic energy of the secondary electron (eV)
 *
 ***********************************************************************/
UREAL Csbrems(UREAL e, int32 m, Tcutoff *cutoff, Txbre *xbre, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Creates a new secondary bremsstrahlung photon and
 *		stores its state in the secondary stack
 *	Input:
 *      elost -> energy of the secondary photon (eV)
 *
 ***********************************************************************/
void Cputbre(UREAL elost, Tdpmpart *dpmpart, Trdpmstck *dpmstck, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Transports an electron until it either escapes from the
 *		universe or its energy drops below Eabs
 *	Input:
 *      electron initial state
 *    Output:
 *      deposits energy in counters and updates secondary stack
 *
 ***********************************************************************/
void Celectr(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmmat *dpmmat,
			 Tdpmjmp *dpmjmp, Tdpmsrc *dpmsrc, Tdpmsub *dpmsub, dTdpmsim *dpmsim,
			 Trdpmscsr *dpmscsr, Trdpmlam *dpmlam, Tcgo *cgo, Tcgeom3 *cgeom3,
			 Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw, Trdpmlbr *dpmlbr,
			 Tdpmesc_new *dpmesc_new, Tdpmrpt *dpmrpt, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,
			 Trdpmstck *dpmstck, Tcutoff *cutoff, Txbre *xbre, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Creates secondaries for annhilation photons, stores
 *		states in the secondary stack
 *	Input:
 *
 ***********************************************************************/
void Cputann(Tdpmpart *dpmpart, Trdpmstck *dpmstck, Tdpmsrc *dpmsrc, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Mean free path prepared to play the Woodcock trick
 *	Input:
 *		e -> kinetic energy in eV
 *	Output:
 *		Minimum mean free path in cm
 *
 ***********************************************************************/
UREAL Clamwck(UREAL e, Trdpmwck *dpmwck);

/**********************************************************************
 *
 *  Function: Inverse Compton mean free path --3spline interpolation
 *	Input:
 *		matid -> material id#
 *		e -> kinetic energy in eV
 *	Output:
 *		Inverse total mean free path in cm^2/g
 *
 ***********************************************************************/
UREAL Cicptip(int32 Lmatid, UREAL e, Trdpmcmp *dpmcmp);

/**********************************************************************
 *
 *  Function: Photon total inverse mean free path --3spline
 *		interpolation
 *	Input:
 *		matid -> material id#
 *		e -> kinetic energy in eV
 *	Output:
 *		Total inverse mean free path in cm^2/g
 *
 ***********************************************************************/
UREAL Citphip(int32 Lmatid, UREAL e, Trdpmlph *dpmlph);

/**********************************************************************
 *
 *  Function: Compton angular deviation of the secondary electron.
 *	Input:
 *		energy -> photon energy in eV
 *		efrac -> fraction of initial energy kept by 2nd photon
 *		costhe-> photon scattering angle
 *  Output:                                                      *
 *  	-> cos(theta) of the 2nd electron
 *
 ***********************************************************************/
UREAL Ccomele(UREAL Lenergy, UREAL efrac, UREAL costhe);

/**********************************************************************
 *
 *  Function: Creates a secondary electron from a Compton interaction
 *		and stores its state in the secondary stack.
 *	Input:
 *		elost -> energy of the secondary electron being created
 *		efrac -> fraction of initial energy kept by 2nd photon
 *		phi -> athimutal angle
 *		costhe -> photon scattering angle
 *
 ***********************************************************************/
void Cputcom(UREAL elost, UREAL efrac, UREAL phi, UREAL costhe,
			 Tdpmpart *dpmpart, Trdpmstck *dpmstck);

/**********************************************************************
 *
 *  Function: Inverse pair production mean free path --3spline interpol
 *	Input:
 *		matid -> material id#
 *		e -> kinetic energy in eV
 *  Output:
 *  	Inverse total mean free path in cm^2/g
 *
 ***********************************************************************/
UREAL Cipapip(int32 Lmatid, UREAL e, Trdpmpap *dpmpap);

/**********************************************************************
 *
 *  Function: Samples a Compton event following Klein-Nishina DCS
 *	Input:
 *		energy -> photon energy in eV
 *  Output:
 *  	efrac -> fraction of initial energy kept by 2nd photon
 *		costhe -> cos(theta) of the 2nd photon
 *
 ***********************************************************************/
void Ccomsam(UREAL Lenergy, UREAL *efrac, UREAL *costhe, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Transports a photon until it either escapes from the
 *		universe or its energy drops below EabsPhoton
 *	Input:
 *		photon initial state
 *	Output:
 *		deposits energy in counters and updates secondary stack
 *
 ***********************************************************************/
void Cphoton(Tdpmpart *dpmpart, Trdpmstck *dpmstck, Tdpmvox *dpmvox,
			 Tdpmsrc *dpmsrc, Tdpmesc_new *dpmesc_new, Trdpmwck *dpmwck, Trdpmcmp *dpmcmp,
			 Tdpmrpt *dpmrpt, Tcgeom3 *cgeom3, Trdpmlph *dpmlph, Trdpmpap *dpmpap, Trseed *rseed);

/**********************************************************************
 *
 *  Function: Retrieves a particle from the secondary stack and fills
 *		up the current particle common with its state
 *  Output:
 *  	-> 0 if there was no particle in the stack, 1 else
 *
 ***********************************************************************/
int32 Cscndry(Tdpmpart *dpmpart, Trdpmstck *dpmstck);

/**********************************************************************
 *
 *  Function: Dumps all tmp counters into the corresponding final
 *		counters
 *  Output:
 *
 ***********************************************************************/
void Cdumpe(Tdpmesc *dpmesc);


void Cdumpe_new(Trdpmesc *dpmesc,Tdpmesc_new *dpmesc_new);
/**********************************************************************
 *
 *  pzc_simulate 初始化
 *  启动设备，并将所有的全局数据传输到设备端
 *
 **********************************************************************/
void pzc_simulate(Tdpmsim *dpmsim, Tdpmvox *dpmvox, Tdpmsrc *dpmsrc,
                  dTdpmesc *dpmesc, Trdpmwck *dpmwck, Trdpmcmp *dpmcmp,
                  Tdpmrpt *dpmrpt, Trdpmlph *dpmlph, Trdpmpap *dpmpap,
                  Tdpmmat *dpmmat, Tdpmsub *dpmsub, Trdpmscsr *dpmscsr,
                  Trdpmlam *dpmlam, Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw,
                  Trdpmlbr *dpmlbr, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,
                  Tcutoff *cutoff, Txbre *xbre, Tsd *sd0, Trseed *rseed0,
                  Tmpi_para *mpi_para, double *tlaps);

/**********************************************************************
 *
 *  Function: 计算每个进程分配的粒子数
 *
 *  Output:
 *
 ***********************************************************************/
int32 Ccomp_np(int32 nprtcl_all, int32 gloabl_id, int32 gloabl_num);

/**********************************************************************
 *
 *   计时函数
 *
 ***********************************************************************/
void timers(double *et);
#endif
//_CFUN_H