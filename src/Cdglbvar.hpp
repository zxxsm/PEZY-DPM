/*
 * File: Cdglbvar.h
 * Author: Qinglin Wang
 * Email: wangqinglin.thu@gmail.com
 * Description: Parallization of double 
 */
#include "Cglbvar.hpp"
/**********************************************************************
 *	
 *				Redefinition of Data type
 *
 **********************************************************************/
// typedef double DUREAL;    //放在 "Cglbvar.hpp" 中了  
/* Type for random 32 or 64-bit real, e.g. double or float */
							   
/**********************************************************************
 * 
 *  Shared for all particles and threads.
 *  共享的数据类型
 *  进程级并行 ： 每个进程单独所有，需要每个进程独立初始化/读取后BECAST广播
 *  线程级并行 ： 所有线程共享所有；
 *  SCx异构并行： 所有线程共享所有；需要进行CPU->SCx的数据拷贝；
 *  s1-20：DPM 原始结构体
 *  s21：MPI参数，新增结构体
 **********************************************************************/
/*
 *  s1-/dpmsim/
 * 	dpmsim used in 
 * 	Function electr()
 * 	PROPERTY：Shared
 *  
 *  此结构体都是从dpm.in中读取
 *  maxhis：最大历史次数
 *  atime：只在控制函数command函数中调用
 */
typedef struct
{
	DUREAL atime;
	int32 maxhis;
}dTdpmsim;

/*	
 *  s2-/dpmsrc/
 * 	dpmsrc used in 
 * 	Function source()
 * 	PROPERTY：Shared
 * 
 *  此结构体都是从dpm.in中读取
 *  esrc:源能量（eV）-------->Source energy (eV)
 *  eabs:电子吸收能量-------->Electron absorption energy
 *  eabsph:光子吸收能量-------->Photon absorption energy
 *  param: 方形源束尺寸；用0表示笔形射束（厘米）-------->Square source beam size; use 0 for pencil beam (cm):
 *  intype:光子/电子类型-------->Particle Type (-1 for e-, 0 for photons)
 */	
typedef struct
{
	DUREAL esrc,eabs,eabsph,param;
	int32 intype;
}dTdpmsrc;

/*
 *  s3-/dpmmat/
 * 	dpmmat used in 
 * 	Function flight()
 * 	PROPERTY：shared
 * 
 *  在rmater中完成读取和计算计算
 */
typedef struct
{
	DUREAL matden[maxmat],zmass[maxmat],z2mass[maxmat];
	DUREAL wcc,wcb;
	int32 nmat;
}dTdpmmat;

/*
 *  s4-/dpmscsr/
 * 	dpmscsr used in 
 * 	Function stepip()
 * 	PROPERTY：shared
 * 
 *  在读取函数rstep中就已经算好了
 */
typedef struct
{
	DUREAL escsr[nscsr],scssp[nscsr],scsspa[nscsr];
	DUREAL scsspb[nscsr],scsspc[nscsr],scsspd[nscsr];
	DUREAL idless;
}dTdpmscsr;
/*
 * 	s5-/dpmscpw/ 
 *  dpmscpw used in 
 * 	Function scpwip(matid,e)
 * 	PROPERTY：Shared
 * 
 *  从文件中读取第一个 TMFP 并设置插值矩阵，rscpw函数读取
 */
typedef struct
{
	DUREAL escp[nscp], scpsp[nscp];
	DUREAL scpspa[maxmat*nscp],scpspb[maxmat*nscp];
	DUREAL idlesc;
}dTdpmscpw;
/*
 * 	s6-/xdpmbw/ 
 * 	xdpmbw used in 
 * 	xbwip()
 * 	PROPERTY：Shared
 * 
 *  从文件中读取筛选参数数据并设置插值矩阵，rbw函数读取
 */
typedef struct
{
	DUREAL ebw[nbw],bwsp[nbw],bwspa[nbw];
	DUREAL bwspb[nbw],bwspc[nbw],bwspd[nbw];
	DUREAL idlebw;
}dTxdpmbw;


/*
 *  s8-/dpmlam/
 * 	dpmlam used in 
 * 	lamoip()
 * 	PROPERTY：shared
 * 
 *  在读取函数rlammo中就已经算好了
 *  在lamoip中计算穆勒平均自由程
 */
typedef struct
{
	DUREAL elam[nlam],lamsp[nlam],lamspa[nlam];
	DUREAL lamspb[nlam],lamspc[nlam],lamspd[nlam];
	DUREAL idlela;
}dTdpmlam;
/*	
 *  s9-/dpmlbr/
 * 	dpmlbr used in 
 * 	Function ilabip(matid,ie)
 * 	PROPERTY：Shared
 * 
 *  从文件中读取 lambda_radiative 数据并设置插值矩阵，rlabre函数读取
 */
typedef struct
{
	DUREAL elra[nlra], lrasp[nlra];
	DUREAL lraspa[maxmat*nlra],lraspb[maxmat*nlra];
	DUREAL idlelr;
}dTdpmlbr;
/*
 *  s10-/dpmstpw/
 * 	dpmstpw used in 
 * 	Function rstpip(matid,e)
 * 	PROPERTY：Shared
 * 
 *  从文件中读取 StopPower_restr数据，并设置插值矩阵，rrstpw函数读取
 *  maxmat表示材料数
 */
typedef struct
{
	DUREAL est[nst], stsp[nst];
	DUREAL stspa[maxmat*nst],stspb[maxmat*nst];
	DUREAL idlest;
}dTdpmstpw;
/*
 *  s11-/xbre/
 * 	xbre used in 
 * 	dpmscsr SBREMS()
 * 	PROPERTY：Shared
 * 
 *  使用 PENELOPE 例程读取 brems 常数以进行采样，rbcon函数读取
 */
typedef struct
{
	DUREAL f0[maxmat*3];
	DUREAL bcb[maxmat];
}dTxbre;

/*
 *  s12-/dpmlph/
 * 	dpmlph used in 
 * 	itphip()
 * 	PROPERTY：shared
 * 
 *  从文件中读取光子总反向平均自由路径数据并设置插值矩阵，rlamph函数读取
 */
typedef struct
{
	DUREAL elaph[nlaph],lamph[nlaph];
	DUREAL lampha[maxmat*nlaph], lamphb[maxmat*nlaph];
	DUREAL lamphc[maxmat*nlaph], lamphd[maxmat*nlaph];
	DUREAL idleph;
}dTdpmlph;
/*
 *  s13-/dpmcmp/
 * 	dpmcmp used in 
 * 	icptip()
 * 	PROPERTY：shared
 * 
 *  从文件中读取康普顿反平均自由路径数据并设置插值矩阵，rcompt函数读取
 */
typedef struct
{
	DUREAL ecmpt[ncmpt],compt[ncmpt];
	DUREAL compta[maxmat*ncmpt], comptb[maxmat*ncmpt];
	DUREAL comptc[maxmat*ncmpt], comptd[maxmat*ncmpt];
	DUREAL idlecp;
}dTdpmcmp;
/*
 *  s14-/dpmpte/
 * 	dpmpte used in 
 * 	ipheip()
 * 	dpm计算中没有调用后续相关函数ipheip
 */

/*
 *  s15-/dpmpap/
 * 	dpmpap used in 
 * 	ipapip()
 * 	PROPERTY：shared
 * 
 *  从输入文件读取体素几何图形，rpairp函数读取
 */
typedef struct
{
	DUREAL epair[npair],pairp[npair];
	DUREAL pairpa[maxmat*npair], pairpb[maxmat*npair];
	DUREAL pairpc[maxmat*npair], pairpd[maxmat*npair];
	DUREAL idlepp;
}dTdpmpap;
/*
 *  s16-/dpmvox/
 * 	dpmvox used in 
 * 	Function source()
 * 	PROPERTY：Shared
 * 
 *  dx,dy,dz：x,y,z的间隔距离
 *  xmid,ymid,zmid：整体的中心点
 *  dens：材料密度
 *  mat:材料编号
 *  Unxvox,Unyvox,Unzvox：体素个数
 */	
typedef struct
{
	DUREAL ddx,ddy,ddz,xmid,ymid,zmid;
	DUREAL dens[nxyz];
	int32 mat[nxyz];
	int32 Unxvox,Unyvox,Unzvox;
}dTdpmvox;

/*
 * 	s18-/cutoff/
 *  cutoff used in 
 * 	sammo()
 * 	PROPERTY：Shared
 * 
 *  在iniion 和 inibre中生成
 */
typedef struct
{
	DUREAL wcion;
	DUREAL wcbre;
}dTcutoff;

/*
 *  s19-/dpmwck/
 * 	dpmwck used in 
 * 	lamwck()
 * 	PROPERTY：private
 * 
 *  在iniwck中生成
 */
typedef struct
{
	DUREAL a0wck[nwck],a1wck[nwck];
	DUREAL idlewk, wcke0;
}dTdpmwck;
/*  
 *  s20-/dpmsub/
 * 	dpmsub used in 
 * 	Function subabs()
 * 	PROPERTY：shared
 * 
 *  在inisub中生成
 */
typedef struct
{
	DUREAL subden,subfac,substp;
}dTdpmsub;
/*
 *  t2-/dpmesc/
 * 	dpmesc used in 
 * 	Function score()
 * 	PROPERTY：Private but nead to reduction (Maybe).
 * 
 *  记录1个原生粒子所派生的 所有粒子的计数累加和；
 *  在init中初始化；得分信息/需要规约;
 *  本身就是512位对其 
 * 
 *  SCx异构并行：xxxxxxx
 */
typedef struct
{
	//614400*8*3+614400*4= 16.4MB
	//如果是4096个线程都要有需要65.6GB的内存大小；
	//超过了pezy的栈空间
	//只能用memcopy开辟
	DUREAL escore[nxyz],escor2[nxyz],etmp[nxyz];
	int32 lasthi[nxyz];
}dTdpmesc;