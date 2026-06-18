/*
 * File: Cglbvar.hpp
 * Author: Zhuxiaoxong
 * Email: zhuxiaoxiong_cn@163.com
 */
#ifndef _GBL_H 
#define _GBL_H

#include <math.h>


/**********************************************************************
 *	
 *				Redefinition of Data type
 *
 **********************************************************************/
typedef int  int32;
#ifdef DPM_USE_DOUBLE
typedef double UREAL;
#else
typedef float UREAL;
#endif
typedef double DUREAL; 							   
/**********************************************************************
 *	
 *				MACRO Definition 
 *
 **********************************************************************/
//定义SCx与CPU的计算性能比
//1表示只使用SCx
#define SCx_ratio 1  

//总共使用的最大种子数
//目前使用的是10e5的rseed文件
#define nsd 100000

//使用OPENMP
#ifdef OMP_P
//THREAD_NUM可以调整
#define THREAD_NUM 64
#endif

#define Cmax(i,j) ((i)>(j)?(i):(j))

#define nscsr 512 // nscsr=2**9
#define nlam 256 //nlam=2**8
#define nbw 512 //nbw=2**9
#define nxvox 64  //2**6
#define nyvox 64  //2**6
#define nzvox 150
#define nxyz 614400 //nxyz=nxvox*nyvox*nzvox
#define nyzvox 9600 //nyzvox=nyvox*nzvox
#define maxmat 5
#define nst 1024 //nst=2**10
#define nscp 2048 //nscp=2**11
#define nlra 512 //nlra=2**9
#define nuq 128 //nuq=2**7
#define neq 256 //neq=2**8
#define ncmpt 256 //ncmpt=2**8
#define maxns 128 //maxns=2**7
#define nwck 4096 //nwck=2**12
#define nlaph 2048 //nlaph=2**11
#define npair 1024 //npair=2**10

// #define zero 1.0E-30
// #define inf 1.0E+30
// #define pi 3.1415926535897932E0 //pi=3.1415926535897932d0
// #define SZERO 1.0E-14 //SZERO=1.0d-14
// #define ZERO 1.0E-90 //ZERO=1.0d-90
// #define mc2 510.9991E+03 //mc2=510.9991d3
// #define kcmax 0.4999999999E0 //kcmax=0.4999999999d0
// #define rev 5.1099906E+05 //REV=5.1099906D5  */

// #define twopi 6.2831853071795864E0 //twopi=2.0d0*pi
// #define mc2sq2 3.613309286948971E+05 //mc2sq2=7.07106781d-1*mc2
// #define imc2 1.956950609110662E-06 //imc2=1.d0/mc2
// #define twomc2 1.0219982E+06//twomc2=2.0d0*mc2 
// #define USCALE 4.656612873077393E-10 //1.0E0/2.0E0**31 //USCALE=1.0D0/2.0D0**31 

// float 自动转换
const UREAL zero=   1.0E-14; //1.0E-30
const UREAL inf=    1.0E+14;   //1.0E+30
const UREAL pi=     3.1415927E0; //pi=3.1415926535897932d0
const UREAL SZERO=  1.0E-7; //SZERO=1.0d-14
const UREAL ZERO=   1.0E-38; //ZERO=1.0d-90
const UREAL mc2=    510.9991E+03; //mc2=510.9991d3
const UREAL kcmax=  0.4999999E0; //kcmax=0.4999999999d0
const UREAL rev =   5.1099906E+05; //REV=5.1099906D5  */
const UREAL twopi=  6.2831853E0; //twopi=2.0d0*pi
const UREAL mc2sq2= 3.6133093E+05; //mc2sq2=7.07106781d-1*mc2
const UREAL imc2=   1.9569506E-06; //imc2=1.d0/mc2
const UREAL twomc2= 1.0219982E+06;//twomc2=2.0d0*mc2 
const UREAL USCALE= 4.6566129E-10; //1.0E0/2.0E0**31 //USCALE=1.0D0/2.0D0**31 
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
	UREAL atime;
	int32 maxhis;
}Tdpmsim;

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
	UREAL esrc,eabs,eabsph,param;
	int32 intype;
}Tdpmsrc;

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
	UREAL matden[maxmat],zmass[maxmat],z2mass[maxmat];
	UREAL wcc,wcb;
	int32 nmat;
}Tdpmmat;

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
	UREAL escsr[nscsr],scssp[nscsr],scsspa[nscsr];
	UREAL scsspb[nscsr],scsspc[nscsr],scsspd[nscsr];
	UREAL idless;
}Tdpmscsr;
typedef struct
{
	UREAL scsspx[4*nscsr];
	UREAL escsr, idless;
}Trdpmscsr;
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
	UREAL escp[nscp], scpsp[nscp];
	UREAL scpspa[maxmat*nscp],scpspb[maxmat*nscp];
	UREAL idlesc;
}Tdpmscpw;
/* Re-struct */
typedef struct
{
	UREAL scpspx[2*maxmat*nscp];
	UREAL escp, idlesc;
}Trdpmscpw; 
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
	UREAL ebw[nbw],bwsp[nbw],bwspa[nbw];
	UREAL bwspb[nbw],bwspc[nbw],bwspd[nbw];
	UREAL idlebw;
}Txdpmbw;
/* Re-struct */
typedef struct
{
	UREAL bwspx[4*nbw];
	UREAL ebw, idlebw;
}Trxdpmbw;
/*
 *  s7-/xdpmq2D/
 * 	xdpmq2D used in 
 * 	xq2Dip()
 * 	PROPERTY：Shared

 *  从文件中读取 q 表面数据并设置插值矩阵，rqsurf函数读取
 */
typedef struct
{
	UREAL q[nuq*neq];
	UREAL le0q,idleq,iduq;
}Txdpmq2D;

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
	UREAL elam[nlam],lamsp[nlam],lamspa[nlam];
	UREAL lamspb[nlam],lamspc[nlam],lamspd[nlam];
	UREAL idlela;
}Tdpmlam;
/* Re-struct */
typedef struct
{
	UREAL lamspx[4*nlam];
	UREAL elam, idlela;
}Trdpmlam;
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
	UREAL elra[nlra], lrasp[nlra];
	UREAL lraspa[maxmat*nlra],lraspb[maxmat*nlra];
	UREAL idlelr;
}Tdpmlbr;
/* Re-struct */
typedef struct
{
	UREAL lraspx[2*maxmat*nlra];
	UREAL elra, idlelr;
}Trdpmlbr; 

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
	UREAL est[nst], stsp[nst];
	UREAL stspa[maxmat*nst],stspb[maxmat*nst];
	UREAL idlest;
}Tdpmstpw;
/* Re-struct */
typedef struct
{
	UREAL stspx[2*maxmat*nst];
	UREAL est, idlest;
}Trdpmstpw; 

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
	UREAL f0[maxmat*3];
	UREAL bcb[maxmat];
}Txbre;

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
	UREAL elaph[nlaph],lamph[nlaph];
	UREAL lampha[maxmat*nlaph], lamphb[maxmat*nlaph];
	UREAL lamphc[maxmat*nlaph], lamphd[maxmat*nlaph];
	UREAL idleph;
}Tdpmlph;
/* Re-struct */
typedef struct
{
	UREAL lamphx[4*maxmat*nlaph];
	UREAL elaph, idleph;
}Trdpmlph; 

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
	UREAL ecmpt[ncmpt],compt[ncmpt];
	UREAL compta[maxmat*ncmpt], comptb[maxmat*ncmpt];
	UREAL comptc[maxmat*ncmpt], comptd[maxmat*ncmpt];
	UREAL idlecp;
}Tdpmcmp;
/* Re-struct */
typedef struct
{
	UREAL comptx[4*maxmat*ncmpt];
	UREAL ecmpt, idlecp;
}Trdpmcmp; 
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
	UREAL epair[npair],pairp[npair];
	UREAL pairpa[maxmat*npair], pairpb[maxmat*npair];
	UREAL pairpc[maxmat*npair], pairpd[maxmat*npair];
	UREAL idlepp;
}Tdpmpap;
/* Re-struct */
typedef struct
{
	UREAL pairpx[4*maxmat*npair];
	UREAL epair, idlepp;
}Trdpmpap;

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
	UREAL ddx,ddy,ddz,xmid,ymid,zmid;
	UREAL dens[nxyz];
	int32 mat[nxyz];
	int32 Unxvox,Unyvox,Unzvox;
}Tdpmvox;

/*
 *  s17-/dpmrpt/
 * 	dpmrpt used in 
 * 	Function score()
 * 	PROPERTY：shared
 * 
 *  记录的是Rectangular Region of Interest的像素
 *  记录的是  计分的长方形区域像素  
 * 	PROPERTY：shared
 */
typedef struct
{
	int32 nxini,nxfin,nyini,nyfin,nzini,nzfin;
}Tdpmrpt;

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
	UREAL wcion;
	UREAL wcbre;
}Tcutoff;

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
	UREAL a0wck[nwck],a1wck[nwck];
	UREAL idlewk, wcke0;
}Tdpmwck;
/* Re-struct */
typedef struct
{
	UREAL awck[2*nwck];
	UREAL idlewk, wcke0;
}Trdpmwck; 

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
	UREAL subden,subfac,substp;
}Tdpmsub;

/*
 *  s21-/mpi_para/
 * 	mpi_para used in 
 * 	Function simulate()/sim()
 * 	PROPERTY：Shared
 */	
typedef struct
{
	int32 myrank,pro_num,ierr;
}Tmpi_para;
/**********************************************************************
 * 
 *  Private for all threads.每个线程单独所有；
 *  进程/线程私有数据类型
 *  进程级并行 ： 每个进程单独所有；
 *  线程级并行 ： 每个线程单独所有；
 *  SCx异构并行： 每个结构体单独设计；
 **********************************************************************/
/*  
 *  t1-/dpmstck/
 * 	dpmstck used in 
 * 	putmol()：通过莫勒相互作用产生一个新的二次电子，并将其状态存储在二次堆栈中
 *  putbre()：创建一个新的二次轫致辐射光子，并将其状态存储在二次堆栈中
 *  photon()：传输光子，直到它逃离宇宙或其能量降到 EabsPhoton 以下为止，电子对效应中创建二级电子
 *  putcom()：通过康普顿相互作用产生一个次级电子，并将其状态存储在二级堆栈中
 *  putann()：为湮灭光子创建二级粒子，将其状态存储在二级堆栈中
 *  scndry()：从二级堆栈中检索一个粒子，并将其状态填充到当前粒子公共堆栈中
 * 	PROPERTY：private 
 *  
 *  记录1个原生粒子所派生的所有粒子的基本信息
 *  sptype 派生粒子类型 -1电子/0光子
 *  senerg 派生粒子能量
 *  svx/svy/svz
 *  maxns表示每个 初始粒子 能够派生的最大个数。
 *  
 *  SCx异构并行：每个线程单独所有
 */

// typedef struct
// {	
// 	//10760B
// 	//64byte 对齐需要15个int
// 	//超过了pezy的栈空间
// 	//只能用memcopy开辟
// 	UREAL swght[maxns], senerg[maxns], svx[maxns], svy[maxns];
// 	UREAL svz[maxns], sx[maxns], sy[maxns], sz[maxns];
// 	int32 sptype[maxns], sxvox[maxns], syvox[maxns], szvox[maxns];
// 	int32 sabsvx[maxns]; 
// 	int32 nnsec; 

// 	int32 isup[15];//64byte对齐
// }Tdpmstck;

/* Re-struct */
typedef struct
{
	UREAL fltdata[8*maxns];
	int32 intdata[5*maxns];
	int32 nsec;

	int32 isup[15];//64byte对齐 
}Trdpmstck; 

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
	UREAL escore[nxyz],escor2[nxyz],etmp[nxyz];
	int32 lasthi[nxyz];
}Tdpmesc;

/* Re-struct */
typedef struct
{
	DUREAL escorex[nxyz*2];
}Trdpmesc;
/*
 * 	dpmesc_new used in 
 * 	Function score_new()
 * 	PROPERTY：Private but nead to reduction (Maybe).
 *  
 *  记录1个原生粒子所派生的 所有粒子的计数累加和；
 *  在init中初始化；得分信息/需要规约;
 *  本身就是512位对其 
 * 
 *  SCx异构并行：xxxxxxx
 *  由zxx新定义
 *  estmp用于记录etmp值
 *  esidx用于记录对应的absvox
 *  tail用于记录尾部
 *  单个总大小为1.5KB,可以放在SCx的栈中（把栈调到3KB）
 */

typedef struct
{
	DUREAL estmp[maxns];
	int32 esidx[maxns];
	int32 tail;
}Tdpmesc_new;
/*
 *  t3-/RSEED/
 * 	rseed used in 
 * 	rng()
 * 	PROPERTY：single
 * 
 *  程序原始初始的随机种子， 在并行版本中被 分段随机种子替代；
 *  同时用于记录随机种子结构体，SCx异构/多线程并行中，
 *  从随机种子池中读取种子
 */
typedef struct
{
	int32 IISEED1,IISEED2;
	int32 IICNT;
}Trseed;

/*
 *  t4-/sd/
 * 	rseed used in 
 * 	rng()
 * 	PROPERTY：single
 * 
 *  新的分段随机种子数组，nsd表示种子分段数量
 *  可以理解为是一个种子池
 */
typedef struct
{
	int32 Iseed[2*nsd];
}Tsd; //1-dB

/**********************************************************************
 * 
 *  Private for all particles in each threads.每个粒子单独所有
 *  进程/线程私有数据类型
 *  进程级并行 ： 每个进程单独所有；
 *  线程级并行 ： 每个线程单独所有；
 *  SCx异构并行： 每个结构体单独设计； 
 **********************************************************************/
/*
 * 	p1-/dpmpart/
 *  dpmpart used in 
 * 	Function source()
 * 	PROPERTY：Private for each particle
 * 
 *  记录当前计算粒子的信息
 *  ptype -> -1（如果是电子），0（如果是光子）
 *  energy -> kinetic energy 能量 -> 动能
 *  {vx,vy,vz} -> 飞行方向
 *  {x,y,z} -> 位置
 *  xvox,yvox,zvox, -> 体素坐标索引
 *  absvox 将xvox,yvox,zvox转换成1维数组索引
 *  nhist 第n个粒子
 *  
 *  SCx异构并行：每个线程单独所有
 */
typedef struct
{
	UREAL energy,vx,vy,vz,x,y,z; //7
	int32 xvox,yvox,zvox,absvox,nnhist; //5
	int32 pptype; //1
}Tdpmpart;

/*
 * 	p2-/cgeom3/
 * 	cgeom3 used in 
 * 	Function where()
 * 	PROPERTY：Private for each particle
 * 
 *  Maybe can't shared owe to that setv() must be called 
 *	previously every time v changes.
 *  在inigeo中初始化，dx/dy/dz与idx/idy/idz，idx/idy/idz 记录的是dx/dy/dz的倒数
 *  后续计算中，只有ivx,ivy,ivz发生变化
 *  ivx/ivy/ivz 记录的是vx/vy/vz的倒数
 *  
 *  SCx异构并行：每个线程单独所有
 */	
typedef struct
{
	UREAL dxx,dyy,dzz,idxx,idyy,idzz,ivx,ivy,ivz; //9
}Tcgeom3;

/*  
 *  p3-/dpmjmp/
 * 	dpmjmp used in 
 * 	Function flight()
 * 	PROPERTY：Private
 * 
 *  记录电子发生作用的信息
 *  fuelel 弹性燃料
 *  fuelxt 发生作用带来的燃料损失
 *  fuelmo 穆勒燃料
 *  fuelbr 轫致辐射燃料 
 *  burnel 弹性燃烧率
 *  burnmo 穆勒燃烧率
 *  burnbr 轫致辐射燃烧率
 *  smax   将要相交/穿过体素面的距离（厘米） 
 *  voxden 体素密度
 *  ebefor 初始能量
 *  matid 材料比编号
 *  event 事件
 *  modeel 模式
 * 
 *  SCx异构并行：每个线程单独所有
 */
typedef struct
{
	UREAL fuelel,fuelmo,fuelbr,fuelxt; //4
	UREAL burnel,burnmo,burnbr,smax,voxden,ebefor; //6
	int32 matid,event,modeel; //3
}Tdpmjmp;

/*
 *  p4-/cgo/
 * 	cgo used in 
 * 	Function inters() / chvox()
 * 	PROPERTY：Private
 * 
 *  记录电子发生作用的时移动穿越哪个截面
 *  index=1/2/3 表示xyz上个方向
 *  dvox=1/-1 表示向前还是向后移动
 */ 
typedef struct
{
	int32 index,dvox; //2
}Tcgo;

/**********************************************************************
 * 
 *  Varibles for Debugging.
 *
 **********************************************************************/
/*
 * 	callcnt used in 
 * 	
 * 	PROPERTY：Shared
 */
typedef struct
{
	int32 cnt1, cnt2, cnt3;
}Tcallcnt;

/**********************************************************************
 *	
 *				Data Declaration
 *
 **********************************************************************/
/*
 * Shared Data among all threads
 */
/* dpmpart */
#define energy ((*dpmpart).energy)
#define vx ((*dpmpart).vx)
#define vy ((*dpmpart).vy)
#define vz ((*dpmpart).vz)
#define x ((*dpmpart).x)
#define y ((*dpmpart).y)
#define z ((*dpmpart).z)
#define xvox ((*dpmpart).xvox)
#define yvox ((*dpmpart).yvox)
#define zvox ((*dpmpart).zvox)
#define absvox ((*dpmpart).absvox)
#define nhist ((*dpmpart).nnhist)
#define ptype ((*dpmpart).pptype)

/* dpmsrc */
#define _esrc ((*dpmsrc).esrc)
#define _eabs ((*dpmsrc).eabs)
#define _eabsph ((*dpmsrc).eabsph)
#define _param ((*dpmsrc).param)
#define _intype ((*dpmsrc).intype)

/* dpmvox */
// /dpmvox/ ddx,ddy,ddz,xmid,ymid,zmid
#define _ddx ((*dpmvox).ddx)
#define _ddy ((*dpmvox).ddy)
#define _ddz ((*dpmvox).ddz)
#define _xmid ((*dpmvox).xmid)
#define _ymid ((*dpmvox).ymid)
#define _zmid ((*dpmvox).zmid)
#define _dens(i) ((*dpmvox).dens[(i)-1])
#define _mat(i) ((*dpmvox).mat[(i)-1])
#define _Unxvox ((*dpmvox).Unxvox)
#define _Unyvox ((*dpmvox).Unyvox)
#define _Unzvox ((*dpmvox).Unzvox)

/* cgeom3 */
#define dx ((*cgeom3).dxx)
#define dy ((*cgeom3).dyy)
#define dz ((*cgeom3).dzz)
#define idx ((*cgeom3).idxx)
#define idy ((*cgeom3).idyy)
#define idz ((*cgeom3).idzz)
#define ivx ((*cgeom3).ivx)
#define ivy ((*cgeom3).ivy)
#define ivz ((*cgeom3).ivz)

/* rdpmstpw */
#define _est ((*dpmstpw).est)
#define _stspa(i,j) ((*dpmstpw).stspx[((j)-1)*maxmat*2+((i)-1)*2])
#define _stspb(i,j) ((*dpmstpw).stspx[((j)-1)*maxmat*2+((i)-1)*2+1])
#define _idlest ((*dpmstpw).idlest)

/* rdpmscpw */
#define _escp ((*dpmscpw).escp)
#define _scpspa(i,j) ((*dpmscpw).scpspx[((j)-1)*maxmat*2+((i)-1)*2])
#define _scpspb(i,j) ((*dpmscpw).scpspx[((j)-1)*maxmat*2+((i)-1)*2+1])
#define _idlesc ((*dpmscpw).idlesc)

/* rdpmlbr */
#define _elra ((*dpmlbr).elra)
#define _lraspa(i,j) ((*dpmlbr).lraspx[((j)-1)*maxmat*2+((i)-1)*2])
#define _lraspb(i,j) ((*dpmlbr).lraspx[((j)-1)*maxmat*2+((i)-1)*2+1])
#define _idlelr ((*dpmlbr).idlelr)

/* dpmrpt */
#define nxini ((*dpmrpt).nxini)
#define nxfin ((*dpmrpt).nxfin)
#define nyini ((*dpmrpt).nyini)
#define nyfin ((*dpmrpt).nyfin)
#define nzini ((*dpmrpt).nzini)
#define nzfin ((*dpmrpt).nzfin)

/* dpmesc */
#define escore(i) ((*dpmesc).escore[(i)-1])
#define escor2(i) ((*dpmesc).escor2[(i)-1])
#define etmp(i) ((*dpmesc).etmp[(i)-1])
#define lasthi(i) ((*dpmesc).lasthi[(i)-1])

/* rdpmesc */ //p-02
#define _escore(i) ((*dpmesc).escorex[((i)-1)*2])
#define _escor2(i) ((*dpmesc).escorex[((i)-1)*2+1])

/* dpmesc_new */
//本索引从0开始，因为tail从0开始
#define estmp(i) ((*dpmesc_new).estmp[i])
#define esidx(i) ((*dpmesc_new).esidx[i])
#define tail ((*dpmesc_new).tail)

/* dpmsub */
#define _subden ((*dpmsub).subden)
#define _subfac ((*dpmsub).subfac)
#define _substp ((*dpmsub).substp)

/* dpmmat */
#define _matden(i) ((*dpmmat).matden[(i)-1])
#define _zmass(i) ((*dpmmat).zmass[(i)-1])
#define _z2mass(i) ((*dpmmat).z2mass[(i)-1])
#define _wcc ((*dpmmat).wcc)
#define _wcb ((*dpmmat).wcb)
#define _nmat ((*dpmmat).nmat)

/* dpmjmp */
#define fuelel ((*dpmjmp).fuelel)
#define fuelmo ((*dpmjmp).fuelmo)
#define fuelbr ((*dpmjmp).fuelbr)
#define fuelxt ((*dpmjmp).fuelxt)
#define burnel ((*dpmjmp).burnel)
#define burnmo ((*dpmjmp).burnmo)
#define burnbr ((*dpmjmp).burnbr)
#define smax ((*dpmjmp).smax)
#define voxden ((*dpmjmp).voxden)
#define ebefor ((*dpmjmp).ebefor)
#define matid ((*dpmjmp).matid)
#define event ((*dpmjmp).event)
#define modeel ((*dpmjmp).modeel)

/* dpmsim */
#define _atime ((*dpmsim).atime)
#define _maxhis ((*dpmsim).maxhis)

/* rdpmscsr */
#define _escsr ((*dpmscsr).escsr)
#define _scsspa(i) ((*dpmscsr).scsspx[((i)-1)*4])
#define _scsspb(i) ((*dpmscsr).scsspx[((i)-1)*4+1])
#define _scsspc(i) ((*dpmscsr).scsspx[((i)-1)*4+2])
#define _scsspd(i) ((*dpmscsr).scsspx[((i)-1)*4+3])
#define _idless ((*dpmscsr).idless)

/* rdpmlam */
#define _elam ((*dpmlam).elam)
#define _lamspa(i) ((*dpmlam).lamspx[((i)-1)*4])
#define _lamspb(i) ((*dpmlam).lamspx[((i)-1)*4+1])
#define _lamspc(i) ((*dpmlam).lamspx[((i)-1)*4+2])
#define _lamspd(i) ((*dpmlam).lamspx[((i)-1)*4+3])
#define _idlela ((*dpmlam).idlela)

/* rxdpmbw */
#define _ebw ((*xdpmbw).ebw)
#define _bwspa(i) ((*xdpmbw).bwspx[((i)-1)*4])
#define _bwspb(i) ((*xdpmbw).bwspx[((i)-1)*4+1])
#define _bwspc(i) ((*xdpmbw).bwspx[((i)-1)*4+2])
#define _bwspd(i) ((*xdpmbw).bwspx[((i)-1)*4+3])
#define _idlebw ((*xdpmbw).idlebw)

/* xdpmq2D */
#define q(i,j) ((*xdpmq2D).q[(i)-1+((j)-1)*nuq])
#define le0q ((*xdpmq2D).le0q)
#define idleq ((*xdpmq2D).idleq)
#define iduq ((*xdpmq2D).iduq)

/* cutoff */
#define _wcion ((*cutoff).wcion)
#define _wcbre ((*cutoff).wcbre)

/* rdpmstck */
#define _swght(i) ((*dpmstck).fltdata[((i)-1)*8])
#define _senerg(i) ((*dpmstck).fltdata[((i)-1)*8+1])
#define _svx(i) ((*dpmstck).fltdata[((i)-1)*8+2])
#define _svy(i) ((*dpmstck).fltdata[((i)-1)*8+3])
#define _svz(i) ((*dpmstck).fltdata[((i)-1)*8+4])
#define _sx(i) ((*dpmstck).fltdata[((i)-1)*8+5])
#define _sy(i) ((*dpmstck).fltdata[((i)-1)*8+6])
#define _sz(i) ((*dpmstck).fltdata[((i)-1)*8+7])

#define _sptype(i) ((*dpmstck).intdata[((i)-1)*5])
#define _sxvox(i) ((*dpmstck).intdata[((i)-1)*5+1])
#define _syvox(i) ((*dpmstck).intdata[((i)-1)*5+2])
#define _szvox(i) ((*dpmstck).intdata[((i)-1)*5+3])
#define _sabsvx(i) ((*dpmstck).intdata[((i)-1)*5+4])
#define _nsec ((*dpmstck).nsec)

/* xbre */
#define f0(i,j) ((*xbre).f0[(i)-1+((j)-1)*maxmat])
#define bcb(i) ((*xbre).bcb[(i)-1])

/* rdpmwck */
#define _a0wck(i) ((*dpmwck).awck[((i)-1)*2])
#define _a1wck(i) ((*dpmwck).awck[((i)-1)*2+1])
#define _idlewk ((*dpmwck).idlewk)
#define _wcke0 ((*dpmwck).wcke0)

/* rdpmlph */
#define _elaph ((*dpmlph).elaph)
#define _lampha(i,j) ((*dpmlph).lamphx[((j)-1)*maxmat*4+((i)-1)*4])
#define _lamphb(i,j) ((*dpmlph).lamphx[((j)-1)*maxmat*4+((i)-1)*4+1])
#define _lamphc(i,j) ((*dpmlph).lamphx[((j)-1)*maxmat*4+((i)-1)*4+2])
#define _lamphd(i,j) ((*dpmlph).lamphx[((j)-1)*maxmat*4+((i)-1)*4+3])
#define _idleph ((*dpmlph).idleph)

/* rdpmcmp */
#define _ecmpt ((*dpmcmp).ecmpt)
#define _compta(i,j) ((*dpmcmp).comptx[((j)-1)*maxmat*4+((i)-1)*4])
#define _comptb(i,j) ((*dpmcmp).comptx[((j)-1)*maxmat*4+((i)-1)*4+1])
#define _comptc(i,j) ((*dpmcmp).comptx[((j)-1)*maxmat*4+((i)-1)*4+2])
#define _comptd(i,j) ((*dpmcmp).comptx[((j)-1)*maxmat*4+((i)-1)*4+3])
#define _idlecp ((*dpmcmp).idlecp)

/* rdpmpap */
#define _epair ((*dpmpap).epair)
#define _pairpa(i,j) ((*dpmpap).pairpx[((j)-1)*maxmat*4+((i)-1)*4])
#define _pairpb(i,j) ((*dpmpap).pairpx[((j)-1)*maxmat*4+((i)-1)*4+1])
#define _pairpc(i,j) ((*dpmpap).pairpx[((j)-1)*maxmat*4+((i)-1)*4+2])
#define _pairpd(i,j) ((*dpmpap).pairpx[((j)-1)*maxmat*4+((i)-1)*4+3])
#define _idlepp ((*dpmpap).idlepp)

/* rseed */
#define ISEED1 ((*rseed).IISEED1)
#define ISEED2 ((*rseed).IISEED2)
#define ICNT ((*rseed).IICNT)

/* cgo */
#define index ((*cgo).index)
#define dvox ((*cgo).dvox)

/* callcnt */
#define cnt1 ((*callcnt).cnt1)
#define cnt2 ((*callcnt).cnt2)
#define cnt3 ((*callcnt).cnt3)

#endif
